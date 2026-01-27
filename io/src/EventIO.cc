/* -- C++ -- */
/**
 *  @file  io/src/EventIO.cc
 *
 *  @brief Implementation of event-level IO.
 */

#include "EventIO.hh"

#include <cctype>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <system_error>
#include <unistd.h>
#include <vector>

#include <Compression.h>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RSnapshotOptions.hxx>

#include <TFile.h>
#include <TFileMerger.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>

namespace nuxsec
{
namespace event
{
std::string EventIO::sanitise_root_key(std::string s)
{
    for (char &c : s)
    {
        const unsigned char u = static_cast<unsigned char>(c);
        if (!(std::isalnum(u) || c == '_'))
            c = '_';
    }
    if (s.empty())
        s = "sample";
    return s;
}

void EventIO::init(const std::string &out_path,
                   const Header &header,
                   const std::vector<SampleInfo> &sample_refs,
                   const std::string &event_schema_tsv,
                   const std::string &schema_tag)
{
    const std::filesystem::path output_path(out_path);
    if (!output_path.parent_path().empty())
    {
        std::error_code ec;
        std::filesystem::create_directories(output_path.parent_path(), ec);
        if (ec)
        {
            throw std::runtime_error(
                "EventIO::init: failed to create output directory: "
                + output_path.parent_path().string()
                + " (" + ec.message() + ")");
        }
    }

    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie())
        throw std::runtime_error("EventIO::init: failed to create output file: " + out_path);

    TObjString(header.analysis_name.c_str()).Write("analysis_name");
    TObjString(header.provenance_tree.c_str()).Write("provenance_tree");
    TObjString(header.event_tree.c_str()).Write("event_tree");
    TObjString(header.sample_list_source.c_str()).Write("sample_list_source");

    if (!event_schema_tsv.empty())
    {
        const std::string key = schema_tag.empty() ? "event_schema"
                                                   : ("event_schema_" + sanitise_root_key(schema_tag));
        TObjString(event_schema_tsv.c_str()).Write(key.c_str());
    }

    TTree tref("sample_refs", "Sample references (source + POT/triggers)");
    std::string sample_name;
    std::string sample_rootio_path;
    int sample_origin = -1;
    int beam_mode = -1;
    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    tref.Branch("sample_name", &sample_name);
    tref.Branch("sample_rootio_path", &sample_rootio_path);
    tref.Branch("sample_origin", &sample_origin);
    tref.Branch("beam_mode", &beam_mode);
    tref.Branch("subrun_pot_sum", &subrun_pot_sum);
    tref.Branch("db_tortgt_pot_sum", &db_tortgt_pot_sum);
    tref.Branch("db_tor101_pot_sum", &db_tor101_pot_sum);

    for (const auto &r : sample_refs)
    {
        sample_name = r.sample_name;
        sample_rootio_path = r.sample_rootio_path;
        sample_origin = r.sample_origin;
        beam_mode = r.beam_mode;
        subrun_pot_sum = r.subrun_pot_sum;
        db_tortgt_pot_sum = r.db_tortgt_pot_sum;
        db_tor101_pot_sum = r.db_tor101_pot_sum;
        tref.Fill();
    }

    tref.Write();
    fout->Close();
}

EventIO::EventIO(std::string out_path, OpenMode mode)
    : m_path(std::move(out_path)), m_mode(mode)
{
    const char *opt = (m_mode == OpenMode::kRead) ? "READ" : "UPDATE";
    std::unique_ptr<TFile> f(TFile::Open(m_path.c_str(), opt));
    if (!f || f->IsZombie())
        throw std::runtime_error("EventIO: failed to open " + m_path);
    f->Close();
}

std::string EventIO::sample_tree_name(const std::string &sample_name,
                                      const std::string &tree_prefix) const
{
    const std::string p = tree_prefix.empty() ? "events" : tree_prefix;
    return sanitise_root_key(p) + "_" + sanitise_root_key(sample_name);
}

ULong64_t EventIO::snapshot_event_list(ROOT::RDF::RNode node,
                                       const std::string &sample_name,
                                       const std::vector<std::string> &columns,
                                       const std::string &selection,
                                       const std::string &tree_prefix,
                                       bool overwrite_if_exists) const
{
    ROOT::RDF::RNode filtered = std::move(node);
    if (!selection.empty() && selection != "true")
    {
        filtered = filtered.Filter(selection, "eventio_selection");
    }

    const std::string tree_name = sample_tree_name(sample_name, tree_prefix);
    if (!overwrite_if_exists)
    {
        std::unique_ptr<TFile> check_file(TFile::Open(m_path.c_str(), "READ"));
        if (check_file && check_file->Get(tree_name.c_str()))
        {
            std::cerr << "[EventIO] stage=snapshot_skip_existing"
                      << " sample=" << sample_name
                      << " tree=" << tree_name
                      << " output=" << m_path
                      << "\n";
            return 0;
        }
    }

    // Snapshot to a fresh temp file (RECREATE) to avoid Snapshot+UPDATE corner-cases.
    // Then merge the temp file into the real output file with TFileMerger.
    std::filesystem::path tmp_dir;
    if (const char *p = std::getenv("NUXSEC_TMPDIR"); p && *p)
        tmp_dir = p;
    else if (const char *p = std::getenv("TMPDIR"); p && *p)
        tmp_dir = p;
    else
        tmp_dir = std::filesystem::temp_directory_path();

    const std::string tmp_file =
        (tmp_dir / ("nuxsec_snapshot_" + tree_name + "_" + std::to_string(::getpid()) + ".root")).string();

    ROOT::RDF::RSnapshotOptions options;
    options.fMode = "RECREATE";
    options.fOverwriteIfExists = false; // irrelevant for RECREATE
    options.fLazy = true;
    // Make snapshot writing less bursty and less CPU-heavy:
    //  - LZ4 is typically much faster than ZLIB for analysis ntuples
    //  - negative AutoFlush => flush by buffer size (bytes), smoothing long stalls
    options.fCompressionAlgorithm = ROOT::kLZ4;
    options.fCompressionLevel = 1;
    options.fAutoFlush = -50LL * 1024 * 1024; // ~50 MB
    // Avoid deep splitting (default is 99) which can explode branch count/cost for complex types.
    options.fSplitLevel = 0;

    auto count = filtered.Count();
    constexpr ULong64_t progress_every = 1000;
    const auto start_time = std::chrono::steady_clock::now();
    count.OnPartialResult(progress_every,
                          [sample_name, start_time](ULong64_t processed)
                          {
                              const auto now = std::chrono::steady_clock::now();
                              const double elapsed_seconds =
                                  std::chrono::duration_cast<std::chrono::duration<double>>(now - start_time).count();
                              std::cerr << "[EventIO] stage=snapshot_progress"
                                        << " sample=" << sample_name
                                        << " processed=" << processed
                                        << " elapsed_seconds=" << elapsed_seconds
                                        << "\n";
                          });
    auto snapshot = filtered.Snapshot(tree_name, tmp_file, columns, options);
    std::cerr << "[EventIO] stage=snapshot_run"
              << " sample=" << sample_name
              << " tmp_file=" << tmp_file
              << "\n";
    ROOT::RDF::RunGraphs({count, snapshot});

    std::cerr << "[EventIO] stage=snapshot_merge_begin"
              << " sample=" << sample_name
              << " tmp_file=" << tmp_file
              << " out_file=" << m_path
              << "\n";
    {
        TFileMerger merger(kTRUE, kTRUE);
        merger.SetFastMethod(kTRUE);
        if (!merger.OutputFile(m_path.c_str(), "UPDATE"))
            throw std::runtime_error("EventIO: failed to open output for merge: " + m_path);
        if (!merger.AddFile(tmp_file.c_str()))
            throw std::runtime_error("EventIO: failed to add temp file to merger: " + tmp_file);
        if (!merger.Merge())
            throw std::runtime_error("EventIO: merge failed for temp file: " + tmp_file);
    }
    std::cerr << "[EventIO] stage=snapshot_merge_done"
              << " sample=" << sample_name
              << "\n";

    {
        std::error_code ec;
        std::filesystem::remove(tmp_file, ec);
        if (ec)
            std::cerr << "[EventIO] warning=failed_to_remove_tmp_file path=" << tmp_file
                      << " err=" << ec.message() << "\n";
    }

    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    std::cerr << "[EventIO] stage=snapshot_complete"
              << " sample=" << sample_name
              << " processed=" << count.GetValue()
              << " elapsed_seconds=" << elapsed_seconds
              << "\n";
    return count.GetValue();
}

} // namespace event
} // namespace nuxsec
