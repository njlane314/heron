/* -- C++ -- */
/**
 *  @file  io/src/EventIO.cc
 *
 *  @brief Implementation of event-level IO.
 */

#include "EventIO.hh"

#include <cctype>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <system_error>
#include <vector>

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RSnapshotOptions.hxx>

#include <TFile.h>
#include <TObjString.h>
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

    ROOT::RDF::RSnapshotOptions options;
    options.fMode = "UPDATE";
    options.fOverwriteIfExists = overwrite_if_exists;
    options.fLazy = true;
    options.fAutoFlush = 100000;

    auto count = filtered.Count();
    constexpr ULong64_t progress_every = 1000;
    const auto start_time = std::chrono::steady_clock::now();
    count.OnPartialResult(progress_every,
                          [sample_name](ULong64_t processed)
                          {
                              std::cerr << "[EventIO] stage=snapshot_progress"
                                        << " sample=" << sample_name
                                        << " processed=" << processed
                                        << "\n";
                          });
    auto snapshot = filtered.Snapshot(tree_name, m_path, columns, options);
    ROOT::RDF::RunGraphs({count, snapshot});
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
