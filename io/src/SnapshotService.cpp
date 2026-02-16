/* -- C++ -- */
/**
 *  @file  io/src/SnapshotService.cpp
 *
 *  @brief Implementation of snapshot helpers for event-level output.
 */

#include "SnapshotService.hh"

#include <algorithm>
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
#include <TObject.h>
#include <TTree.h>


std::string SnapshotService::sanitise_root_key(std::string s)
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

namespace
{
void append_tree_fast(const std::string &out_path,
                      const std::string &scratch_file,
                      const std::string &tree_name)
{
    std::unique_ptr<TFile> fin(TFile::Open(scratch_file.c_str(), "READ"));
    if (!fin || fin->IsZombie())
        throw std::runtime_error("SnapshotService: failed to open scratch snapshot file: " + scratch_file);

    TTree *tin = dynamic_cast<TTree *>(fin->Get(tree_name.c_str()));
    if (!tin)
        throw std::runtime_error("SnapshotService: scratch snapshot missing tree: " + tree_name + " in " + scratch_file);

    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
        throw std::runtime_error("SnapshotService: failed to open output for append: " + out_path);

    TTree *tout = dynamic_cast<TTree *>(fout->Get(tree_name.c_str()));
    fout->cd();

    if (!tout)
    {
        std::unique_ptr<TTree> cloned(tin->CloneTree(-1, "fast"));
        cloned->SetName(tree_name.c_str());
        cloned->Write(tree_name.c_str(), TObject::kOverwrite);
    }
    else
    {
        tout->SetDirectory(fout.get());
        const Long64_t n = tout->CopyEntries(tin, -1, "fast");
        (void)n;
        tout->Write("", TObject::kOverwrite);
    }

    fout->Close();
    fin->Close();
}

std::filesystem::path snapshot_scratch_dir()
{
    const char *user = std::getenv("USER");
    if (!user || !*user)
        throw std::runtime_error("SnapshotService: USER must be set to resolve snapshot staging directory");

    return std::filesystem::path("/exp/uboone/data/users") / user / "staging";
}
} // namespace

ULong64_t SnapshotService::snapshot_event_list_merged(ROOT::RDF::RNode node,
                                                      const std::string &out_path,
                                                      int sample_id,
                                                      const std::string &sample_name,
                                                      const std::vector<std::string> &columns,
                                                      const std::string &selection,
                                                      const std::string &tree_name_in)
{
    ROOT::RDF::RNode filtered = std::move(node);
    if (!selection.empty() && selection != "true")
        filtered = filtered.Filter(selection, "eventio_selection");

    const std::string tree_name = sanitise_root_key(tree_name_in.empty() ? "events" : tree_name_in);

    filtered = filtered.Define("sample_id", [sample_id]() { return sample_id; });

    std::vector<std::string> snapshot_cols = columns;
    if (std::find(snapshot_cols.begin(), snapshot_cols.end(), "sample_id") == snapshot_cols.end())
        snapshot_cols.push_back("sample_id");

    std::filesystem::path scratch_dir = snapshot_scratch_dir();
    {
        std::error_code ec;
        std::filesystem::create_directories(scratch_dir, ec);
        if (ec)
        {
            throw std::runtime_error(
                "SnapshotService: failed to create scratch directory: "
                + scratch_dir.string()
                + " (" + ec.message() + ")");
        }
    }

    const std::string scratch_file =
        (scratch_dir / ("heron_snapshot_" + tree_name + "_" + sanitise_root_key(sample_name) + "_"
                        + std::to_string(::getpid()) + ".root"))
            .string();

    ROOT::RDF::RSnapshotOptions options;
    options.fMode = "RECREATE";
    options.fOverwriteIfExists = false;
    options.fLazy = true;
    options.fCompressionAlgorithm = ROOT::kLZ4;
    options.fCompressionLevel = 1;
    options.fAutoFlush = -50LL * 1024 * 1024;
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
                              std::cerr << "[SnapshotService] stage=snapshot_progress"
                                        << " sample=" << sample_name
                                        << " processed=" << processed
                                        << " elapsed_seconds=" << elapsed_seconds
                                        << "\n";
                          });

    auto snapshot = filtered.Snapshot(tree_name, scratch_file, snapshot_cols, options);
    std::cerr << "[SnapshotService] stage=snapshot_run"
              << " sample=" << sample_name
              << " scratch_file=" << scratch_file
              << "\n";
    (void)snapshot.GetValue();

    std::cerr << "[SnapshotService] stage=append_begin"
              << " sample=" << sample_name
              << " scratch_file=" << scratch_file
              << " out_file=" << out_path
              << " tree=" << tree_name
              << "\n";
    append_tree_fast(out_path, scratch_file, tree_name);
    std::cerr << "[SnapshotService] stage=append_done sample=" << sample_name << "\n";

    {
        std::error_code ec;
        std::filesystem::remove(scratch_file, ec);
        if (ec)
            std::cerr << "[SnapshotService] warning=failed_to_remove_scratch_file path=" << scratch_file
                      << " err=" << ec.message() << "\n";
    }

    return count.GetValue();
}

ULong64_t SnapshotService::snapshot_event_list(ROOT::RDF::RNode node,
                                               const std::string &out_path,
                                               const std::string &sample_name,
                                               const std::vector<std::string> &columns,
                                               const std::string &selection,
                                               const std::string &tree_prefix,
                                               bool overwrite_if_exists)
{
    ROOT::RDF::RNode filtered = std::move(node);
    if (!selection.empty() && selection != "true")
    {
        filtered = filtered.Filter(selection, "eventio_selection");
    }

    const std::string tree_name =
        sanitise_root_key(tree_prefix.empty() ? "events" : tree_prefix) + "_" + sanitise_root_key(sample_name);
    if (!overwrite_if_exists)
    {
        std::unique_ptr<TFile> check_file(TFile::Open(out_path.c_str(), "READ"));
        if (check_file && check_file->Get(tree_name.c_str()))
        {
            std::cerr << "[SnapshotService] stage=snapshot_skip_existing"
                      << " sample=" << sample_name
                      << " tree=" << tree_name
                      << " output=" << out_path
                      << "\n";
            return 0;
        }
    }

    std::filesystem::path scratch_dir = snapshot_scratch_dir();
    {
        std::error_code ec;
        std::filesystem::create_directories(scratch_dir, ec);
        if (ec)
        {
            throw std::runtime_error(
                "SnapshotService: failed to create scratch directory: "
                + scratch_dir.string()
                + " (" + ec.message() + ")");
        }
    }

    const std::string scratch_file =
        (scratch_dir / ("heron_snapshot_" + tree_name + "_" + std::to_string(::getpid()) + ".root")).string();

    ROOT::RDF::RSnapshotOptions options;
    options.fMode = "RECREATE";
    options.fOverwriteIfExists = false;
    options.fLazy = true;
    options.fCompressionAlgorithm = ROOT::kLZ4;
    options.fCompressionLevel = 1;
    options.fAutoFlush = -50LL * 1024 * 1024;
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
                              std::cerr << "[SnapshotService] stage=snapshot_progress"
                                        << " sample=" << sample_name
                                        << " processed=" << processed
                                        << " elapsed_seconds=" << elapsed_seconds
                                        << "\n";
                          });
    auto snapshot = filtered.Snapshot(tree_name, scratch_file, columns, options);
    std::cerr << "[SnapshotService] stage=snapshot_run"
              << " sample=" << sample_name
              << " scratch_file=" << scratch_file
              << "\n";
    ROOT::RDF::RunGraphs({count, snapshot});
    (void)snapshot.GetValue();

    std::cerr << "[SnapshotService] stage=snapshot_merge_begin"
              << " sample=" << sample_name
              << " scratch_file=" << scratch_file
              << " out_file=" << out_path
              << "\n";
    {
        TFileMerger merger(kTRUE, kFALSE);
        merger.SetFastMethod(kFALSE);
        if (!merger.OutputFile(out_path.c_str(), "UPDATE"))
            throw std::runtime_error("SnapshotService: failed to open output for merge: " + out_path);
        if (!merger.AddFile(scratch_file.c_str()))
            throw std::runtime_error("SnapshotService: failed to add scratch file to merger: " + scratch_file);
        if (!merger.Merge())
            throw std::runtime_error("SnapshotService: merge failed for scratch file: " + scratch_file);
    }
    std::cerr << "[SnapshotService] stage=snapshot_merge_done"
              << " sample=" << sample_name
              << "\n";

    {
        std::error_code ec;
        std::filesystem::remove(scratch_file, ec);
        if (ec)
            std::cerr << "[SnapshotService] warning=failed_to_remove_scratch_file path=" << scratch_file
                      << " err=" << ec.message() << "\n";
    }

    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    std::cerr << "[SnapshotService] stage=snapshot_complete"
              << " sample=" << sample_name
              << " processed=" << count.GetValue()
              << " elapsed_seconds=" << elapsed_seconds
              << "\n";
    return count.GetValue();
}
