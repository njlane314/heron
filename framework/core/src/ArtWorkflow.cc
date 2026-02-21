/* -- C++ -- */
/**
 *  @file  framework/core/src/ArtWorkflow.cc
 *
 *  @brief Provenance generation workflow (invoked by the unified heron CLI).
 */

#include <chrono>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx> // ROOT::EnableImplicitMT

#include "AppUtils.hh"
#include "ArtCLI.hh"
#include "StatusMonitor.hh"

int run(const ArtArgs &art_args, const std::string &log_prefix)
{
    ROOT::EnableImplicitMT();

    std::filesystem::path out_path(art_args.art_path);
    if (!out_path.parent_path().empty())
    {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const auto files = read_paths(art_args.input.filelist_path);

    Provenance rec;
    rec.input = art_args.input;
    rec.input_files = files;
    rec.kind = art_args.sample_origin;
    rec.beam = art_args.beam_mode;

    if (rec.kind == SampleIO::SampleOrigin::kUnknown &&
        is_selection_data_file(files.front()))
    {
        rec.kind = SampleIO::SampleOrigin::kData;
    }

    const auto start_time = std::chrono::steady_clock::now();
    log_scan_start(log_prefix);

    StatusMonitor status_monitor(
        log_prefix,
        "action=art_scan status=running message=scan_in_progress");
    rec.summary = SubRunInventoryService::scan_subruns(files);
    status_monitor.stop();

    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    log_scan_finish(log_prefix, rec.summary.n_entries, elapsed_seconds);


    std::ostringstream log_message;
    log_message << "action=input_register status=complete input=" << rec.input.input_name
                << " files=" << rec.input_files.size()
                << " pairs=" << rec.summary.unique_pairs.size()
                << " pot_sum=" << rec.summary.pot_sum;
    log_success(log_prefix, log_message.str());

    ArtFileProvenanceIO::write(rec, art_args.art_path);

    return 0;
}
