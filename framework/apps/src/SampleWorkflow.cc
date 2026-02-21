/* -- C++ -- */
/**
 *  @file  framework/apps/src/SampleWorkflow.cc
 *
 *  @brief Sample aggregation workflow (invoked by the unified heron CLI).
 */

#include <chrono>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "AppUtils.hh"
#include "SampleCLI.hh"
#include "StatusMonitor.hh"

int run(const SampleArgs &sample_args, const std::string &log_prefix)
{
    const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    const auto files = read_paths(sample_args.filelist_path);

    std::filesystem::path output_path(sample_args.output_path);
    if (!output_path.parent_path().empty())
    {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::filesystem::path sample_list_path(sample_args.sample_list_path);
    if (!sample_list_path.parent_path().empty())
    {
        std::filesystem::create_directories(sample_list_path.parent_path());
    }

    const auto start_time = std::chrono::steady_clock::now();
    log_sample_start(log_prefix, files.size());

    StatusMonitor status_monitor(
        log_prefix,
        "action=sample_build status=running message=processing");
    SampleIO::Sample sample =
        NormalisationService::build_sample(sample_args.sample_name,
                                           files,
                                           db_path);

    status_monitor.stop();

    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    log_sample_finish(log_prefix, sample.inputs.size(), elapsed_seconds);

    SampleIO::write(sample, sample_args.output_path);
    update_sample_list(sample_args.sample_list_path, sample, sample_args.output_path);

    std::ostringstream log_message;
    log_message << "action=sample_write status=complete sample=" << sample.sample_name
                << " inputs=" << sample.inputs.size()
                << " pot_sum=" << sample.subrun_pot_sum
                << " db_tortgt_pot_sum=" << sample.db_tortgt_pot_sum
                << " normalisation=" << sample.normalisation
                << " normalised_pot_sum=" << sample.normalised_pot_sum
                << " output=" << sample_args.output_path
                << " sample_list=" << sample_args.sample_list_path;
    log_success(log_prefix, log_message.str());

    return 0;
}
