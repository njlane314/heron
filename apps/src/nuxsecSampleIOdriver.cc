/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecSampleIOdriver.cc
 *
 *  @brief Main entrypoint for Sample aggregation.
 */

#include "SampleCLI.hh"

#include <sstream>
#include <string>
#include <vector>

namespace nuxsec
{

namespace app
{

int run(const sample::Args &sample_args, const std::string &log_prefix)
{
    const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    const auto files = nuxsec::app::read_paths(sample_args.filelist_path);

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
    nuxsec::app::sample::log_sample_start(log_prefix, files.size());
    
    nuxsec::sample::SampleIO::Sample sample = nuxsec::NormalisationService::build_sample(sample_args.sample_name, files, db_path);
    
    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    nuxsec::app::sample::log_sample_finish(log_prefix, sample.inputs.size(), elapsed_seconds);
    
    nuxsec::sample::SampleIO::write(sample, sample_args.output_path);
    nuxsec::app::sample::update_sample_list(sample_args.sample_list_path, sample, sample_args.output_path);

    std::ostringstream log_message;
    log_message << "action=sample_write status=complete sample=" << sample.sample_name
                << " inputs=" << sample.inputs.size()
                << " pot_sum=" << sample.subrun_pot_sum
                << " db_tortgt_pot_sum=" << sample.db_tortgt_pot_sum
                << " normalisation=" << sample.normalisation
                << " normalised_pot_sum=" << sample.normalised_pot_sum
                << " output=" << sample_args.output_path
                << " sample_list=" << sample_args.sample_list_path;
    nuxsec::app::log::log_success(log_prefix, log_message.str());

    return 0;
}

}

}

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        "nuxsecSampleIOdriver",
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::sample::Args sample_args =
                nuxsec::app::sample::parse_args(
                    args, "Usage: nuxsecSampleIOdriver NAME:FILELIST");
            return nuxsec::app::run(sample_args, "nuxsecSampleIOdriver");
        });
}
