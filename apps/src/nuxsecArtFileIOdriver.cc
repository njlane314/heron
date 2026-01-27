/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecArtFileIOdriver.cc
 *
 *  @brief Main entrypoint for Art file provenance generation.
 */

#include "ArtCLI.hh"

#include <sstream>
#include <string>
#include <vector>

namespace nuxsec
{

namespace app
{

int run(const art::Args &art_args, const std::string &log_prefix)
{
    ROOT::EnableImplicitMT();
    
    std::filesystem::path out_path(art_args.art_path);
    if (!out_path.parent_path().empty())
    {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const auto files = nuxsec::app::read_paths(art_args.input.filelist_path);

    nuxsec::art::Provenance rec;
    rec.input = art_args.input;
    rec.input_files = files;
    rec.kind = art_args.sample_origin;
    rec.beam = art_args.beam_mode;

    if (rec.kind == nuxsec::sample::SampleIO::SampleOrigin::kUnknown &&
        nuxsec::app::art::is_selection_data_file(files.front()))
    {
        rec.kind = nuxsec::sample::SampleIO::SampleOrigin::kData;
    }

    const auto start_time = std::chrono::steady_clock::now();
    nuxsec::app::art::log_scan_start(log_prefix);
    
    rec.summary = nuxsec::SubRunInventoryService::scan_subruns(files);
    
    const auto end_time = std::chrono::steady_clock::now();
    const double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    nuxsec::app::art::log_scan_finish(log_prefix, rec.summary.n_entries, elapsed_seconds);

    rec.summary.pot_sum *= 1;
    rec.scale = 1;

    std::ostringstream log_message;
    log_message << "action=input_register status=complete input=" << rec.input.input_name
                << " files=" << rec.input_files.size()
                << " pairs=" << rec.summary.unique_pairs.size()
                << " pot_sum=" << rec.summary.pot_sum;
    nuxsec::app::log::log_success(log_prefix, log_message.str());

    nuxsec::ArtFileProvenanceIO::write(rec, art_args.art_path);

    return 0;
}

}

}

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        "nuxsecArtFileIOdriver",
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::art::Args art_args =
                nuxsec::app::art::parse_args(
                    args,
                    "Usage: nuxsecArtFileIOdriver INPUT_NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]");
            return nuxsec::app::run(art_args, "nuxsecArtFileIOdriver");
        });
}
