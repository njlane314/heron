/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecArtFileIOdriver.cc
 *
 *  @brief Main entrypoint for Art file provenance generation.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ArtCommand.hh"

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::art::Args art_args =
                nuxsec::app::art::parse_args(
                    args,
                    "Usage: nuxsecArtFileIOdriver INPUT_NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]");
            return nuxsec::app::art::run(art_args, "nuxsecArtFileIOdriver");
        });
}
