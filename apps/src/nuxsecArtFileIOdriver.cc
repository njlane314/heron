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
    return nuxsec::app::run_with_exceptions(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::ArtArgs art_args =
                nuxsec::app::parse_art_args(args,
                                           "Usage: nuxsecArtFileIOdriver NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]");
            return nuxsec::app::run_art(art_args, "nuxsecArtFileIOdriver");
        });
}
