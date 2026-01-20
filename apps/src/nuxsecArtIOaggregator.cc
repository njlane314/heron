/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecArtIOaggregator.cc
 *
 *  @brief Main entrypoint for Art file provenance generation.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppCommandHelpers.hh"

int main(int argc, char **argv)
{
    try
    {
        std::vector<std::string> args;
        args.reserve(static_cast<size_t>(argc > 0 ? argc - 1 : 0));
        for (int i = 1; i < argc; ++i)
        {
            args.emplace_back(argv[i]);
        }

        const nuxsec::app::ArtArgs art_args =
            nuxsec::app::parse_art_args(args, "Usage: nuxsecArtIOaggregator NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]");
        return nuxsec::app::run_artio(art_args, "nuxsecArtIOaggregator");
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
