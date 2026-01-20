/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecSampleIOaggregator.cc
 *
 *  @brief Main entrypoint for Sample aggregation.
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

        const nuxsec::app::SampleArgs sample_args =
            nuxsec::app::parse_sample_args(args, "Usage: nuxsecSampleIOaggregator NAME:FILELIST");
        return nuxsec::app::run_sample(sample_args, "nuxsecSampleIOaggregator");
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
