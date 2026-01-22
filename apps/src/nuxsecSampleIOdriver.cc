/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecSampleIOdriver.cc
 *
 *  @brief Main entrypoint for Sample aggregation.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "SampleCommand.hh"

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::sample::Args sample_args =
                nuxsec::app::sample::parse_args(
                    args, "Usage: nuxsecSampleIOdriver NAME:FILELIST");
            return nuxsec::app::sample::run(sample_args, "nuxsecSampleIOdriver");
        });
}
