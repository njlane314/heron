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
    return nuxsec::app::run_with_exceptions(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::SampleArgs sample_args =
                nuxsec::app::parse_sample_args(args, "Usage: nuxsecSampleIOdriver NAME:FILELIST");
            return nuxsec::app::run_sample(sample_args, "nuxsecSampleIOdriver");
        });
}
