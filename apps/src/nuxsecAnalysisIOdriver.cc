/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecAnalysisIOdriver.cc
 *
 *  @brief Build analysis output from aggregated samples.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AnalysisCommand.hh"

int main(int argc, char **argv)
{
    return nuxsec::app::run_with_exceptions(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::AnalysisArgs analysis_args = nuxsec::app::parse_analysis_args(
                args, "Usage: nuxsecAnalysisIOdriver SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]");
            return nuxsec::app::run_analysis(analysis_args, "nuxsecAnalysisIOdriver");
        });
}
