/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecEventIOdriver.cc
 *
 *  @brief Build event-level output from aggregated samples.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "EventCommand.hh"

int main(int argc, char **argv)
{
    return nuxsec::app::run_with_exceptions(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::EventIOArgs event_args = nuxsec::app::parse_eventio_args(
                args, "Usage: nuxsecEventIOdriver SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]");
            return nuxsec::app::run_eventio(event_args, "nuxsecEventIOdriver");
        });
}
