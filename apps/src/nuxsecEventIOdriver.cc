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

#include "EventDriver.hh"

int main(int argc, char **argv)
{
    return nuxsec::app::run_guarded(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::event::Args event_args =
                nuxsec::app::event::parse_args(
                    args, "Usage: nuxsecEventIOdriver SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]");
            return nuxsec::app::event::run(event_args, "nuxsecEventIOdriver");
        });
}
