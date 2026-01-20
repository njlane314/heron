/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecTemplateMaker.cc
 *
 *  @brief Build binned template histograms from aggregated samples.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppTemplateCommand.hh"

int main(int argc, char **argv)
{
    return nuxsec::app::run_with_exceptions(
        [argc, argv]()
        {
            const std::vector<std::string> args = nuxsec::app::collect_args(argc, argv);
            const nuxsec::app::TemplateArgs tpl_args = nuxsec::app::parse_template_args(
                args, "Usage: nuxsecTemplateMaker SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]");
            return nuxsec::app::run_template(tpl_args, "nuxsecTemplateMaker");
        });
}
