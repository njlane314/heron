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

        const nuxsec::app::TemplateArgs tpl_args = nuxsec::app::parse_template_args(
            args, "Usage: nuxsecTemplateMaker SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]");
        return nuxsec::app::run_template(tpl_args, "nuxsecTemplateMaker");
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
