/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecSystMaker.cc
 *
 *  @brief Build systematic template variations for existing templates.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AppUtils.hh"
#include "SystematicsBuilder.hh"
#include "SystematicsSpec.hh"
#include "TemplateSpec.hh"

namespace
{

std::vector<nuxsec::SampleListEntry> read_sample_list(const std::string &list_path)
{
    const auto app_entries = nuxsec::app::read_sample_list(list_path);
    std::vector<nuxsec::SampleListEntry> entries;
    entries.reserve(app_entries.size());
    for (const auto &entry : app_entries)
    {
        nuxsec::SampleListEntry out;
        out.sample_name = entry.sample_name;
        out.sample_kind = entry.sample_kind;
        out.beam_mode = entry.beam_mode;
        out.output_path = entry.output_path;
        entries.push_back(std::move(out));
    }

    return entries;
}

struct Args
{
    std::string list_path;
    std::string tree_name;
    std::string template_root;
    int nthreads = 1;
};

Args parse_args(int argc, char **argv)
{
    if (argc < 4 || argc > 5)
    {
        throw std::runtime_error("Usage: nuxsecSystMaker SAMPLE_LIST.tsv TREE_NAME TEMPLATES.root [NTHREADS]");
    }

    Args args;
    args.list_path = nuxsec::app::trim(argv[1]);
    args.tree_name = nuxsec::app::trim(argv[2]);
    args.template_root = nuxsec::app::trim(argv[3]);
    if (argc == 5)
    {
        args.nthreads = std::stoi(nuxsec::app::trim(argv[4]));
    }

    if (args.list_path.empty() || args.tree_name.empty() || args.template_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path)");
    }
    if (args.nthreads < 1)
    {
        throw std::runtime_error("Invalid thread count (must be >= 1)");
    }

    return args;
}

} // namespace

int main(int argc, char **argv)
{
    return nuxsec::app::run_with_exceptions(
        [argc, argv]()
        {
            const Args args = parse_args(argc, argv);

            const std::vector<nuxsec::TemplateSpec1D> fit_specs =
                nuxsec::read_template_spec_1d_tsv("inputs/fit_templates.tsv");

            const auto entries = read_sample_list(args.list_path);

            nuxsec::SystematicsBuilder::Options opt;
            opt.nthreads = args.nthreads;

            nuxsec::SystematicsBuilder::BuildAll(entries,
                                                 args.tree_name,
                                                 fit_specs,
                                                 args.template_root,
                                                 nuxsec::DefaultSystematics(),
                                                 opt);

            std::cerr << "[nuxsecSystMaker] wrote systematics into " << args.template_root << "\n";
            return 0;
        });
}
