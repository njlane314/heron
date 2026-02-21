/* -- C++ -- */
/**
 *  @file  framework/core/src/EventCLI.cc
 *
 *  @brief Non-inline EventCLI helpers for argument parsing and ROOT input checks.
 */

#include "EventCLI.hh"

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "AppUtils.hh"

void ensure_tree_present(const SampleIO::Sample &sample,
                         const std::string &tree_name)
{
    if (sample.inputs.empty())
    {
        throw std::runtime_error("Event inputs missing ROOT files for sample: " + sample.sample_name);
    }

    std::vector<std::string> files = SampleIO::resolve_root_files(sample);
    if (files.empty())
    {
        throw std::runtime_error("Event inputs missing ROOT files for sample: " + sample.sample_name);
    }

    for (const auto &path : files)
    {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie())
        {
            throw std::runtime_error("Event input failed to open ROOT file: " + path);
        }

        TTree *tree = nullptr;
        f->GetObject(tree_name.c_str(), tree);
        if (!tree)
        {
            throw std::runtime_error(
                "Event input missing tree '" + tree_name + "' in " + path);
        }
    }
}

EventArgs parse_event_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() < 2 || args.size() > 4)
    {
        throw std::runtime_error(usage);
    }

    EventArgs out;
    out.list_path = trim(args.at(0));
    out.output_root = trim(args.at(1));
    if (args.size() > 2)
    {
        out.selection = trim(args.at(2));
    }
    if (args.size() > 3)
    {
        out.columns_tsv_path = trim(args.at(3));
    }

    if (out.list_path.empty() || out.output_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path)");
    }

    std::filesystem::path output_root(out.output_root);
    if (output_root.is_relative() && output_root.parent_path().empty())
    {
        const std::filesystem::path event_dir =
            stage_output_dir("HERON_EVENT_DIR", "event");
        out.output_root = (event_dir / output_root).string();
    }

    return out;
}
