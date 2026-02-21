/* -- C++ -- */
/**
 *  @file  apps/include/EventCLI.hh
 *
 *  @brief CLI helpers that drive event-level workflows, including configuration,
 *         input selection, and summary output for analysis-ready datasets.
 */
#ifndef HERON_APPS_EVENTCLI_H
#define HERON_APPS_EVENTCLI_H

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>

#include "AppLog.hh"
#include "AppUtils.hh"
#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"
#include "EventListIO.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"




inline void log_event_start(const std::string &log_prefix, const size_t sample_count)
{
    log_info(
        log_prefix,
        "action=event_build status=start samples=" +
            format_count(static_cast<long long>(sample_count)));
}

inline void log_event_finish(const std::string &log_prefix,
                             const size_t sample_count,
                             const double elapsed_seconds)
{
    std::ostringstream out;
    out << "action=event_build status=complete samples="
        << format_count(static_cast<long long>(sample_count))
        << " elapsed_s=" << std::fixed << std::setprecision(1)
        << elapsed_seconds;
    log_success(log_prefix, out.str());
}

inline void ensure_tree_present(const SampleIO::Sample &sample,
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

struct EventArgs
{
    std::string list_path;
    std::string output_root;
    std::string selection;
    std::string columns_tsv_path;
};

struct EventInput
{
    SampleListEntry entry;
    SampleIO::Sample sample;
};

inline EventArgs parse_event_args(const std::vector<std::string> &args, const std::string &usage)
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

int run(const EventArgs &event_args, const std::string &log_prefix);




#endif // HERON_APPS_EVENTCLI_H
