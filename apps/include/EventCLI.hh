/* -- C++ -- */
#ifndef NUXSEC_APPS_EVENTCLI_H
#define NUXSEC_APPS_EVENTCLI_H

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
#include "EventIO.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"

namespace nuxsec
{

namespace app
{

namespace event
{

inline void log_event_start(const std::string &log_prefix, const size_t sample_count)
{
    nuxsec::app::log::log_info(
        log_prefix,
        "action=event_build status=start samples=" +
            nuxsec::app::log::format_count(static_cast<long long>(sample_count)));
}

inline void log_event_finish(const std::string &log_prefix,
                             const size_t sample_count,
                             const double elapsed_seconds)
{
    std::ostringstream out;
    out << "action=event_build status=complete samples="
        << nuxsec::app::log::format_count(static_cast<long long>(sample_count))
        << " elapsed_s=" << std::fixed << std::setprecision(1)
        << elapsed_seconds;
    nuxsec::app::log::log_success(log_prefix, out.str());
}

inline void ensure_tree_present(const nuxsec::sample::SampleIO::Sample &sample,
                                const std::string &tree_name)
{
    if (sample.inputs.empty())
    {
        throw std::runtime_error("Event inputs missing ROOT files for sample: " + sample.sample_name);
    }

    std::vector<std::string> files = nuxsec::sample::SampleIO::resolve_root_files(sample);
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

struct Args
{
    std::string list_path;
    std::string output_root;
    std::string selection;
    std::string columns_tsv_path;
};

struct Input
{
    nuxsec::app::sample::SampleListEntry entry;
    nuxsec::sample::SampleIO::Sample sample;
};

inline Args parse_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() < 2 || args.size() > 4)
    {
        throw std::runtime_error(usage);
    }

    Args out;
    out.list_path = nuxsec::app::trim(args.at(0));
    out.output_root = nuxsec::app::trim(args.at(1));
    if (args.size() > 2)
    {
        out.selection = nuxsec::app::trim(args.at(2));
    }
    if (args.size() > 3)
    {
        out.columns_tsv_path = nuxsec::app::trim(args.at(3));
    }

    if (out.list_path.empty() || out.output_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path)");
    }

    std::filesystem::path output_root(out.output_root);
    if (output_root.is_relative() && output_root.parent_path().empty())
    {
        const std::filesystem::path event_dir =
            nuxsec::app::stage_output_dir("NUXSEC_EVENT_DIR", "event");
        out.output_root = (event_dir / output_root).string();
    }

    return out;
}

int run(const Args &event_args, const std::string &log_prefix);

} // namespace event

} // namespace app

} // namespace nuxsec

#endif // NUXSEC_APPS_EVENTCLI_H
