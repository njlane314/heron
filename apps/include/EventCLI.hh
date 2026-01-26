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

#include "AppUtils.hh"
#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"
#include "EventIO.hh"
#include "RDataFrameFactory.hh"
#include "SampleIO.hh"
#include "SampleUtils.hh"

namespace nuxsec
{

namespace app
{

namespace event
{

inline std::string format_count(const size_t count)
{
    std::ostringstream out;
    if (count >= 1000000)
    {
        out << std::fixed << std::setprecision(1)
            << (static_cast<double>(count) / 1000000.0) << "M";
    }
    else if (count >= 1000)
    {
        out << (count / 1000) << "k";
    }
    else
    {
        out << count;
    }
    return out.str();
}

inline void log_event_start(const std::string &log_prefix, const size_t sample_count)
{
    std::cerr << "[" << log_prefix << "] "
              << "Building events for " << format_count(sample_count) << " samples\n";
}

inline void log_event_finish(const std::string &log_prefix,
                             const size_t sample_count,
                             const double elapsed_seconds)
{
    std::cerr << "[" << log_prefix << "] "
              << "Completed event build for " << format_count(sample_count) << " samples in "
              << std::fixed << std::setprecision(1)
              << elapsed_seconds << "s\n";
}

inline void ensure_tree_present(const nuxsec::sample::SampleIO::Sample &sample,
                                const std::string &tree_name)
{
    if (sample.inputs.empty())
    {
        throw std::runtime_error("Event inputs missing ROOT files for sample: " + sample.sample_name);
    }

    std::vector<std::string> files = nuxsec::RDataFrameFactory::collect_files(sample);
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
};

struct Input
{
    nuxsec::app::SampleListEntry entry;
    nuxsec::sample::SampleIO::Sample sample;
};

inline Args parse_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 2)
    {
        throw std::runtime_error(usage);
    }

    Args out;
    out.list_path = nuxsec::app::trim(args.at(0));
    out.output_root = nuxsec::app::trim(args.at(1));

    if (out.list_path.empty() || out.output_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path)");
    }

    return out;
}

int run(const Args &event_args, const std::string &log_prefix);

} // namespace event

} // namespace app

} // namespace nuxsec

#endif // NUXSEC_APPS_EVENTCLI_H
