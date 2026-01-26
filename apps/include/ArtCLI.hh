/* -- C++ -- */
#ifndef NUXSEC_APPS_ARTCLI_H
#define NUXSEC_APPS_ARTCLI_H

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AppUtils.hh"
#include "ArtFileProvenanceIO.hh"
#include "SampleIO.hh"
#include "SubRunInventoryService.hh"

namespace nuxsec
{

namespace app
{

namespace art
{

inline std::string format_count(const long long count)
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

inline void log_scan_start(const std::string &log_prefix)
{
    std::cerr << "[" << log_prefix << "] "
              << "Scanning SubRun entries\n";
}

inline void log_scan_finish(const std::string &log_prefix,
                            const long long total,
                            const double elapsed_seconds)
{
    std::cerr << "[" << log_prefix << "] "
              << "Completed scan of " << format_count(total) << " entries in "
              << std::fixed << std::setprecision(1)
              << elapsed_seconds << "s\n";
}

inline bool is_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct Args
{
    std::string art_path;
    nuxsec::art::Input input;
    nuxsec::sample::SampleIO::SampleOrigin sample_origin =
        nuxsec::sample::SampleIO::SampleOrigin::kUnknown;
    nuxsec::sample::SampleIO::BeamMode beam_mode = nuxsec::sample::SampleIO::BeamMode::kUnknown;
};

inline Args parse_input(const std::string &input)
{
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= input.size())
    {
        const size_t pos = input.find(':', start);
        if (pos == std::string::npos)
        {
            fields.push_back(nuxsec::app::trim(input.substr(start)));
            break;
        }
        fields.push_back(nuxsec::app::trim(input.substr(start, pos - start)));
        start = pos + 1;
    }

    if (fields.size() < 2)
    {
        throw std::runtime_error("Bad input definition (expected NAME:FILELIST): " + input);
    }

    Args out;
    out.input.input_name = fields[0];
    out.input.filelist_path = fields[1];

    if (out.input.input_name.empty() || out.input.filelist_path.empty())
    {
        throw std::runtime_error("Bad input definition: " + input);
    }

    if (fields.size() >= 4)
    {
        out.sample_origin = nuxsec::sample::SampleIO::parse_sample_origin(fields[2]);
        out.beam_mode = nuxsec::sample::SampleIO::parse_beam_mode(fields[3]);
        if (out.sample_origin == nuxsec::sample::SampleIO::SampleOrigin::kUnknown)
        {
            throw std::runtime_error("Bad input sample kind: " + fields[2]);
        }
        if (out.beam_mode == nuxsec::sample::SampleIO::BeamMode::kUnknown)
        {
            throw std::runtime_error("Bad input beam mode: " + fields[3]);
        }
    }
    else if (fields.size() != 2)
    {
        throw std::runtime_error("Bad input definition (expected NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]): " + input);
    }

    out.art_path = "build/out/art/art_prov_" + out.input.input_name + ".root";

    return out;
}

inline Args parse_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(usage);
    }

    return parse_input(args[0]);
}

int run(const Args &art_args, const std::string &log_prefix);

}

}

}

#endif
