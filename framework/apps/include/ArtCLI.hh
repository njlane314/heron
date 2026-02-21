/* -- C++ -- */
/**
 *  @file  apps/include/ArtCLI.hh
 *
 *  @brief CLI helpers that orchestrate art file processing and related
 *         provenance-driven workflows, coordinating input discovery and
 *         reporting.
 */
#ifndef HERON_APPS_ARTCLI_H
#define HERON_APPS_ARTCLI_H

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AppLog.hh"
#include "AppUtils.hh"
#include "ArtFileProvenanceIO.hh"
#include "SampleIO.hh"
#include "SubRunInventoryService.hh"

inline void log_scan_start(const std::string &log_prefix)
{
    log_info(log_prefix, "action=subrun_scan status=start");
}

inline void log_scan_finish(const std::string &log_prefix,
                            const long long total,
                            const double elapsed_seconds)
{
    std::ostringstream out;
    out << "action=subrun_scan status=complete entries="
        << format_count(total)
        << " elapsed_s=" << std::fixed << std::setprecision(1)
        << elapsed_seconds;
    log_success(log_prefix, out.str());
}

inline bool is_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct ArtArgs
{
    std::string art_path;
    Input input;
    SampleIO::SampleOrigin sample_origin =
        SampleIO::SampleOrigin::kUnknown;
    SampleIO::BeamMode beam_mode = SampleIO::BeamMode::kUnknown;
};

inline ArtArgs parse_art_input(const std::string &input)
{
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= input.size())
    {
        const size_t pos = input.find(':', start);
        if (pos == std::string::npos)
        {
            fields.push_back(trim(input.substr(start)));
            break;
        }
        fields.push_back(trim(input.substr(start, pos - start)));
        start = pos + 1;
    }

    if (fields.size() < 2)
    {
        throw std::runtime_error("Bad input definition (expected NAME:FILELIST): " + input);
    }

    ArtArgs out;
    out.input.input_name = fields[0];
    out.input.filelist_path = fields[1];

    if (out.input.input_name.empty() || out.input.filelist_path.empty())
    {
        throw std::runtime_error("Bad input definition: " + input);
    }

    if (fields.size() >= 4)
    {
        out.sample_origin = SampleIO::parse_sample_origin(fields[2]);
        out.beam_mode = SampleIO::parse_beam_mode(fields[3]);
        if (out.sample_origin == SampleIO::SampleOrigin::kUnknown)
        {
            throw std::runtime_error("Bad input sample kind: " + fields[2]);
        }
        if (out.beam_mode == SampleIO::BeamMode::kUnknown)
        {
            throw std::runtime_error("Bad input beam mode: " + fields[3]);
        }
    }
    else if (fields.size() != 2)
    {
        throw std::runtime_error("Bad input definition (expected NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]): " + input);
    }

    const char *output_dir = getenv_cstr("HERON_OUTPUT_DIR");
    if (!output_dir)
    {
        throw std::runtime_error("HERON_OUTPUT_DIR must be set for the art command");
    }

    const std::filesystem::path art_dir =
        std::filesystem::path(output_dir) / "art";
    out.art_path = (art_dir / ("art_prov_" + out.input.input_name + ".root")).string();

    return out;
}

inline ArtArgs parse_art_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(usage);
    }

    return parse_art_input(args[0]);
}

int run(const ArtArgs &art_args, const std::string &log_prefix);

#endif
