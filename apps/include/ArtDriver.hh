/* -- C++ -- */
#ifndef NUXSEC_APPS_ARTDRIVER_H
#define NUXSEC_APPS_ARTDRIVER_H

#include <filesystem>
#include <iostream>
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

inline int run(const Args &art_args, const std::string &log_prefix)
{
    const double pot_scale = 1e12;

    std::filesystem::path out_path(art_args.art_path);
    if (!out_path.parent_path().empty())
    {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const auto files = nuxsec::app::read_paths(art_args.input.filelist_path);

    nuxsec::art::Provenance rec;
    rec.input = art_args.input;
    rec.input_files = files;
    rec.kind = art_args.sample_origin;
    rec.beam = art_args.beam_mode;

    if (rec.kind == nuxsec::sample::SampleIO::SampleOrigin::kUnknown &&
        is_selection_data_file(files.front()))
    {
        rec.kind = nuxsec::sample::SampleIO::SampleOrigin::kData;
    }

    rec.summary = nuxsec::SubRunInventoryService::scan_subruns(files);

    rec.summary.pot_sum *= pot_scale;
    rec.scale = pot_scale;

    std::cerr << "[" << log_prefix << "] add input=" << rec.input.input_name
              << " files=" << rec.input_files.size()
              << " pairs=" << rec.summary.unique_pairs.size()
              << " pot_sum=" << rec.summary.pot_sum
              << "\n";

    nuxsec::ArtFileProvenanceIO::write(rec, art_args.art_path);

    return 0;
}

}

}

}

#endif
