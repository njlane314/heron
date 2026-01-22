/* -- C++ -- */
#ifndef NUXSEC_APPS_ARTCOMMAND_H
#define NUXSEC_APPS_ARTCOMMAND_H

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Utils.hh"
#include "ArtFileProvenanceIO.hh"
#include "SampleIO.hh"
#include "SubRunInventoryService.hh"

namespace nuxsec
{

namespace app
{

inline bool is_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct ArtArgs
{
    std::string artio_path;
    nuxsec::artio::Stage stage_cfg;
    sample::SampleIO::SampleOrigin sample_kind = sample::SampleIO::SampleOrigin::kUnknown;
    sample::SampleIO::BeamMode beam_mode = sample::SampleIO::BeamMode::kUnknown;
};

inline ArtArgs parse_art_spec(const std::string &spec)
{
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= spec.size())
    {
        const size_t pos = spec.find(':', start);
        if (pos == std::string::npos)
        {
            fields.push_back(nuxsec::app::trim(spec.substr(start)));
            break;
        }
        fields.push_back(nuxsec::app::trim(spec.substr(start, pos - start)));
        start = pos + 1;
    }

    if (fields.size() < 2)
    {
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST): " + spec);
    }

    ArtArgs out;
    out.stage_cfg.stage_name = fields[0];
    out.stage_cfg.filelist_path = fields[1];

    if (out.stage_cfg.stage_name.empty() || out.stage_cfg.filelist_path.empty())
    {
        throw std::runtime_error("Bad stage spec: " + spec);
    }

    if (fields.size() >= 4)
    {
        out.sample_kind = sample::SampleIO::parse_sample_kind(fields[2]);
        out.beam_mode = sample::SampleIO::parse_beam_mode(fields[3]);
        if (out.sample_kind == sample::SampleIO::SampleOrigin::kUnknown)
        {
            throw std::runtime_error("Bad stage sample kind: " + fields[2]);
        }
        if (out.beam_mode == sample::SampleIO::BeamMode::kUnknown)
        {
            throw std::runtime_error("Bad stage beam mode: " + fields[3]);
        }
    }
    else if (fields.size() != 2)
    {
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]): " + spec);
    }

    out.artio_path = "build/out/art/art_prov_" + out.stage_cfg.stage_name + ".root";

    return out;
}

inline ArtArgs parse_art_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(usage);
    }

    return parse_art_spec(args[0]);
}

inline int run_artio(const ArtArgs &art_args, const std::string &log_prefix)
{
    const double pot_scale = 1e12;

    std::filesystem::path out_path(art_args.artio_path);
    if (!out_path.parent_path().empty())
    {
        std::filesystem::create_directories(out_path.parent_path());
    }

    const auto files = nuxsec::app::read_file_list(art_args.stage_cfg.filelist_path);

    nuxsec::artio::Provenance rec;
    rec.cfg = art_args.stage_cfg;
    rec.input_files = files;
    rec.kind = art_args.sample_kind;
    rec.beam = art_args.beam_mode;

    if (rec.kind == sample::SampleIO::SampleOrigin::kUnknown && is_selection_data_file(files.front()))
    {
        rec.kind = sample::SampleIO::SampleOrigin::kData;
    }

    rec.subrun = nuxsec::SubRunInventoryService::scan_subrun_tree(files);

    rec.subrun.pot_sum *= pot_scale;
    rec.scale = pot_scale;

    std::cerr << "[" << log_prefix << "] add stage=" << rec.cfg.stage_name
              << " files=" << rec.input_files.size()
              << " pairs=" << rec.subrun.unique_pairs.size()
              << " pot_sum=" << rec.subrun.pot_sum
              << "\n";

    nuxsec::ArtFileProvenanceIO::write(rec, art_args.artio_path);

    return 0;
}

}

}

#endif
