/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecArtIOaggregator.cc
 *
 *  @brief Main entrypoint for Art file provenance generation.
 */

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AppUtils.hh"
#include "ArtFileProvenanceRootIO.hh"
#include "RunInfoSqliteReader.hh"
#include "SubrunTreeScanner.hh"

namespace
{

bool is_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct Args
{
    std::string artio_path;
    nuxsec::StageCfg stage_cfg;
    nuxsec::SampleKind sample_kind = nuxsec::SampleKind::kUnknown;
    nuxsec::BeamMode beam_mode = nuxsec::BeamMode::kUnknown;
};

Args parse_args(int argc, char **argv)
{
    if (argc != 2)
    {
        throw std::runtime_error("Usage: nuxsecArtIOaggregator NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]");
    }

    const std::string spec = argv[1];

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

    Args args;
    args.stage_cfg.stage_name = fields[0];
    args.stage_cfg.filelist_path = fields[1];

    if (args.stage_cfg.stage_name.empty() || args.stage_cfg.filelist_path.empty())
    {
        throw std::runtime_error("Bad stage spec: " + spec);
    }

    if (fields.size() >= 4)
    {
        args.sample_kind = nuxsec::parse_sample_kind(fields[2]);
        args.beam_mode = nuxsec::parse_beam_mode(fields[3]);
        if (args.sample_kind == nuxsec::SampleKind::kUnknown)
        {
            throw std::runtime_error("Bad stage sample kind: " + fields[2]);
        }
        if (args.beam_mode == nuxsec::BeamMode::kUnknown)
        {
            throw std::runtime_error("Bad stage beam mode: " + fields[3]);
        }
    }
    else if (fields.size() != 2)
    {
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]): " + spec);
    }

    args.artio_path = "build/out/art/art_prov_" + args.stage_cfg.stage_name + ".root";

    return args;
}

}

int main(int argc, char **argv)
{
    using namespace nuxsec;

    try
    {
        const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
        const double pot_scale = 1e12;

        const Args args = parse_args(argc, argv);
        std::filesystem::path out_path(args.artio_path);
        if (!out_path.parent_path().empty())
        {
            std::filesystem::create_directories(out_path.parent_path());
        }

        RunInfoSqliteReader db(db_path);

        const auto files = nuxsec::app::read_file_list(args.stage_cfg.filelist_path);

        ArtFileProvenance rec;
        rec.cfg = args.stage_cfg;
        rec.input_files = files;
        rec.kind = args.sample_kind;
        rec.beam = args.beam_mode;

        if (rec.kind == SampleKind::kUnknown && is_selection_data_file(files.front()))
        {
            rec.kind = SampleKind::kData;
        }

        rec.subrun = scanSubrunTree(files);
        rec.runinfo = db.sumRunInfo(rec.subrun.unique_pairs);

        rec.subrun.pot_sum *= pot_scale;
        rec.runinfo.tortgt_sum *= pot_scale;
        rec.runinfo.tor101_sum *= pot_scale;
        rec.runinfo.tor860_sum *= pot_scale;
        rec.runinfo.tor875_sum *= pot_scale;
        rec.db_tortgt_pot = rec.runinfo.tortgt_sum;
        rec.db_tor101_pot = rec.runinfo.tor101_sum;
        rec.scale = pot_scale;

        std::cerr << "[nuxsecArtIOaggregator] add stage=" << rec.cfg.stage_name
                  << " files=" << rec.input_files.size()
                  << " pairs=" << rec.subrun.unique_pairs.size()
                  << " pot_sum=" << rec.subrun.pot_sum
                  << " tortgt=" << rec.runinfo.tortgt_sum
                  << "\n";

        ArtFileProvenanceRootIO::write(rec, args.artio_path);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
