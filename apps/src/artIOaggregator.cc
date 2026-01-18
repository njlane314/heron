/* -- C++ -- */
/**
 *  @file  apps/src/artIOaggregator.cc
 *
 *  @brief Main entrypoint for ArtIO provenance generation.
 */

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "AppUtils.hh"
#include "ArtProvenanceIO.hh"
#include "RunInfoDB.hh"
#include "SubRunScanner.hh"

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
};

Args parse_args(int argc, char **argv)
{
    if (argc != 2)
    {
        throw std::runtime_error("Usage: artIOaggregator NAME:FILELIST");
    }

    const std::string spec = argv[1];
    const auto pos = spec.find(':');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST): " + spec);
    }

    Args args;
    args.stage_cfg.stage_name = nuxsec::app::trim(spec.substr(0, pos));
    args.stage_cfg.filelist_path = nuxsec::app::trim(spec.substr(pos + 1));

    if (args.stage_cfg.stage_name.empty() || args.stage_cfg.filelist_path.empty())
    {
        throw std::runtime_error("Bad stage spec: " + spec);
    }

    args.artio_path = "./ArtIO_" + args.stage_cfg.stage_name + ".root";

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

        RunInfoDB db(db_path);

        const auto files = nuxsec::app::read_file_list(args.stage_cfg.filelist_path);

        ArtProvenance rec;
        rec.cfg = args.stage_cfg;
        rec.input_files = files;

        if (is_selection_data_file(files.front()))
        {
            rec.kind = SampleKind::kData;
        }

        rec.subrun = scan_subrun_tree(files);
        rec.runinfo = db.sum_runinfo(rec.subrun.unique_pairs);

        rec.subrun.pot_sum *= pot_scale;
        rec.runinfo.tortgt_sum *= pot_scale;
        rec.runinfo.tor101_sum *= pot_scale;
        rec.runinfo.tor860_sum *= pot_scale;
        rec.runinfo.tor875_sum *= pot_scale;
        rec.db_tortgt_pot = rec.runinfo.tortgt_sum;
        rec.db_tor101_pot = rec.runinfo.tor101_sum;
        rec.scale = pot_scale;

        std::cerr << "[artIOaggregator] add stage=" << rec.cfg.stage_name
                  << " files=" << rec.input_files.size()
                  << " pairs=" << rec.subrun.unique_pairs.size()
                  << " pot_sum=" << rec.subrun.pot_sum
                  << " tortgt=" << rec.runinfo.tortgt_sum
                  << "\n";

        ArtProvenanceIO::write(rec, args.artio_path);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
