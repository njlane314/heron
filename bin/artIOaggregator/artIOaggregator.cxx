/**
 *  @file  bin/artIOaggregator/artIOaggregator.cxx
 *
 *  @brief Main entrypoint for ArtIO provenance generation
 */

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "NuIO/ArtProvenanceIO.h"
#include "NuIO/RunInfoDB.h"
#include "NuIO/SubRunScanner.h"

namespace
{

std::string trim(std::string s)
{
    auto notspace = [](unsigned char c)
    {
        return std::isspace(c) == 0;
    };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

std::vector<std::string> read_file_list(const std::string &filelist_path)
{
    std::ifstream fin(filelist_path);
    if (!fin)
    {
        throw std::runtime_error("Failed to open filelist: " + filelist_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    std::vector<std::string> files;
    std::string line;
    while (std::getline(fin, line))
    {
        line = trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        files.push_back(line);
    }
    if (files.empty())
    {
        throw std::runtime_error("Filelist is empty: " + filelist_path);
    }
    return files;
}

bool is_nu_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct Args
{
    std::string artio_path = "./ArtIO.root";
    nuio::StageCfg stage_cfg;
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
    args.stage_cfg.stage_name = trim(spec.substr(0, pos));
    args.stage_cfg.filelist_path = trim(spec.substr(pos + 1));

    if (args.stage_cfg.stage_name.empty() || args.stage_cfg.filelist_path.empty())
    {
        throw std::runtime_error("Bad stage spec: " + spec);
    }

    return args;
}

}

int main(int argc, char **argv)
{
    using namespace nuio;

    try
    {
        const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
        const double pot_scale = 1e12;

        const Args args = parse_args(argc, argv);

        RunInfoDB db(db_path);

        const auto files = read_file_list(args.stage_cfg.filelist_path);

        ArtProvenance rec;
        rec.cfg = args.stage_cfg;
        rec.input_files = files;

        if (is_nu_selection_data_file(files.front()))
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
