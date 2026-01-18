/**
 *  @file  bin/artIOaggregator/artIOaggregator.cxx
 *
 *  @brief Main entrypoint for ArtIO manifest generation
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

#include "NuIO/ArtIOManifestIO.h"
#include "NuIO/RunInfoDB.h"
#include "NuIO/SubRunScanner.h"

namespace
{

std::string Trim(std::string s)
{
    auto notspace = [](unsigned char c)
    {
        return std::isspace(c) == 0;
    };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

std::vector<std::string> ReadFileList(const std::string &filelistPath)
{
    std::ifstream fin(filelistPath);
    if (!fin)
    {
        throw std::runtime_error("Failed to open filelist: " + filelistPath +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    std::vector<std::string> files;
    std::string line;
    while (std::getline(fin, line))
    {
        line = Trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        files.push_back(line);
    }
    if (files.empty())
    {
        throw std::runtime_error("Filelist is empty: " + filelistPath);
    }
    return files;
}

bool IsNuSelectionDataFile(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct CLI
{
    std::string artio_path = "./ArtIO.root";
    nuio::StageCfg stage_cfg;
};

CLI ParseArgs(int argc, char **argv)
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

    CLI cli;
    cli.stage_cfg.stage_name = Trim(spec.substr(0, pos));
    cli.stage_cfg.filelist_path = Trim(spec.substr(pos + 1));

    if (cli.stage_cfg.stage_name.empty() || cli.stage_cfg.filelist_path.empty())
    {
        throw std::runtime_error("Bad stage spec: " + spec);
    }

    return cli;
}

}

int main(int argc, char **argv)
{
    using namespace nuio;

    try
    {
        const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
        const double pot_scale = 1e12;

        const CLI cli = ParseArgs(argc, argv);

        std::vector<std::string> existing = ArtIOManifestIO::ListStages(cli.artio_path);
        std::sort(existing.begin(), existing.end());

        RunInfoDB db(db_path);

        if (std::binary_search(existing.begin(), existing.end(), cli.stage_cfg.stage_name))
        {
            std::cerr << "[artIOaggregator] exists stage=" << cli.stage_cfg.stage_name << "\n";
            return 0;
        }

        const auto files = ReadFileList(cli.stage_cfg.filelist_path);

        ArtProvenance rec;
        rec.cfg = cli.stage_cfg;
        rec.input_files = files;

        if (IsNuSelectionDataFile(files.front()))
        {
            rec.kind = SampleKind::kData;
        }

        rec.subrun = ScanSubRunTree(files);
        rec.runinfo = db.SumRuninfo(rec.subrun.unique_pairs);

        std::cerr << "[artIOaggregator] add stage=" << rec.cfg.stage_name
                  << " files=" << rec.input_files.size()
                  << " pairs=" << rec.subrun.unique_pairs.size()
                  << " pot_sum=" << rec.subrun.pot_sum
                  << " tortgt=" << rec.runinfo.tortgt_sum * pot_scale
                  << "\n";

        ArtIOManifestIO::WriteStage(cli.artio_path, db_path, pot_scale, rec);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
