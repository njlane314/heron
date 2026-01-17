/**
 *  @file  bin/artIOpartition/artIOpartition.cxx
 *
 *  @brief Main entrypoint for ArtIO partition manifest generation
 */

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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
            continue;
        files.push_back(line);
    }
    if (files.empty())
        throw std::runtime_error("Filelist is empty: " + filelistPath);
    return files;
}

std::vector<std::string> ReadArgsFile(const std::string &argsPath)
{
    std::ifstream fin(argsPath);
    if (!fin)
    {
        throw std::runtime_error("Failed to open args file: " + argsPath +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    std::vector<std::string> args;
    std::string line;
    while (std::getline(fin, line))
    {
        line = Trim(line);
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        std::string tok;
        while (iss >> tok)
            args.push_back(tok);
    }
    return args;
}

bool IsNuSelectionDataFile(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

std::vector<std::string> ExpandArgs(int argc, char **argv)
{
    std::vector<std::string> expanded;
    expanded.reserve(static_cast<size_t>(argc));
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0 && std::string(argv[i]) == "--args")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("Missing value for --args");
            const auto fileArgs = ReadArgsFile(argv[++i]);
            expanded.insert(expanded.end(), fileArgs.begin(), fileArgs.end());
            continue;
        }
        expanded.emplace_back(argv[i]);
    }
    return expanded;
}

struct StageSpec
{
    nuio::StageCfg cfg;
};

struct CLI
{
    std::string artio_path = "./ArtIO.root";
    std::vector<StageSpec> stages;
};

CLI ParseArgs(const std::vector<std::string> &args)
{
    CLI cli;

    for (size_t i = 1; i < args.size(); ++i)
    {
        const std::string &a = args[i];

        auto need = [&](const char *opt) -> std::string
        {
            if (i + 1 >= args.size())
                throw std::runtime_error(std::string("Missing value for ") + opt);
            return args[++i];
        };

        if (a == "--args")
        {
            throw std::runtime_error("Unexpected --args");
        }
        else if (a == "--artio")
        {
            cli.artio_path = need("--artio");
        }
        else if (a == "--stage")
        {
            const std::string spec = need("--stage");
            const auto pos = spec.find(':');
            if (pos == std::string::npos)
                throw std::runtime_error("Bad --stage spec (expected NAME:FILELIST): " + spec);

            StageSpec st;
            st.cfg.stage_name = Trim(spec.substr(0, pos));
            st.cfg.filelist_path = Trim(spec.substr(pos + 1));

            if (st.cfg.stage_name.empty() || st.cfg.filelist_path.empty())
                throw std::runtime_error("Bad --stage spec: " + spec);

            cli.stages.push_back(std::move(st));
        }
        else
        {
            throw std::runtime_error("Unknown argument: " + a);
        }
    }

    if (cli.stages.empty())
        throw std::runtime_error("No --stage specified.");

    return cli;
}

bool HasStageName(const std::vector<std::string> &sorted_names, const std::string &name)
{
    return std::binary_search(sorted_names.begin(), sorted_names.end(), name);
}

}

int main(int argc, char **argv)
{
    using namespace nuio;

    try
    {
        const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
        const double pot_scale = 1e12;

        const auto args = ExpandArgs(argc, argv);
        const CLI cli = ParseArgs(args);

        std::vector<std::string> existing = ArtIOManifestIO::ListStages(cli.artio_path);
        std::sort(existing.begin(), existing.end());

        RunInfoDB db(db_path);

        std::vector<ArtIOStage> out;
        out.reserve(cli.stages.size());

        for (const auto &st : cli.stages)
        {
            if (HasStageName(existing, st.cfg.stage_name))
            {
                std::cerr << "[artIOpartition] exists stage=" << st.cfg.stage_name << "\n";
                continue;
            }

            const auto files = ReadFileList(st.cfg.filelist_path);

            ArtIOStage rec;
            rec.cfg = st.cfg;
            rec.n_input_files = static_cast<long long>(files.size());

            if (IsNuSelectionDataFile(files.front()))
                rec.kind = SampleKind::kData;

            rec.subrun = ScanSubRunTree(files);
            rec.runinfo = db.SumRuninfoForSelection(rec.subrun.unique_pairs);

            std::cerr << "[artIOpartition] add stage=" << rec.cfg.stage_name
                      << " files=" << rec.n_input_files
                      << " pairs=" << rec.subrun.unique_pairs.size()
                      << " pot_sum=" << rec.subrun.pot_sum
                      << " tortgt=" << rec.runinfo.tortgt_sum * pot_scale
                      << "\n";

            out.push_back(std::move(rec));
        }

        ArtIOManifestIO::AppendStages(cli.artio_path, db_path, pot_scale, out);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
