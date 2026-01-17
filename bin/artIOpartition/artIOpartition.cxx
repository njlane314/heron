/**
 *  @file  bin/artIOpartition/artIOpartition.cxx
 *
 *  @brief Main entrypoint for ArtIO partition manifest generation
 */

#include <TFile.h>
#include <TTree.h>

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

struct FilePeek
{
    std::optional<bool> isData;
    std::optional<bool> isNuMI;
};

FilePeek PeekEventFlags(const std::string &rootFile)
{
    FilePeek out;

    std::unique_ptr<TFile> f(TFile::Open(rootFile.c_str(), "READ"));
    if (!f || f->IsZombie())
        return out;

    TTree *t = dynamic_cast<TTree *>(f->Get("EventSelectionFilter"));
    if (!t)
        return out;

    const bool has_is_data = (t->GetBranch("is_data") != nullptr);
    const bool has_is_numi = (t->GetBranch("is_numi") != nullptr);

    if (!has_is_data && !has_is_numi)
        return out;
    if (t->GetEntries() <= 0)
        return out;

    Bool_t is_data = false;
    Bool_t is_numi = false;

    if (has_is_data)
        t->SetBranchAddress("is_data", &is_data);
    if (has_is_numi)
        t->SetBranchAddress("is_numi", &is_numi);

    t->GetEntry(0);

    if (has_is_data)
        out.isData = static_cast<bool>(is_data);
    if (has_is_numi)
        out.isNuMI = static_cast<bool>(is_numi);

    return out;
}

struct StageSpec
{
    nuio::StageCfg cfg;
};

StageSpec ParseStageSpec(const std::string &spec)
{
    const auto pos = spec.find(':');
    if (pos == std::string::npos)
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST): " + spec);

    StageSpec st;
    st.cfg.stage_name = Trim(spec.substr(0, pos));
    st.cfg.filelist_path = Trim(spec.substr(pos + 1));

    if (st.cfg.stage_name.empty() || st.cfg.filelist_path.empty())
        throw std::runtime_error("Bad stage spec: " + spec);

    return st;
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

        if (argc != 2)
            throw std::runtime_error("Usage: artIOpartition STAGE:FILELIST");

        const StageSpec stage = ParseStageSpec(argv[1]);
        const std::string artio_path = stage.cfg.stage_name + ".root";

        std::vector<std::string> existing = ArtIOManifestIO::ListStages(artio_path);
        std::sort(existing.begin(), existing.end());

        RunInfoDB db(db_path);

        std::vector<ArtIOStage> out;
        out.reserve(1);

        if (HasStageName(existing, stage.cfg.stage_name))
        {
            std::cerr << "[artIOpartition] exists stage=" << stage.cfg.stage_name << "\n";
            return 0;
        }

        const auto files = ReadFileList(stage.cfg.filelist_path);

        ArtIOStage rec;
        rec.cfg = stage.cfg;
        rec.n_input_files = static_cast<long long>(files.size());

        const auto peek = PeekEventFlags(files.front());
        if (peek.isNuMI.has_value())
            rec.beam = (*peek.isNuMI ? BeamMode::kNuMI : BeamMode::kBNB);
        if (peek.isData.has_value())
            rec.kind = (*peek.isData ? SampleKind::kData : SampleKind::kUnknown);

        rec.subrun = ScanSubRunTree(files);
        rec.runinfo = db.SumRuninfoForSelection(rec.subrun.unique_pairs);

        std::cerr << "[artIOpartition] add stage=" << rec.cfg.stage_name
                  << " files=" << rec.n_input_files
                  << " pairs=" << rec.subrun.unique_pairs.size()
                  << " pot_sum=" << rec.subrun.pot_sum
                  << " tortgt=" << rec.runinfo.tortgt_sum * pot_scale
                  << "\n";

        out.push_back(std::move(rec));

        ArtIOManifestIO::AppendStages(artio_path, db_path, pot_scale, out);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
