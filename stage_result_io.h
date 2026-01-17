/**
 *  @file  stage_result_io.h
 *
 *  @brief Data structures and IO helpers for NuCondenser stage results
 */

#ifndef NUCOND_STAGE_RESULT_IO_H
#define NUCOND_STAGE_RESULT_IO_H

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nucond
{

static inline std::string ToLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    return s;
}

enum class SampleKind
{
    kUnknown = 0,
    kData,
    kEXT,
    kMCOverlay,
    kMCDirt,
    kMCStrangeness
};

static inline const char* SampleKindName(SampleKind k)
{
    switch (k)
    {
    case SampleKind::kData:
        return "data";
    case SampleKind::kEXT:
        return "ext";
    case SampleKind::kMCOverlay:
        return "mc_overlay";
    case SampleKind::kMCDirt:
        return "mc_dirt";
    case SampleKind::kMCStrangeness:
        return "mc_strangeness";
    default:
        return "unknown";
    }
}

enum class BeamMode
{
    kUnknown = 0,
    kNuMI,
    kBNB
};

static inline const char* BeamModeName(BeamMode b)
{
    switch (b)
    {
    case BeamMode::kNuMI:
        return "numi";
    case BeamMode::kBNB:
        return "bnb";
    default:
        return "unknown";
    }
}

struct RunSubrun
{
    int run = 0;
    int subrun = 0;
};

struct SubRunSummary
{
    double pot_sum = 0.0;
    long long n_entries = 0;
    std::vector<RunSubrun> unique_pairs;
};

struct DBSums
{
    double tortgt_sum = 0.0;
    double tor101_sum = 0.0;
    double tor860_sum = 0.0;
    double tor875_sum = 0.0;

    long long EA9CNT_sum = 0;
    long long E1DCNT_sum = 0;
    long long EXTTrig_sum = 0;
    long long Gate1Trig_sum = 0;
    long long Gate2Trig_sum = 0;

    long long n_pairs_loaded = 0;
};

struct StageConfig
{
    std::string stage_name;
    std::string filelist_path;
};

struct StageResult
{
    StageConfig cfg;
    SampleKind kind = SampleKind::kUnknown;
    BeamMode beam = BeamMode::kUnknown;

    std::vector<std::string> input_files;

    SubRunSummary subrun;
    DBSums dbsums;

    double scale = 1.0;

    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;
};

static inline SampleKind SampleKindFromName(const std::string& name)
{
    const std::string v = ToLower(name);
    if (v == "data")
        return SampleKind::kData;
    if (v == "ext")
        return SampleKind::kEXT;
    if (v == "mc_overlay")
        return SampleKind::kMCOverlay;
    if (v == "mc_dirt")
        return SampleKind::kMCDirt;
    if (v == "mc_strangeness")
        return SampleKind::kMCStrangeness;
    return SampleKind::kUnknown;
}

static inline BeamMode BeamModeFromName(const std::string& name)
{
    const std::string v = ToLower(name);
    if (v == "numi")
        return BeamMode::kNuMI;
    if (v == "bnb")
        return BeamMode::kBNB;
    return BeamMode::kUnknown;
}

class StageResultIO
{
  public:
    static void Write(const StageResult& r, const std::string& outFile)
    {
        std::unique_ptr<TFile> f(TFile::Open(outFile.c_str(), "UPDATE"));
        if (!f || f->IsZombie())
        {
            throw std::runtime_error("Failed to open merged output file for UPDATE: " + outFile);
        }

        TDirectory* d = f->GetDirectory("NuCondenser");
        if (!d)
            d = f->mkdir("NuCondenser");
        d->cd();

        TNamed("stage_name", r.cfg.stage_name.c_str()).Write("stage_name", TObject::kOverwrite);
        TNamed("sample_kind", SampleKindName(r.kind)).Write("sample_kind", TObject::kOverwrite);
        TNamed("beam_mode", BeamModeName(r.beam)).Write("beam_mode", TObject::kOverwrite);

        TParameter<double>("subrun_pot_sum", r.subrun.pot_sum).Write("subrun_pot_sum", TObject::kOverwrite);
        TParameter<long long>("subrun_entries", r.subrun.n_entries).Write("subrun_entries", TObject::kOverwrite);
        TParameter<long long>("unique_run_subrun_pairs", static_cast<long long>(r.subrun.unique_pairs.size()))
            .Write("unique_run_subrun_pairs", TObject::kOverwrite);

        TParameter<double>("db_tortgt_sum_raw", r.dbsums.tortgt_sum).Write("db_tortgt_sum_raw", TObject::kOverwrite);
        TParameter<double>("db_tor101_sum_raw", r.dbsums.tor101_sum).Write("db_tor101_sum_raw", TObject::kOverwrite);

        TParameter<double>("db_tortgt_pot", r.db_tortgt_pot).Write("db_tortgt_pot", TObject::kOverwrite);
        TParameter<double>("db_tor101_pot", r.db_tor101_pot).Write("db_tor101_pot", TObject::kOverwrite);

        TParameter<long long>("db_ea9cnt_sum", r.dbsums.EA9CNT_sum).Write("db_ea9cnt_sum", TObject::kOverwrite);
        TParameter<long long>("db_exttrig_sum", r.dbsums.EXTTrig_sum).Write("db_exttrig_sum", TObject::kOverwrite);
        TParameter<long long>("db_gate1trig_sum", r.dbsums.Gate1Trig_sum).Write("db_gate1trig_sum", TObject::kOverwrite);
        TParameter<long long>("db_gate2trig_sum", r.dbsums.Gate2Trig_sum).Write("db_gate2trig_sum", TObject::kOverwrite);

        TParameter<double>("scale_factor", r.scale).Write("scale_factor", TObject::kOverwrite);

        {
            TObjArray arr;
            arr.SetOwner(true);
            for (const auto& in : r.input_files)
            {
                arr.Add(new TObjString(in.c_str()));
            }
            arr.Write("input_files", TObject::kSingleKey | TObject::kOverwrite);
        }

        {
            TTree rs("run_subrun", "Unique (run,subrun) pairs used for DB sums");
            Int_t run = 0;
            Int_t subrun = 0;
            rs.Branch("run", &run, "run/I");
            rs.Branch("subrun", &subrun, "subrun/I");
            for (const auto& p : r.subrun.unique_pairs)
            {
                run = p.run;
                subrun = p.subrun;
                rs.Fill();
            }
            rs.Write("run_subrun", TObject::kOverwrite);
        }

        f->Write();
        f->Close();
    }

    static StageResult Read(const std::string& inFile)
    {
        std::unique_ptr<TFile> f(TFile::Open(inFile.c_str(), "READ"));
        if (!f || f->IsZombie())
        {
            throw std::runtime_error("Failed to open merged output file for READ: " + inFile);
        }

        TDirectory* d = f->GetDirectory("NuCondenser");
        if (!d)
        {
            throw std::runtime_error("Missing NuCondenser directory in file: " + inFile);
        }
        d->cd();

        StageResult r;
        r.cfg.stage_name = ReadNamedString(d, "stage_name");
        r.kind = SampleKindFromName(ReadNamedString(d, "sample_kind"));
        r.beam = BeamModeFromName(ReadNamedString(d, "beam_mode"));

        r.subrun.pot_sum = ReadParam<double>(d, "subrun_pot_sum");
        r.subrun.n_entries = ReadParam<long long>(d, "subrun_entries");

        r.dbsums.tortgt_sum = ReadParam<double>(d, "db_tortgt_sum_raw");
        r.dbsums.tor101_sum = ReadParam<double>(d, "db_tor101_sum_raw");
        r.db_tortgt_pot = ReadParam<double>(d, "db_tortgt_pot");
        r.db_tor101_pot = ReadParam<double>(d, "db_tor101_pot");

        r.dbsums.EA9CNT_sum = ReadParam<long long>(d, "db_ea9cnt_sum");
        r.dbsums.EXTTrig_sum = ReadParam<long long>(d, "db_exttrig_sum");
        r.dbsums.Gate1Trig_sum = ReadParam<long long>(d, "db_gate1trig_sum");
        r.dbsums.Gate2Trig_sum = ReadParam<long long>(d, "db_gate2trig_sum");
        r.scale = ReadParam<double>(d, "scale_factor");

        r.input_files = ReadInputFiles(d);
        r.subrun.unique_pairs = ReadRunSubrunPairs(d);
        return r;
    }

  private:
    static std::string ReadNamedString(TDirectory* d, const char* key)
    {
        TObject* obj = d->Get(key);
        auto* named = dynamic_cast<TNamed*>(obj);
        if (!named)
        {
            throw std::runtime_error("Missing TNamed for key: " + std::string(key));
        }
        return std::string(named->GetTitle());
    }

    template <typename T>
    static T ReadParam(TDirectory* d, const char* key)
    {
        TObject* obj = d->Get(key);
        auto* param = dynamic_cast<TParameter<T>*>(obj);
        if (!param)
        {
            throw std::runtime_error("Missing TParameter for key: " + std::string(key));
        }
        return param->GetVal();
    }

    static std::vector<std::string> ReadInputFiles(TDirectory* d)
    {
        std::vector<std::string> files;
        TObject* obj = d->Get("input_files");
        auto* arr = dynamic_cast<TObjArray*>(obj);
        if (!arr)
        {
            throw std::runtime_error("Missing input_files array");
        }
        const int n = arr->GetEntries();
        files.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i)
        {
            auto* entry = dynamic_cast<TObjString*>(arr->At(i));
            if (!entry)
            {
                throw std::runtime_error("Invalid entry in input_files array");
            }
            files.emplace_back(entry->GetString().Data());
        }
        return files;
    }

    static std::vector<RunSubrun> ReadRunSubrunPairs(TDirectory* d)
    {
        TObject* obj = d->Get("run_subrun");
        auto* tree = dynamic_cast<TTree*>(obj);
        if (!tree)
        {
            throw std::runtime_error("Missing run_subrun tree");
        }

        Int_t run = 0;
        Int_t subrun = 0;
        tree->SetBranchAddress("run", &run);
        tree->SetBranchAddress("subrun", &subrun);

        const Long64_t n = tree->GetEntries();
        std::vector<RunSubrun> pairs;
        pairs.reserve(static_cast<size_t>(n));
        for (Long64_t i = 0; i < n; ++i)
        {
            tree->GetEntry(i);
            pairs.push_back(RunSubrun{static_cast<int>(run), static_cast<int>(subrun)});
        }
        return pairs;
    }
};

} // namespace nucond

#endif // NUCOND_STAGE_RESULT_IO_H
