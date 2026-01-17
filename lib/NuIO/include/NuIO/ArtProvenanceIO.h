/**
 *  @file  lib/NuIO/include/NuIO/ArtProvenanceIO.h
 *
 *  @brief Data structures and IO helpers for NuIO aggregator stage provenance
 */

#ifndef NUIO_ART_PROVENANCE_IO_H
#define NUIO_ART_PROVENANCE_IO_H

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

namespace nuio
{

inline std::string ToLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c)
                   { return static_cast<char>(std::tolower(c)); });
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

const char* SampleKindName(SampleKind k);

enum class BeamMode
{
    kUnknown = 0,
    kNuMI,
    kBNB
};

const char* BeamModeName(BeamMode b);

struct RunSubrun
{
    int run = 0;
    int subrun = 0;
};

struct SubRunInfo
{
    double pot_sum = 0.0;
    long long n_entries = 0;
    std::vector<RunSubrun> unique_pairs;
};

struct RunInfoSums
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

struct StageCfg
{
    std::string stage_name;
    std::string filelist_path;
};

struct ArtProvenance
{
    StageCfg cfg;
    SampleKind kind = SampleKind::kUnknown;
    BeamMode beam = BeamMode::kUnknown;

    std::vector<std::string> input_files;

    SubRunInfo subrun;
    RunInfoSums runinfo;

    double scale = 1.0;

    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;
};

SampleKind SampleKindFromName(const std::string& name);
BeamMode BeamModeFromName(const std::string& name);

class ArtProvenanceIO
{
  public:
    static void Write(const ArtProvenance& r, const std::string& outFile);
    static ArtProvenance Read(const std::string& inFile);

  private:
    static std::string ReadNamedString(TDirectory* d, const char* key);

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

    static std::vector<std::string> ReadInputFiles(TDirectory* d);
    static std::vector<RunSubrun> ReadRunSubrunPairs(TDirectory* d);
};

} // namespace nuio

#endif // NUIO_ART_PROVENANCE_IO_H
