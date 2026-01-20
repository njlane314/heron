/* -- C++ -- */
/**
 *  @file  io/include/ArtFileProvenanceIO.hh
 *
 *  @brief Data structures and IO helpers for Art file provenance ROOT IO.
 */

#ifndef Nuxsec_IO_ART_FILE_PROVENANCE_IO_H_INCLUDED
#define Nuxsec_IO_ART_FILE_PROVENANCE_IO_H_INCLUDED

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TTree.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SampleTypes.hh"

namespace nuxsec
{

struct RunSubrun
{
    int run = 0;
    int subrun = 0;
};

struct SubrunTreeSummary
{
    double pot_sum = 0.0;
    long long n_entries = 0;
    std::vector<RunSubrun> unique_pairs;
};

struct StageCfg
{
    std::string stage_name;
    std::string filelist_path;
};

struct ArtFileProvenance
{
    StageCfg cfg;
    SampleKind kind = SampleKind::kUnknown;
    BeamMode beam = BeamMode::kUnknown;

    std::vector<std::string> input_files;

    SubrunTreeSummary subrun;

    double scale = 1.0;
};

class ArtFileProvenanceIO
{
  public:
    static void write(const ArtFileProvenance &r, const std::string &out_file);
    static ArtFileProvenance read(const std::string &in_file);
    static ArtFileProvenance read(const std::string &in_file, SampleKind kind, BeamMode beam);

  private:
    static std::string read_named_string(TDirectory *d, const char *key);

    template <typename T>
    static T read_param(TDirectory *d, const char *key)
    {
        TObject *obj = d->Get(key);
        auto *param = dynamic_cast<TParameter<T> *>(obj);
        if (!param)
        {
            throw std::runtime_error("Missing TParameter for key: " + std::string(key));
        }
        return param->GetVal();
    }

    static std::vector<std::string> read_input_files(TDirectory *d);
    static std::vector<RunSubrun> read_run_subrun_pairs(TDirectory *d);
    static ArtFileProvenance read_directory(TDirectory *d, SampleKind kind, BeamMode beam);
};

} // namespace nuxsec

#endif // Nuxsec_IO_ART_FILE_PROVENANCE_IO_H_INCLUDED
