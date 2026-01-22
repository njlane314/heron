/* -- C++ -- */
/**
 *  @file  io/include/ArtFileProvenanceIO.hh
 *
 *  @brief Data structures and IO helpers for Art file provenance ROOT IO.
 */

#ifndef NUXSEC_IO_ART_FILE_PROVENANCE_IO_H
#define NUXSEC_IO_ART_FILE_PROVENANCE_IO_H

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

#include "SampleIO.hh"

namespace nuxsec
{

using SampleIO = sample::SampleIO;

namespace artio
{

struct RunSubrunPair
{
    int run = 0;
    int subrun = 0;
};

struct SubrunSummary
{
    double pot_sum = 0.0;
    long long n_entries = 0;
    std::vector<RunSubrunPair> unique_pairs;
};

struct Stage
{
    std::string stage_name;
    std::string filelist_path;
};

struct Provenance
{
    Stage cfg;
    SampleIO::SampleOrigin kind = SampleIO::SampleOrigin::kUnknown;
    SampleIO::BeamMode beam = SampleIO::BeamMode::kUnknown;

    std::vector<std::string> input_files;

    SubrunSummary subrun;

    double scale = 1.0;
};

} // namespace artio

class ArtFileProvenanceIO
{
  public:
    static void write(const artio::Provenance &r, const std::string &out_file);
    static artio::Provenance read(const std::string &in_file);
    static artio::Provenance read(const std::string &in_file,
                                  SampleIO::SampleOrigin kind,
                                  SampleIO::BeamMode beam);

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
    static std::vector<artio::RunSubrunPair> read_run_subrun_pairs(TDirectory *d);
    static artio::Provenance read_directory(TDirectory *d,
                                            SampleIO::SampleOrigin kind,
                                            SampleIO::BeamMode beam);
};

} // namespace nuxsec

#endif // NUXSEC_IO_ART_FILE_PROVENANCE_IO_H
