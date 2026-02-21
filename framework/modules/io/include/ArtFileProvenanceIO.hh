/* -- C++ -- */
/**
 *  @file  framework/io/include/ArtFileProvenanceIO.hh
 *
 *  @brief Data structures and IO helpers for Art file provenance ROOT IO,
 *         including provenance records, mapping, and persistence support.
 */

#ifndef HERON_IO_ART_FILE_PROVENANCE_IO_H
#define HERON_IO_ART_FILE_PROVENANCE_IO_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TTree.h>

#include "SampleIO.hh"



struct Subrun
{
    int run = 0;
    int subrun = 0;
};

struct Summary
{
    double pot_sum = 0.0;
    long long n_entries = 0;
    std::vector<Subrun> unique_pairs;
};

struct Input
{
    std::string input_name;
    std::string filelist_path;
};

struct Provenance
{
    Input input;
    SampleIO::SampleOrigin kind = SampleIO::SampleOrigin::kUnknown;
    SampleIO::BeamMode beam = SampleIO::BeamMode::kUnknown;

    std::vector<std::string> input_files;

    Summary summary;

    double scale = 1.0;
};


class ArtFileProvenanceIO
{
  public:
    static void write(const Provenance &r, const std::string &out_file);
    static Provenance read(const std::string &in_file);
    static Provenance read(const std::string &in_file,
                                SampleIO::SampleOrigin kind,
                                SampleIO::BeamMode beam);

  private:
    static std::string read_named_string(TDirectory *d, const char *key);
    static SampleIO::SampleOrigin read_sample_origin(TDirectory *d);

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
    static std::vector<Subrun> read_run_subrun_pairs(TDirectory *d);
    static Provenance read_directory(TDirectory *d,
                                          SampleIO::SampleOrigin kind,
                                          SampleIO::BeamMode beam);
};


#endif // HERON_IO_ART_FILE_PROVENANCE_IO_H
