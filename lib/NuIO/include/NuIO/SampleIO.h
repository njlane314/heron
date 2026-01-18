/**
 *  @file  lib/NuIO/include/NuIO/SampleIO.h
 *
 *  @brief Data structures and IO helpers for sample aggregation
 */

#ifndef NUIO_SAMPLE_IO_H
#define NUIO_SAMPLE_IO_H

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <string>
#include <vector>

#include "NuIO/ArtProvenanceIO.h"

namespace nuio
{

struct SampleStage
{
    std::string stage_name;
    std::string artio_path;

    double subrun_pot_sum = 0.0;
    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;

    double normalization = 1.0;
    double normalized_pot_sum = 0.0;
};

struct Sample
{
    std::string sample_name;
    SampleKind kind = SampleKind::kUnknown;
    BeamMode beam = BeamMode::kUnknown;

    std::vector<SampleStage> stages;

    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    double normalization = 1.0;
    double normalized_pot_sum = 0.0;
};

class SampleIO
{
  public:
    static Sample aggregate(const std::string &sample_name, const std::vector<std::string> &artio_files);
    static void write(const Sample &sample, const std::string &out_file);
    static Sample read(const std::string &in_file);

  private:
    static double compute_normalization(double subrun_pot_sum, double db_tortgt_pot);
    static SampleStage make_stage(const ArtProvenance &prov, const std::string &artio_path);
};

} // namespace nuio

#endif // NUIO_SAMPLE_IO_H
