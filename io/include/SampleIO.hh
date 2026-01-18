/* -- C++ -- */
/**
 *  @file  io/include/SampleIO.hh
 *
 *  @brief Data structures and IO helpers for sample aggregation.
 */

#ifndef Nuxsec_IO_SAMPLE_IO_H_INCLUDED
#define Nuxsec_IO_SAMPLE_IO_H_INCLUDED

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <string>
#include <vector>

#include "ArtProvenanceIO.hh"

namespace nuxsec
{

struct SampleStage
{
    std::string stage_name;
    std::string artio_path;

    double subrun_pot_sum = 0.0;
    double db_tortgt_pot = 0.0;
    double db_tor101_pot = 0.0;

    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;
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

    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;
};

class SampleIO
{
  public:
    static Sample aggregate(const std::string &sample_name, const std::vector<std::string> &artio_files);
    static void write(const Sample &sample, const std::string &out_file);
    static Sample read(const std::string &in_file);

  private:
    static double compute_normalisation(double subrun_pot_sum, double db_tortgt_pot);
    static SampleStage make_stage(const ArtProvenance &prov, const std::string &artio_path);
};

} // namespace nuxsec

#endif // Nuxsec_IO_SAMPLE_IO_H_INCLUDED
