/* -- C++ -- */
/**
 *  @file  io/include/Sample.hh
 *
 *  @brief Data structures for sample aggregation.
 */

#ifndef NUXSEC_SAMPLE_SAMPLE_H
#define NUXSEC_SAMPLE_SAMPLE_H

#include <string>
#include <vector>

#include "SampleTypes.hh"

namespace nuxsec
{

struct SampleFragment
{
    std::string fragment_name;
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

    std::vector<SampleFragment> fragments;

    double subrun_pot_sum = 0.0;
    double db_tortgt_pot_sum = 0.0;
    double db_tor101_pot_sum = 0.0;

    double normalisation = 1.0;
    double normalised_pot_sum = 0.0;
};

} // namespace nuxsec

#endif // NUXSEC_SAMPLE_SAMPLE_H
