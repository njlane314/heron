/* -- C++ -- */
/**
 *  @file  sample/include/Sample.hh
 *
 *  @brief Sample structures for aggregated art outputs.
 */

#ifndef Nuxsec_SAMPLE_SAMPLE_H_INCLUDED
#define Nuxsec_SAMPLE_SAMPLE_H_INCLUDED

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

    double normalization = 1.0;
    double normalized_pot_sum = 0.0;
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

    double normalization = 1.0;
    double normalized_pot_sum = 0.0;
};

} // namespace nuxsec

#endif // Nuxsec_SAMPLE_SAMPLE_H_INCLUDED
