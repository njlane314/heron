/* -- C++ -- */
/**
 *  @file  sample/include/SampleAggregator.hh
 *
 *  @brief Aggregation helpers for building samples from art provenance.
 */

#ifndef Nuxsec_SAMPLE_SAMPLE_AGGREGATOR_H_INCLUDED
#define Nuxsec_SAMPLE_SAMPLE_AGGREGATOR_H_INCLUDED

#include <string>
#include <vector>

#include "ArtFileProvenanceRootIO.hh"
#include "Sample.hh"

namespace nuxsec
{

class SampleAggregator
{
  public:
    static Sample Aggregate(const std::string &sample_name, const std::vector<std::string> &artio_files);

  private:
    static double ComputeNormalization(double subrun_pot_sum, double db_tortgt_pot);
    static SampleFragment MakeFragment(const ArtFileProvenance &prov, const std::string &artio_path);
};

} // namespace nuxsec

#endif // Nuxsec_SAMPLE_SAMPLE_AGGREGATOR_H_INCLUDED
