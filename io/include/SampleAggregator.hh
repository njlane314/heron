/* -- C++ -- */
/**
 *  @file  io/include/SampleAggregator.hh
 *
 *  @brief Sample aggregation helpers.
 */

#ifndef NUXSEC_SAMPLE_SAMPLE_AGGREGATOR_H
#define NUXSEC_SAMPLE_SAMPLE_AGGREGATOR_H

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"
#include "Sample.hh"

namespace nuxsec
{

class SampleAggregator
{
  public:
    static Sample aggregate(const std::string &sample_name,
                            const std::vector<std::string> &artio_files,
                            const std::string &db_path);

  private:
    static double compute_normalisation(double subrun_pot_sum, double db_tortgt_pot);
    static SampleFragment make_fragment(const ArtFileProvenance &prov,
                                        const std::string &artio_path,
                                        double db_tortgt_pot,
                                        double db_tor101_pot);
};

} // namespace nuxsec

#endif // NUXSEC_SAMPLE_SAMPLE_AGGREGATOR_H
