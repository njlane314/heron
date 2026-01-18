/* -- C++ -- */
/**
 *  @file  io/include/SampleAggregator.hh
 *
 *  @brief Sample aggregation helpers.
 */

#ifndef Nuxsec_IO_SAMPLE_AGGREGATOR_H_INCLUDED
#define Nuxsec_IO_SAMPLE_AGGREGATOR_H_INCLUDED

#include <string>
#include <vector>

#include "ArtFileProvenanceRootIO.hh"
#include "Sample.hh"

namespace nuxsec
{

class SampleAggregator
{
  public:
    static Sample aggregate(const std::string &sample_name, const std::vector<std::string> &artio_files);

  private:
    static double compute_normalisation(double subrun_pot_sum, double db_tortgt_pot);
    static SampleFragment make_fragment(const ArtFileProvenance &prov, const std::string &artio_path);
};

} // namespace nuxsec

#endif // Nuxsec_IO_SAMPLE_AGGREGATOR_H_INCLUDED
