/* -- C++ -- */
/**
 *  @file  io/include/SampleNormalisationService.hh
 *
 *  @brief Sample normalisation service helpers.
 */

#ifndef NUXSEC_SAMPLE_SAMPLE_NORMALISATION_SERVICE_H
#define NUXSEC_SAMPLE_SAMPLE_NORMALISATION_SERVICE_H

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"
#include "SampleIO.hh"

namespace nuxsec
{

using SampleIO = sample::SampleIO;

class SampleNormalisationService
{
  public:
    static SampleIO::Sample aggregate(const std::string &sample_name,
                                      const std::vector<std::string> &artio_files,
                                      const std::string &db_path);

  private:
    static double compute_normalisation(double subrun_pot_sum, double db_tortgt_pot);
    static SampleIO::ProvenanceInput make_fragment(const artio::Provenance &prov,
                                                  const std::string &artio_path,
                                                  double db_tortgt_pot,
                                                  double db_tor101_pot);
};

} // namespace nuxsec

#endif // NUXSEC_SAMPLE_SAMPLE_NORMALISATION_SERVICE_H
