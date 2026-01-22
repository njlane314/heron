/* -- C++ -- */
/**
 *  @file  io/include/NormalisationService.hh
 *
 *  @brief Sample normalisation service helpers.
 */

#ifndef NUXSEC_IO_NORMALISATION_SERVICE_H
#define NUXSEC_IO_NORMALISATION_SERVICE_H

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"
#include "SampleIO.hh"

namespace nuxsec
{

using SampleIO = sample::SampleIO;

class NormalisationService
{
  public:
    static SampleIO::Sample aggregate(const std::string &sample_name,
                                      const std::vector<std::string> &art_files,
                                      const std::string &db_path);

  private:
    static double compute_normalisation(double subrun_pot_sum, double db_tortgt_pot);
    static SampleIO::ProvenanceInput make_entry(const art::Provenance &prov,
                                                const std::string &art_path,
                                                double db_tortgt_pot,
                                                double db_tor101_pot);
};

} // namespace nuxsec

#endif // NUXSEC_IO_NORMALISATION_SERVICE_H
