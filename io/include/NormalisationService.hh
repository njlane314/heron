/* -- C++ -- */
/**
 *  @file  io/include/NormalisationService.hh
 *
 *  @brief Sample normalisation service helpers that define scaling and weight
 *         adjustments for consistent analysis bookkeeping.
 */

#ifndef HERON_IO_NORMALISATION_SERVICE_H
#define HERON_IO_NORMALISATION_SERVICE_H

#include <string>
#include <vector>

#include "ArtFileProvenanceIO.hh"
#include "SampleIO.hh"


class NormalisationService
{
  public:
    static SampleIO::Sample build_sample(const std::string &sample_name,
                                                 const std::vector<std::string> &art_files,
                                                 const std::string &db_path);

  private:
    static double compute_normalisation(double subrun_pot_sum, double db_tortgt_pot);
    static SampleIO::ProvenanceInput make_entry(const Provenance &prov,
                                                        const std::string &art_path,
                                                        double db_tortgt_pot,
                                                        double db_tor101_pot);
};


#endif // HERON_IO_NORMALISATION_SERVICE_H
