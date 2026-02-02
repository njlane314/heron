/* -- C++ -- */
/**
 *  @file  ana/include/ColumnDerivationService.hh
 *
 *  @brief Variable definitions for analysis RDataFrame processing
 */

#ifndef NUXSEC_ANA_COLUMN_DERIVATION_SERVICE_H
#define NUXSEC_ANA_COLUMN_DERIVATION_SERVICE_H

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

namespace nuxsec
{

enum class Type
{
    kUnknown = 0, ///< Unknown or unset source type.
    kData,        ///< On-beam data sample.
    kExt,         ///< Off-beam external/background data sample.
    kMC           ///< Simulated Monte Carlo sample.
};

enum class Channel
{
    Unknown = 0,      ///< Unclassified channel (fallback for unknown categorisation).
    OutFV,            ///< Interaction outside the truth fiducial volume.
    External,         ///< External (non-neutrino) or out-of-volume background.
    NC,               ///< Neutral-current interaction in fiducial volume.
    CCS1,             ///< Charged-current interaction with exactly one strange hadron.
    CCSgt1,           ///< Charged-current interaction with multiple strange hadrons.
    ECCC,             ///< Electron-neutrino charged-current interaction.
    MuCC0pi_ge1p,     ///< Muon-neutrino charged-current with 0 pions and ≥1 proton.
    MuCC1pi,          ///< Muon-neutrino charged-current with exactly one charged pion.
    MuCCPi0OrGamma,   ///< Muon-neutrino charged-current with π0 or photon activity.
    MuCCNpi,          ///< Muon-neutrino charged-current with more than one pion.
    MuCCOther,        ///< Other muon-neutrino charged-current topologies.
    DataInclusive     ///< Inclusive data channel (non-MC).
};

struct ProcessorEntry
{
    Type source = Type::kUnknown;
    double pot_nom = 0.0;
    double pot_eqv = 0.0;
    double trig_nom = 0.0;
    double trig_eqv = 0.0;
};

class ColumnDerivationService
{
  public:
    ROOT::RDF::RNode define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const;
    static const ColumnDerivationService &instance();

  private:
    static const double kRecognisedPurityMin;
    static const double kRecognisedCompletenessMin;

};

} // namespace nuxsec

#endif // NUXSEC_ANA_COLUMN_DERIVATION_SERVICE_H
