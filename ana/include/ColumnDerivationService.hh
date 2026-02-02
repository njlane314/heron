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
    kUnknown = 0,
    kData,
    kExt,
    kMC
};

enum class Channel
{
    Unknown = 0,
    OutFV,
    External,
    NC,
    CCS1,
    CCSgt1,
    ECCC,
    MuCC0pi_ge1p,
    MuCC1pi,
    MuCCPi0OrGamma,
    MuCCNpi,
    MuCCOther,
    DataInclusive
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
