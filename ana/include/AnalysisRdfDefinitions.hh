/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisRdfDefinitions.hh
 *
 *  @brief Variable definitions for analysis RDataFrame processing
 */

#ifndef Nuxsec_ANA_ANALYSIS_RDF_DEFINITIONS_H_INCLUDED
#define Nuxsec_ANA_ANALYSIS_RDF_DEFINITIONS_H_INCLUDED

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

namespace nuxsec
{

enum class SourceKind
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
    SourceKind source = SourceKind::kUnknown;
    double pot_nom = 0.0;
    double pot_eqv = 0.0;
    double trig_nom = 0.0;
    double trig_eqv = 0.0;
};

/** \brief Apply analysis variable definitions to an RDataFrame. */
class AnalysisRdfDefinitions
{
  public:
    ROOT::RDF::RNode Define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const;
    static const AnalysisRdfDefinitions &Instance();

  private:
    static constexpr double kRecognisedPurityMin;
    static constexpr double kRecognisedCompletenessMin;

    static constexpr float kTrainingFraction;
    static constexpr bool kTrainingIncludeExt;

    static bool IsInTruthVolume(float x, float y, float z) noexcept;
    static bool IsInRecoVolume(float x, float y, float z) noexcept;
};

} // namespace nuxsec

#endif // Nuxsec_ANA_ANALYSIS_RDF_DEFINITIONS_H_INCLUDED
