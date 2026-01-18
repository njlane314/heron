/* -- C++ -- */
/**
 *  @file  lib/NuAna/include/NuAna/NuAnalysisProcessor.hh
 *
 *  @brief Variable definitions for analysis RDataFrame processing
 */

#ifndef Nu_ANALYSIS_PROCESSOR_H_INCLUDED
#define Nu_ANALYSIS_PROCESSOR_H_INCLUDED

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

namespace nuana
{

enum class NuSource
{
    kUnknown = 0,
    kData,
    kExt,
    kMC
};

enum class NuChannel
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

struct NuProcessorEntry
{
    NuSource source = NuSource::kUnknown;
    double pot_nom = 0.0;
    double pot_eqv = 0.0;
    double trig_nom = 0.0;
    double trig_eqv = 0.0;
};

/** \brief Apply analysis variable definitions to an RDataFrame. */
class NuAnalysisProcessor
{
  public:
    ROOT::RDF::RNode Run(ROOT::RDF::RNode node, const NuProcessorEntry &rec) const;
    static const NuAnalysisProcessor &Processor();

  private:
    static constexpr double kRecognisedPurityMin;
    static constexpr double kRecognisedCompletenessMin;

    static constexpr float kTrainingFraction;
    static constexpr bool kTrainingIncludeExt;

    static bool IsInTruthVolume(float x, float y, float z) noexcept;
    static bool IsInRecoVolume(float x, float y, float z) noexcept;
};

} // namespace nuana

#endif // Nu_ANALYSIS_PROCESSOR_H_INCLUDED
