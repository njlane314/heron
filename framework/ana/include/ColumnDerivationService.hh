/* -- C++ -- */
/**
 *  @file  ana/include/ColumnDerivationService.hh
 *
 *  @brief Variable definitions for analysis RDataFrame processing, capturing
 *         derived columns and related transformation helpers.
 */

#ifndef HERON_ANA_COLUMN_DERIVATION_SERVICE_H
#define HERON_ANA_COLUMN_DERIVATION_SERVICE_H

#include <ROOT/RDataFrame.hxx>

#include "AnalysisChannels.hh"

enum class Type
{
    kUnknown = 0, ///< Unknown or unset source type.
    kData,        ///< On-beam data sample.
    kExt,         ///< Off-beam external/background data sample.
    kMC           ///< Simulated Monte Carlo sample.
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
    virtual ~ColumnDerivationService() = default;

    static const ColumnDerivationService &instance();
    virtual ROOT::RDF::RNode define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const = 0;
};


#endif // HERON_ANA_COLUMN_DERIVATION_SERVICE_H
