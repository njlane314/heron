/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisConfigService.hh
 *
 *  @brief Compiled analysis configuration service that aggregates selections,
 *         columns, and configuration metadata for processing.
 */

#ifndef HERON_ANA_ANALYSIS_CONFIG_SERVICE_H
#define HERON_ANA_ANALYSIS_CONFIG_SERVICE_H

#include <string>

#include "ColumnDerivationService.hh"
#include "SampleIO.hh"


class AnalysisConfigService
{
  public:
    virtual ~AnalysisConfigService() = default;

    virtual const std::string &name() const noexcept = 0;
    virtual const std::string &tree_name() const noexcept = 0;
    virtual ProcessorEntry make_processor(const SampleIO::Sample &sample) const noexcept
    {
        ProcessorEntry proc_entry;

        switch (sample.origin)
        {
        case SampleIO::SampleOrigin::kData:
            proc_entry.source = Type::kData;
            break;
        case SampleIO::SampleOrigin::kEXT:
            proc_entry.source = Type::kExt;
            proc_entry.trig_nom = sample.db_tor101_pot_sum;
            proc_entry.trig_eqv = sample.subrun_pot_sum;
            break;
        case SampleIO::SampleOrigin::kOverlay:
        case SampleIO::SampleOrigin::kDirt:
        case SampleIO::SampleOrigin::kStrangeness:
            proc_entry.source = Type::kMC;
            proc_entry.pot_nom = sample.db_tortgt_pot_sum;
            proc_entry.pot_eqv = sample.subrun_pot_sum;
            break;
        default:
            proc_entry.source = Type::kUnknown;
            break;
        }

        return proc_entry;
    }
};


#endif // HERON_ANA_ANALYSIS_CONFIG_SERVICE_H
