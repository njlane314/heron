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
    virtual ProcessorEntry make_processor(const SampleIO::Sample &sample) const noexcept = 0;
};


#endif // HERON_ANA_ANALYSIS_CONFIG_SERVICE_H
