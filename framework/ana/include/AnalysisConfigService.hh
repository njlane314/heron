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
#include <vector>

#include "ColumnDerivationService.hh"
#include "SampleIO.hh"


class AnalysisConfigService final
{
  public:
    static const AnalysisConfigService &instance();

    const std::string &name() const noexcept { return m_name; }
    const std::string &tree_name() const noexcept { return m_tree_name; }
    ProcessorEntry make_processor(const SampleIO::Sample &sample) const noexcept;

  private:
    AnalysisConfigService();

    std::string m_name;
    std::string m_tree_name;
};


#endif // HERON_ANA_ANALYSIS_CONFIG_SERVICE_H
