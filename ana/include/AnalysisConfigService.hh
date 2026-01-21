/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisConfigService.hh
 *
 *  @brief Compiled analysis configuration service.
 */

#ifndef NUXSEC_ANA_ANALYSIS_CONFIG_SERVICE_H
#define NUXSEC_ANA_ANALYSIS_CONFIG_SERVICE_H

#include <string>
#include <vector>

#include "SampleIO.hh"

namespace nuxsec
{

struct ProcessorEntry;
using SampleIO = sample::SampleIO;

/** \brief Compiled analysis configuration. */
class AnalysisConfigService final
{
  public:
    static const AnalysisConfigService &instance();

    const std::string &name() const noexcept { return m_name; }
    const std::string &tree_name() const noexcept { return m_tree_name; }
    ProcessorEntry make_processor_entry(const SampleIO::Sample &sample) const noexcept;

  private:
    AnalysisConfigService();

    std::string m_name;
    std::string m_tree_name;
};

} // namespace nuxsec

#endif // NUXSEC_ANA_ANALYSIS_CONFIG_SERVICE_H
