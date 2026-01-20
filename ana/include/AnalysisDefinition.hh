/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisDefinition.hh
 *
 *  @brief Compiled analysis definition for template production.
 */

#ifndef NUXSEC_ANA_ANALYSIS_DEFINITION_H
#define NUXSEC_ANA_ANALYSIS_DEFINITION_H

#include <string>
#include <vector>

#include "TemplateSpec.hh"
#include "SampleIO.hh"

namespace nuxsec
{

struct ProcessorEntry;
using SampleIO = sample::SampleIO;

/** \brief Compiled analysis configuration for template production. */
class AnalysisDefinition final
{
  public:
    static const AnalysisDefinition &instance();

    const std::string &name() const noexcept { return m_name; }
    const std::string &tree_name() const noexcept { return m_tree_name; }
    const std::vector<TemplateSpec1D> &templates_1d() const noexcept { return m_templates_1d; }
    std::string templates_1d_to_tsv() const;

    ProcessorEntry make_processor_entry(const SampleIO::Sample &sample) const noexcept;

  private:
    AnalysisDefinition();

    std::string m_name;
    std::string m_tree_name;
    std::vector<TemplateSpec1D> m_templates_1d;
};

} // namespace nuxsec

#endif // NUXSEC_ANA_ANALYSIS_DEFINITION_H
