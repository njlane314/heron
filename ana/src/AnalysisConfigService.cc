/* -- C++ -- */
/**
 *  @file  ana/src/AnalysisConfigService.cc
 *
 *  @brief Compiled analysis configuration service.
 */

#include "AnalysisConfigService.hh"

#include "ColumnDerivationService.hh"
#include "SampleIO.hh"

namespace nuxsec
{

const AnalysisConfigService &AnalysisConfigService::instance()
{
    static const AnalysisConfigService analysis{};
    return analysis;
}

AnalysisConfigService::AnalysisConfigService()
{
    m_name = "nuxsec_default_v1";
    m_tree_name = "MyTree";

}

ProcessorEntry AnalysisConfigService::make_processor_entry(const SampleIO::Sample &sample) const noexcept
{
    ProcessorEntry proc_entry;

    switch (sample.kind)
    {
    case SampleIO::SampleKind::kData:
        proc_entry.source = SourceKind::kData;
        break;
    case SampleIO::SampleKind::kEXT:
        proc_entry.source = SourceKind::kExt;
        proc_entry.trig_nom = sample.db_tor101_pot_sum;
        proc_entry.trig_eqv = sample.subrun_pot_sum;
        break;
    case SampleIO::SampleKind::kOverlay:
    case SampleIO::SampleKind::kDirt:
    case SampleIO::SampleKind::kStrangeness:
        proc_entry.source = SourceKind::kMC;
        proc_entry.pot_nom = sample.db_tortgt_pot_sum;
        proc_entry.pot_eqv = sample.subrun_pot_sum;
        break;
    default:
        proc_entry.source = SourceKind::kUnknown;
        break;
    }

    return proc_entry;
}

} // namespace nuxsec
