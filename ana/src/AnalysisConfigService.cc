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

ProcessorEntry AnalysisConfigService::make_processor(const SampleIO::Sample &sample) const noexcept
{
    ProcessorEntry proc_entry;

    switch (sample.kind)
    {
    case SampleIO::SampleOrigin::kData:
        proc_entry.source = SourceKind::kData;
        break;
    case SampleIO::SampleOrigin::kEXT:
        proc_entry.source = SourceKind::kExt;
        proc_entry.trig_nom = sample.db_tor101_pot_sum;
        proc_entry.trig_eqv = sample.subrun_pot_sum;
        break;
    case SampleIO::SampleOrigin::kOverlay:
    case SampleIO::SampleOrigin::kDirt:
    case SampleIO::SampleOrigin::kStrangeness:
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
