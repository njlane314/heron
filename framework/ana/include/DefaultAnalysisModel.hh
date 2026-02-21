/* -- C++ -- */
/**
 *  @file  ana/include/DefaultAnalysisModel.hh
 *
 *  @brief Default analysis model that provides concrete service logic.
 */

#ifndef HERON_ANA_DEFAULT_ANALYSIS_MODEL_H
#define HERON_ANA_DEFAULT_ANALYSIS_MODEL_H

#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "AnalysisConfigService.hh"
#include "AnalysisModel.hh"
#include "ColumnDerivationService.hh"
#include "EventSampleFilterService.hh"
#include "RDataFrameService.hh"
#include "SelectionService.hh"

class DefaultAnalysisModel : public heron::AnalysisModel,
                             public AnalysisConfigService,
                             public ColumnDerivationService,
                             public EventSampleFilterService,
                             public RDataFrameService,
                             public SelectionService
{
  public:
    static const DefaultAnalysisModel &instance();

    const std::string &name() const noexcept override;
    const std::string &tree_name() const noexcept override;
    ProcessorEntry make_processor(const SampleIO::Sample &sample) const noexcept override;

    ROOT::RDF::RNode define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const override;

    const char *filter_stage(SampleIO::SampleOrigin origin) const override;
    ROOT::RDF::RNode apply(ROOT::RDF::RNode node, SampleIO::SampleOrigin origin) const override;

    ROOT::RDataFrame load_sample(const SampleIO::Sample &sample,
                                 const std::string &tree_name) const override;
    ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                      const std::vector<Column> &definitions) const override;

    int slice_required_count() const noexcept override;
    float slice_min_topology_score() const noexcept override;
    float muon_min_track_score() const noexcept override;
    float muon_min_track_length() const noexcept override;
    float muon_max_track_distance() const noexcept override;
    unsigned muon_required_generation() const noexcept override;

    ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, SelectionEntry *selection = NULL) const override;
    ROOT::RDF::RNode decorate(ROOT::RDF::RNode node) const override;
    std::string selection_label(Preset p) const override;
    bool is_in_truth_volume(float x, float y, float z) const noexcept override;
    bool is_in_reco_volume(float x, float y, float z) const noexcept override;

  private:
    DefaultAnalysisModel();

    std::string m_name;
    std::string m_tree_name;
};

#endif // HERON_ANA_DEFAULT_ANALYSIS_MODEL_H
