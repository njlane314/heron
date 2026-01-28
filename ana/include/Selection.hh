/* -- C++ -- */
/**
 *  @file  ana/include/Selection.hh
 *
 *  @brief Selection helpers for analysis filters and summaries.
 */

#ifndef NUXSEC_ANA_SELECTION_H
#define NUXSEC_ANA_SELECTION_H

#include <cstddef>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "ColumnDerivationService.hh"

namespace nuxsec
{
namespace selection
{

struct FrameView
{
    ROOT::RDF::RNode node;

    ROOT::RDF::RNode rnode() const { return node; }
};

struct Entry
{
    SourceKind source = SourceKind::kUnknown;
    FrameView nominal;
};

enum class Preset
{
    Empty,
    Trigger,
    Slice,
    Fiducial,
    Topology,
    Muon,
    InclusiveMuCC
};

/** \brief Selection helpers for analysis filters and summaries. */
class SelectionService
{
  public:
    static const float trigger_min_beam_pe;
    static const float trigger_max_veto_pe;

    static const int slice_required_count;
    static const float slice_min_topology_score;

    static const float topology_min_contained_fraction;
    static const float topology_min_cluster_fraction;

    static const float muon_min_track_score;
    static const float muon_min_track_length;
    static const float muon_max_track_distance;
    static const unsigned muon_required_generation;

    static ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, const Entry &rec);
    static bool is_in_truth_volume(float x, float y, float z) noexcept;
    static bool is_in_reco_volume(float x, float y, float z) noexcept;
};

inline ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, const Entry &rec)
{
    return SelectionService::apply(node, p, rec);
}

} // namespace selection
} // namespace nuxsec

#endif // NUXSEC_ANA_SELECTION_H
