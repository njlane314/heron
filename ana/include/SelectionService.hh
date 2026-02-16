/* -- C++ -- */
/**
 *  @file  ana/include/SelectionService.hh
 *
 *  @brief Selection helpers for analysis filters and summaries, describing
 *         selection names, cuts, and bookkeeping outputs.
 */

#ifndef HERON_ANA_SELECTION_SERVICE_H
#define HERON_ANA_SELECTION_SERVICE_H

#include <cstddef>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "ColumnDerivationService.hh"



struct Frame
{
    ROOT::RDF::RNode node;
    ROOT::RDF::RNode rnode() const { return node; }
};

struct SelectionEntry
{
    Type source = Type::kUnknown;
    Frame nominal;
};

enum class Preset
{
    Empty,
    Trigger,
    Slice,
    Fiducial,
    Topology,
    Muon
};

class SelectionService
{
  public:
    static const float trigger_min_beam_pe;
    static const float trigger_max_veto_pe;

    static const int slice_required_count;
    static const float slice_min_topology_score;

    static const float muon_min_track_score;
    static const float muon_min_track_length;
    static const float muon_max_track_distance;
    static const float muon_min_mipness_median_plane_score;
    static const unsigned muon_required_generation;

    static ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, const SelectionEntry &rec);
    static ROOT::RDF::RNode decorate(ROOT::RDF::RNode node, Preset p, const SelectionEntry &rec);
    static ROOT::RDF::RNode decorate(ROOT::RDF::RNode node, const SelectionEntry &rec);
    static std::string selection_label(Preset p);
    static bool is_in_truth_volume(float x, float y, float z) noexcept;
    static bool is_in_reco_volume(float x, float y, float z) noexcept;
};

inline ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, const SelectionEntry &rec)
{
    return SelectionService::apply(node, p, rec);
}



#endif // HERON_ANA_SELECTION_SERVICE_H
