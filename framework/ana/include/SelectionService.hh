/* -- C++ -- */
/**
 *  @file  ana/include/SelectionService.hh
 *
 *  @brief Selection helpers for analysis filters and summaries, describing
 *         selection names, cuts, and bookkeeping outputs.
 */

#ifndef HERON_ANA_SELECTION_SERVICE_H
#define HERON_ANA_SELECTION_SERVICE_H

#include <string>

#include <ROOT/RDataFrame.hxx>

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
    Muon
};

class SelectionService
{
  public:
    virtual ~SelectionService() = default;

    virtual int slice_required_count() const noexcept = 0;
    virtual float slice_min_topology_score() const noexcept = 0;
    virtual float muon_min_track_score() const noexcept = 0;
    virtual float muon_min_track_length() const noexcept = 0;
    virtual float muon_max_track_distance() const noexcept = 0;
    virtual unsigned muon_required_generation() const noexcept = 0;

    virtual ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, SelectionEntry *selection = NULL) const = 0;
    virtual ROOT::RDF::RNode decorate(ROOT::RDF::RNode node) const = 0;
    virtual std::string selection_label(Preset p) const = 0;
    virtual bool is_in_truth_volume(float x, float y, float z) const noexcept = 0;
    virtual bool is_in_reco_volume(float x, float y, float z) const noexcept = 0;
};


#endif // HERON_ANA_SELECTION_SERVICE_H
