/* -- C++ -- */
/**
 *  @file  ana/src/SelectionService.cc
 *
 *  @brief Selection helpers for analysis filters and summaries.
 */

#include "Selection.hh"

namespace nuxsec
{
namespace selection
{

const float SelectionService::trigger_min_beam_pe = 0.f;
const float SelectionService::trigger_max_veto_pe = 20.f;

const int SelectionService::slice_required_count = 1;
const float SelectionService::slice_min_topology_score = 0.06f;

const float SelectionService::topology_min_contained_fraction = 0.0f;
const float SelectionService::topology_min_cluster_fraction = 0.5f;

const float SelectionService::muon_min_track_score = 0.5f;
const float SelectionService::muon_min_track_length = 10.0f;
const float SelectionService::muon_max_track_distance = 4.0f;
const unsigned SelectionService::muon_required_generation = 2u;

namespace
{

constexpr float min_x = 5.f;
constexpr float max_x = 251.f;
constexpr float min_y = -110.f;
constexpr float max_y = 110.f;
constexpr float min_z = 20.f;
constexpr float max_z = 986.f;

constexpr float reco_gap_min_z = 675.f;
constexpr float reco_gap_max_z = 775.f;

template <typename T>
bool is_within(const T &value, float low, float high)
{
    return value > low && value < high;
}

template <typename X, typename Y, typename Z>
bool is_in_active_volume(const X &x, const Y &y, const Z &z)
{
    return is_within(x, min_x, max_x) &&
           is_within(y, min_y, max_y) &&
           is_within(z, min_z, max_z);
}

}

ROOT::RDF::RNode SelectionService::apply(ROOT::RDF::RNode node, Preset p, const Entry &rec)
{
    switch (p)
    {
    case Preset::Empty:
        return node;
    case Preset::Trigger:
        return node.Filter(
            [src = rec.source](float pe_beam, float pe_veto, int sw) {
                const bool requires_dataset_gate = (src == SourceKind::kMC);
                const bool dataset_gate = requires_dataset_gate
                                              ? (pe_beam > trigger_min_beam_pe &&
                                                 pe_veto < trigger_max_veto_pe && sw > 0)
                                              : true;
                return dataset_gate;
            },
            {"optical_filter_pe_beam", "optical_filter_pe_veto", "software_trigger"});
    case Preset::Slice:
        return node.Filter(
            [](int ns, float topo) {
                return ns == slice_required_count && topo > slice_min_topology_score;
            },
            {"num_slices", "topological_score"});
    case Preset::Fiducial:
        return node.Filter([](bool fv) { return fv; }, {"in_reco_fiducial"});
    case Preset::Topology:
        return node.Filter(
            [](float cf, float cl) {
                return cf >= topology_min_contained_fraction &&
                       cl >= topology_min_cluster_fraction;
            },
            {"contained_fraction", "slice_cluster_fraction"});
    case Preset::Muon:
        return node.Filter(
            [](const ROOT::RVec<float> &scores,
               const ROOT::RVec<float> &lengths,
               const ROOT::RVec<float> &distances,
               const ROOT::RVec<unsigned> &generations) {
                const auto n = scores.size();
                for (std::size_t i = 0; i < n; ++i)
                {
                    const bool passes = scores[i] > muon_min_track_score &&
                                        lengths[i] > muon_min_track_length &&
                                        distances[i] < muon_max_track_distance &&
                                        generations[i] == muon_required_generation;
                    if (passes)
                    {
                        return true;
                    }
                }
                return false;
            },
            {"track_shower_scores",
             "track_length",
             "track_distance_to_vertex",
             "pfp_generations"});
    case Preset::InclusiveMuCC:
    default:
    {
        auto filtered = apply(node, Preset::Trigger, rec);
        filtered = apply(filtered, Preset::Slice, rec);
        filtered = apply(filtered, Preset::Fiducial, rec);
        filtered = apply(filtered, Preset::Topology, rec);
        return apply(filtered, Preset::Muon, rec);
    }
    }
}

bool SelectionService::is_in_truth_volume(float x, float y, float z) noexcept
{
    return is_in_active_volume(x, y, z);
}

bool SelectionService::is_in_reco_volume(float x, float y, float z) noexcept
{
    return is_in_active_volume(x, y, z) && (z < reco_gap_min_z || z > reco_gap_max_z);
}

} // namespace selection
} // namespace nuxsec
