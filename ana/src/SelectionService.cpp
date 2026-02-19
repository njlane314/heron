/* -- C++ -- */
/**
 *  @file  ana/src/SelectionService.cpp
 *
 *  @brief Selection helpers for analysis filters and summaries.
 */

#include "SelectionService.hh"

#include <algorithm>
#include <string>
#include <vector>


const float SelectionService::trigger_min_beam_pe = 0.f;
const float SelectionService::trigger_max_veto_pe = 20.f;

const int SelectionService::slice_required_count = 1;
const float SelectionService::slice_min_topology_score = 0.06f;

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

inline bool passes_trigger(Type src, float pe_beam, float pe_veto, int sw)
{
    const bool requires_dataset_gate = (src == Type::kMC);
    (void)pe_beam;
    (void)pe_veto;
    return requires_dataset_gate ? (sw > 0) : true;
}

inline bool passes_slice(int ns, float topo)
{
    return ns == SelectionService::slice_required_count &&
           topo > SelectionService::slice_min_topology_score;
}

inline bool passes_muon(const ROOT::RVec<float> &scores,
                        const ROOT::RVec<float> &lengths,
                        const ROOT::RVec<float> &distances,
                        const ROOT::RVec<unsigned> &generations)
{
    const auto n = scores.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i >= lengths.size() || i >= distances.size() || i >= generations.size())
            return false;
        const bool ok = scores[i] > SelectionService::muon_min_track_score &&
                        lengths[i] > SelectionService::muon_min_track_length &&
                        distances[i] < SelectionService::muon_max_track_distance &&
                        generations[i] == SelectionService::muon_required_generation;
        if (ok)
            return true;
    }
    return false;
}

inline ROOT::RDF::RNode filter_on(ROOT::RDF::RNode node, const char *col)
{
    return node.Filter([](bool pass) { return pass; }, {col});
}

} // namespace

ROOT::RDF::RNode SelectionService::apply(ROOT::RDF::RNode node, Preset p, const SelectionEntry &rec)
{
    node = decorate(node, p, rec);
    switch (p)
    {
    case Preset::Empty:
        return node;
    case Preset::Trigger:
        return filter_on(node, "sel_trigger");
    case Preset::Slice:
        return filter_on(node, "sel_slice");
    case Preset::Fiducial:
        return filter_on(node, "sel_fiducial");
    case Preset::Topology:
        return filter_on(node, "sel_topology");
    case Preset::Muon:
        return filter_on(node, "sel_muon");
    default:
        return node;
    }
}

ROOT::RDF::RNode SelectionService::decorate(ROOT::RDF::RNode node, Preset p, const SelectionEntry &rec)
{
    std::vector<std::string> names = node.GetColumnNames();
    auto has = [&](const std::string &name) {
        return std::find(names.begin(), names.end(), name) != names.end();
    };
    auto define_if_missing = [&](const char *name, auto &&f, std::initializer_list<const char *> deps) {
        if (has(name))
            return;
        const std::vector<std::string> columns{deps.begin(), deps.end()};
        node = node.Define(name, std::forward<decltype(f)>(f), columns);
        names.emplace_back(name);
    };

    switch (p)
    {
    case Preset::Empty:
        return node;
    case Preset::Trigger:
        define_if_missing(
            "sel_trigger",
            [src = rec.source](float pe_beam, float pe_veto, int sw) {
                return passes_trigger(src, pe_beam, pe_veto, sw);
            },
            {"optical_filter_pe_beam", "optical_filter_pe_veto", "software_trigger"});
        return node;
    case Preset::Slice:
        define_if_missing(
            "sel_slice",
            [](int ns, float topo) { return passes_slice(ns, topo); },
            {"num_slices", "topological_score"});
        return node;
    case Preset::Fiducial:
    case Preset::Topology:
    case Preset::Muon:
        define_if_missing(
            "sel_slice",
            [](int ns, float topo) { return passes_slice(ns, topo); },
            {"num_slices", "topological_score"});
        define_if_missing(
            "sel_fiducial",
            [](bool slice, bool fv) { return slice && fv; },
            {"sel_slice", "in_reco_fiducial"});
        define_if_missing(
            "sel_topology",
            [](bool fid) { return fid; },
            {"sel_fiducial"});
        if (p == Preset::Muon)
        {
            define_if_missing(
                "sel_muon",
                [](bool topo,
                   const ROOT::RVec<float> &scores,
                   const ROOT::RVec<float> &lengths,
                   const ROOT::RVec<float> &distances,
                   const ROOT::RVec<unsigned> &generations) {
                    if (!topo)
                        return false;
                    return passes_muon(scores, lengths, distances, generations);
                },
                {"sel_topology",
                 "track_shower_scores",
                 "track_length",
                 "track_distance_to_vertex",
                 "pfp_generations"});
        }
        return node;
    default:
        return node;
    }
}

ROOT::RDF::RNode SelectionService::decorate(ROOT::RDF::RNode node, const SelectionEntry &rec)
{
    node = decorate(node, Preset::Trigger, rec);
    node = decorate(node, Preset::Muon, rec);

    std::vector<std::string> names = node.GetColumnNames();
    auto has = [&](const std::string &name) {
        return std::find(names.begin(), names.end(), name) != names.end();
    };
    auto define_if_missing = [&](const char *name, auto &&f, std::initializer_list<const char *> deps) {
        if (has(name))
            return;
        const std::vector<std::string> columns{deps.begin(), deps.end()};
        node = node.Define(name, std::forward<decltype(f)>(f), columns);
        names.emplace_back(name);
    };

    define_if_missing(
        "sel_inclusive_mu_cc",
        [](bool mu) { return mu; },
        {"sel_muon"});
    define_if_missing("sel_reco_fv", [](bool fv) { return fv; }, {"in_reco_fiducial"});
    define_if_missing("sel_triggered_slice", [](bool t, bool s) { return t && s; }, {"sel_trigger", "sel_slice"});
    define_if_missing("sel_triggered_muon", [](bool t, bool m) { return t && m; }, {"sel_trigger", "sel_muon"});
    return node;
}

std::string SelectionService::selection_label(Preset p)
{
    switch (p)
    {
    case Preset::Trigger:
        return "Trigger Selection";
    case Preset::Slice:
        return "Slice Selection";
    case Preset::Fiducial:
        return "Fiducial Selection";
    case Preset::Topology:
        return "Topology Selection";
    case Preset::Muon:
        return "Muon Selection";
    case Preset::Empty:
    default:
        return "Empty Selection";
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
