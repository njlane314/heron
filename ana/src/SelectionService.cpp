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

#include "SampleIO.hh"


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

    if (p == Preset::Empty)
        return node;

    if (p == Preset::Trigger || p == Preset::Slice || p == Preset::Fiducial || p == Preset::Muon)
    {
        if (has("beam_mode") && has("run") && has("software_trigger") && has("software_trigger_pre") && has("software_trigger_post"))
        {
            define_if_missing(
                "sel_trigger",
                [](int beam_mode, int run, int sw, int sw_pre, int sw_post) {
                    const int numi_beam_mode = static_cast<int>(SampleIO::BeamMode::kNuMI);
                    if (beam_mode == numi_beam_mode)
                    {
                        constexpr int numi_run_boundary = 16880;
                        if (run < numi_run_boundary)
                            return sw_pre > 0;
                        return sw_post > 0;
                    }
                    return sw > 0;
                },
                {"beam_mode", "run", "software_trigger", "software_trigger_pre", "software_trigger_post"});
        }
        else
        {
            define_if_missing("sel_trigger", [](int sw) { return sw > 0; }, {"software_trigger"});
        }
    }

    if (p == Preset::Slice || p == Preset::Fiducial || p == Preset::Muon)
    {
        define_if_missing(
            "sel_slice",
            [](int ns, float topo) { return passes_slice(ns, topo); },
            {"num_slices", "topological_score"});
    }

    if (p == Preset::Fiducial || p == Preset::Muon)
    {
        define_if_missing(
            "sel_fiducial",
            [](bool slice, bool fv) { return slice && fv; },
            {"sel_slice", "in_reco_fiducial"});
        define_if_missing(
            "sel_topology",
            [](bool fid) { return fid; },
            {"sel_fiducial"});
    }

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
