// plot/macro/plotParticleLevelMipnessSelectionStacked.C
//
// Particle-level stacked histogram of track MIPness, using the cumulative
// selection up to (and including) the reco-fiducial stage by default.
//
// Run with:
//   ./nuxsec macro plotParticleLevelMipnessSelectionStacked.C
//   ./nuxsec macro plotParticleLevelMipnessSelectionStacked.C \
//     'plotParticleLevelMipnessSelectionStacked("./scratch/out/event_list_myana.root")'
//   ./nuxsec macro plotParticleLevelMipnessSelectionStacked.C \
//     'plotParticleLevelMipnessSelectionStacked("", "w_nominal", "", false, true, "sel_trigger && sel_triggered_slice && sel_reco_fv && sel_triggered_muon")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TFile.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &p)
{
    const auto n = p.size();
    if (n < 5 || p.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
    if (!f || f->IsZombie())
        return false;

    const bool has_refs = (f->Get("sample_refs") != nullptr);
    const bool has_events_tree = (f->Get("events") != nullptr);
    const bool has_event_tree_key = (f->Get("event_tree") != nullptr);

    return has_refs && (has_events_tree || has_event_tree_key);
}

bool has_column(const std::vector<std::string> &cols, const std::string &name)
{
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

ROOT::RVec<float> derive_track_mipness_median_plane_score(const ROOT::RVec<float> &T70u,
                                                           const ROOT::RVec<float> &Q50u,
                                                           const ROOT::RVec<float> &Q90u,
                                                           const ROOT::RVec<float> &T70v,
                                                           const ROOT::RVec<float> &Q50v,
                                                           const ROOT::RVec<float> &Q90v,
                                                           const ROOT::RVec<float> &T70y,
                                                           const ROOT::RVec<float> &Q50y,
                                                           const ROOT::RVec<float> &Q90y,
                                                           double Mref_mev_per_cm = 2.12,
                                                           double tail_weight_w = 1.0)
{
    const std::size_t n = std::min({T70u.size(), Q50u.size(), Q90u.size(),
                                    T70v.size(), Q50v.size(), Q90v.size(),
                                    T70y.size(), Q50y.size(), Q90y.size()});

    ROOT::RVec<float> out(n, std::numeric_limits<float>::quiet_NaN());

    auto D = [Mref_mev_per_cm, tail_weight_w](float T70, float Q50, float Q90) -> double {
        if (!(std::isfinite(T70) && std::isfinite(Q50) && std::isfinite(Q90)))
            return std::numeric_limits<double>::quiet_NaN();
        if (!(T70 > 0.0f && Q50 > 0.0f && Q90 > 0.0f))
            return std::numeric_limits<double>::quiet_NaN();
        return std::log(static_cast<double>(T70) / Mref_mev_per_cm) +
               tail_weight_w * std::log(static_cast<double>(Q90) / static_cast<double>(Q50));
    };

    for (std::size_t i = 0; i < n; ++i)
    {
        std::vector<double> ds;
        ds.reserve(3);

        const double Du = D(T70u[i], Q50u[i], Q90u[i]);
        const double Dv = D(T70v[i], Q50v[i], Q90v[i]);
        const double Dy = D(T70y[i], Q50y[i], Q90y[i]);

        if (std::isfinite(Du))
            ds.push_back(Du);
        if (std::isfinite(Dv))
            ds.push_back(Dv);
        if (std::isfinite(Dy))
            ds.push_back(Dy);

        if (ds.empty())
            continue;

        std::sort(ds.begin(), ds.end());
        const double med = ds[ds.size() / 2];
        out[i] = static_cast<float>(std::exp(-med));
    }

    return out;
}

} // namespace

int plotParticleLevelMipnessSelectionStacked(const std::string &event_list_path = "",
                                             const std::string &mc_weight = "w_nominal",
                                             const std::string &extra_sel = "",
                                             bool use_logy = false,
                                             bool include_data = true,
                                             const std::string &selection_upto = "sel_trigger && sel_triggered_slice && sel_reco_fv",
                                             const std::string &particle_pdg_branch = "backtracked_pdg_codes")
{
    ROOT::EnableImplicitMT();

    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;

    if (!looks_like_event_list_root(input_path))
    {
        std::cerr << "[plotParticleLevelMipnessSelectionStacked] input is not an event-list root file: " << input_path << "\n";
        return 1;
    }

    EventListIO el(input_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) {
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf;
    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<size_t>(sid)]);
                                   },
                                   {"sample_id"});
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

    const auto cols = rdf.GetColumnNames();
    const bool has_mipness_scores = has_column(cols, "track_mipness_median_plane_score");
    const bool has_mipness_stats =
        has_column(cols, "track_dedx_T70_u") && has_column(cols, "track_dedx_Q50_u") && has_column(cols, "track_dedx_Q90_u") &&
        has_column(cols, "track_dedx_T70_v") && has_column(cols, "track_dedx_Q50_v") && has_column(cols, "track_dedx_Q90_v") &&
        has_column(cols, "track_dedx_T70_y") && has_column(cols, "track_dedx_Q50_y") && has_column(cols, "track_dedx_Q90_y");

    if (!has_mipness_scores && has_mipness_stats)
    {
        std::cout << "[plotParticleLevelMipnessSelectionStacked] deriving track_mipness_median_plane_score from dE/dx quantiles\n";
        node_mc = node_mc.Define(
            "track_mipness_median_plane_score",
            derive_track_mipness_median_plane_score,
            {"track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
             "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
             "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});
        node_ext = node_ext.Define(
            "track_mipness_median_plane_score",
            derive_track_mipness_median_plane_score,
            {"track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
             "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
             "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});
        node_data = node_data.Define(
            "track_mipness_median_plane_score",
            derive_track_mipness_median_plane_score,
            {"track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
             "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
             "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});
    }
    else if (!has_mipness_scores)
    {
        std::cerr << "[plotParticleLevelMipnessSelectionStacked] missing required MIPness columns. Need either:\n"
                  << "  - track_mipness_median_plane_score\n"
                  << "  - or track_dedx_{T70,Q50,Q90}_{u,v,y}\n";
        return 1;
    }

    std::vector<Entry> entries;
    entries.reserve(include_data ? 3 : 2);

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;
    entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
    Entry &e_mc = entries.back();
    mc.push_back(&e_mc);

    ProcessorEntry rec_ext;
    rec_ext.source = Type::kExt;
    entries.emplace_back(make_entry(std::move(node_ext), rec_ext));
    Entry &e_ext = entries.back();
    mc.push_back(&e_ext);

    Entry *p_data = nullptr;
    if (include_data)
    {
        ProcessorEntry rec_data;
        rec_data.source = Type::kData;
        entries.emplace_back(make_entry(std::move(node_data), rec_data));
        p_data = &entries.back();
        data.push_back(p_data);
    }

    const std::string stage_sel = selection_upto.empty() ? std::string("true") : selection_upto;
    const std::string combined_sel = extra_sel.empty()
                                         ? "(" + stage_sel + ")"
                                         : "(" + stage_sel + ") && (" + extra_sel + ")";

    e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(combined_sel);
    e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(combined_sel);
    if (p_data != nullptr)
        p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(combined_sel);

    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = false;
    opt.show_ratio = include_data;
    opt.show_ratio_band = include_data;
    opt.signal_channels = Channels::signal_keys();
    opt.x_title = "Track MIPness score";
    opt.y_title = "Particles";
    opt.analysis_region_label = "Selection up to MIPness stage";
    opt.image_format = "pdf";
    opt.particle_level = true;
    opt.particle_pdg_branch = particle_pdg_branch;

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();
    opt.run_numbers = {"1"};

    TH1DModel spec;
    spec.id = "particle_level_track_mipness_selection_upto";
    spec.name = "particle_level_track_mipness_selection_upto";
    spec.title = "Track MIPness (particle-level stack)";
    spec.expr = "track_mipness_median_plane_score";
    spec.weight = mc_weight;
    spec.nbins = 50;
    spec.xmin = 0.0;
    spec.xmax = 2.0;
    spec.sel = Preset::Empty;

    std::cout << "[plotParticleLevelMipnessSelectionStacked] selection=" << combined_sel
              << ", include_data=" << (include_data ? "true" : "false")
              << ", use_logy=" << (use_logy ? "true" : "false")
              << ", particle_pdg_branch=" << opt.particle_pdg_branch << "\n";

    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    std::cout << "[plotParticleLevelMipnessSelectionStacked] wrote "
              << plotter.options().out_dir << "/" << spec.id << "."
              << plotter.options().image_format << "\n";

    return 0;
}
