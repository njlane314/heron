// plot/macro/plotInclusiveMuCCSelectionStages.C
//
// Suite of stage-by-stage diagnostic plots for the “inclusive νμ CC” selection:
//
//   Stage 0: pre-trigger (no selection)
//   Stage 1: trigger
//   Stage 2: trigger + slice
//   Stage 3: trigger + slice + fiducial (reco FV)
//   Stage 4: trigger + muon (final)
//
// For each stage we plot “pre” and “post” distributions of the variables that define
// the next selection step, plus muon-candidate diagnostics and MIPness.
//
// MIPness (from your Fig. 63 / Eqs 7.1–7.3):
//   D_p        = ln(T70_p / Mref) + w * ln(Q90_p / Q50_p)
//   MIPness_p  = exp(-D_p)
//   MIPness_med uses the median of {D_u, D_v, D_y} (over finite planes)
//
// This macro computes candidate-track scalars from the per-event track vectors:
//   - “longest track” (trk_longest_*) for pre-muon diagnostics
//   - “muon-candidate track” (mu_cand_*) using the same thresholds as SelectionService
//
// Run with:
//   ./heron macro plotInclusiveMuCCSelectionStages.C
//   ./heron macro plotInclusiveMuCCSelectionStages.C \
//       'plotInclusiveMuCCSelectionStages("./scratch/out/event_list_myana.root")'
//
//   // With data overlay + log-y
//   ./heron macro plotInclusiveMuCCSelectionStages.C \
//       'plotInclusiveMuCCSelectionStages("./scratch/out/event_list_myana.root","true",true,true)'
//
// Notes / requirements:
//   - Input must be an event-list ROOT file (EventListIO format).
//   - Your event list MUST contain the raw cut variables and TrackAnalysis dE/dx stats used here.
//     See the “required columns” list printed on startup if something is missing.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "MacroGuard.hh"

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

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

std::string getenv_or(const char *key, const std::string &fallback)
{
    const char *v = std::getenv(key);
    if (!v || std::string(v).empty())
        return fallback;
    return std::string(v);
}

std::string sanitize_for_filename(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (const char c : s)
    {
        const bool ok =
            (c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') ||
            (c == '_');
        out.push_back(ok ? c : '_');
    }
    return out;
}

bool has_column(const std::vector<std::string> &cols, const std::string &name)
{
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

ROOT::RDF::RNode define_track_pick_columns(ROOT::RDF::RNode node,
                                          float mu_score_cut,
                                          float mu_len_cut_cm,
                                          float mu_dist_cut_cm,
                                          unsigned mu_req_gen)
{
    using ROOT::RVec;

    // Longest-track index (pre-muon diagnostics): pick max length among finite lengths.
    node = node.Define(
        "trk_longest_idx",
        [](const RVec<float> &lengths) -> int {
            int best = -1;
            float best_len = -std::numeric_limits<float>::infinity();
            for (std::size_t i = 0; i < lengths.size(); ++i)
            {
                const float L = lengths[i];
                if (!std::isfinite(L))
                    continue;
                if (L > best_len)
                {
                    best_len = L;
                    best = static_cast<int>(i);
                }
            }
            return best;
        },
        {"track_length"});

    // Muon-candidate index: same logic as SelectionService::passes_muon, but return the (chosen) track index.
    // If multiple pass, pick the longest.
    node = node.Define(
        "mu_cand_idx",
        [mu_score_cut, mu_len_cut_cm, mu_dist_cut_cm, mu_req_gen](const RVec<float> &scores,
                                                                  const RVec<float> &lengths,
                                                                  const RVec<float> &distances,
                                                                  const RVec<unsigned> &generations) -> int {
            const std::size_t n = std::min(std::min(scores.size(), lengths.size()),
                                           std::min(distances.size(), generations.size()));
            int best = -1;
            float best_len = -std::numeric_limits<float>::infinity();

            for (std::size_t i = 0; i < n; ++i)
            {
                const float s = scores[i];
                const float L = lengths[i];
                const float d = distances[i];

                if (!(std::isfinite(s) && std::isfinite(L) && std::isfinite(d)))
                    continue;

                const bool ok = (s > mu_score_cut) &&
                                (L > mu_len_cut_cm) &&
                                (d < mu_dist_cut_cm) &&
                                (generations[i] == mu_req_gen);

                if (!ok)
                    continue;

                if (L > best_len)
                {
                    best_len = L;
                    best = static_cast<int>(i);
                }
            }
            return best;
        },
        {"track_shower_scores", "track_length", "track_distance_to_vertex", "pfp_generations"});

    // Index -> scalar pickers (float)
    node = node.Define(
        "trk_longest_score",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"trk_longest_idx", "track_shower_scores"});

    node = node.Define(
        "trk_longest_length",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"trk_longest_idx", "track_length"});

    node = node.Define(
        "trk_longest_distance",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"trk_longest_idx", "track_distance_to_vertex"});

    node = node.Define(
        "mu_cand_score",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"mu_cand_idx", "track_shower_scores"});

    node = node.Define(
        "mu_cand_length",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"mu_cand_idx", "track_length"});

    node = node.Define(
        "mu_cand_distance",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"mu_cand_idx", "track_distance_to_vertex"});

    // Index -> scalar pickers (int) for generation diagnostics
    node = node.Define(
        "trk_longest_generation",
        [](int idx, const RVec<unsigned> &v) -> int {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return -1;
            return static_cast<int>(v[static_cast<std::size_t>(idx)]);
        },
        {"trk_longest_idx", "pfp_generations"});

    node = node.Define(
        "mu_cand_generation",
        [](int idx, const RVec<unsigned> &v) -> int {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return -1;
            return static_cast<int>(v[static_cast<std::size_t>(idx)]);
        },
        {"mu_cand_idx", "pfp_generations"});

    return node;
}

ROOT::RDF::RNode define_mipness_columns(ROOT::RDF::RNode node,
                                       double Mref_mev_per_cm,
                                       double tail_weight_w)
{
    using ROOT::RVec;

    auto mipness_from_stats = [Mref_mev_per_cm, tail_weight_w](float T70, float Q50, float Q90) -> float {
        if (!(std::isfinite(T70) && std::isfinite(Q50) && std::isfinite(Q90)))
            return std::numeric_limits<float>::quiet_NaN();
        if (!(T70 > 0.0f && Q50 > 0.0f && Q90 > 0.0f))
            return std::numeric_limits<float>::quiet_NaN();

        const double D = std::log(static_cast<double>(T70) / Mref_mev_per_cm) +
                         tail_weight_w * std::log(static_cast<double>(Q90) / static_cast<double>(Q50));
        return static_cast<float>(std::exp(-D));
    };

    auto D_from_stats = [Mref_mev_per_cm, tail_weight_w](float T70, float Q50, float Q90) -> double {
        if (!(std::isfinite(T70) && std::isfinite(Q50) && std::isfinite(Q90)))
            return std::numeric_limits<double>::quiet_NaN();
        if (!(T70 > 0.0f && Q50 > 0.0f && Q90 > 0.0f))
            return std::numeric_limits<double>::quiet_NaN();
        return std::log(static_cast<double>(T70) / Mref_mev_per_cm) +
               tail_weight_w * std::log(static_cast<double>(Q90) / static_cast<double>(Q50));
    };

    auto mipness_med_at = [mipness_from_stats, D_from_stats](int idx,
                                                             const RVec<float> &T70u, const RVec<float> &Q50u, const RVec<float> &Q90u,
                                                             const RVec<float> &T70v, const RVec<float> &Q50v, const RVec<float> &Q90v,
                                                             const RVec<float> &T70y, const RVec<float> &Q50y, const RVec<float> &Q90y) -> float {
        auto get = [](int i, const RVec<float> &v) -> float {
            if (i < 0 || static_cast<std::size_t>(i) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(i)];
        };

        std::vector<double> Ds;
        Ds.reserve(3);

        const double Du = D_from_stats(get(idx, T70u), get(idx, Q50u), get(idx, Q90u));
        const double Dv = D_from_stats(get(idx, T70v), get(idx, Q50v), get(idx, Q90v));
        const double Dy = D_from_stats(get(idx, T70y), get(idx, Q50y), get(idx, Q90y));

        if (std::isfinite(Du))
            Ds.push_back(Du);
        if (std::isfinite(Dv))
            Ds.push_back(Dv);
        if (std::isfinite(Dy))
            Ds.push_back(Dy);

        if (Ds.empty())
            return std::numeric_limits<float>::quiet_NaN();

        std::sort(Ds.begin(), Ds.end());
        const double Dmed = Ds[Ds.size() / 2];
        return static_cast<float>(std::exp(-Dmed));
    };

    auto mipness_y_at = [mipness_from_stats](int idx,
                                             const RVec<float> &T70y,
                                             const RVec<float> &Q50y,
                                             const RVec<float> &Q90y) -> float {
        if (idx < 0)
            return std::numeric_limits<float>::quiet_NaN();
        if (static_cast<std::size_t>(idx) >= T70y.size() ||
            static_cast<std::size_t>(idx) >= Q50y.size() ||
            static_cast<std::size_t>(idx) >= Q90y.size())
            return std::numeric_limits<float>::quiet_NaN();
        return mipness_from_stats(T70y[static_cast<std::size_t>(idx)],
                                  Q50y[static_cast<std::size_t>(idx)],
                                  Q90y[static_cast<std::size_t>(idx)]);
    };

    // Longest-track MIPness
    node = node.Define(
        "trk_longest_mipness_med",
        mipness_med_at,
        {"trk_longest_idx",
         "track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
         "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
         "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});

    node = node.Define(
        "trk_longest_mipness_y",
        mipness_y_at,
        {"trk_longest_idx", "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});

    // Muon-candidate-track MIPness
    node = node.Define(
        "mu_cand_mipness_med",
        mipness_med_at,
        {"mu_cand_idx",
         "track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
         "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
         "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});

    node = node.Define(
        "mu_cand_mipness_y",
        mipness_y_at,
        {"mu_cand_idx", "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});

    return node;
}

void save_2d(const ROOT::RDF::RNode &node,
             const std::string &x,
             const std::string &y,
             const std::string &w,
             const std::string &sel,
             const std::string &tag,
             int nbx, double xmin, double xmax,
             int nby, double ymin, double ymax)
{
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    const std::string out_dir = getenv_or("HERON_PLOT_DIR", "./scratch/plots");
    const std::string out_fmt = getenv_or("HERON_PLOT_FORMAT", "pdf");
    gSystem->mkdir(out_dir.c_str(), /*recursive=*/true);

    ROOT::RDF::RNode n = node;
    if (!sel.empty())
        n = n.Filter(sel);

    const std::string hname = "h2_" + sanitize_for_filename(tag);
    auto h2 = n.Histo2D({hname.c_str(), "", nbx, xmin, xmax, nby, ymin, ymax}, x, y, w);

    TCanvas c(("c2_" + sanitize_for_filename(tag)).c_str(), "", 900, 800);
    c.SetRightMargin(0.14);
    c.SetLogz(true);

    h2->SetTitle((";" + x + ";" + y).c_str());
    h2->GetZaxis()->SetTitle("Events");
    h2->Draw("COLZ");

    const std::string out_path = out_dir + "/" + sanitize_for_filename(tag) + "." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plotInclusiveMuCCSelectionStages] wrote " << out_path << "\n";
}

} // namespace

int plotInclusiveMuCCSelectionStages(const std::string &samples_tsv = "",
                                     const std::string &extra_sel = "true",
                                     bool use_logy = false,
                                     bool include_data = false,
                                     double Mref_mev_per_cm = 2.10,
                                     double mip_tail_weight = 1.0,
                                     float mu_score_cut = 0.5f,
                                     float mu_len_cut_cm = 10.0f,
                                     float mu_dist_cut_cm = 4.0f,
                                     unsigned mu_req_gen = 2u)
{
  return heron::macro::run_with_guard("plotInclusiveMuCCSelectionStages", [&]() -> int {
    if (implicit_mt_enabled())
        ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotInclusiveMuCCSelectionStages] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotInclusiveMuCCSelectionStages] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    // Check required columns early (fail-fast with a helpful message).
    const std::vector<std::string> cols = rdf.GetColumnNames();
    const std::vector<std::string> required = {
        // selection flags
        "sel_trigger",
        "sel_triggered_slice",
        "sel_reco_fv",
        "sel_triggered_muon",

        // trigger vars
        "optical_filter_pe_beam",
        "optical_filter_pe_veto",
        "software_trigger",

        // slice vars
        "num_slices",
        "topological_score",

        // fiducial diagnostic
        "reco_neutrino_vertex_sce_x",
        "reco_neutrino_vertex_sce_y",
        "reco_neutrino_vertex_sce_z",

        // track vectors used in muon selection
        "track_shower_scores",
        "track_length",
        "track_distance_to_vertex",
        "pfp_generations",

        // TrackAnalysis robust dE/dx stats (for MIPness)
        "track_dedx_T70_u",
        "track_dedx_T70_v",
        "track_dedx_T70_y",
        "track_dedx_Q50_u",
        "track_dedx_Q50_v",
        "track_dedx_Q50_y",
        "track_dedx_Q90_u",
        "track_dedx_Q90_v",
        "track_dedx_Q90_y",

        // plotting / weights
        "analysis_channels",
        "w_nominal"};

    std::vector<std::string> missing;
    for (const auto &c : required)
        if (!has_column(cols, c))
            missing.push_back(c);

    if (!missing.empty())
    {
        std::cerr << "[plotInclusiveMuCCSelectionStages] missing required event-list columns:\n";
        for (const auto &c : missing)
            std::cerr << "  - " << c << "\n";
        std::cerr << "Rebuild the event list including these branches (default schema: framework/core/config/event_columns.tsv).\n";
        return 2;
    }

    // Split into sources (same pattern as plotStackedHistTrueVertex)
    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) {
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<std::size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf;
    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<std::size_t>(sid)]);
                                   },
                                   {"sample_id"});
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

    if (!extra_sel.empty())
    {
        node_mc = node_mc.Filter(extra_sel);
        node_ext = node_ext.Filter(extra_sel);
        node_data = node_data.Filter(extra_sel);
    }

    // Add candidate-track scalars and MIPness columns.
    node_mc = define_track_pick_columns(node_mc, mu_score_cut, mu_len_cut_cm, mu_dist_cut_cm, mu_req_gen);
    node_ext = define_track_pick_columns(node_ext, mu_score_cut, mu_len_cut_cm, mu_dist_cut_cm, mu_req_gen);
    node_data = define_track_pick_columns(node_data, mu_score_cut, mu_len_cut_cm, mu_dist_cut_cm, mu_req_gen);

    node_mc = define_mipness_columns(node_mc, Mref_mev_per_cm, mip_tail_weight);
    node_ext = define_mipness_columns(node_ext, Mref_mev_per_cm, mip_tail_weight);
    node_data = define_mipness_columns(node_data, Mref_mev_per_cm, mip_tail_weight);

    // Configure plotter
    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = true;
    opt.show_ratio = include_data;
    opt.show_ratio_band = include_data;
    opt.signal_channels = Channels::signal_keys();
    opt.y_title = "Events";
    opt.image_format = getenv_or("HERON_PLOT_FORMAT", "pdf");

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();
    opt.run_numbers = {"1"};

    const std::string weight = "w_nominal";

    // Draw helper: build temporary Entries with the requested selection and draw a stack.
    const auto draw_one = [&](const std::string &plot_tag,
                              const std::string &sel,
                              const std::string &expr,
                              int nbins,
                              double xmin,
                              double xmax,
                              const std::string &x_title,
                              bool force_logy = false) {
        ROOT::RDF::RNode mc_sel = node_mc;
        ROOT::RDF::RNode ext_sel = node_ext;
        ROOT::RDF::RNode data_sel = node_data;

        if (!sel.empty())
        {
            mc_sel = mc_sel.Filter(sel);
            ext_sel = ext_sel.Filter(sel);
            if (include_data)
                data_sel = data_sel.Filter(sel);
        }

        // Per-plot log-y override (keep global setting for everything else).
        const bool old_logy = opt.use_log_y;
        opt.use_log_y = (use_logy || force_logy);

        // Build stage-local Entry list
        std::vector<Entry> entries;
        entries.reserve(include_data ? 3 : 2);

        std::vector<const Entry *> mc;
        std::vector<const Entry *> data;

        ProcessorEntry rec_mc;
        rec_mc.source = Type::kMC;

        ProcessorEntry rec_ext;
        rec_ext.source = Type::kExt;

        entries.emplace_back(make_entry(std::move(mc_sel), rec_mc));
        mc.push_back(&entries.back());

        entries.emplace_back(make_entry(std::move(ext_sel), rec_ext));
        mc.push_back(&entries.back());

        if (include_data)
        {
            ProcessorEntry rec_data;
            rec_data.source = Type::kData;
            entries.emplace_back(make_entry(std::move(data_sel), rec_data));
            data.push_back(&entries.back());
        }

        opt.x_title = x_title.empty() ? expr : x_title;

        // Stage-specific filename prefix via the "name" passed to make_spec:
        // we keep spec.expr as the real draw expression.
        const std::string spec_name = plot_tag + "__" + expr;
        TH1DModel spec = make_spec(spec_name, nbins, xmin, xmax, weight);
        spec.expr = expr;
        spec.sel = Preset::Empty;

        if (include_data)
            plotter.draw_stack(spec, mc, data);
        else
            plotter.draw_stack(spec, mc);

        opt.use_log_y = old_logy;
    };

    // -------------------------------------------------------------------------
    // Stage selections (cumulative)
    // -------------------------------------------------------------------------
    const std::string sel_pre_trigger = "true";
    const std::string sel_post_trigger = "sel_trigger";

    const std::string sel_pre_slice = "sel_trigger";
    const std::string sel_post_slice = "sel_triggered_slice";

    const std::string sel_pre_fv = "sel_triggered_slice";
    const std::string sel_post_fv = "sel_triggered_slice && sel_reco_fv";

    const std::string sel_pre_mu = "sel_triggered_slice && sel_reco_fv";
    const std::string sel_post_mu = "sel_triggered_muon";

    // -------------------------------------------------------------------------
    // Trigger diagnostics
    // -------------------------------------------------------------------------
    draw_one("trigger_pre", sel_pre_trigger, "optical_filter_pe_beam", 50, -10.0, 500.0, "Beam optical PE");
    draw_one("trigger_post", sel_post_trigger, "optical_filter_pe_beam", 50, -10.0, 500.0, "Beam optical PE");

    draw_one("trigger_pre", sel_pre_trigger, "optical_filter_pe_veto", 50, -10.0, 200.0, "Veto optical PE");
    draw_one("trigger_post", sel_post_trigger, "optical_filter_pe_veto", 50, -10.0, 200.0, "Veto optical PE");

    draw_one("trigger_pre", sel_pre_trigger, "software_trigger", 6, -0.5, 5.5, "Software trigger");
    draw_one("trigger_post", sel_post_trigger, "software_trigger", 6, -0.5, 5.5, "Software trigger");

    // -------------------------------------------------------------------------
    // Slice diagnostics
    // -------------------------------------------------------------------------
    draw_one("slice_pre", sel_pre_slice, "num_slices", 8, -0.5, 7.5, "Number of Pandora slices");
    draw_one("slice_post", sel_post_slice, "num_slices", 8, -0.5, 7.5, "Number of Pandora slices");

    draw_one("slice_pre", sel_pre_slice, "topological_score", 50, 0.0, 1.0, "Topological score", /*force_logy=*/true);
    draw_one("slice_post", sel_post_slice, "topological_score", 50, 0.0, 1.0, "Topological score", /*force_logy=*/true);

    // -------------------------------------------------------------------------
    // Fiducial diagnostics (reco SCE vertex)
    // -------------------------------------------------------------------------
    draw_one("fv_pre", sel_pre_fv, "reco_neutrino_vertex_sce_z", 50, -50.0, 1100.0, "Reco ν vertex z [cm]");
    draw_one("fv_post", sel_post_fv, "reco_neutrino_vertex_sce_z", 50, -50.0, 1100.0, "Reco ν vertex z [cm]");

    draw_one("fv_pre", sel_pre_fv, "reco_neutrino_vertex_sce_x", 50, -50.0, 300.0, "Reco ν vertex x [cm]");
    draw_one("fv_post", sel_post_fv, "reco_neutrino_vertex_sce_x", 50, -50.0, 300.0, "Reco ν vertex x [cm]");

    draw_one("fv_pre", sel_pre_fv, "reco_neutrino_vertex_sce_y", 50, -180.0, 180.0, "Reco ν vertex y [cm]");
    draw_one("fv_post", sel_post_fv, "reco_neutrino_vertex_sce_y", 50, -180.0, 180.0, "Reco ν vertex y [cm]");

    // -------------------------------------------------------------------------
    // Muon-candidate diagnostics
    //   Pre: use “longest track” scalars (exists even when event fails muon selection)
    //   Post: use “mu_cand_*” scalars (track that passes the muon cuts)
    // -------------------------------------------------------------------------

    const std::string sel_pre_mu_with_track = sel_pre_mu + " && (trk_longest_idx >= 0)";
    const std::string sel_post_mu_with_track = sel_post_mu + " && (mu_cand_idx >= 0)";

    // Track-likeness proxy (your shower-vs-track score)
    draw_one("mu_pre", sel_pre_mu_with_track, "trk_longest_score", 50, 0.0, 1.0, "Longest-track shower→track score");
    draw_one("mu_post", sel_post_mu_with_track, "mu_cand_score", 50, 0.0, 1.0, "Muon-candidate shower→track score");

    // Geometry/proximity
    draw_one("mu_pre", sel_pre_mu_with_track, "trk_longest_length", 50, 0.0, 800.0, "Longest-track length [cm]");
    draw_one("mu_post", sel_post_mu_with_track, "mu_cand_length", 50, 0.0, 800.0, "Muon-candidate length [cm]");

    draw_one("mu_pre", sel_pre_mu_with_track, "trk_longest_distance", 50, 0.0, 50.0, "Longest-track start distance to ν vtx [cm]");
    draw_one("mu_post", sel_post_mu_with_track, "mu_cand_distance", 50, 0.0, 50.0, "Muon-candidate start distance to ν vtx [cm]");

    // PFP generation (explicit muon-selection requirement; previously missing)
    draw_one("mu_pre", sel_pre_mu_with_track, "trk_longest_generation", 6, -0.5, 5.5, "Longest-track PFP generation");
    draw_one("mu_post", sel_post_mu_with_track, "mu_cand_generation", 6, -0.5, 5.5, "Muon-candidate PFP generation");

    // MIPness (Y plane and median-over-planes)
    draw_one("mu_pre", sel_pre_mu_with_track + " && (trk_longest_mipness_y == trk_longest_mipness_y)",
             "trk_longest_mipness_y", 50, 0.0, 2.0, "Longest-track MIPness (Y plane)");

    draw_one("mu_post", sel_post_mu_with_track + " && (mu_cand_mipness_y == mu_cand_mipness_y)",
             "mu_cand_mipness_y", 50, 0.0, 2.0, "Muon-candidate MIPness (Y plane)");

    draw_one("mu_pre", sel_pre_mu_with_track + " && (trk_longest_mipness_med == trk_longest_mipness_med)",
             "trk_longest_mipness_med", 50, 0.0, 2.0, "Longest-track MIPness (median planes)");

    draw_one("mu_post", sel_post_mu_with_track + " && (mu_cand_mipness_med == mu_cand_mipness_med)",
             "mu_cand_mipness_med", 50, 0.0, 2.0, "Muon-candidate MIPness (median planes)");

    // Optional: 2D correlations (MC only), written directly to HERON_PLOT_DIR
    save_2d(node_mc,
            "trk_longest_score",
            "trk_longest_mipness_med",
            "w_nominal",
            sel_pre_mu_with_track + " && (trk_longest_mipness_med == trk_longest_mipness_med)",
            "mipness2d_pre_longest_score_vs_mipnessmed",
            50, 0.0, 1.0,
            50, 0.0, 2.0);

    save_2d(node_mc,
            "mu_cand_score",
            "mu_cand_mipness_med",
            "w_nominal",
            sel_post_mu_with_track + " && (mu_cand_mipness_med == mu_cand_mipness_med)",
            "mipness2d_post_mucand_score_vs_mipnessmed",
            50, 0.0, 1.0,
            50, 0.0, 2.0);

    return 0;

  });
}
