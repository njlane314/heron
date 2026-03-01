#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "SampleCLI.hh"
#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "PlotEnv.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"

using namespace nu;

namespace {
constexpr const char *kEventWeight = "w_nominal";

int env_int(const char *name, int fallback) {
    const char *v = gSystem->Getenv(name);
    if (!v || !*v)
        return fallback;

    char *end = nullptr;
    const long out = std::strtol(v, &end, 10);
    if (end == v || (end != nullptr && *end != '\0'))
        return fallback;
    return static_cast<int>(out);
}

bool debug_enabled() {
    return env_int("HERON_DEBUG_ADAPTIVE_BINNING", 0) != 0;
}

void stage_log(const std::string &msg) {
    std::cout << "[plot_model_logit_adaptive_binning] " << msg << "\n";
    std::cout.flush();
}

void debug_log(const std::string &msg) {
    if (!debug_enabled())
        return;
    stage_log(msg);
}

std::shared_ptr<const std::vector<char>>
subtract_masks(std::shared_ptr<const std::vector<char>> include,
               std::shared_ptr<const std::vector<char>> exclude) {
    const size_t n = std::max(include ? include->size() : size_t{0},
                              exclude ? exclude->size() : size_t{0});
    auto out = std::make_shared<std::vector<char>>(n, 0);

    for (size_t i = 0; i < n; ++i) {
        const bool in =
            include && i < include->size() && (*include)[i];
        const bool ex =
            exclude && i < exclude->size() && (*exclude)[i];
        (*out)[i] = (in && !ex) ? 1 : 0;
    }

    return out;
}

struct ScoreRange {
    double lo = 0.0;
    double hi = 0.0;
    bool valid = false;
};

ScoreRange find_score_range(ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext) {
    ScoreRange out;

    debug_log("find_score_range: booking min/max/count actions");

    auto min_mc = node_mc.Min<float>("inf_score_0");
    auto max_mc = node_mc.Max<float>("inf_score_0");
    auto min_ext = node_ext.Min<float>("inf_score_0");
    auto max_ext = node_ext.Max<float>("inf_score_0");
    auto n_mc = node_mc.Count();
    auto n_ext = node_ext.Count();

    const auto nsel_mc = n_mc.GetValue();
    const auto nsel_ext = n_ext.GetValue();
    const bool has_mc = (nsel_mc > 0);
    const bool has_ext = (nsel_ext > 0);
    if (!has_mc && !has_ext)
        return out;

    out.valid = true;
    if (has_mc && has_ext) {
        out.lo = std::min(static_cast<double>(min_mc.GetValue()),
                          static_cast<double>(min_ext.GetValue()));
        out.hi = std::max(static_cast<double>(max_mc.GetValue()),
                          static_cast<double>(max_ext.GetValue()));
    } else if (has_mc) {
        out.lo = min_mc.GetValue();
        out.hi = max_mc.GetValue();
    } else {
        out.lo = min_ext.GetValue();
        out.hi = max_ext.GetValue();
    }

    if (!std::isfinite(out.lo) || !std::isfinite(out.hi))
        out.valid = false;
    if (out.valid && out.hi <= out.lo) {
        const double span = (std::abs(out.lo) > 0.0 ? 0.05 * std::abs(out.lo)
                                                    : 1.0);
        out.lo -= span;
        out.hi += span;
    }

    debug_log("find_score_range: selected entries mc=" +
              std::to_string(nsel_mc) + ", ext=" + std::to_string(nsel_ext) +
              ", range=[" + std::to_string(out.lo) + ", " +
              std::to_string(out.hi) + "]");

    return out;
}

void draw_stack_plots(Plotter &plotter, std::vector<const Entry *> &mc,
                      std::vector<const Entry *> &data, bool include_data,
                      const std::vector<double> &adaptive_edges) {
    auto &opt = plotter.options();

    TH1DModel spec = make_spec("inf_score_0", adaptive_edges, kEventWeight);
    spec.sel = Preset::Empty;
    opt.x_title = "Inference score [0]";
    debug_log("draw_stack_plots: drawing with " +
              std::to_string(std::max(0,
                                      static_cast<int>(adaptive_edges.size()) -
                                          1)) +
              " adaptive bins");
    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);
}

std::vector<double> make_adaptive_score_bins(
    ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext,
    const std::string &signal_sel, int n_fine_bins,
    double nmin_signal, double nmin_background, int max_bins,
    double score_lo, double score_hi) {
    std::vector<double> edges;
    if (n_fine_bins < 10)
        n_fine_bins = 10;
    if (max_bins < 1)
        max_bins = 1;
    edges.reserve(static_cast<size_t>(max_bins) + 1);
    stage_log("adaptive binning: begin"
              " (n_fine_bins=" + std::to_string(n_fine_bins) +
              ", nmin_signal=" + std::to_string(nmin_signal) +
              ", nmin_background=" + std::to_string(nmin_background) +
              ", max_bins=" + std::to_string(max_bins) + ")");

    ROOT::RDF::TH1DModel model_sig_w("h_sig_w", "", n_fine_bins, score_lo,
                                     score_hi);
    ROOT::RDF::TH1DModel model_bkg_mc_w("h_bkg_mc_w", "", n_fine_bins,
                                        score_lo, score_hi);
    ROOT::RDF::TH1DModel model_ext_w("h_ext_w", "", n_fine_bins, score_lo,
                                     score_hi);

    auto node_mc_tagged = node_mc.Define("__is_sig_bin__", signal_sel);

    auto h_sig_w = node_mc_tagged.Filter(
                           [](bool is_sig) { return is_sig; },
                           {"__is_sig_bin__"})
                       .Histo1D(model_sig_w, "inf_score_0", "__w__");
    auto h_bkg_mc_w =
        node_mc_tagged.Filter(
                         [](bool is_sig) { return !is_sig; },
                         {"__is_sig_bin__"})
            .Histo1D(model_bkg_mc_w, "inf_score_0", "__w__");
    auto h_ext_w = node_ext.Histo1D(model_ext_w, "inf_score_0", "__w__");

    debug_log("adaptive binning: forcing MC fine histograms");
    (void)h_sig_w->GetEntries();
    (void)h_bkg_mc_w->GetEntries();
    debug_log("adaptive binning: forcing EXT fine histogram");
    (void)h_ext_w->GetEntries();

    const double width =
        (score_hi - score_lo) / static_cast<double>(n_fine_bins);
    edges.push_back(score_lo);

    double acc_sw = 0.0;
    double acc_bw = 0.0;

    for (int bin = n_fine_bins; bin >= 1; --bin) {
        acc_sw += h_sig_w->GetBinContent(bin);
        acc_bw += h_bkg_mc_w->GetBinContent(bin) + h_ext_w->GetBinContent(bin);

        // Apply adaptive thresholds to the expected yields in the current
        // coarse bin. Background includes both non-signal MC and EXT.
        const bool pass = (acc_sw >= nmin_signal && acc_bw >= nmin_background);
        const int n_splits = static_cast<int>(edges.size()) - 1;
        const bool can_add_split = (n_splits < max_bins - 1);

        // We scan from high score to low score, so when a bin closes the new
        // boundary is the low edge of the accumulated interval, not the high
        // edge of the current fine bin. Once we already have max_bins - 1
        // interior splits, keep merging the remaining low-score region into
        // the final bin.
        if (bin == 1 || (pass && can_add_split)) {
            const double edge =
                score_lo + static_cast<double>(bin - 1) * width;
            edges.push_back(edge);
            acc_sw = 0.0;
            acc_bw = 0.0;
        }
    }

    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    if (edges.empty() || edges.front() > score_lo)
        edges.insert(edges.begin(), score_lo);
    if (edges.back() < score_hi)
        edges.push_back(score_hi);

    stage_log("adaptive binning: built " +
              std::to_string(std::max(0, static_cast<int>(edges.size()) - 1)) +
              " bins");

    return edges;
}

} // namespace

int plot_model_logit_adaptive_binning(
    const std::string &samples_tsv = "",
    const std::string &extra_sel = "sel_muon", bool use_logy = true,
    bool include_data = false, const std::string &signal_sel = "is_signal",
    const std::string &mc_weight = "w_nominal",
    double raw_threshold_min = -15.0, double raw_threshold_max = 15.0,
    int n_fine_bins = 200, double nmin_signal = 10.0,
    double nmin_background = 30.0, int max_bins = 15) {
    (void)samples_tsv;
    (void)mc_weight;
    const std::string list_path =
        "/exp/uboone/data/users/nlane/heron/out/event/events_muon.root";
    stage_log("input=" + list_path);

    if (!looks_like_event_list_root(list_path)) {
        std::cerr
            << "[plot_model_logit_adaptive_binning] input is not an event "
               "list ROOT file: "
            << list_path << "\n";
        return 2;
    }

    const int n_rdf_threads = std::max(1, env_int("HERON_RDF_THREADS", 1));
    if (n_rdf_threads > 1) {
        ROOT::EnableImplicitMT(static_cast<unsigned>(n_rdf_threads));
        stage_log("RDataFrame threads=" + std::to_string(n_rdf_threads));
    } else {
        stage_log("RDataFrame threads=1 (ImplicitMT disabled)");
    }
    gErrorIgnoreLevel = kWarning;

    stage_log("opening event list");
    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();
    auto mask_pure_mc = subtract_masks(mask_mc, mask_ext);

    const bool has_analysis_channels = rdf.HasColumn("analysis_channels");
    const bool has_is_signal = rdf.HasColumn("is_signal");

    ROOT::RDF::RNode base =
        rdf.Define("inf_score_0",
                   [](const ROOT::RVec<float> &scores) {
                       if (scores.empty())
                           return -1.0f;
                       return scores[0];
                   },
                   {"inf_scores"});

    if (has_analysis_channels) {
        base = base.Define(
            "is_signal_label",
            [](int analysis_channel) {
                return (analysis_channel == 15 || analysis_channel == 16) ? 1
                                                                          : 0;
            },
            {"analysis_channels"});
    } else if (has_is_signal) {
        base = base.Define("is_signal_label",
                           [](bool is_signal) { return is_signal ? 1 : 0; },
                           {"is_signal"});
    } else {
        stage_log("neither analysis_channels nor is_signal exists; defaulting "
                  "is_signal_label=false");
        base = base.Define("is_signal_label", []() { return 0; });
    }

    ROOT::RDF::RNode selected = base;
    if (!extra_sel.empty()) {
        const bool named_column = rdf.HasColumn(extra_sel);
        stage_log("applying selection '" + extra_sel +
                  (named_column ? "' [column]" : "' [expression]"));
        if (named_column) {
            selected = selected.Filter(
                [](bool pass) { return pass; }, {extra_sel});
        } else {
            selected = selected.Filter(extra_sel);
        }
    }

    ROOT::RDF::RNode node_ext =
        filter_by_sample_mask(selected, mask_ext).Define("__w__", kEventWeight);
    ROOT::RDF::RNode node_mc =
        filter_by_sample_mask(selected, mask_pure_mc).Define("__w__", kEventWeight);
    ROOT::RDF::RNode node_data =
        filter_by_sample_mask(selected, mask_data);

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
    opt.run_numbers = {"1"};
    opt.image_format = "pdf";

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    stage_log("computing score range");
    ScoreRange score_range =
        find_score_range(node_mc, node_ext);
    if (!score_range.valid) {
        score_range.lo = raw_threshold_min;
        score_range.hi = raw_threshold_max;
    }

    stage_log("score range: [" + std::to_string(score_range.lo) + ", " +
              std::to_string(score_range.hi) + "]");

    const std::string adaptive_signal_sel =
        (signal_sel == "is_signal" && !has_is_signal)
            ? "is_signal_label"
            : signal_sel;
    stage_log("computing adaptive edges");
    const std::vector<double> adaptive_edges = make_adaptive_score_bins(
        node_mc, node_ext,
        adaptive_signal_sel, n_fine_bins, nmin_signal, nmin_background,
        max_bins, score_range.lo, score_range.hi);
    std::cout << "[plot_model_logit_adaptive_binning] adaptive edges:";
    for (double e : adaptive_edges)
        std::cout << " " << e;
    std::cout << "\n";
    std::cout.flush();

    stage_log("building plot entries");
    std::vector<Entry> entries;
    entries.reserve(include_data ? 3 : 2);

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;
    entries.emplace_back(make_entry(node_mc, rec_mc));
    mc.push_back(&entries.back());

    ProcessorEntry rec_ext;
    rec_ext.source = Type::kExt;
    entries.emplace_back(make_entry(node_ext, rec_ext));
    mc.push_back(&entries.back());

    if (include_data) {
        ProcessorEntry rec_data;
        rec_data.source = Type::kData;
        entries.emplace_back(make_entry(node_data, rec_data));
        data.push_back(&entries.back());
    }

    stage_log("drawing stacked plot");
    draw_stack_plots(plotter, mc, data, include_data, adaptive_edges);
    stage_log("done");
    return 0;
}
