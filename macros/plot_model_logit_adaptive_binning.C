#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
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
struct FineBinStats {
    double lo = 0.0;
    double hi = 0.0;
    double s_w = 0.0;
    double s_w2 = 0.0;
    double b_w = 0.0;
    double b_w2 = 0.0;
};

struct ScoreRange {
    double lo = 0.0;
    double hi = 0.0;
    bool valid = false;
};

ScoreRange find_score_range(ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext) {
    ScoreRange out;

    auto min_mc = node_mc.Min<float>("inf_score_0");
    auto max_mc = node_mc.Max<float>("inf_score_0");
    auto min_ext = node_ext.Min<float>("inf_score_0");
    auto max_ext = node_ext.Max<float>("inf_score_0");

    const bool has_mc = (node_mc.Count().GetValue() > 0);
    const bool has_ext = (node_ext.Count().GetValue() > 0);
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

    return out;
}

void draw_stack_plots(Plotter &plotter, std::vector<const Entry *> &mc,
                      std::vector<const Entry *> &data, bool include_data,
                      const std::vector<double> &adaptive_edges) {
    auto &opt = plotter.options();

    TH1DModel spec = make_spec("inf_score_0", adaptive_edges, "w_nominal");
    spec.sel = Preset::Empty;
    opt.x_title = "Inference score [0]";
    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);
}

double effective_count(double sum_w, double sum_w2) {
    if (sum_w2 <= 0.0)
        return 0.0;
    return (sum_w * sum_w) / sum_w2;
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

    std::vector<FineBinStats> fine;
    fine.reserve(static_cast<size_t>(n_fine_bins));

    ROOT::RDF::TH1DModel model_sig_w("h_sig_w", "", n_fine_bins, score_lo,
                                     score_hi);
    ROOT::RDF::TH1DModel model_sig_w2("h_sig_w2", "", n_fine_bins,
                                      score_lo, score_hi);
    ROOT::RDF::TH1DModel model_bkg_mc_w("h_bkg_mc_w", "", n_fine_bins,
                                        score_lo, score_hi);
    ROOT::RDF::TH1DModel model_bkg_mc_w2("h_bkg_mc_w2", "", n_fine_bins,
                                         score_lo, score_hi);
    ROOT::RDF::TH1DModel model_ext_w("h_ext_w", "", n_fine_bins, score_lo,
                                     score_hi);
    ROOT::RDF::TH1DModel model_ext_w2("h_ext_w2", "", n_fine_bins,
                                      score_lo, score_hi);

    auto node_mc_tagged = node_mc.Define("__is_sig_bin__", signal_sel);

    auto h_sig_w = node_mc_tagged.Filter("__is_sig_bin__")
                       .Histo1D(model_sig_w, "inf_score_0", "__w__");
    auto h_sig_w2 = node_mc_tagged.Filter("__is_sig_bin__")
                        .Histo1D(model_sig_w2, "inf_score_0", "__w2__");
    auto h_bkg_mc_w =
        node_mc_tagged.Filter("!__is_sig_bin__")
            .Histo1D(model_bkg_mc_w, "inf_score_0", "__w__");
    auto h_bkg_mc_w2 =
        node_mc_tagged.Filter("!__is_sig_bin__")
            .Histo1D(model_bkg_mc_w2, "inf_score_0", "__w2__");
    auto h_ext_w = node_ext.Histo1D(model_ext_w, "inf_score_0", "__w__");
    auto h_ext_w2 = node_ext.Histo1D(model_ext_w2, "inf_score_0", "__w2__");

    const double width =
        (score_hi - score_lo) / static_cast<double>(n_fine_bins);
    for (int i = 0; i < n_fine_bins; ++i) {
        const double lo = score_lo + i * width;
        const double hi = score_lo + (i + 1) * width;
        const int bin = i + 1;

        FineBinStats stat;
        stat.lo = lo;
        stat.hi = hi;
        stat.s_w = h_sig_w->GetBinContent(bin);
        stat.s_w2 = h_sig_w2->GetBinContent(bin);
        stat.b_w =
            h_bkg_mc_w->GetBinContent(bin) + h_ext_w->GetBinContent(bin);
        stat.b_w2 =
            h_bkg_mc_w2->GetBinContent(bin) + h_ext_w2->GetBinContent(bin);
        fine.push_back(stat);
    }

    edges.push_back(score_lo);

    double acc_sw = 0.0;
    double acc_sw2 = 0.0;
    double acc_bw = 0.0;
    double acc_bw2 = 0.0;

    for (int i = n_fine_bins - 1; i >= 0; --i) {
        acc_sw += fine[static_cast<size_t>(i)].s_w;
        acc_sw2 += fine[static_cast<size_t>(i)].s_w2;
        acc_bw += fine[static_cast<size_t>(i)].b_w;
        acc_bw2 += fine[static_cast<size_t>(i)].b_w2;

        const double neff_s = effective_count(acc_sw, acc_sw2);
        const double neff_b = effective_count(acc_bw, acc_bw2);
        (void)neff_s;
        (void)neff_b;

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
        if (i == 0 || (pass && can_add_split)) {
            edges.push_back(fine[static_cast<size_t>(i)].lo);
            acc_sw = 0.0;
            acc_sw2 = 0.0;
            acc_bw = 0.0;
            acc_bw2 = 0.0;
        }
    }

    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    if (edges.empty() || edges.front() > score_lo)
        edges.insert(edges.begin(), score_lo);
    if (edges.back() < score_hi)
        edges.push_back(score_hi);

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
    const std::string list_path =
        "/exp/uboone/data/users/nlane/heron/out/event/events_muon.root";
    std::cout << "[plot_model_logit_adaptive_binning] input=" << list_path
              << "\n";

    if (!looks_like_event_list_root(list_path)) {
        std::cerr
            << "[plot_model_logit_adaptive_binning] input is not an event "
               "list ROOT file: "
            << list_path << "\n";
        return 2;
    }

    ROOT::EnableImplicitMT();

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    auto filter_by_mask = [](ROOT::RDF::RNode node,
                             std::shared_ptr<const std::vector<char>> mask) {
        return node.Filter(
            [mask](int sid) {
                return sid >= 0 && sid < static_cast<int>(mask->size()) &&
                       (*mask)[static_cast<std::size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base =
        rdf.Define("inf_score_0",
                   [](const ROOT::RVec<float> &scores) {
                       if (scores.empty())
                           return -1.0f;
                       return scores[0];
                   },
                   {"inf_scores"});

    if (rdf.HasColumn("analysis_channels")) {
        base = base.Define(
            "is_signal_label",
            [](int analysis_channel) {
                return (analysis_channel == 15 || analysis_channel == 16) ? 1
                                                                          : 0;
            },
            {"analysis_channels"});
    } else {
        base = base.Define("is_signal_label",
                           [](bool is_signal) { return is_signal ? 1 : 0; },
                           {"is_signal"});
    }

    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext)
                                    .Define("__w__", mc_weight)
                                    .Define("__w2__", "__w__*__w__");
    ROOT::RDF::RNode node_mc =
        filter_by_mask(base, mask_mc)
            .Filter(
                [mask_ext](int sid) {
                    return !(sid >= 0 &&
                             sid < static_cast<int>(mask_ext->size()) &&
                             (*mask_ext)[static_cast<std::size_t>(sid)]);
                },
                {"sample_id"})
            .Define("__w__", mc_weight)
            .Define("__w2__", "__w__*__w__");
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

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
    if (include_data) {
        ProcessorEntry rec_data;
        rec_data.source = Type::kData;
        entries.emplace_back(make_entry(std::move(node_data), rec_data));
        p_data = &entries.back();
        data.push_back(p_data);
    }

    if (!extra_sel.empty()) {
        const bool named_column = rdf.HasColumn(extra_sel);

        if (named_column) {
            e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(
                [](bool pass) { return pass; }, {extra_sel});
            e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(
                [](bool pass) { return pass; }, {extra_sel});
            if (p_data != nullptr)
                p_data->selection.nominal.node =
                    p_data->selection.nominal.node.Filter(
                        [](bool pass) { return pass; }, {extra_sel});
        } else {
            e_mc.selection.nominal.node =
                e_mc.selection.nominal.node.Filter(extra_sel);
            e_ext.selection.nominal.node =
                e_ext.selection.nominal.node.Filter(extra_sel);
            if (p_data != nullptr)
                p_data->selection.nominal.node =
                    p_data->selection.nominal.node.Filter(extra_sel);
        }
    }

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

    ScoreRange score_range =
        find_score_range(e_mc.selection.nominal.node, e_ext.selection.nominal.node);
    if (!score_range.valid) {
        score_range.lo = raw_threshold_min;
        score_range.hi = raw_threshold_max;
    }

    std::cout << "[plot_model_logit_adaptive_binning] score range: ["
              << score_range.lo << ", " << score_range.hi << "]\n";

    const std::string adaptive_signal_sel =
        (signal_sel == "is_signal" && !rdf.HasColumn("is_signal"))
            ? "is_signal_label"
            : signal_sel;
    const std::vector<double> adaptive_edges = make_adaptive_score_bins(
        e_mc.selection.nominal.node, e_ext.selection.nominal.node,
        adaptive_signal_sel, n_fine_bins, nmin_signal, nmin_background,
        max_bins, score_range.lo, score_range.hi);
    draw_stack_plots(plotter, mc, data, include_data, adaptive_edges);
    std::cout << "[plot_model_logit_adaptive_binning] adaptive edges:";
    for (double e : adaptive_edges)
        std::cout << " " << e;
    std::cout << "\n";
    std::cout << "[plot_model_logit_adaptive_binning] done\n";
    return 0;
}
