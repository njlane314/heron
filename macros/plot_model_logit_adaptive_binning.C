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
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
struct RocCurve {
    std::vector<double> fpr;
    std::vector<double> tpr;
    double auc = 0.0;
};

struct MetricScan {
    std::vector<double> x;
    std::vector<double> eff;
    std::vector<double> pur;
    std::vector<double> effpur;
    double best_thr = 0.0;
    double best_effpur = -1.0;
};

struct FineBinStats {
    double lo = 0.0;
    double hi = 0.0;
    double s_w = 0.0;
    double s_w2 = 0.0;
    double b_w = 0.0;
    double b_w2 = 0.0;
};

bool implicit_mt_enabled() {
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

RocCurve make_roc_curve(const std::vector<float> &scores,
                        const std::vector<int> &labels) {
    RocCurve out;
    if (scores.size() != labels.size() || scores.empty())
        return out;

    size_t n_pos = 0;
    size_t n_neg = 0;
    std::vector<std::pair<float, int>> ranked;
    ranked.reserve(scores.size());

    for (size_t i = 0; i < scores.size(); ++i) {
        const int y = labels[i] ? 1 : 0;
        ranked.emplace_back(scores[i], y);
        if (y)
            ++n_pos;
        else
            ++n_neg;
    }

    if (n_pos == 0 || n_neg == 0)
        return out;

    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<float, int> &a,
                 const std::pair<float, int> &b) { return a.first > b.first; });

    double tp = 0.0;
    double fp = 0.0;
    out.fpr.push_back(0.0);
    out.tpr.push_back(0.0);

    size_t i = 0;
    while (i < ranked.size()) {
        const float thr = ranked[i].first;
        while (i < ranked.size() && ranked[i].first == thr) {
            if (ranked[i].second)
                tp += 1.0;
            else
                fp += 1.0;
            ++i;
        }
        out.tpr.push_back(tp / static_cast<double>(n_pos));
        out.fpr.push_back(fp / static_cast<double>(n_neg));
    }

    for (size_t k = 1; k < out.fpr.size(); ++k) {
        const double dx = out.fpr[k] - out.fpr[k - 1];
        const double yavg = 0.5 * (out.tpr[k] + out.tpr[k - 1]);
        out.auc += dx * yavg;
    }

    return out;
}

MetricScan scan_thresholds(ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext,
                           const std::string &signal_sel, int n_thresholds,
                           const std::string &score_expr, double thr_min,
                           double thr_max) {
    MetricScan out;
    if (n_thresholds < 2)
        n_thresholds = 2;
    if (thr_max < thr_min)
        std::swap(thr_min, thr_max);

    const double signal_total =
        *(node_mc.Filter(signal_sel).Sum<double>("__w__"));
    if (signal_total <= 0.0)
        return out;

    out.x.reserve(static_cast<size_t>(n_thresholds));
    out.eff.reserve(static_cast<size_t>(n_thresholds));
    out.pur.reserve(static_cast<size_t>(n_thresholds));
    out.effpur.reserve(static_cast<size_t>(n_thresholds));

    for (int i = 0; i < n_thresholds; ++i) {
        const double frac =
            static_cast<double>(i) / static_cast<double>(n_thresholds - 1);
        const double thr = thr_min + frac * (thr_max - thr_min);
        const std::string pass_expr =
            "(" + score_expr + " >= " + std::to_string(thr) + ")";

        const double signal_pass =
            *(node_mc.Filter("(" + signal_sel + ") && " + pass_expr)
                  .Sum<double>("__w__"));
        const double selected_mc =
            *(node_mc.Filter(pass_expr).Sum<double>("__w__"));
        const double selected_ext =
            *(node_ext.Filter(pass_expr).Sum<double>("__w__"));

        const double selected_all = selected_mc + selected_ext;
        const double efficiency = signal_pass / signal_total;
        const double purity =
            (selected_all > 0.0) ? signal_pass / selected_all : 0.0;
        const double effpur = efficiency * purity;

        out.x.push_back(thr);
        out.eff.push_back(efficiency);
        out.pur.push_back(purity);
        out.effpur.push_back(effpur);

        if (effpur > out.best_effpur) {
            out.best_effpur = effpur;
            out.best_thr = thr;
        }
    }

    return out;
}

void draw_metric_plot(const MetricScan &scan, const std::string &canvas_id,
                      const std::string &title, const std::string &x_title,
                      const std::string &output_stem) {
    if (scan.x.empty())
        return;

    gStyle->SetOptStat(0);
    TCanvas c(canvas_id.c_str(), title.c_str(), 900, 700);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.07);
    c.SetBottomMargin(0.12);

    TGraph g_eff(static_cast<int>(scan.x.size()), scan.x.data(),
                 scan.eff.data());
    TGraph g_pur(static_cast<int>(scan.x.size()), scan.x.data(),
                 scan.pur.data());
    TGraph g_effpur(static_cast<int>(scan.x.size()), scan.x.data(),
                    scan.effpur.data());

    g_eff.SetTitle((";" + x_title + ";metric value").c_str());
    g_eff.SetLineColor(kBlue + 1);
    g_eff.SetMarkerColor(kBlue + 1);
    g_eff.SetLineWidth(3);
    g_eff.SetMarkerStyle(20);

    g_pur.SetLineColor(kRed + 1);
    g_pur.SetMarkerColor(kRed + 1);
    g_pur.SetLineWidth(3);
    g_pur.SetMarkerStyle(21);

    g_effpur.SetLineColor(kGreen + 2);
    g_effpur.SetMarkerColor(kGreen + 2);
    g_effpur.SetLineWidth(3);
    g_effpur.SetMarkerStyle(22);

    g_eff.Draw("ALP");
    g_pur.Draw("LP SAME");
    g_effpur.Draw("LP SAME");
    g_eff.GetYaxis()->SetRangeUser(0.0, 1.05);

    TLegend leg(0.17, 0.70, 0.56, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&g_eff, "efficiency", "lp");
    leg.AddEntry(&g_pur, "purity", "lp");
    leg.AddEntry(&g_effpur, "efficiency #times purity", "lp");
    leg.Draw();

    c.RedrawAxis();
    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());
    std::cout << "[plot_model_logit] saved plot: " << out << "\n";
}

void draw_stack_plots(Plotter &plotter, std::vector<const Entry *> &mc,
                      std::vector<const Entry *> &data, bool include_data) {
    auto &opt = plotter.options();

    TH1DModel spec = make_spec("inf_score_0", 50, -15.0, 15.0, "w_nominal");
    spec.sel = Preset::Empty;
    opt.x_title = "Inference score [0]";
    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    TH1DModel spec_sigmoid =
        make_spec("inf_score_0_sigmoid", 50, 0.0, 1.0, "w_nominal");
    spec_sigmoid.sel = Preset::Empty;
    opt.x_title = "Sigmoid(inference score [0])";
    if (include_data)
        plotter.draw_stack(spec_sigmoid, mc, data);
    else
        plotter.draw_stack(spec_sigmoid, mc);
}

double effective_count(double sum_w, double sum_w2) {
    if (sum_w2 <= 0.0)
        return 0.0;
    return (sum_w * sum_w) / sum_w2;
}

std::vector<double> make_adaptive_sigmoid_bins(
    ROOT::RDF::RNode node_mc, const std::string &signal_sel, int n_fine_bins,
    double nmin_signal, double nmin_background, int max_bins) {
    std::vector<double> edges;
    if (n_fine_bins < 10)
        n_fine_bins = 10;

    std::vector<FineBinStats> fine;
    fine.reserve(static_cast<size_t>(n_fine_bins));

    const double width = 1.0 / static_cast<double>(n_fine_bins);
    for (int i = 0; i < n_fine_bins; ++i) {
        const double lo = i * width;
        const double hi = (i + 1) * width;

        const std::string in_bin =
            "(inf_score_0_sigmoid >= " + std::to_string(lo) +
            ") && (inf_score_0_sigmoid < " + std::to_string(hi) + ")";
        const std::string sig_bin = "(" + signal_sel + ") && " + in_bin;
        const std::string bkg_bin = "(!(" + signal_sel + ")) && " + in_bin;

        FineBinStats stat;
        stat.lo = lo;
        stat.hi = hi;
        stat.s_w = *(node_mc.Filter(sig_bin).Sum<double>("__w__"));
        stat.s_w2 = *(node_mc.Filter(sig_bin).Sum<double>("__w2__"));
        stat.b_w = *(node_mc.Filter(bkg_bin).Sum<double>("__w__"));
        stat.b_w2 = *(node_mc.Filter(bkg_bin).Sum<double>("__w2__"));
        fine.push_back(stat);
    }

    edges.push_back(0.0);

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
        const bool pass = (neff_s >= nmin_signal && neff_b >= nmin_background);

        if (pass || i == 0) {
            edges.push_back(fine[static_cast<size_t>(i)].hi);
            acc_sw = 0.0;
            acc_sw2 = 0.0;
            acc_bw = 0.0;
            acc_bw2 = 0.0;
            if (static_cast<int>(edges.size()) - 1 >= max_bins && i > 0)
                break;
        }
    }

    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    if (edges.empty() || edges.front() > 0.0)
        edges.insert(edges.begin(), 0.0);
    if (edges.back() < 1.0)
        edges.push_back(1.0);

    return edges;
}

void draw_roc_plot(ROOT::RDF::RNode auc_node) {
    auto v_logit = auc_node.Take<float>("inf_score_0");
    auto v_sigmoid = auc_node.Take<float>("inf_score_0_sigmoid");
    auto v_label = auc_node.Take<int>("is_signal_label");

    const RocCurve roc_logit = make_roc_curve(*v_logit, *v_label);
    const RocCurve roc_sigmoid = make_roc_curve(*v_sigmoid, *v_label);

    if (roc_logit.fpr.empty() || roc_sigmoid.fpr.empty()) {
        std::cerr << "[plot_model_logit] unable to compute ROC/AUC (missing "
                     "signal/background entries).\n";
        return;
    }

    gStyle->SetOptStat(0);
    TCanvas c("c_infscore_auc", "c_infscore_auc", 800, 700);

    TGraph g_logit(static_cast<int>(roc_logit.fpr.size()), roc_logit.fpr.data(),
                   roc_logit.tpr.data());
    g_logit.SetTitle(";False positive rate;True positive rate");
    g_logit.SetLineColor(kBlue + 1);
    g_logit.SetLineWidth(3);
    g_logit.Draw("AL");

    TGraph g_sigmoid(static_cast<int>(roc_sigmoid.fpr.size()),
                     roc_sigmoid.fpr.data(), roc_sigmoid.tpr.data());
    g_sigmoid.SetLineColor(kRed + 1);
    g_sigmoid.SetLineStyle(2);
    g_sigmoid.SetLineWidth(3);
    g_sigmoid.Draw("L SAME");

    TLine diagonal(0.0, 0.0, 1.0, 1.0);
    diagonal.SetLineStyle(3);
    diagonal.SetLineColor(kGray + 2);
    diagonal.Draw("SAME");

    TLegend leg(0.16, 0.70, 0.60, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    char label_logit[64];
    char label_sigmoid[64];
    std::snprintf(label_logit, sizeof(label_logit), "Raw logit (AUC = %.4f)",
                  roc_logit.auc);
    std::snprintf(label_sigmoid, sizeof(label_sigmoid),
                  "Sigmoid(logit) (AUC = %.4f)", roc_sigmoid.auc);
    leg.AddEntry(&g_logit, label_logit, "l");
    leg.AddEntry(&g_sigmoid, label_sigmoid, "l");
    leg.Draw();

    const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
    const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
    gSystem->mkdir(out_dir.c_str(), true);
    const std::string out_path =
        out_dir + "/first_inference_score_roc_auc." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plot_model_logit] wrote " << out_path << "\n";
}

void draw_effpur_plots(ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext,
                       const std::string &signal_sel, int n_thresholds,
                       double raw_threshold_min, double raw_threshold_max,
                       const std::string &output_stem) {
    const MetricScan scan_sigmoid =
        scan_thresholds(node_mc, node_ext, signal_sel, n_thresholds,
                        "inf_score_0_sigmoid", 0.0, 1.0);
    if (scan_sigmoid.x.empty()) {
        std::cerr
            << "[plot_model_logit] signal denominator is <= 0 for signal_sel='"
            << signal_sel << "'.\n";
        return;
    }

    std::cout << "[plot_model_logit] best sigmoid threshold="
              << scan_sigmoid.best_thr
              << " with efficiency x purity=" << scan_sigmoid.best_effpur
              << "\n";
    draw_metric_plot(scan_sigmoid, "c_first_inf_score_effpur",
                     "First inference-score threshold scan",
                     "Sigmoid(inference score [0]) threshold", output_stem);

    const MetricScan scan_raw =
        scan_thresholds(node_mc, node_ext, signal_sel, n_thresholds,
                        "inf_score_0", raw_threshold_min, raw_threshold_max);
    std::cout << "[plot_model_logit] best raw threshold=" << scan_raw.best_thr
              << " with efficiency x purity=" << scan_raw.best_effpur << "\n";
    draw_metric_plot(scan_raw, "c_first_inf_score_effpur_raw",
                     "First inference-score raw-threshold scan",
                     "Inference score [0] threshold", output_stem + "_raw");
}
} // namespace

int plot_model_logit_adaptive_binning(
    const std::string &samples_tsv = "",
    const std::string &extra_sel = "sel_muon", bool use_logy = true,
    bool include_data = false, const std::string &signal_sel = "is_signal",
    const std::string &mc_weight = "w_nominal", int n_thresholds = 101,
    double raw_threshold_min = -15.0, double raw_threshold_max = 15.0,
    int n_fine_bins = 200, double nmin_signal = 10.0,
    double nmin_background = 30.0, int max_bins = 15,
    const std::string &output_stem = "first_inference_score_eff_purity") {
    const std::string list_path =
        samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plot_model_logit_adaptive_binning] input=" << list_path
              << "\n";

    if (!looks_like_event_list_root(list_path)) {
        std::cerr
            << "[plot_model_logit_adaptive_binning] input is not an event "
               "list ROOT file: "
            << list_path << "\n";
        return 2;
    }

    if (implicit_mt_enabled())
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
                   {"inf_scores"})
            .Define("inf_score_0_sigmoid",
                    [](float x) { return 1.0f / (1.0f + std::exp(-x)); },
                    {"inf_score_0"});

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

    ROOT::RDF::RNode auc_node = filter_by_mask(base, mask_mc);

    if (!extra_sel.empty()) {
        const bool named_column = rdf.HasColumn(extra_sel);

        if (named_column) {
            e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(
                [](bool pass) { return pass; }, {extra_sel});
            e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(
                [](bool pass) { return pass; }, {extra_sel});
            auc_node =
                auc_node.Filter([](bool pass) { return pass; }, {extra_sel});
            if (p_data != nullptr)
                p_data->selection.nominal.node =
                    p_data->selection.nominal.node.Filter(
                        [](bool pass) { return pass; }, {extra_sel});
        } else {
            e_mc.selection.nominal.node =
                e_mc.selection.nominal.node.Filter(extra_sel);
            e_ext.selection.nominal.node =
                e_ext.selection.nominal.node.Filter(extra_sel);
            auc_node = auc_node.Filter(extra_sel);
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

    draw_stack_plots(plotter, mc, data, include_data);
    const std::vector<double> adaptive_edges = make_adaptive_sigmoid_bins(
        e_mc.selection.nominal.node, signal_sel, n_fine_bins, nmin_signal,
        nmin_background, max_bins);
    std::cout << "[plot_model_logit_adaptive_binning] adaptive edges:";
    for (double e : adaptive_edges)
        std::cout << " " << e;
    std::cout << "\n";
    draw_roc_plot(auc_node);
    draw_effpur_plots(e_mc.selection.nominal.node, e_ext.selection.nominal.node,
                      signal_sel, n_thresholds, raw_threshold_min,
                      raw_threshold_max, output_stem);

    std::cout << "[plot_model_logit_adaptive_binning] done\n";
    return 0;
}
