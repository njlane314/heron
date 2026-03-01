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
#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
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
#include "PlotEnv.hh"
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

    const int n_bins = n_thresholds - 1;
    ROOT::RDF::TH1DModel model_signal("h_signal_scan", "", n_bins, thr_min,
                                      thr_max);
    ROOT::RDF::TH1DModel model_mc("h_mc_scan", "", n_bins, thr_min, thr_max);
    ROOT::RDF::TH1DModel model_ext("h_ext_scan", "", n_bins, thr_min,
                                   thr_max);

    auto h_signal =
        node_mc.Filter(signal_sel).Histo1D(model_signal, score_expr, "__w__");
    auto h_mc = node_mc.Histo1D(model_mc, score_expr, "__w__");
    auto h_ext = node_ext.Histo1D(model_ext, score_expr, "__w__");

    const TH1D &h_sig_ref = *h_signal;
    const TH1D &h_mc_ref = *h_mc;
    const TH1D &h_ext_ref = *h_ext;

    const double signal_total = h_sig_ref.Integral(0, n_bins + 1);
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
        const int bin = h_sig_ref.FindFixBin(thr);
        const int first_bin = std::max(0, std::min(bin, n_bins + 1));

        const double signal_pass = h_sig_ref.Integral(first_bin, n_bins + 1);
        const double selected_mc = h_mc_ref.Integral(first_bin, n_bins + 1);
        const double selected_ext = h_ext_ref.Integral(first_bin, n_bins + 1);

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

void draw_adaptive_score_plot(ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext,
                              ROOT::RDF::RNode node_data, bool include_data,
                              const std::vector<double> &edges) {
    if (edges.size() < 2)
        return;

    const int nbins = static_cast<int>(edges.size()) - 1;
    ROOT::RDF::TH1DModel model_mc("h_score_mc", ";Inference score [0];Events", nbins,
                                  edges.data());
    ROOT::RDF::TH1DModel model_ext("h_score_ext", ";Inference score [0];Events", nbins,
                                   edges.data());
    ROOT::RDF::TH1DModel model_data("h_score_data", ";Inference score [0];Events", nbins,
                                    edges.data());

    auto h_mc = node_mc.Histo1D(model_mc, "inf_score_0", "__w__");
    auto h_ext = node_ext.Histo1D(model_ext, "inf_score_0", "__w__");

    std::unique_ptr<ROOT::RDF::RResultPtr<TH1D>> h_data;
    if (include_data)
        h_data = std::make_unique<ROOT::RDF::RResultPtr<TH1D>>(
            node_data.Histo1D(model_data, "inf_score_0"));

    gStyle->SetOptStat(0);
    TCanvas c("c_infscore_adaptive", "c_infscore_adaptive", 900, 700);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.07);
    c.SetBottomMargin(0.12);

    TH1D *p_mc = h_mc.GetPtr();
    TH1D *p_ext = h_ext.GetPtr();
    p_mc->SetFillColor(kAzure + 1);
    p_mc->SetLineColor(kAzure + 3);
    p_ext->SetFillColor(kOrange + 7);
    p_ext->SetLineColor(kOrange + 3);

    THStack stack("hs_infscore_adaptive", ";Inference score [0];Events");
    stack.Add(p_ext);
    stack.Add(p_mc);
    stack.Draw("HIST");

    TH1D *p_data = nullptr;
    if (include_data) {
        p_data = (*h_data).GetPtr();
        p_data->SetLineColor(kBlack);
        p_data->SetMarkerColor(kBlack);
        p_data->SetMarkerStyle(20);
        p_data->Draw("E1 SAME");
    }

    TLegend leg(0.62, 0.68, 0.89, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(p_mc, "MC", "f");
    leg.AddEntry(p_ext, "EXT", "f");
    if (p_data != nullptr)
        leg.AddEntry(p_data, "Data", "ep");
    leg.Draw();

    c.RedrawAxis();
    const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
    const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
    gSystem->mkdir(out_dir.c_str(), true);
    const std::string out_path =
        out_dir + "/first_inference_score_adaptive_binning." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plot_model_logit] wrote " << out_path << "\n";
}

double effective_count(double sum_w, double sum_w2) {
    if (sum_w2 <= 0.0)
        return 0.0;
    return (sum_w * sum_w) / sum_w2;
}

std::vector<double> make_adaptive_score_bins(
    ROOT::RDF::RNode node_mc, const std::string &signal_sel, int n_fine_bins,
    double nmin_signal, double nmin_background, int max_bins,
    double score_lo, double score_hi) {
    std::vector<double> edges;
    if (n_fine_bins < 10)
        n_fine_bins = 10;

    std::vector<FineBinStats> fine;
    fine.reserve(static_cast<size_t>(n_fine_bins));

    ROOT::RDF::TH1DModel model_sig_w("h_sig_w", "", n_fine_bins, score_lo,
                                     score_hi);
    ROOT::RDF::TH1DModel model_sig_w2("h_sig_w2", "", n_fine_bins,
                                      score_lo, score_hi);
    ROOT::RDF::TH1DModel model_bkg_w("h_bkg_w", "", n_fine_bins, score_lo,
                                     score_hi);
    ROOT::RDF::TH1DModel model_bkg_w2("h_bkg_w2", "", n_fine_bins,
                                      score_lo, score_hi);

    auto h_sig_w = node_mc.Filter(signal_sel)
                       .Histo1D(model_sig_w, "inf_score_0", "__w__");
    auto h_sig_w2 = node_mc.Filter(signal_sel)
                        .Histo1D(model_sig_w2, "inf_score_0", "__w2__");
    auto h_bkg_w =
        node_mc.Filter("!(" + signal_sel + ")")
            .Histo1D(model_bkg_w, "inf_score_0", "__w__");
    auto h_bkg_w2 =
        node_mc.Filter("!(" + signal_sel + ")")
            .Histo1D(model_bkg_w2, "inf_score_0", "__w2__");

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
        stat.b_w = h_bkg_w->GetBinContent(bin);
        stat.b_w2 = h_bkg_w2->GetBinContent(bin);
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

    if (edges.empty() || edges.front() > score_lo)
        edges.insert(edges.begin(), score_lo);
    if (edges.back() < score_hi)
        edges.push_back(score_hi);

    return edges;
}

void draw_roc_plot(ROOT::RDF::RNode auc_node) {
    auto v_logit = auc_node.Take<float>("inf_score_0");
    auto v_label = auc_node.Take<int>("is_signal_label");

    const RocCurve roc_logit = make_roc_curve(*v_logit, *v_label);

    if (roc_logit.fpr.empty()) {
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

    TLine diagonal(0.0, 0.0, 1.0, 1.0);
    diagonal.SetLineStyle(3);
    diagonal.SetLineColor(kGray + 2);
    diagonal.Draw("SAME");

    TLegend leg(0.16, 0.70, 0.60, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    char label_logit[64];
    std::snprintf(label_logit, sizeof(label_logit), "Raw logit (AUC = %.4f)",
                  roc_logit.auc);
    leg.AddEntry(&g_logit, label_logit, "l");
    leg.Draw();

    const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
    const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
    gSystem->mkdir(out_dir.c_str(), true);
    const std::string out_path =
        out_dir + "/first_inference_score_roc_auc." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plot_model_logit] wrote " << out_path << "\n";
}

void draw_effpur_plot(ROOT::RDF::RNode node_mc, ROOT::RDF::RNode node_ext,
                      const std::string &signal_sel, int n_thresholds,
                      double raw_threshold_min, double raw_threshold_max,
                      const std::string &output_stem) {
    const MetricScan scan_raw =
        scan_thresholds(node_mc, node_ext, signal_sel, n_thresholds,
                        "inf_score_0", raw_threshold_min, raw_threshold_max);
    if (scan_raw.x.empty()) {
        std::cerr
            << "[plot_model_logit] signal denominator is <= 0 for signal_sel='"
            << signal_sel << "'.\n";
        return;
    }

    std::cout << "[plot_model_logit] best raw threshold=" << scan_raw.best_thr
              << " with efficiency x purity=" << scan_raw.best_effpur << "\n";
    draw_metric_plot(scan_raw, "c_first_inf_score_effpur_raw",
                     "First inference-score threshold scan",
                     "Inference score [0] threshold", output_stem);
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
            ;

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

    ROOT::RDF::RNode auc_node = filter_by_mask(base, mask_mc);

    if (!extra_sel.empty()) {
        const bool named_column = rdf.HasColumn(extra_sel);

        if (named_column) {
            node_mc = node_mc.Filter([](bool pass) { return pass; }, {extra_sel});
            node_ext = node_ext.Filter([](bool pass) { return pass; }, {extra_sel});
            auc_node =
                auc_node.Filter([](bool pass) { return pass; }, {extra_sel});
            if (include_data)
                node_data =
                    node_data.Filter([](bool pass) { return pass; }, {extra_sel});
        } else {
            node_mc = node_mc.Filter(extra_sel);
            node_ext = node_ext.Filter(extra_sel);
            auc_node = auc_node.Filter(extra_sel);
            if (include_data)
                node_data = node_data.Filter(extra_sel);
        }
    }

    (void)use_logy;

    ScoreRange score_range = find_score_range(node_mc, node_ext);
    if (!score_range.valid) {
        score_range.lo = raw_threshold_min;
        score_range.hi = raw_threshold_max;
    }

    std::cout << "[plot_model_logit_adaptive_binning] score range: ["
              << score_range.lo << ", " << score_range.hi << "]\n";

    const std::vector<double> adaptive_edges = make_adaptive_score_bins(
        node_mc, signal_sel, n_fine_bins, nmin_signal, nmin_background,
        max_bins, score_range.lo, score_range.hi);
    std::cout << "[plot_model_logit_adaptive_binning] adaptive edges:";
    for (double e : adaptive_edges)
        std::cout << " " << e;
    std::cout << "\n";
    draw_adaptive_score_plot(node_mc, node_ext, node_data, include_data,
                             adaptive_edges);
    draw_roc_plot(auc_node);
    draw_effpur_plot(node_mc, node_ext, signal_sel, n_thresholds,
                     score_range.lo, score_range.hi, output_stem);

    std::cout << "[plot_model_logit_adaptive_binning] done\n";
    return 0;
}
