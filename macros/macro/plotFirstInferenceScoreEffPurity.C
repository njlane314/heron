// plot/macro/plotFirstInferenceScoreEffPurity.C
//
// Scan a threshold on the first inference-score entry and plot:
//   - efficiency
//   - purity (MC + EXT in denominator)
//   - efficiency x purity
//
// Run with:
//   ./heron macro plotFirstInferenceScoreEffPurity.C
//   ./heron macro plotFirstInferenceScoreEffPurity.C \
//     'plotFirstInferenceScoreEffPurity("./scratch/out/event_list_myana.root")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "MacroGuard.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &path)
{
    const auto n = path.size();
    if (n < 5 || path.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> input_file(TFile::Open(path.c_str(), "READ"));
    if (!input_file || input_file->IsZombie())
        return false;

    const bool has_refs = (input_file->Get("sample_refs") != nullptr);
    const bool has_events_tree = (input_file->Get("events") != nullptr);
    const bool has_event_tree_key = (input_file->Get("event_tree") != nullptr);
    return has_refs && (has_events_tree || has_event_tree_key);
}

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}
}

int plotFirstInferenceScoreEffPurity(const std::string &event_list_path = "",
                                     const std::string &base_sel = "sel_muon",
                                     const std::string &signal_sel = "is_signal",
                                     const std::string &mc_weight = "w_nominal",
                                     int n_thresholds = 101,
                                     double raw_threshold_min = -15.0,
                                     double raw_threshold_max = 15.0,
                                     const std::string &output_stem = "first_inference_score_eff_purity")
{
  return heron::macro::run_with_guard("plotFirstInferenceScoreEffPurity", [&]() -> int {
    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
    std::cout << "[plotFirstInferenceScoreEffPurity] input=" << input_path << "\n";

    if (!looks_like_event_list_root(input_path))
    {
        std::cerr << "[plotFirstInferenceScoreEffPurity] input is not an event-list root file: "
                  << input_path << "\n";
        return 1;
    }

    if (n_thresholds < 2)
        n_thresholds = 2;

    if (implicit_mt_enabled())
        ROOT::EnableImplicitMT();

    EventListIO el(input_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();

    auto filter_by_mask = [](ROOT::RDF::RNode node, std::shared_ptr<const std::vector<char>> mask) {
        return node.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<std::size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf.Define("inf_score_0",
                                       [](const ROOT::RVec<float> &scores) {
                                           if (scores.empty())
                                               return -1.0f;
                                           return scores[0];
                                       },
                                       {"inf_scores"})
                               .Define("inf_score_0_sigmoid",
                                       [](float x) { return 1.0f / (1.0f + std::exp(-x)); },
                                       {"inf_score_0"});

    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<std::size_t>(sid)]);
                                   },
                                   {"sample_id"})
                                   .Define("__w__", mc_weight);

    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext)
                                    .Define("__w__", mc_weight);

    if (!base_sel.empty())
    {
        if (rdf.HasColumn(base_sel))
        {
            node_mc = node_mc.Filter([](bool pass) { return pass; }, {base_sel});
            node_ext = node_ext.Filter([](bool pass) { return pass; }, {base_sel});
        }
        else
        {
            node_mc = node_mc.Filter(base_sel);
            node_ext = node_ext.Filter(base_sel);
        }
    }

    const double signal_total = *(node_mc.Filter(signal_sel).Sum<double>("__w__"));
    if (signal_total <= 0.0)
    {
        std::cerr << "[plotFirstInferenceScoreEffPurity] signal denominator is <= 0 for signal_sel='"
                  << signal_sel << "'.\n";
        return 1;
    }

    std::vector<double> x;
    std::vector<double> y_eff;
    std::vector<double> y_pur;
    std::vector<double> y_effpur;
    x.reserve(static_cast<std::size_t>(n_thresholds));
    y_eff.reserve(static_cast<std::size_t>(n_thresholds));
    y_pur.reserve(static_cast<std::size_t>(n_thresholds));
    y_effpur.reserve(static_cast<std::size_t>(n_thresholds));

    double best_thr = 0.0;
    double best_effpur = -1.0;

    for (int i = 0; i < n_thresholds; ++i)
    {
        const double thr = static_cast<double>(i) / static_cast<double>(n_thresholds - 1);
        const std::string pass_expr = "(inf_score_0_sigmoid >= " + std::to_string(thr) + ")";

        const double signal_pass = *(node_mc.Filter("(" + signal_sel + ") && " + pass_expr).Sum<double>("__w__"));
        const double selected_mc = *(node_mc.Filter(pass_expr).Sum<double>("__w__"));
        const double selected_ext = *(node_ext.Filter(pass_expr).Sum<double>("__w__"));

        const double selected_all = selected_mc + selected_ext;
        const double efficiency = signal_pass / signal_total;
        const double purity = (selected_all > 0.0) ? signal_pass / selected_all : 0.0;
        const double effpur = efficiency * purity;

        x.push_back(thr);
        y_eff.push_back(efficiency);
        y_pur.push_back(purity);
        y_effpur.push_back(effpur);

        if (effpur > best_effpur)
        {
            best_effpur = effpur;
            best_thr = thr;
        }
    }

    std::cout << "[plotFirstInferenceScoreEffPurity] best threshold=" << best_thr
              << " with efficiency x purity=" << best_effpur << "\n";

    if (raw_threshold_max < raw_threshold_min)
        std::swap(raw_threshold_min, raw_threshold_max);

    std::vector<double> x_raw;
    std::vector<double> y_raw_eff;
    std::vector<double> y_raw_pur;
    std::vector<double> y_raw_effpur;
    x_raw.reserve(static_cast<std::size_t>(n_thresholds));
    y_raw_eff.reserve(static_cast<std::size_t>(n_thresholds));
    y_raw_pur.reserve(static_cast<std::size_t>(n_thresholds));
    y_raw_effpur.reserve(static_cast<std::size_t>(n_thresholds));

    double best_raw_thr = raw_threshold_min;
    double best_raw_effpur = -1.0;

    for (int i = 0; i < n_thresholds; ++i)
    {
        const double frac = static_cast<double>(i) / static_cast<double>(n_thresholds - 1);
        const double thr = raw_threshold_min + frac * (raw_threshold_max - raw_threshold_min);
        const std::string pass_expr = "(inf_score_0 >= " + std::to_string(thr) + ")";

        const double signal_pass = *(node_mc.Filter("(" + signal_sel + ") && " + pass_expr).Sum<double>("__w__"));
        const double selected_mc = *(node_mc.Filter(pass_expr).Sum<double>("__w__"));
        const double selected_ext = *(node_ext.Filter(pass_expr).Sum<double>("__w__"));

        const double selected_all = selected_mc + selected_ext;
        const double efficiency = signal_pass / signal_total;
        const double purity = (selected_all > 0.0) ? signal_pass / selected_all : 0.0;
        const double effpur = efficiency * purity;

        x_raw.push_back(thr);
        y_raw_eff.push_back(efficiency);
        y_raw_pur.push_back(purity);
        y_raw_effpur.push_back(effpur);

        if (effpur > best_raw_effpur)
        {
            best_raw_effpur = effpur;
            best_raw_thr = thr;
        }
    }

    std::cout << "[plotFirstInferenceScoreEffPurity] best raw threshold=" << best_raw_thr
              << " with efficiency x purity=" << best_raw_effpur << "\n";

    gStyle->SetOptStat(0);
    TCanvas c("c_first_inf_score_effpur", "First inference-score threshold scan", 900, 700);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.07);
    c.SetBottomMargin(0.12);

    TGraph g_eff(static_cast<int>(x.size()), x.data(), y_eff.data());
    TGraph g_pur(static_cast<int>(x.size()), x.data(), y_pur.data());
    TGraph g_effpur(static_cast<int>(x.size()), x.data(), y_effpur.data());

    g_eff.SetTitle(";Sigmoid(inference score [0]) threshold;metric value");
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

    TCanvas c_raw("c_first_inf_score_effpur_raw", "First inference-score raw-threshold scan", 900, 700);
    c_raw.SetLeftMargin(0.11);
    c_raw.SetRightMargin(0.07);
    c_raw.SetBottomMargin(0.12);

    TGraph g_raw_eff(static_cast<int>(x_raw.size()), x_raw.data(), y_raw_eff.data());
    TGraph g_raw_pur(static_cast<int>(x_raw.size()), x_raw.data(), y_raw_pur.data());
    TGraph g_raw_effpur(static_cast<int>(x_raw.size()), x_raw.data(), y_raw_effpur.data());

    g_raw_eff.SetTitle(";Inference score [0] threshold;metric value");
    g_raw_eff.SetLineColor(kBlue + 1);
    g_raw_eff.SetMarkerColor(kBlue + 1);
    g_raw_eff.SetLineWidth(3);
    g_raw_eff.SetMarkerStyle(20);

    g_raw_pur.SetLineColor(kRed + 1);
    g_raw_pur.SetMarkerColor(kRed + 1);
    g_raw_pur.SetLineWidth(3);
    g_raw_pur.SetMarkerStyle(21);

    g_raw_effpur.SetLineColor(kGreen + 2);
    g_raw_effpur.SetMarkerColor(kGreen + 2);
    g_raw_effpur.SetLineWidth(3);
    g_raw_effpur.SetMarkerStyle(22);

    g_raw_eff.Draw("ALP");
    g_raw_pur.Draw("LP SAME");
    g_raw_effpur.Draw("LP SAME");

    g_raw_eff.GetYaxis()->SetRangeUser(0.0, 1.05);

    TLegend leg_raw(0.17, 0.70, 0.56, 0.88);
    leg_raw.SetBorderSize(0);
    leg_raw.SetFillStyle(0);
    leg_raw.AddEntry(&g_raw_eff, "efficiency", "lp");
    leg_raw.AddEntry(&g_raw_pur, "purity", "lp");
    leg_raw.AddEntry(&g_raw_effpur, "efficiency #times purity", "lp");
    leg_raw.Draw();

    c_raw.RedrawAxis();

    const auto raw_out = plot_output_file(output_stem + "_raw").string();
    c_raw.SaveAs(raw_out.c_str());

    std::cout << "[plotFirstInferenceScoreEffPurity] saved plot: " << out << "\n";
    std::cout << "[plotFirstInferenceScoreEffPurity] saved raw-threshold plot: " << raw_out << "\n";
    std::cout << "[plotFirstInferenceScoreEffPurity] done\n";

    return 0;

  });
}
