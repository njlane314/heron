// macros/plot_inference_score_effpur_scan.C
//
// Scan thresholds on the first inference score entry and plot:
//   - signal efficiency
//   - purity (MC + EXT in denominator)
//   - efficiency x purity
//
// In the style of the modern Heron macros:
//   - Uses EventListIO + SelectionService::decorate
//   - Uses a single set of histograms and cumulative integrals (fast)
//   - Adds uncertainty bars using exact Clopper-Pearson intervals
//
// Definitions (per threshold t):
//   N_S,tot      = raw count of signal events after base_sel (MC only; EXT excluded)
//   N_S,pass(t)  = raw count of signal events with score >= t
//   N_all,pass(t)= raw count of all non-data events (MC + EXT) with score >= t
//
//   efficiency(t) = N_S,pass(t) / N_S,tot
//   purity(t)     = N_S,pass(t) / N_all,pass(t)
//   effpur(t)     = efficiency(t) * purity(t)
//
// Uncertainties:
//   - efficiency: Clopper-Pearson on (N_S,pass | N_S,tot)
//   - purity:     Clopper-Pearson on (N_S,pass | N_all,pass)
//   - effpur:     interval obtained by multiplying the above CP intervals:
//                 [eff_lo*pur_lo, eff_hi*pur_hi]
//
// Notes:
//   - The scan is cumulative; adjacent points are highly correlated.
//   - The CP intervals use raw counts (weights are not binomial). Weighted
//     efficiencies/purities are also computed and reported for reference.
//
// Run with:
//   ./heron macro plot_inference_score_effpur_scan.C
//   ./heron macro plot_inference_score_effpur_scan.C \
//     'plot_inference_score_effpur_scan("./scratch/out/event_list_myana.root")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TEfficiency.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

struct MetricStyle
{
    const char *label = "";
    int color = kBlack;
    int marker = 20;
};

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

// numerically-stable sigmoid
double sigmoid(double x)
{
    if (!std::isfinite(x))
        return 0.0;
    if (x >= 0.0)
    {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
}

double cumulative_hist_sum(const TH1 &h, int first_bin)
{
    const int last_bin = h.GetNbinsX() + 1; // include overflow
    double out = 0.0;
    for (int b = first_bin; b <= last_bin; ++b)
        out += h.GetBinContent(b);
    return out;
}

double pass_fraction_lo_cp(ULong64_t n_total, ULong64_t n_pass, double cl)
{
    if (n_total == 0)
        return 0.0;
    return TEfficiency::ClopperPearson(static_cast<unsigned int>(n_total),
                                       static_cast<unsigned int>(n_pass),
                                       cl,
                                       false);
}

double pass_fraction_hi_cp(ULong64_t n_total, ULong64_t n_pass, double cl)
{
    if (n_total == 0)
        return 0.0;
    return TEfficiency::ClopperPearson(static_cast<unsigned int>(n_total),
                                       static_cast<unsigned int>(n_pass),
                                       cl,
                                       true);
}

void set_graph_style(TGraphAsymmErrors &g, const MetricStyle &s)
{
    g.SetLineColor(s.color);
    g.SetMarkerColor(s.color);
    g.SetLineWidth(2);
    g.SetMarkerStyle(s.marker);
    g.SetMarkerSize(1.0);
}

// Book and compute curves from two histograms:
//   - signal histogram (MC only, filtered by signal_sel)
//   - all histogram (MC+EXT, includes signal)
// Both histograms must have the same binning, matching the threshold grid.
struct ScanCurves
{
    std::vector<double> x;
    std::vector<double> eff;
    std::vector<double> pur;
    std::vector<double> effpur;

    std::vector<double> eff_eyl;
    std::vector<double> eff_eyh;
    std::vector<double> pur_eyl;
    std::vector<double> pur_eyh;
    std::vector<double> effpur_eyl;
    std::vector<double> effpur_eyh;

    // Weighted central values (no binomial interval).
    std::vector<double> eff_w;
    std::vector<double> pur_w;
    std::vector<double> effpur_w;

    double best_x_raw = 0.0;
    double best_effpur_raw = -1.0;
    double best_x_w = 0.0;
    double best_effpur_w = -1.0;

    ULong64_t n_sig_total_raw = 0;
    double w_sig_total = 0.0;
};

ScanCurves compute_scan(const std::vector<double> &thresholds,
                        const TH1D &h_sig_raw,
                        const TH1D &h_all_raw,
                        const TH1D &h_sig_w,
                        const TH1D &h_all_w,
                        double cl)
{
    const int n = static_cast<int>(thresholds.size());
    ScanCurves out;
    out.x = thresholds;
    out.eff.resize(static_cast<std::size_t>(n), 0.0);
    out.pur.resize(static_cast<std::size_t>(n), 0.0);
    out.effpur.resize(static_cast<std::size_t>(n), 0.0);

    out.eff_eyl.resize(static_cast<std::size_t>(n), 0.0);
    out.eff_eyh.resize(static_cast<std::size_t>(n), 0.0);
    out.pur_eyl.resize(static_cast<std::size_t>(n), 0.0);
    out.pur_eyh.resize(static_cast<std::size_t>(n), 0.0);
    out.effpur_eyl.resize(static_cast<std::size_t>(n), 0.0);
    out.effpur_eyh.resize(static_cast<std::size_t>(n), 0.0);

    out.eff_w.resize(static_cast<std::size_t>(n), 0.0);
    out.pur_w.resize(static_cast<std::size_t>(n), 0.0);
    out.effpur_w.resize(static_cast<std::size_t>(n), 0.0);

    // Totals are the histogram integrals including underflow + overflow.
    // (Events outside the plotted score range must still count in the
    //  efficiency denominator; they simply never pass the threshold grid.)
    const ULong64_t n_sig_total_raw =
        static_cast<ULong64_t>(std::llround(cumulative_hist_sum(h_sig_raw, 0)));
    const double w_sig_total = cumulative_hist_sum(h_sig_w, 0);
    out.n_sig_total_raw = n_sig_total_raw;
    out.w_sig_total = w_sig_total;

    for (int i = 0; i < n; ++i)
    {
        const int first_bin = i + 1;

        const ULong64_t n_sig_pass_raw = static_cast<ULong64_t>(std::llround(cumulative_hist_sum(h_sig_raw, first_bin)));
        const ULong64_t n_all_pass_raw = static_cast<ULong64_t>(std::llround(cumulative_hist_sum(h_all_raw, first_bin)));

        const double w_sig_pass = cumulative_hist_sum(h_sig_w, first_bin);
        const double w_all_pass = cumulative_hist_sum(h_all_w, first_bin);

        // Raw central values.
        const double eff = (n_sig_total_raw > 0) ? (static_cast<double>(n_sig_pass_raw) / static_cast<double>(n_sig_total_raw)) : 0.0;
        const double pur = (n_all_pass_raw > 0) ? (static_cast<double>(n_sig_pass_raw) / static_cast<double>(n_all_pass_raw)) : 0.0;
        const double effpur = eff * pur;

        out.eff[static_cast<std::size_t>(i)] = eff;
        out.pur[static_cast<std::size_t>(i)] = pur;
        out.effpur[static_cast<std::size_t>(i)] = effpur;

        // CP intervals on raw fractions.
        const double eff_lo = (n_sig_total_raw > 0) ? pass_fraction_lo_cp(n_sig_total_raw, n_sig_pass_raw, cl) : 0.0;
        const double eff_hi = (n_sig_total_raw > 0) ? pass_fraction_hi_cp(n_sig_total_raw, n_sig_pass_raw, cl) : 0.0;

        const double pur_lo = (n_all_pass_raw > 0) ? pass_fraction_lo_cp(n_all_pass_raw, n_sig_pass_raw, cl) : 0.0;
        const double pur_hi = (n_all_pass_raw > 0) ? pass_fraction_hi_cp(n_all_pass_raw, n_sig_pass_raw, cl) : 0.0;

        out.eff_eyl[static_cast<std::size_t>(i)] = std::max(0.0, eff - eff_lo);
        out.eff_eyh[static_cast<std::size_t>(i)] = std::max(0.0, eff_hi - eff);

        out.pur_eyl[static_cast<std::size_t>(i)] = std::max(0.0, pur - pur_lo);
        out.pur_eyh[static_cast<std::size_t>(i)] = std::max(0.0, pur_hi - pur);

        const double effpur_lo = eff_lo * pur_lo;
        const double effpur_hi = eff_hi * pur_hi;

        out.effpur_eyl[static_cast<std::size_t>(i)] = std::max(0.0, effpur - effpur_lo);
        out.effpur_eyh[static_cast<std::size_t>(i)] = std::max(0.0, effpur_hi - effpur);

        // Weighted central values for reference.
        const double eff_w = (w_sig_total > 0.0) ? (w_sig_pass / w_sig_total) : 0.0;
        const double pur_w = (w_all_pass > 0.0) ? (w_sig_pass / w_all_pass) : 0.0;
        const double effpur_w = eff_w * pur_w;

        out.eff_w[static_cast<std::size_t>(i)] = eff_w;
        out.pur_w[static_cast<std::size_t>(i)] = pur_w;
        out.effpur_w[static_cast<std::size_t>(i)] = effpur_w;

        // Track best thresholds.
        if (effpur > out.best_effpur_raw)
        {
            out.best_effpur_raw = effpur;
            out.best_x_raw = thresholds[static_cast<std::size_t>(i)];
        }
        if (effpur_w > out.best_effpur_w)
        {
            out.best_effpur_w = effpur_w;
            out.best_x_w = thresholds[static_cast<std::size_t>(i)];
        }
    }

    return out;
}

} // namespace

int plot_inference_score_effpur_scan(const std::string &event_list_path = "",
                                     const std::string &base_sel = "sel_muon",
                                     const std::string &signal_sel = "is_signal",
                                     const std::string &mc_weight = "w_nominal",
                                     int n_thresholds = 101,
                                     double raw_threshold_min = -15.0,
                                     double raw_threshold_max = 15.0,
                                     const std::string &output_stem = "inference_score_effpur_scan",
                                     double cl = 0.682689492137086)
{
    return heron::macro::run_with_guard("plot_inference_score_effpur_scan", [&]() -> int {
        ROOT::EnableImplicitMT();
        TH1::SetDefaultSumw2();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_inference_score_effpur_scan] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inference_score_effpur_scan] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        if (n_thresholds < 2)
            n_thresholds = 2;

        if (raw_threshold_max < raw_threshold_min)
            std::swap(raw_threshold_min, raw_threshold_max);
        if (raw_threshold_max == raw_threshold_min)
            raw_threshold_max = raw_threshold_min + 1.0;

        // --- Prepare threshold grids and histogram bin edges ---
        // Sigmoid scan in [0,1].
        std::vector<double> thr_sig(static_cast<std::size_t>(n_thresholds));
        std::vector<double> edges_sig(static_cast<std::size_t>(n_thresholds + 1));
        const double step_sig = 1.0 / static_cast<double>(n_thresholds - 1);
        for (int i = 0; i < n_thresholds; ++i)
        {
            const double t = static_cast<double>(i) * step_sig;
            thr_sig[static_cast<std::size_t>(i)] = t;
            edges_sig[static_cast<std::size_t>(i)] = t;
        }
        edges_sig.back() = 1.0 + step_sig;

        // Raw scan.
        std::vector<double> thr_raw(static_cast<std::size_t>(n_thresholds));
        std::vector<double> edges_raw(static_cast<std::size_t>(n_thresholds + 1));
        const double step_raw = (raw_threshold_max - raw_threshold_min) / static_cast<double>(n_thresholds - 1);
        for (int i = 0; i < n_thresholds; ++i)
        {
            const double t = raw_threshold_min + static_cast<double>(i) * step_raw;
            thr_raw[static_cast<std::size_t>(i)] = t;
            edges_raw[static_cast<std::size_t>(i)] = t;
        }
        edges_raw.back() = raw_threshold_max + step_raw;

        // --- Load event list and define score columns ---
        EventListIO el(input_path);

        ROOT::RDF::RNode base = SelectionService::decorate(el.rdf())
                                   .Define("inf_score_0",
                                           [](const ROOT::RVec<float> &scores) {
                                               // Missing score -> very negative (fails any reasonable cut).
                                               return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                           },
                                           {"inf_scores"})
                                   .Define("inf_score_0_sigmoid",
                                           [](double x) { return sigmoid(x); },
                                           {"inf_score_0"});

        // --- Sample masks ---
        auto mask_ext = el.mask_for_ext();
        auto mask_mc_like = el.mask_for_mc_like(); // non-data (includes EXT)

        // Node with ALL non-data (MC + EXT). This is the purity denominator.
        ROOT::RDF::RNode node_all = filter_by_sample_mask(base, mask_mc_like, "sample_id")
                                        .Define("__w__", mc_weight);

        // MC-only (exclude EXT). This is used for signal efficiency.
        ROOT::RDF::RNode node_mc = filter_not_sample_mask(node_all, mask_ext, "sample_id");

        // Apply base selection.
        if (!base_sel.empty())
        {
            if (has_column(node_all, base_sel))
            {
                node_all = node_all.Filter([](bool pass) { return pass; }, {base_sel});
                node_mc = node_mc.Filter([](bool pass) { return pass; }, {base_sel});
            }
            else
            {
                node_all = node_all.Filter(base_sel);
                node_mc = node_mc.Filter(base_sel);
            }
        }

        // Signal node (MC only).
        ROOT::RDF::RNode node_sig = node_mc.Filter(signal_sel);

        // Quick sanity check.
        const ULong64_t n_sig_total_raw_check = *node_sig.Count();
        if (n_sig_total_raw_check == 0)
        {
            std::cerr << "[plot_inference_score_effpur_scan] signal denominator is 0 for signal_sel='"
                      << signal_sel << "' after base_sel='" << base_sel << "'.\n";
            return 1;
        }

        // --- Book histograms once (fast scans) ---
        // Use unique histogram names to avoid ROOT object-name collisions.
        ROOT::RDF::TH1DModel h_sig_sigmoid_raw_model("h_sig_sigmoid_raw", "", n_thresholds, edges_sig.data());
        ROOT::RDF::TH1DModel h_all_sigmoid_raw_model("h_all_sigmoid_raw", "", n_thresholds, edges_sig.data());
        ROOT::RDF::TH1DModel h_sig_sigmoid_w_model("h_sig_sigmoid_w", "", n_thresholds, edges_sig.data());
        ROOT::RDF::TH1DModel h_all_sigmoid_w_model("h_all_sigmoid_w", "", n_thresholds, edges_sig.data());

        auto h_sig_sigmoid_raw = node_sig.Histo1D(h_sig_sigmoid_raw_model, "inf_score_0_sigmoid");
        auto h_all_sigmoid_raw = node_all.Histo1D(h_all_sigmoid_raw_model, "inf_score_0_sigmoid");

        auto h_sig_sigmoid_w = node_sig.Histo1D(h_sig_sigmoid_w_model, "inf_score_0_sigmoid", "__w__");
        auto h_all_sigmoid_w = node_all.Histo1D(h_all_sigmoid_w_model, "inf_score_0_sigmoid", "__w__");

        ROOT::RDF::TH1DModel h_sig_raw_raw_model("h_sig_raw_raw", "", n_thresholds, edges_raw.data());
        ROOT::RDF::TH1DModel h_all_raw_raw_model("h_all_raw_raw", "", n_thresholds, edges_raw.data());
        ROOT::RDF::TH1DModel h_sig_raw_w_model("h_sig_raw_w", "", n_thresholds, edges_raw.data());
        ROOT::RDF::TH1DModel h_all_raw_w_model("h_all_raw_w", "", n_thresholds, edges_raw.data());

        auto h_sig_raw_raw = node_sig.Histo1D(h_sig_raw_raw_model, "inf_score_0");
        auto h_all_raw_raw = node_all.Histo1D(h_all_raw_raw_model, "inf_score_0");

        auto h_sig_raw_w = node_sig.Histo1D(h_sig_raw_w_model, "inf_score_0", "__w__");
        auto h_all_raw_w = node_all.Histo1D(h_all_raw_w_model, "inf_score_0", "__w__");

        // Force evaluation now.
        const TH1D &hs_sigmoid_raw = *h_sig_sigmoid_raw;
        const TH1D &ha_sigmoid_raw = *h_all_sigmoid_raw;
        const TH1D &hs_sigmoid_w = *h_sig_sigmoid_w;
        const TH1D &ha_sigmoid_w = *h_all_sigmoid_w;

        const TH1D &hs_raw_raw = *h_sig_raw_raw;
        const TH1D &ha_raw_raw = *h_all_raw_raw;
        const TH1D &hs_raw_w = *h_sig_raw_w;
        const TH1D &ha_raw_w = *h_all_raw_w;

        // --- Compute curves ---
        const ScanCurves scan_sig = compute_scan(thr_sig, hs_sigmoid_raw, ha_sigmoid_raw, hs_sigmoid_w, ha_sigmoid_w, cl);
        const ScanCurves scan_raw = compute_scan(thr_raw, hs_raw_raw, ha_raw_raw, hs_raw_w, ha_raw_w, cl);

        std::cout << "[plot_inference_score_effpur_scan] sigmoid scan best (raw eff*pur) thr="
                  << scan_sig.best_x_raw << " value=" << scan_sig.best_effpur_raw
                  << " | best (weighted eff*pur) thr=" << scan_sig.best_x_w
                  << " value=" << scan_sig.best_effpur_w << "\n";

        std::cout << "[plot_inference_score_effpur_scan] raw scan best (raw eff*pur) thr="
                  << scan_raw.best_x_raw << " value=" << scan_raw.best_effpur_raw
                  << " | best (weighted eff*pur) thr=" << scan_raw.best_x_w
                  << " value=" << scan_raw.best_effpur_w << "\n";

        // --- Plot: sigmoid threshold scan ---
        {
            Plotter plotter;
            plotter.set_global_style();
            gStyle->SetOptStat(0);

            TCanvas c("c_inf_score_effpur_sigmoid", "Inference-score threshold scan (sigmoid): efficiency, purity, and efficiency#timespurity", 1150, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            TH1D h_frame("h_frame_effpur_sigmoid",
                         ";Sigmoid(inference score [0]) threshold;metric value",
                         100,
                         0.0,
                         1.0);
            h_frame.SetMinimum(0.0);
            h_frame.SetMaximum(1.05);
            h_frame.Draw("AXIS");

            std::vector<double> ex0(scan_sig.x.size(), 0.0);

            TGraphAsymmErrors g_eff(static_cast<int>(scan_sig.x.size()),
                                    scan_sig.x.data(),
                                    scan_sig.eff.data(),
                                    ex0.data(),
                                    ex0.data(),
                                    scan_sig.eff_eyl.data(),
                                    scan_sig.eff_eyh.data());

            TGraphAsymmErrors g_pur(static_cast<int>(scan_sig.x.size()),
                                    scan_sig.x.data(),
                                    scan_sig.pur.data(),
                                    ex0.data(),
                                    ex0.data(),
                                    scan_sig.pur_eyl.data(),
                                    scan_sig.pur_eyh.data());

            TGraphAsymmErrors g_effpur(static_cast<int>(scan_sig.x.size()),
                                       scan_sig.x.data(),
                                       scan_sig.effpur.data(),
                                       ex0.data(),
                                       ex0.data(),
                                       scan_sig.effpur_eyl.data(),
                                       scan_sig.effpur_eyh.data());

            const MetricStyle s_eff{"efficiency", kBlue + 1, 20};
            const MetricStyle s_pur{"purity", kRed + 1, 21};
            const MetricStyle s_effpur{"efficiency #times purity", kGreen + 2, 22};

            set_graph_style(g_eff, s_eff);
            set_graph_style(g_pur, s_pur);
            set_graph_style(g_effpur, s_effpur);

            g_eff.Draw("PZ SAME");
            g_eff.Draw("L SAME");

            g_pur.Draw("PZ SAME");
            g_pur.Draw("L SAME");

            g_effpur.Draw("PZ SAME");
            g_effpur.Draw("L SAME");

            // Mark best weighted threshold (commonly used for optimisation).
            TLine best_line(scan_sig.best_x_w, 0.0, scan_sig.best_x_w, 1.05);
            best_line.SetLineStyle(2);
            best_line.SetLineWidth(2);
            best_line.Draw();

            TLegend leg(0.14, 0.70, 0.55, 0.88);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(&g_eff, "efficiency (raw, CP errors)", "lp");
            leg.AddEntry(&g_pur, "purity (raw, CP errors)", "lp");
            leg.AddEntry(&g_effpur, "efficiency #times purity", "lp");

            std::ostringstream best_lbl;
            best_lbl << "best weighted eff#timespur @ t=" << std::fixed << std::setprecision(3) << scan_sig.best_x_w;
            leg.AddEntry(&best_line, best_lbl.str().c_str(), "l");
            leg.Draw();

            c.RedrawAxis();

            const auto out = plot_output_file(output_stem).string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_inference_score_effpur_scan] saved plot: " << out << "\n";
        }

        // --- Plot: raw-threshold scan ---
        {
            Plotter plotter;
            plotter.set_global_style();
            gStyle->SetOptStat(0);

            TCanvas c("c_inf_score_effpur_raw", "Inference-score threshold scan (raw): efficiency, purity, and efficiency#timespurity", 1150, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            TH1D h_frame("h_frame_effpur_raw",
                         ";Inference score [0] threshold;metric value",
                         100,
                         raw_threshold_min,
                         raw_threshold_max);
            h_frame.SetMinimum(0.0);
            h_frame.SetMaximum(1.05);
            h_frame.Draw("AXIS");

            std::vector<double> ex0(scan_raw.x.size(), 0.0);

            TGraphAsymmErrors g_eff(static_cast<int>(scan_raw.x.size()),
                                    scan_raw.x.data(),
                                    scan_raw.eff.data(),
                                    ex0.data(),
                                    ex0.data(),
                                    scan_raw.eff_eyl.data(),
                                    scan_raw.eff_eyh.data());

            TGraphAsymmErrors g_pur(static_cast<int>(scan_raw.x.size()),
                                    scan_raw.x.data(),
                                    scan_raw.pur.data(),
                                    ex0.data(),
                                    ex0.data(),
                                    scan_raw.pur_eyl.data(),
                                    scan_raw.pur_eyh.data());

            TGraphAsymmErrors g_effpur(static_cast<int>(scan_raw.x.size()),
                                       scan_raw.x.data(),
                                       scan_raw.effpur.data(),
                                       ex0.data(),
                                       ex0.data(),
                                       scan_raw.effpur_eyl.data(),
                                       scan_raw.effpur_eyh.data());

            const MetricStyle s_eff{"efficiency", kBlue + 1, 20};
            const MetricStyle s_pur{"purity", kRed + 1, 21};
            const MetricStyle s_effpur{"efficiency #times purity", kGreen + 2, 22};

            set_graph_style(g_eff, s_eff);
            set_graph_style(g_pur, s_pur);
            set_graph_style(g_effpur, s_effpur);

            g_eff.Draw("PZ SAME");
            g_eff.Draw("L SAME");

            g_pur.Draw("PZ SAME");
            g_pur.Draw("L SAME");

            g_effpur.Draw("PZ SAME");
            g_effpur.Draw("L SAME");

            TLine best_line(scan_raw.best_x_w, 0.0, scan_raw.best_x_w, 1.05);
            best_line.SetLineStyle(2);
            best_line.SetLineWidth(2);
            best_line.Draw();

            TLegend leg(0.14, 0.70, 0.55, 0.88);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(&g_eff, "efficiency (raw, CP errors)", "lp");
            leg.AddEntry(&g_pur, "purity (raw, CP errors)", "lp");
            leg.AddEntry(&g_effpur, "efficiency #times purity", "lp");

            std::ostringstream best_lbl;
            best_lbl << "best weighted eff#timespur @ t=" << std::fixed << std::setprecision(3) << scan_raw.best_x_w;
            leg.AddEntry(&best_line, best_lbl.str().c_str(), "l");
            leg.Draw();

            c.RedrawAxis();

            const auto out = plot_output_file(output_stem + "_raw").string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_inference_score_effpur_scan] saved raw-threshold plot: " << out << "\n";
        }

        // --- Report some totals for context ---
        std::cout << "\n[plot_inference_score_effpur_scan] totals after base_sel='" << base_sel << "'\n";
        std::cout << "  signal_sel='" << signal_sel << "'\n";
        std::cout << "  signal total raw  = " << scan_raw.n_sig_total_raw << "\n";
        std::cout << "  signal total w    = " << std::fixed << std::setprecision(3) << scan_raw.w_sig_total << "\n";
        std::cout << "  cl (CP intervals) = " << cl << "\n";
        std::cout << "[plot_inference_score_effpur_scan] done\n";

        return 0;
    });
}

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif
