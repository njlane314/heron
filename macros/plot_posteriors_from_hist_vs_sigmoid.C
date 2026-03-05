// plot/macro/plot_posteriors_from_hist_vs_sigmoid.C
//
// Compare the posterior signal probability inferred from score histograms to
// the network's sigmoid(score) output.
//
// In bins of the first inference score s = inf_scores[0], define:
//   p_hist_raw(s_bin) = N_S,raw(bin) / N_total,raw(bin)
// and compute exact Clopper-Pearson intervals (binomial) on p_hist_raw.
//
// Also compute the (optionally analysis-weighted) posterior and mean predicted
// probability in the same score bins:
//   p_hist_w(bin)   = SumW_S(bin) / SumW_total(bin)
//   <sigmoid>_w(bin)= SumW(sigmoid(score)) / SumW_total(bin)
//
// If the score is an (approximate) log-likelihood ratio, then the posterior
// under the sample priors should satisfy:
//   P(S|s) = sigmoid( s + log(pi_S/pi_B) )
// This macro optionally overlays the *mean* of sigmoid(s + log_prior_odds)
// in each bin, where log_prior_odds is estimated from the weighted totals in
// the plotted score range.
//
// Run with:
//   ./heron macro plot_posteriors_from_hist_vs_sigmoid.C
//   ./heron macro plot_posteriors_from_hist_vs_sigmoid.C \
//     'plot_posteriors_from_hist_vs_sigmoid("./scratch/out/event_list_myana.root")'
//
// Outputs:
//   - a single canvas with:
//       * p_hist_raw (points with CP error bars)
//       * p_hist_w   (open markers, no CP bars)
//       * <sigmoid>_w (line + markers)
//       * <sigmoid(s + log_prior_odds)>_w (dashed)

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
#include <TGraph.h>
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
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace {

bool has_column(ROOT::RDF::RNode node, const std::string& name) {
  const auto cols = node.GetColumnNames();
  return std::find(cols.begin(), cols.end(), name) != cols.end();
}

double cp_lo(unsigned int n_total, unsigned int n_pass, double cl) {
  return TEfficiency::ClopperPearson(n_total, n_pass, cl, false);
}

double cp_hi(unsigned int n_total, unsigned int n_pass, double cl) {
  return TEfficiency::ClopperPearson(n_total, n_pass, cl, true);
}

void set_style_raw(TGraphAsymmErrors& g) {
  g.SetLineColor(kBlack);
  g.SetMarkerColor(kBlack);
  g.SetLineWidth(2);
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.0);
}

void set_style_open(TGraph& g) {
  g.SetLineColor(kBlack);
  g.SetMarkerColor(kBlack);
  g.SetLineWidth(2);
  g.SetMarkerStyle(24);
  g.SetMarkerSize(1.0);
}

void set_style_pred(TGraph& g, int ls) {
  g.SetLineColor(kRed + 1);
  g.SetMarkerColor(kRed + 1);
  g.SetLineWidth(2);
  g.SetLineStyle(ls);
  g.SetMarkerStyle(21);
  g.SetMarkerSize(1.0);
}

}  // namespace

int plot_posteriors_from_hist_vs_sigmoid(const std::string& event_list_path = "",
                                         const std::string& base_sel = "sel_muon",
                                         const std::string& signal_sel = "is_signal",
                                         const std::string& mc_weight = "w_nominal",
                                         int nbins = 30,
                                         double raw_xmin = -15.0,
                                         double raw_xmax = 15.0,
                                         double cl = 0.682689492137086,
                                         bool overlay_prior_shift = true,
                                         const std::string& output_stem =
                                             "plot_posteriors_from_hist_vs_sigmoid") {
  return heron::macro::run_with_guard("plot_posteriors_from_hist_vs_sigmoid", [&]() -> int {
    ROOT::EnableImplicitMT();
    TH1::SetDefaultSumw2();

    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
    std::cout << "[plot_posteriors_from_hist_vs_sigmoid] input=" << input_path << "\n";

    if (!looks_like_event_list_root(input_path)) {
      std::cerr << "[plot_posteriors_from_hist_vs_sigmoid] input is not an event-list root file: "
                << input_path << "\n";
      return 1;
    }

    if (nbins < 1) nbins = 1;
    if (raw_xmax < raw_xmin) std::swap(raw_xmin, raw_xmax);
    if (raw_xmax == raw_xmin) raw_xmax = raw_xmin + 1.0;

    EventListIO el(input_path);

    ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf())
                               .Define("inf_score_0",
                                       [](const ROOT::RVec<float>& scores) {
                                         return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                       },
                                       {"inf_scores"})
                               .Define("pred_sigmoid",
                                       [](double x) {
                                         if (!std::isfinite(x)) return 0.0;
                                         // Stable-ish sigmoid for large |x|.
                                         if (x >= 0) {
                                           const double z = std::exp(-x);
                                           return 1.0 / (1.0 + z);
                                         }
                                         const double z = std::exp(x);
                                         return z / (1.0 + z);
                                       },
                                       {"inf_score_0"});

    auto mask_mc_like = el.mask_for_mc_like();
    ROOT::RDF::RNode node =
        filter_by_sample_mask(rdf, mask_mc_like, "sample_id").Define("__w__", mc_weight).Define(
            "__wpred__", "__w__ * pred_sigmoid");

    if (!base_sel.empty()) {
      if (has_column(rdf, base_sel))
        node = node.Filter([](bool pass) { return pass; }, {base_sel});
      else
        node = node.Filter(base_sel);
    }

    ROOT::RDF::RNode node_sig = node.Filter(signal_sel);

    // Histograms in score bins.
    ROOT::RDF::TH1DModel h_raw_total_m("h_raw_total", "", nbins, raw_xmin, raw_xmax);
    ROOT::RDF::TH1DModel h_raw_sig_m("h_raw_sig", "", nbins, raw_xmin, raw_xmax);
    ROOT::RDF::TH1DModel h_raw_predsum_m("h_raw_predsum", "", nbins, raw_xmin, raw_xmax);

    ROOT::RDF::TH1DModel h_w_total_m("h_w_total", "", nbins, raw_xmin, raw_xmax);
    ROOT::RDF::TH1DModel h_w_sig_m("h_w_sig", "", nbins, raw_xmin, raw_xmax);
    ROOT::RDF::TH1DModel h_w_predsum_m("h_w_predsum", "", nbins, raw_xmin, raw_xmax);

    auto h_raw_total = node.Histo1D(h_raw_total_m, "inf_score_0");
    auto h_raw_sig = node_sig.Histo1D(h_raw_sig_m, "inf_score_0");
    auto h_raw_predsum = node.Histo1D(h_raw_predsum_m, "inf_score_0", "pred_sigmoid");

    auto h_w_total = node.Histo1D(h_w_total_m, "inf_score_0", "__w__");
    auto h_w_sig = node_sig.Histo1D(h_w_sig_m, "inf_score_0", "__w__");
    auto h_w_predsum = node.Histo1D(h_w_predsum_m, "inf_score_0", "__wpred__");

    // Totals in the plotted range for a prior-odds estimate.
    const double sumw_sig_total = *(node_sig.Sum<double>("__w__"));
    const double sumw_all_total = *(node.Sum<double>("__w__"));
    const double sumw_bkg_total = std::max(0.0, sumw_all_total - sumw_sig_total);

    double log_prior_odds = 0.0;
    bool have_prior_odds = false;
    if (sumw_sig_total > 0.0 && sumw_bkg_total > 0.0) {
      log_prior_odds = std::log(sumw_sig_total / sumw_bkg_total);
      have_prior_odds = true;
    }

    // If requested, book shifted-pred sum histogram.
    ROOT::RDF::RResultPtr<TH1D> h_w_predsum_shift;
    if (overlay_prior_shift && have_prior_odds) {
      node = node
                 .Define("pred_shift",
                         [log_prior_odds](double s) {
                           const double x = s + log_prior_odds;
                           if (!std::isfinite(x)) return 0.0;
                           if (x >= 0) {
                             const double z = std::exp(-x);
                             return 1.0 / (1.0 + z);
                           }
                           const double z = std::exp(x);
                           return z / (1.0 + z);
                         },
                         {"inf_score_0"})
                 .Define("__wpred_shift__", "__w__ * pred_shift");

      ROOT::RDF::TH1DModel h_w_predsum_shift_m("h_w_predsum_shift", "", nbins, raw_xmin, raw_xmax);
      h_w_predsum_shift = node.Histo1D(h_w_predsum_shift_m, "inf_score_0", "__wpred_shift__");
    }

    // Build points.
    std::vector<double> x;
    std::vector<double> y_p_raw;
    std::vector<double> y_p_raw_el;
    std::vector<double> y_p_raw_eh;

    std::vector<double> x_w;
    std::vector<double> y_p_w;
    std::vector<double> y_pred_w;
    std::vector<double> y_pred_shift_w;

    x.reserve(static_cast<std::size_t>(nbins));
    y_p_raw.reserve(static_cast<std::size_t>(nbins));
    y_p_raw_el.reserve(static_cast<std::size_t>(nbins));
    y_p_raw_eh.reserve(static_cast<std::size_t>(nbins));

    x_w.reserve(static_cast<std::size_t>(nbins));
    y_p_w.reserve(static_cast<std::size_t>(nbins));
    y_pred_w.reserve(static_cast<std::size_t>(nbins));
    y_pred_shift_w.reserve(static_cast<std::size_t>(nbins));

    const TH1D& Hrt = *h_raw_total;
    const TH1D& Hrs = *h_raw_sig;
    const TH1D& Hrp = *h_raw_predsum;

    const TH1D& Hwt = *h_w_total;
    const TH1D& Hws = *h_w_sig;
    const TH1D& Hwp = *h_w_predsum;

    const TH1D* Hwp_shift = nullptr;
    if (overlay_prior_shift && have_prior_odds) Hwp_shift = &(*h_w_predsum_shift);

    for (int b = 1; b <= nbins; ++b) {
      const double xc = Hrt.GetXaxis()->GetBinCenter(b);

      const double n_tot = Hrt.GetBinContent(b);
      const double n_sig = Hrs.GetBinContent(b);

      if (n_tot > 0.0) {
        const unsigned int Ntot = static_cast<unsigned int>(std::llround(n_tot));
        const unsigned int Nsig = static_cast<unsigned int>(std::llround(n_sig));

        // Defensive: clip.
        const unsigned int Nsig_clip = (Nsig > Ntot) ? Ntot : Nsig;

        const double p = static_cast<double>(Nsig_clip) / static_cast<double>(Ntot);
        const double lo = cp_lo(Ntot, Nsig_clip, cl);
        const double hi = cp_hi(Ntot, Nsig_clip, cl);

        x.push_back(xc);
        y_p_raw.push_back(p);
        y_p_raw_el.push_back(p - lo);
        y_p_raw_eh.push_back(hi - p);
      }

      const double sw_tot = Hwt.GetBinContent(b);
      const double sw_sig = Hws.GetBinContent(b);
      const double sw_pred = Hwp.GetBinContent(b);

      if (sw_tot > 0.0) {
        const double p_w = sw_sig / sw_tot;
        const double mean_pred_w = sw_pred / sw_tot;

        x_w.push_back(xc);
        y_p_w.push_back(std::clamp(p_w, 0.0, 1.0));
        y_pred_w.push_back(std::clamp(mean_pred_w, 0.0, 1.0));

        if (Hwp_shift) {
          const double sw_pred_shift = Hwp_shift->GetBinContent(b);
          const double mean_pred_shift_w = sw_pred_shift / sw_tot;
          y_pred_shift_w.push_back(std::clamp(mean_pred_shift_w, 0.0, 1.0));
        }
      }
    }

    if (x.empty() && x_w.empty()) {
      std::cerr << "[plot_posteriors_from_hist_vs_sigmoid] no populated bins after base selection='"
                << base_sel << "'.\n";
      return 1;
    }

    // Graphs.
    std::vector<double> ex0_raw(x.size(), 0.0);

    TGraphAsymmErrors g_p_raw(static_cast<int>(x.size()),
                              x.data(),
                              y_p_raw.data(),
                              ex0_raw.data(),
                              ex0_raw.data(),
                              y_p_raw_el.data(),
                              y_p_raw_eh.data());
    set_style_raw(g_p_raw);

    TGraph g_p_w(static_cast<int>(x_w.size()), x_w.data(), y_p_w.data());
    set_style_open(g_p_w);

    TGraph g_pred_w(static_cast<int>(x_w.size()), x_w.data(), y_pred_w.data());
    set_style_pred(g_pred_w, 1);

    std::unique_ptr<TGraph> g_pred_shift;
    if (overlay_prior_shift && have_prior_odds && y_pred_shift_w.size() == x_w.size()) {
      g_pred_shift = std::make_unique<TGraph>(static_cast<int>(x_w.size()), x_w.data(), y_pred_shift_w.data());
      set_style_pred(*g_pred_shift, 2);
    }

    // Plot.
    Plotter plotter;
    plotter.set_global_style();
    gStyle->SetOptStat(0);

    TCanvas c("c_posteriors_vs_sigmoid", "Posterior from hist vs sigmoid", 1100, 850);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.12);

    TH1D h_frame("h_frame_post", ";Inference score [0];P(signal | score-bin)", 100, raw_xmin, raw_xmax);
    h_frame.SetMinimum(0.0);
    h_frame.SetMaximum(1.05);
    h_frame.Draw("AXIS");

    // Draw in an order that keeps the CP bars visible.
    g_pred_w.Draw("LP SAME");
    if (g_pred_shift) g_pred_shift->Draw("LP SAME");
    g_p_w.Draw("P SAME");
    g_p_raw.Draw("PZ SAME");

    // Reference line at 0.5
    TLine mid(raw_xmin, 0.5, raw_xmax, 0.5);
    mid.SetLineStyle(3);
    mid.SetLineColor(16);
    mid.Draw();

    TLegend leg(0.12, 0.14, 0.96, 0.34);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetNColumns(2);

    std::ostringstream s1;
    s1 << "p_{hist}^{raw} (CP " << std::fixed << std::setprecision(1) << (cl * 100.0) << "% bars)";
    leg.AddEntry(&g_p_raw, s1.str().c_str(), "lp");
    leg.AddEntry(&g_p_w, "p_{hist}^{weighted}", "p");
    leg.AddEntry(&g_pred_w, "<sigmoid(score)>_{weighted}", "lp");

    if (g_pred_shift) {
      std::ostringstream ss;
      ss << "<sigmoid(score + log #pi_{S}/#pi_{B})>_{weighted}, log #pi_{S}/#pi_{B}=" << std::fixed
         << std::setprecision(3) << log_prior_odds;
      leg.AddEntry(g_pred_shift.get(), ss.str().c_str(), "lp");
    }
    leg.Draw();

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());

    std::cout << "\n[plot_posteriors_from_hist_vs_sigmoid] totals in plotted range (weighted):\n";
    std::cout << "  sumw(signal) = " << sumw_sig_total << "\n";
    std::cout << "  sumw(bkg)    = " << sumw_bkg_total << "\n";
    if (have_prior_odds)
      std::cout << "  log_prior_odds = log(pi_S/pi_B) = " << log_prior_odds << "\n";
    else
      std::cout << "  log_prior_odds not available (need both signal and bkg > 0).\n";

    std::cout << "\n[plot_posteriors_from_hist_vs_sigmoid] saved plot: " << out << "\n";
    std::cout << "[plot_posteriors_from_hist_vs_sigmoid] done\n";

    return 0;
  });
}

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif
