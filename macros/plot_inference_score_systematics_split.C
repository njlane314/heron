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
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"

using namespace nu;

namespace {

struct HistSummary {
  std::vector<double> sumw;
  std::vector<double> sumw2;
};

using UShortRVec = ROOT::RVec<unsigned short>;

struct SplitComponent {
  std::string name;
  TMatrixD total;
  TMatrixD signal;
  TMatrixD background;
  std::size_t n_universes = 0;
  bool is_mc_source = true;
};

bool has_column(ROOT::RDF::RNode node, const std::string& name) {
  const auto cols = node.GetColumnNames();
  return std::find(cols.begin(), cols.end(), name) != cols.end();
}

bool is_simple_identifier(const std::string& expr) {
  if (expr.empty()) return false;
  const unsigned char first = static_cast<unsigned char>(expr.front());
  if (!(std::isalpha(first) || expr.front() == '_')) return false;
  for (size_t i = 1; i < expr.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(expr[i]);
    if (!(std::isalnum(c) || expr[i] == '_')) return false;
  }
  return true;
}

ROOT::RDF::RNode apply_optional_filter(ROOT::RDF::RNode node, const std::string& expr, const std::string& label) {
  if (expr.empty()) return node;

  if (has_column(node, expr)) {
    return node.Filter([](bool pass) { return pass; }, {expr});
  }

  if (is_simple_identifier(expr)) {
    std::cerr << "[plot_inference_score_systematics_split] " << label << " column '" << expr
              << "' is missing; skipping that filter.\n";
    return node;
  }

  return node.Filter(expr);
}

ROOT::RDF::RNode apply_negated_filter(ROOT::RDF::RNode node, const std::string& expr, const std::string& label) {
  if (expr.empty()) return node;

  if (has_column(node, expr)) {
    return node.Filter([](bool pass) { return !pass; }, {expr});
  }

  if (is_simple_identifier(expr)) {
    std::cerr << "[plot_inference_score_systematics_split] " << label << " column '" << expr
              << "' is missing; skipping negated filter.\n";
    return node;
  }

  return node.Filter("!(" + expr + ")");
}

int find_bin(double x, int nbins, double xmin, double xmax, bool fold_overflow) {
  if (nbins <= 0 || !(xmax > xmin)) return -1;
  if (x < xmin) return fold_overflow ? 0 : -1;
  if (x >= xmax) return fold_overflow ? (nbins - 1) : -1;

  const double bw = (xmax - xmin) / static_cast<double>(nbins);
  int bin = static_cast<int>((x - xmin) / bw);
  if (bin < 0) bin = 0;
  if (bin >= nbins) bin = nbins - 1;
  return bin;
}

HistSummary make_empty_hist(int nbins) {
  HistSummary h;
  h.sumw.assign(static_cast<size_t>(nbins), 0.0);
  h.sumw2.assign(static_cast<size_t>(nbins), 0.0);
  return h;
}

TMatrixD make_zero_matrix(int nbins) {
  TMatrixD out(nbins, nbins);
  out.Zero();
  return out;
}

HistSummary build_hist_summary(const std::vector<float>& x, const std::vector<double>& w, int nbins, double xmin,
                               double xmax, bool fold_overflow) {
  HistSummary out = make_empty_hist(nbins);
  const size_t n = std::min(x.size(), w.size());
  for (size_t i = 0; i < n; ++i) {
    const int bin = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (bin < 0) continue;
    out.sumw[static_cast<size_t>(bin)] += w[i];
    out.sumw2[static_cast<size_t>(bin)] += w[i] * w[i];
  }
  return out;
}

HistSummary build_count_hist_summary(const std::vector<float>& x, int nbins, double xmin, double xmax,
                                     bool fold_overflow) {
  HistSummary out = make_empty_hist(nbins);
  for (const float xi : x) {
    const int bin = find_bin(xi, nbins, xmin, xmax, fold_overflow);
    if (bin < 0) continue;
    out.sumw[static_cast<size_t>(bin)] += 1.0;
    out.sumw2[static_cast<size_t>(bin)] += 1.0;
  }
  return out;
}

HistSummary add_hist_summaries(const HistSummary& a, const HistSummary& b) {
  HistSummary out = a;
  const size_t n = std::min(a.sumw.size(), b.sumw.size());
  for (size_t i = 0; i < n; ++i) {
    out.sumw[i] += b.sumw[i];
    out.sumw2[i] += b.sumw2[i];
  }
  return out;
}

TH1D* make_hist1d(const HistSummary& h, const std::string& name, const std::string& title, int nbins, double xmin,
                  double xmax) {
  auto* out = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
  out->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) {
    out->SetBinContent(i + 1, h.sumw[static_cast<size_t>(i)]);
    out->SetBinError(i + 1, std::sqrt(std::max(0.0, h.sumw2[static_cast<size_t>(i)])));
  }
  return out;
}

TH1D* make_band_hist(const HistSummary& h, const TMatrixD& cov, const std::string& name, const std::string& title,
                     int nbins, double xmin, double xmax) {
  auto* out = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
  out->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) {
    out->SetBinContent(i + 1, h.sumw[static_cast<size_t>(i)]);
    out->SetBinError(i + 1, std::sqrt(std::max(0.0, cov(i, i))));
  }
  return out;
}

TH1D* make_fractional_hist(const HistSummary& h, const TMatrixD& cov, const std::string& name,
                           const std::string& title, int nbins, double xmin, double xmax) {
  auto* out = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
  out->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) {
    const double y = h.sumw[static_cast<size_t>(i)];
    const double e = std::sqrt(std::max(0.0, cov(i, i)));
    out->SetBinContent(i + 1, (y > 0.0) ? (e / y) : 0.0);
    out->SetBinError(i + 1, 0.0);
  }
  return out;
}

TH2D* make_matrix_hist(const TMatrixD& m, const std::string& name, const std::string& title, int nbins,
                       double xmin, double xmax) {
  auto* out = new TH2D(name.c_str(), title.c_str(), nbins, xmin, xmax, nbins, xmin, xmax);
  out->SetDirectory(nullptr);
  for (int ix = 0; ix < nbins; ++ix) {
    for (int iy = 0; iy < nbins; ++iy) {
      out->SetBinContent(ix + 1, iy + 1, m(ix, iy));
    }
  }
  return out;
}

TMatrixD make_correlation_matrix(const TMatrixD& cov) {
  const int n = cov.GetNrows();
  TMatrixD corr(n, n);
  corr.Zero();
  for (int i = 0; i < n; ++i) {
    const double sii = std::sqrt(std::max(0.0, cov(i, i)));
    for (int j = 0; j < n; ++j) {
      const double sjj = std::sqrt(std::max(0.0, cov(j, j)));
      const double denom = sii * sjj;
      corr(i, j) = (denom > 0.0) ? cov(i, j) / denom : 0.0;
    }
  }
  return corr;
}

TMatrixD make_diag_covariance(const HistSummary& h) {
  const int n = static_cast<int>(h.sumw.size());
  TMatrixD out(n, n);
  out.Zero();
  for (int i = 0; i < n; ++i) out(i, i) = h.sumw2[static_cast<size_t>(i)];
  return out;
}

TMatrixD make_fullcorr_covariance(const HistSummary& h, double frac) {
  const int n = static_cast<int>(h.sumw.size());
  TMatrixD out(n, n);
  out.Zero();
  const double f2 = frac * frac;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out(i, j) = h.sumw[static_cast<size_t>(i)] * h.sumw[static_cast<size_t>(j)] * f2;
    }
  }
  return out;
}

TMatrixD build_multisim_covariance_from_vectors(const std::vector<float>& x, const std::vector<double>& w_nom,
                                                const std::vector<UShortRVec>& w_univ_ushort,
                                                const std::vector<float>* w_cv_src, int nbins, double xmin,
                                                double xmax, bool fold_overflow, bool average_universes,
                                                std::size_t* n_universes_out) {
  const size_t n_evt = std::min(x.size(), w_nom.size());
  size_t n_universes = 0;
  const size_t n_vec = std::min(n_evt, w_univ_ushort.size());
  for (size_t i = 0; i < n_vec; ++i) {
    n_universes = std::max(n_universes, w_univ_ushort[i].size());
  }

  if (n_universes_out != nullptr) *n_universes_out = n_universes;

  TMatrixD cov(nbins, nbins);
  cov.Zero();
  if (n_universes == 0) return cov;

  std::vector<double> cv_bin(static_cast<size_t>(nbins), 0.0);
  std::vector<std::vector<double>> univ_bin(n_universes, std::vector<double>(static_cast<size_t>(nbins), 0.0));

  for (size_t i = 0; i < n_evt; ++i) {
    const int bin = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (bin < 0) continue;

    const double w0 = w_nom[i];
    cv_bin[static_cast<size_t>(bin)] += w0;

    static const UShortRVec empty;
    const auto& wu = (i < w_univ_ushort.size()) ? w_univ_ushort[i] : empty;
    const double cv_src = (w_cv_src != nullptr && i < w_cv_src->size()) ? static_cast<double>((*w_cv_src)[i]) : 1.0;

    for (size_t u = 0; u < n_universes; ++u) {
      double ratio = 1.0;
      if (u < wu.size()) ratio = 1.0e-3 * static_cast<double>(wu[u]);
      if (w_cv_src != nullptr && cv_src > 0.0) ratio /= cv_src;
      univ_bin[u][static_cast<size_t>(bin)] += w0 * ratio;
    }
  }

  for (size_t u = 0; u < n_universes; ++u) {
    for (int a = 0; a < nbins; ++a) {
      const double da = cv_bin[static_cast<size_t>(a)] - univ_bin[u][static_cast<size_t>(a)];
      for (int b = 0; b < nbins; ++b) {
        const double db = cv_bin[static_cast<size_t>(b)] - univ_bin[u][static_cast<size_t>(b)];
        cov(a, b) += da * db;
      }
    }
  }

  if (average_universes && n_universes > 0) cov *= (1.0 / static_cast<double>(n_universes));
  return cov;
}

TMatrixD build_unisim_covariance_from_vectors(const std::vector<float>& x, const std::vector<double>& w_nom,
                                              const std::vector<UShortRVec>& w_up_ushort,
                                              const std::vector<UShortRVec>& w_dn_ushort,
                                              int nbins, double xmin, double xmax, bool fold_overflow,
                                              std::size_t* n_knobs_out) {
  const size_t n_evt = std::min(x.size(), w_nom.size());
  size_t n_knobs = 0;
  const size_t n_up = std::min(n_evt, w_up_ushort.size());
  const size_t n_dn = std::min(n_evt, w_dn_ushort.size());
  for (size_t i = 0; i < n_up; ++i) n_knobs = std::max(n_knobs, w_up_ushort[i].size());
  for (size_t i = 0; i < n_dn; ++i) n_knobs = std::max(n_knobs, w_dn_ushort[i].size());

  if (n_knobs_out != nullptr) *n_knobs_out = n_knobs;

  TMatrixD cov(nbins, nbins);
  cov.Zero();
  if (n_knobs == 0) return cov;

  std::vector<std::vector<double>> delta(n_knobs, std::vector<double>(static_cast<size_t>(nbins), 0.0));

  for (size_t i = 0; i < n_evt; ++i) {
    const int bin = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (bin < 0) continue;

    static const UShortRVec empty;
    const auto& wup = (i < w_up_ushort.size()) ? w_up_ushort[i] : empty;
    const auto& wdn = (i < w_dn_ushort.size()) ? w_dn_ushort[i] : empty;

    for (size_t k = 0; k < n_knobs; ++k) {
      double rup = 1.0;
      double rdn = 1.0;
      if (k < wup.size()) rup = 1.0e-3 * static_cast<double>(wup[k]);
      if (k < wdn.size()) rdn = 1.0e-3 * static_cast<double>(wdn[k]);
      const double d = 0.5 * w_nom[i] * (rup - rdn);
      delta[k][static_cast<size_t>(bin)] += d;
    }
  }

  for (size_t k = 0; k < n_knobs; ++k) {
    for (int a = 0; a < nbins; ++a) {
      for (int b = 0; b < nbins; ++b) {
        cov(a, b) += delta[k][static_cast<size_t>(a)] * delta[k][static_cast<size_t>(b)];
      }
    }
  }

  return cov;
}

TMatrixD build_multisim_covariance(ROOT::RDF::RNode node, const std::string& score_col, const std::string& weight_col,
                                   const std::string& univ_branch, const std::string& cv_branch, bool divide_by_cv,
                                   int nbins, double xmin, double xmax, bool fold_overflow, bool average_universes,
                                   std::size_t* n_universes_out) {
  auto x_h = node.Take<float>(score_col);
  auto w_h = node.Take<double>(weight_col);
  auto u_h = node.Take<UShortRVec>(univ_branch);

  if (divide_by_cv && !cv_branch.empty()) {
    auto cv_h = node.Take<float>(cv_branch);
    return build_multisim_covariance_from_vectors(*x_h, *w_h, *u_h, &(*cv_h), nbins, xmin, xmax, fold_overflow,
                                                  average_universes, n_universes_out);
  }

  return build_multisim_covariance_from_vectors(*x_h, *w_h, *u_h, nullptr, nbins, xmin, xmax, fold_overflow,
                                                average_universes, n_universes_out);
}

TMatrixD build_unisim_covariance(ROOT::RDF::RNode node, const std::string& score_col, const std::string& weight_col,
                                 const std::string& up_branch, const std::string& dn_branch, int nbins, double xmin,
                                 double xmax, bool fold_overflow, std::size_t* n_knobs_out) {
  auto x_h = node.Take<float>(score_col);
  auto w_h = node.Take<double>(weight_col);
  auto up_h = node.Take<UShortRVec>(up_branch);
  auto dn_h = node.Take<UShortRVec>(dn_branch);

  return build_unisim_covariance_from_vectors(*x_h, *w_h, *up_h, *dn_h, nbins, xmin, xmax, fold_overflow,
                                              n_knobs_out);
}

int colour_for_index(size_t i) {
  static const int colours[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1};
  return colours[i % (sizeof(colours) / sizeof(colours[0]))];
}

void style_fractional_hist(TH1D* h, size_t i) {
  const int c = colour_for_index(i);
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetLineWidth(3);
  h->SetMarkerStyle(20 + static_cast<int>(i % 5));
}

void draw_prediction_breakdown(TH1D* h_total, TH1D* h_band, TH1D* h_sig, TH1D* h_bkg, TH1D* h_ext, TH1D* h_data,
                               const std::string& x_title, const std::string& out_path) {
  if (h_total == nullptr || h_band == nullptr) return;

  gStyle->SetOptStat(0);
  TCanvas c("c_score_prediction_split", "c_score_prediction_split", 980, 760);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetBottomMargin(0.12);

  h_total->SetTitle((";" + x_title + ";Events").c_str());
  h_total->SetLineColor(kBlack);
  h_total->SetLineWidth(3);

  h_band->SetFillColorAlpha(kAzure + 1, 0.35);
  h_band->SetLineColor(kAzure + 2);
  h_band->SetMarkerSize(0.0);

  if (h_sig != nullptr) {
    h_sig->SetLineColor(kRed + 1);
    h_sig->SetLineWidth(3);
  }
  if (h_bkg != nullptr) {
    h_bkg->SetLineColor(kBlue + 1);
    h_bkg->SetLineWidth(3);
  }
  if (h_ext != nullptr) {
    h_ext->SetLineColor(kGreen + 2);
    h_ext->SetLineWidth(3);
  }

  double ymax = std::max(h_total->GetMaximum(), h_band->GetMaximum());
  if (h_sig != nullptr) ymax = std::max(ymax, h_sig->GetMaximum());
  if (h_bkg != nullptr) ymax = std::max(ymax, h_bkg->GetMaximum());
  if (h_ext != nullptr) ymax = std::max(ymax, h_ext->GetMaximum());
  if (h_data != nullptr) ymax = std::max(ymax, h_data->GetMaximum());
  if (ymax <= 0.0) ymax = 1.0;
  h_total->GetYaxis()->SetRangeUser(0.0, 1.35 * ymax);

  h_total->Draw("HIST");
  h_band->Draw("E2 SAME");
  h_total->Draw("HIST SAME");
  if (h_sig != nullptr) h_sig->Draw("HIST SAME");
  if (h_bkg != nullptr) h_bkg->Draw("HIST SAME");
  if (h_ext != nullptr) h_ext->Draw("HIST SAME");

  TLegend leg(0.54, 0.64, 0.90, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(h_total, "Total prediction", "l");
  leg.AddEntry(h_band, "Total uncertainty", "f");
  if (h_sig != nullptr) leg.AddEntry(h_sig, "Truth signal", "l");
  if (h_bkg != nullptr) leg.AddEntry(h_bkg, "Truth background", "l");
  if (h_ext != nullptr) leg.AddEntry(h_ext, "EXT", "l");

  if (h_data != nullptr) {
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetLineWidth(2);
    h_data->Draw("E1 SAME");
    leg.AddEntry(h_data, "Data", "ep");
  }

  leg.Draw();
  c.RedrawAxis();
  c.SaveAs(out_path.c_str());
}

void draw_fractional_uncertainties(const std::vector<TH1D*>& hists, const std::vector<std::string>& labels,
                                   const std::string& x_title, const std::string& out_path) {
  if (hists.empty()) return;

  gStyle->SetOptStat(0);
  TCanvas c("c_score_frac", "c_score_frac", 950, 760);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetBottomMargin(0.12);

  double ymax = 0.0;
  for (auto* h : hists) {
    if (h == nullptr) continue;
    ymax = std::max(ymax, h->GetMaximum());
  }
  if (ymax <= 0.0) ymax = 0.05;

  bool first = true;
  TLegend leg(0.56, 0.56, 0.90, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);

  for (size_t i = 0; i < hists.size(); ++i) {
    auto* h = hists[i];
    if (h == nullptr) continue;
    style_fractional_hist(h, i);
    h->SetTitle((";" + x_title + ";Fractional uncertainty").c_str());
    h->GetYaxis()->SetRangeUser(0.0, 1.25 * ymax);
    h->Draw(first ? "HIST" : "HIST SAME");
    if (i < labels.size()) leg.AddEntry(h, labels[i].c_str(), "l");
    first = false;
  }

  leg.Draw();
  c.RedrawAxis();
  c.SaveAs(out_path.c_str());
}

void draw_matrix(TH2D* h, const std::string& out_path, bool correlation) {
  if (h == nullptr) return;

  gStyle->SetOptStat(0);
  TCanvas c((std::string("c_") + h->GetName()).c_str(), h->GetName(), 900, 800);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.16);
  c.SetBottomMargin(0.12);

  if (correlation) {
    h->SetMinimum(-1.0);
    h->SetMaximum(1.0);
  }

  h->Draw("COLZ");
  c.SaveAs(out_path.c_str());
}

void write_split_component_matrices(TFile& fout, const SplitComponent& comp, int idx, int nbins, double xmin,
                                    double xmax) {
  TH2D* htot = make_matrix_hist(comp.total, "h_cov_component_total_" + std::to_string(idx), comp.name + " total",
                                nbins, xmin, xmax);
  TH2D* hsig = make_matrix_hist(comp.signal, "h_cov_component_signal_" + std::to_string(idx), comp.name + " signal",
                                nbins, xmin, xmax);
  TH2D* hbkg = make_matrix_hist(comp.background, "h_cov_component_background_" + std::to_string(idx),
                                comp.name + " background", nbins, xmin, xmax);

  TMatrixD mixed = comp.total;
  mixed -= comp.signal;
  mixed -= comp.background;
  TH2D* hmix = make_matrix_hist(mixed, "h_cov_component_mixed_" + std::to_string(idx), comp.name + " mixed",
                                nbins, xmin, xmax);

  fout.cd();
  htot->Write();
  hsig->Write();
  hbkg->Write();
  hmix->Write();

  delete htot;
  delete hsig;
  delete hbkg;
  delete hmix;
}

}  // namespace

/*
  This macro propagates event-weight universes into a covariance matrix for the
  first inference score (or another chosen score index), with an explicit truth
  signal/background split driven by `is_signal` (or another boolean expression).

  Observable:
      O_r^u = S_r^u + B_r^u + E_r
  where S_r^u and B_r^u are the truth-signal and truth-background MC yields in
  reco score bin r for universe u, and E_r is any EXT contribution.

  For each enabled multisim source the covariance is
      C_ab = (1/N_u) sum_u (O_a^CV - O_a^u)(O_b^CV - O_b^u)
  if `average_universes=true`, or the raw sum otherwise.

  Because the signal/background partition is built at the event level, the macro
  can also output:
      C^sig, C^bkg, and C^mix = C^MC_total - C^sig - C^bkg
  where C^mix contains the signal-background cross terms for shared universes.

  Important limitation:
      This is still a reco-template covariance. It is not the MCC9 acceptance-only
      signal propagation
          expected_signal_{t,r}^u = S_{t,r}^u N_t^CV
      because that requires an explicit true-bin index t and the true-to-reco
      migration counts N_{t,r}^u, which are not constructed here.
*/
int plot_inference_score_systematics(const std::string& samples_tsv = "",
                                           const char* extra_sel = "sel_muon",
                                           const std::string& nominal_weight = "w_nominal",
                                           int score_index = 0, int nbins = 24, double xmin = -15.0,
                                           double xmax = 15.0, bool fold_overflow = true,
                                           bool average_universes = true, bool include_data = false,
                                           bool add_mcstat = true, bool add_extstat = false,
                                           double mc_fullcorr_frac = 0.0,
                                           double category_fullcorr_frac = 0.0,
                                           const std::string& signal_sel = "is_signal",
                                           const std::string& category_sel = "is_signal",
                                           bool use_genie = true, bool use_reint = true,
                                           bool use_ppfx = true, bool ppfx_is_already_in_nominal = true,
                                           bool use_flux = false, bool use_genie_unisims = false,
                                           const std::string& output_stem = "rareproc_score0_systematics_split") {
  const std::string extra_sel_expr = (extra_sel != nullptr) ? extra_sel : "";
  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plot_inference_score_systematics_split] input=" << list_path << "\n";

  if (!looks_like_event_list_root(list_path)) {
    std::cerr << "[plot_inference_score_systematics_split] input is not an event list ROOT file: " << list_path
              << "\n";
    return 2;
  }

  if (nbins < 1) nbins = 1;
  if (xmax < xmin) std::swap(xmin, xmax);
  if (xmax == xmin) xmax = xmin + 1.0;

  EventListIO el(list_path);
  ROOT::RDataFrame rdf0 = el.rdf();

  auto mask_ext = el.mask_for_ext();
  auto mask_mc = el.mask_for_mc_like();
  auto mask_data = el.mask_for_data();

  auto filter_by_mask = [](ROOT::RDF::RNode node, std::shared_ptr<const std::vector<char>> mask) {
    return node.Filter(
        [mask](int sid) {
          return sid >= 0 && sid < static_cast<int>(mask->size()) && (*mask)[static_cast<std::size_t>(sid)];
        },
        {"sample_id"});
  };

  ROOT::RDF::RNode base = SelectionService::decorate(rdf0).Define(
      "score0",
      [score_index](const ROOT::RVec<float>& scores) {
        if (score_index < 0 || score_index >= static_cast<int>(scores.size())) return -9999.0f;
        return scores[static_cast<size_t>(score_index)];
      },
      {"inf_scores"})
                              .Define("__score_cut_bin__", [] { return 0.5f; });

  ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext).Define("__w__", nominal_weight);
  ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                 .Filter(
                                     [mask_ext](int sid) {
                                       return !(sid >= 0 && sid < static_cast<int>(mask_ext->size()) &&
                                                (*mask_ext)[static_cast<std::size_t>(sid)]);
                                     },
                                     {"sample_id"})
                                 .Define("__w__", nominal_weight);
  ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

  node_mc = apply_optional_filter(node_mc, extra_sel_expr, "extra selection");
  node_ext = apply_optional_filter(node_ext, extra_sel_expr, "extra selection");
  node_data = apply_optional_filter(node_data, extra_sel_expr, "extra selection");

  constexpr double score_cut = 6.9;
  auto pass_score_cut = [](float s) { return s > score_cut; };

  node_mc = node_mc.Filter(pass_score_cut, {"score0"});
  node_ext = node_ext.Filter(pass_score_cut, {"score0"});
  node_data = node_data.Filter(pass_score_cut, {"score0"});

  nbins = 1;
  xmin = 0.0;
  xmax = 1.0;
  fold_overflow = false;

  const std::string hist_col = "__score_cut_bin__";

  if (signal_sel.empty()) {
    std::cerr << "[plot_inference_score_systematics_split] signal_sel must not be empty.\n";
    return 1;
  }

  if (is_simple_identifier(signal_sel) && !has_column(node_mc, signal_sel)) {
    std::cerr << "[plot_inference_score_systematics_split] missing signal selection column '" << signal_sel
              << "'.\n";
    return 1;
  }

  ROOT::RDF::RNode node_sig = apply_optional_filter(node_mc, signal_sel, "signal selection");
  ROOT::RDF::RNode node_bkg = apply_negated_filter(node_mc, signal_sel, "signal selection");

  auto x_mc_h = node_mc.Take<float>(hist_col);
  auto w_mc_h = node_mc.Take<double>("__w__");
  auto x_sig_h = node_sig.Take<float>(hist_col);
  auto w_sig_h = node_sig.Take<double>("__w__");
  auto x_bkg_h = node_bkg.Take<float>(hist_col);
  auto w_bkg_h = node_bkg.Take<double>("__w__");

  const HistSummary h_mc = build_hist_summary(*x_mc_h, *w_mc_h, nbins, xmin, xmax, fold_overflow);
  const HistSummary h_sig = build_hist_summary(*x_sig_h, *w_sig_h, nbins, xmin, xmax, fold_overflow);
  const HistSummary h_bkg = build_hist_summary(*x_bkg_h, *w_bkg_h, nbins, xmin, xmax, fold_overflow);

  HistSummary h_ext = make_empty_hist(nbins);
  auto x_ext_h = node_ext.Take<float>(hist_col);
  auto w_ext_h = node_ext.Take<double>("__w__");
  h_ext = build_hist_summary(*x_ext_h, *w_ext_h, nbins, xmin, xmax, fold_overflow);

  const HistSummary h_pred = add_hist_summaries(h_mc, h_ext);

  std::cout << "[plot_inference_score_systematics_split] total MC yield = "
            << std::accumulate(h_mc.sumw.begin(), h_mc.sumw.end(), 0.0) << "\n";
  std::cout << "[plot_inference_score_systematics_split] signal MC yield = "
            << std::accumulate(h_sig.sumw.begin(), h_sig.sumw.end(), 0.0) << "\n";
  std::cout << "[plot_inference_score_systematics_split] background MC yield = "
            << std::accumulate(h_bkg.sumw.begin(), h_bkg.sumw.end(), 0.0) << "\n";
  std::cout << "[plot_inference_score_systematics_split] EXT yield = "
            << std::accumulate(h_ext.sumw.begin(), h_ext.sumw.end(), 0.0) << "\n";

  std::vector<SplitComponent> components;

  auto add_multisim_source = [&](const std::string& label, const std::string& weight_branch,
                                 const std::string& cv_branch, bool divide_by_cv) {
    if (!has_column(node_mc, weight_branch)) return;

    std::size_t nuniv_total = 0;
    std::size_t nuniv_sig = 0;
    std::size_t nuniv_bkg = 0;

    SplitComponent comp{label, make_zero_matrix(nbins), make_zero_matrix(nbins), make_zero_matrix(nbins), 0, true};
    comp.total = build_multisim_covariance(node_mc, hist_col, "__w__", weight_branch, cv_branch, divide_by_cv,
                                           nbins, xmin, xmax, fold_overflow, average_universes, &nuniv_total);
    comp.signal = build_multisim_covariance(node_sig, hist_col, "__w__", weight_branch, cv_branch, divide_by_cv,
                                            nbins, xmin, xmax, fold_overflow, average_universes, &nuniv_sig);
    comp.background = build_multisim_covariance(node_bkg, hist_col, "__w__", weight_branch, cv_branch, divide_by_cv,
                                                nbins, xmin, xmax, fold_overflow, average_universes, &nuniv_bkg);
    comp.n_universes = std::max(nuniv_total, std::max(nuniv_sig, nuniv_bkg));
    std::cout << "[plot_inference_score_systematics_split] " << label << " universes: " << comp.n_universes
              << (divide_by_cv && !cv_branch.empty() ? " (CV ratio mode)" : "") << "\n";
    components.push_back(std::move(comp));
  };

  auto add_unisim_source = [&](const std::string& label, const std::string& up_branch,
                               const std::string& dn_branch) {
    if (!has_column(node_mc, up_branch) || !has_column(node_mc, dn_branch)) return;

    std::size_t nknob_total = 0;
    std::size_t nknob_sig = 0;
    std::size_t nknob_bkg = 0;

    SplitComponent comp{label, make_zero_matrix(nbins), make_zero_matrix(nbins), make_zero_matrix(nbins), 0, true};
    comp.total = build_unisim_covariance(node_mc, hist_col, "__w__", up_branch, dn_branch, nbins, xmin, xmax,
                                         fold_overflow, &nknob_total);
    comp.signal = build_unisim_covariance(node_sig, hist_col, "__w__", up_branch, dn_branch, nbins, xmin, xmax,
                                          fold_overflow, &nknob_sig);
    comp.background = build_unisim_covariance(node_bkg, hist_col, "__w__", up_branch, dn_branch, nbins, xmin, xmax,
                                              fold_overflow, &nknob_bkg);
    comp.n_universes = std::max(nknob_total, std::max(nknob_sig, nknob_bkg));
    std::cout << "[plot_inference_score_systematics_split] " << label << " knobs: " << comp.n_universes << "\n";
    components.push_back(std::move(comp));
  };

  if (use_genie) add_multisim_source("GENIE multisim", "weightsGenie", "", false);
  if (use_reint) add_multisim_source("Reinteraction multisim", "weightsReint", "", false);
  if (use_ppfx) {
    const bool can_divide = ppfx_is_already_in_nominal && has_column(node_mc, "ppfx_cv");
    add_multisim_source("PPFX multisim", "weightsPPFX", can_divide ? "ppfx_cv" : "", can_divide);
  }
  if (use_flux) add_multisim_source("Flux multisim", "weightsFlux", "", false);
  if (use_genie_unisims) add_unisim_source("GENIE unisim", "weightsGenieUp", "weightsGenieDn");

  if (add_mcstat) {
    SplitComponent comp{"MC stat", make_diag_covariance(h_mc), make_diag_covariance(h_sig), make_diag_covariance(h_bkg),
                        0, true};
    components.push_back(std::move(comp));
  }

  if (add_extstat) {
    SplitComponent comp{"EXT stat", make_diag_covariance(h_ext), make_zero_matrix(nbins), make_zero_matrix(nbins), 0,
                        false};
    components.push_back(std::move(comp));
  }

  if (mc_fullcorr_frac > 0.0) {
    SplitComponent comp{"MC full corr", make_fullcorr_covariance(h_mc, mc_fullcorr_frac),
                        make_fullcorr_covariance(h_sig, mc_fullcorr_frac),
                        make_fullcorr_covariance(h_bkg, mc_fullcorr_frac), 0, true};
    components.push_back(std::move(comp));
  }

  if (category_fullcorr_frac > 0.0) {
    ROOT::RDF::RNode node_cat_total = apply_optional_filter(node_mc, category_sel, "category selection");
    ROOT::RDF::RNode node_cat_sig = apply_optional_filter(node_sig, category_sel, "category selection");
    ROOT::RDF::RNode node_cat_bkg = apply_optional_filter(node_bkg, category_sel, "category selection");

    auto x_cat_tot_h = node_cat_total.Take<float>(hist_col);
    auto w_cat_tot_h = node_cat_total.Take<double>("__w__");
    auto x_cat_sig_h = node_cat_sig.Take<float>(hist_col);
    auto w_cat_sig_h = node_cat_sig.Take<double>("__w__");
    auto x_cat_bkg_h = node_cat_bkg.Take<float>(hist_col);
    auto w_cat_bkg_h = node_cat_bkg.Take<double>("__w__");

    const HistSummary h_cat_total = build_hist_summary(*x_cat_tot_h, *w_cat_tot_h, nbins, xmin, xmax, fold_overflow);
    const HistSummary h_cat_sig = build_hist_summary(*x_cat_sig_h, *w_cat_sig_h, nbins, xmin, xmax, fold_overflow);
    const HistSummary h_cat_bkg = build_hist_summary(*x_cat_bkg_h, *w_cat_bkg_h, nbins, xmin, xmax, fold_overflow);

    SplitComponent comp{"Category full corr", make_fullcorr_covariance(h_cat_total, category_fullcorr_frac),
                        make_fullcorr_covariance(h_cat_sig, category_fullcorr_frac),
                        make_fullcorr_covariance(h_cat_bkg, category_fullcorr_frac), 0, true};
    components.push_back(std::move(comp));
  }

  TMatrixD cov_total = make_zero_matrix(nbins);
  TMatrixD cov_mc_total = make_zero_matrix(nbins);
  TMatrixD cov_sig = make_zero_matrix(nbins);
  TMatrixD cov_bkg = make_zero_matrix(nbins);

  for (const auto& comp : components) {
    cov_total += comp.total;
    if (comp.is_mc_source) cov_mc_total += comp.total;
    cov_sig += comp.signal;
    cov_bkg += comp.background;
  }

  TMatrixD cov_mix = cov_mc_total;
  cov_mix -= cov_sig;
  cov_mix -= cov_bkg;

  const TMatrixD corr_total = make_correlation_matrix(cov_total);
  const TMatrixD corr_sig = make_correlation_matrix(cov_sig);
  const TMatrixD corr_bkg = make_correlation_matrix(cov_bkg);

  TH1D* h_total_hist = make_hist1d(h_pred, "h_score_cv_total", "Total prediction", nbins, xmin, xmax);
  TH1D* h_sig_hist = make_hist1d(h_sig, "h_score_cv_signal", "Truth signal", nbins, xmin, xmax);
  TH1D* h_bkg_hist = make_hist1d(h_bkg, "h_score_cv_background", "Truth background", nbins, xmin, xmax);
  TH1D* h_ext_hist = make_hist1d(h_ext, "h_score_cv_ext", "EXT", nbins, xmin, xmax);
  TH1D* h_band_hist = make_band_hist(h_pred, cov_total, "h_score_total_band", "Total uncertainty band", nbins, xmin,
                                     xmax);

  TH1D* h_frac_total = make_fractional_hist(h_pred, cov_total, "h_score_frac_total", "Total fractional uncertainty",
                                            nbins, xmin, xmax);
  TH1D* h_frac_sig = make_fractional_hist(h_sig, cov_sig, "h_score_frac_signal",
                                          "Signal fractional uncertainty", nbins, xmin, xmax);
  TH1D* h_frac_bkg = make_fractional_hist(h_bkg, cov_bkg, "h_score_frac_background",
                                          "Background fractional uncertainty", nbins, xmin, xmax);

  TH2D* h_cov_total = make_matrix_hist(cov_total, "h_score_cov_total", ";Inference score;Inference score", nbins,
                                       xmin, xmax);
  TH2D* h_corr_total = make_matrix_hist(corr_total, "h_score_corr_total", ";Inference score;Inference score", nbins,
                                        xmin, xmax);

  TH1D* h_data_hist = nullptr;
  if (include_data) {
    auto x_data_h = node_data.Take<float>(hist_col);
    const HistSummary h_data = build_count_hist_summary(*x_data_h, nbins, xmin, xmax, fold_overflow);
    h_data_hist = make_hist1d(h_data, "h_score_data", "Data", nbins, xmin, xmax);
  }

  const std::string cut_label = "score[0] > 6.9";
  for (TH1D* h : {h_total_hist, h_sig_hist, h_bkg_hist, h_ext_hist, h_band_hist, h_frac_total, h_frac_sig, h_frac_bkg,
                  h_data_hist}) {
    if (h != nullptr) h->GetXaxis()->SetBinLabel(1, cut_label.c_str());
  }

  const double y = h_pred.sumw[0];
  const double sigma = std::sqrt(std::max(0.0, cov_total(0, 0)));
  const double frac = (y > 0.0) ? sigma / y : 0.0;
  std::cout << "[plot_inference_score_systematics_split] selected yield = " << y << ", abs unc = " << sigma
            << ", frac unc = " << frac << "\n";

  const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
  const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
  gSystem->mkdir(out_dir.c_str(), true);

  const std::string stem = output_stem.empty() ? "rareproc_score0_systematics_split" : output_stem;
  const std::string score_title = cut_label;
  const std::string pred_path = out_dir + "/" + stem + "_prediction." + out_fmt;
  const std::string frac_path = out_dir + "/" + stem + "_fractional_by_source." + out_fmt;
  const std::string cov_path = out_dir + "/" + stem + "_covariance." + out_fmt;
  const std::string corr_path = out_dir + "/" + stem + "_correlation." + out_fmt;
  const std::string root_path = out_dir + "/" + stem + ".root";

  draw_prediction_breakdown(h_total_hist, h_band_hist, h_sig_hist, h_bkg_hist, h_ext_hist, h_data_hist, score_title,
                            pred_path);
  draw_fractional_uncertainties({h_frac_total}, {"Total"}, score_title, frac_path);
  draw_matrix(h_cov_total, cov_path, false);
  draw_matrix(h_corr_total, corr_path, true);

  TFile fout(root_path.c_str(), "RECREATE");
  h_total_hist->Write();
  h_sig_hist->Write();
  h_bkg_hist->Write();
  h_ext_hist->Write();
  h_band_hist->Write();
  h_frac_total->Write();
  h_frac_sig->Write();
  h_frac_bkg->Write();
  h_cov_total->Write();
  h_corr_total->Write();
  if (h_data_hist != nullptr) h_data_hist->Write();

  TH2D* h_cov_sig = make_matrix_hist(cov_sig, "h_score_cov_signal", ";Inference score;Inference score", nbins, xmin,
                                     xmax);
  TH2D* h_cov_bkg = make_matrix_hist(cov_bkg, "h_score_cov_background", ";Inference score;Inference score", nbins,
                                     xmin, xmax);
  TH2D* h_cov_mix = make_matrix_hist(cov_mix, "h_score_cov_mixed", ";Inference score;Inference score", nbins, xmin,
                                     xmax);
  TH2D* h_corr_sig = make_matrix_hist(corr_sig, "h_score_corr_signal", ";Inference score;Inference score", nbins,
                                      xmin, xmax);
  TH2D* h_corr_bkg = make_matrix_hist(corr_bkg, "h_score_corr_background", ";Inference score;Inference score", nbins,
                                      xmin, xmax);

  h_cov_sig->Write();
  h_cov_bkg->Write();
  h_cov_mix->Write();
  h_corr_sig->Write();
  h_corr_bkg->Write();

  for (size_t i = 0; i < components.size(); ++i) {
    write_split_component_matrices(fout, components[i], static_cast<int>(i), nbins, xmin, xmax);
  }
  fout.Close();

  std::cout << "[plot_inference_score_systematics_split] wrote " << pred_path << "\n";
  std::cout << "[plot_inference_score_systematics_split] wrote " << frac_path << "\n";
  std::cout << "[plot_inference_score_systematics_split] wrote " << cov_path << "\n";
  std::cout << "[plot_inference_score_systematics_split] wrote " << corr_path << "\n";
  std::cout << "[plot_inference_score_systematics_split] wrote " << root_path << "\n";
  std::cout << "[plot_inference_score_systematics_split] done\n";

  delete h_total_hist;
  delete h_sig_hist;
  delete h_bkg_hist;
  delete h_ext_hist;
  delete h_band_hist;
  delete h_frac_total;
  delete h_frac_sig;
  delete h_frac_bkg;
  delete h_cov_total;
  delete h_corr_total;
  delete h_cov_sig;
  delete h_cov_bkg;
  delete h_cov_mix;
  delete h_corr_sig;
  delete h_corr_bkg;
  delete h_data_hist;

  return 0;
}


int plot_inference_score_systematics_split(const std::string& samples_tsv = "",
                                           const char* extra_sel = "sel_muon",
                                           const std::string& nominal_weight = "w_nominal",
                                           int score_index = 0, int nbins = 24, double xmin = -15.0,
                                           double xmax = 15.0, bool fold_overflow = true,
                                           bool average_universes = true, bool include_data = false,
                                           bool add_mcstat = true, bool add_extstat = false,
                                           double mc_fullcorr_frac = 0.0,
                                           double category_fullcorr_frac = 0.0,
                                           const std::string& signal_sel = "is_signal",
                                           const std::string& category_sel = "is_signal",
                                           bool use_genie = true, bool use_reint = true,
                                           bool use_ppfx = true, bool ppfx_is_already_in_nominal = true,
                                           bool use_flux = false, bool use_genie_unisims = false,
                                           const std::string& output_stem = "rareproc_score0_systematics_split") {
  return plot_inference_score_systematics(samples_tsv, extra_sel, nominal_weight, score_index, nbins, xmin, xmax,
                                          fold_overflow, average_universes, include_data, add_mcstat, add_extstat,
                                          mc_fullcorr_frac, category_fullcorr_frac, signal_sel, category_sel,
                                          use_genie, use_reint, use_ppfx, ppfx_is_already_in_nominal, use_flux,
                                          use_genie_unisims, output_stem);
}
