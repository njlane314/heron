#ifndef PLOT_INFERENCE_SCORE_UNIVERSES_CXX
#define PLOT_INFERENCE_SCORE_UNIVERSES_CXX

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TLegend.h>
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
#include <vector>

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"

using namespace nu;

namespace {

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
    std::cerr << "[plot_inference_score_universes] " << label << " column '" << expr
              << "' is missing; skipping that filter.\n";
    return node;
  }

  return node.Filter(expr);
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

std::vector<double> build_cv_hist(const std::vector<float>& x,
                                  const std::vector<double>& w_nom,
                                  int nbins, double xmin, double xmax,
                                  bool fold_overflow) {
  std::vector<double> out(static_cast<size_t>(nbins), 0.0);
  const size_t n = std::min(x.size(), w_nom.size());
  for (size_t i = 0; i < n; ++i) {
    const int b = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (b < 0) continue;
    out[static_cast<size_t>(b)] += w_nom[i];
  }
  return out;
}

std::vector<double> divide_weights_by_cv_source(const std::vector<double>& w_nom,
                                                const std::vector<float>& w_cv_src) {
  std::vector<double> out = w_nom;
  const size_t n = std::min(out.size(), w_cv_src.size());
  for (size_t i = 0; i < n; ++i) {
    const double cv = static_cast<double>(w_cv_src[i]);
    if (std::isfinite(cv) && cv > 0.0) {
      out[i] /= cv;
    }
  }
  return out;
}

std::vector<std::vector<double>> build_universe_hists_from_multisim(
    const std::vector<float>& x,
    const std::vector<double>& w_nom,
    const std::vector<ROOT::RVec<unsigned short>>& w_univ_ushort,
    const std::vector<float>* w_cv_src,
    int nbins, double xmin, double xmax,
    bool fold_overflow) {
  const size_t n_evt = std::min(x.size(), w_nom.size());

  size_t n_universes = 0;
  const size_t n_vec = std::min(n_evt, w_univ_ushort.size());
  for (size_t i = 0; i < n_vec; ++i) n_universes = std::max(n_universes, w_univ_ushort[i].size());

  std::vector<std::vector<double>> out(n_universes, std::vector<double>(static_cast<size_t>(nbins), 0.0));
  if (n_universes == 0) return out;

  for (size_t i = 0; i < n_evt; ++i) {
    const int b = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (b < 0) continue;

    static const ROOT::RVec<unsigned short> empty;
    const auto& wu = (i < w_univ_ushort.size()) ? w_univ_ushort[i] : empty;
    const double cv_src = (w_cv_src != nullptr && i < w_cv_src->size()) ? static_cast<double>((*w_cv_src)[i]) : 1.0;

    for (size_t u = 0; u < n_universes; ++u) {
      double ratio = 1.0;
      if (u < wu.size()) ratio = 1.0e-3 * static_cast<double>(wu[u]);
      if (w_cv_src != nullptr && cv_src > 0.0) ratio /= cv_src;
      out[u][static_cast<size_t>(b)] += w_nom[i] * ratio;
    }
  }

  return out;
}

std::vector<std::vector<double>> build_universe_hists_from_unisim(
    const std::vector<float>& x,
    const std::vector<double>& w_nom,
    const std::vector<ROOT::RVec<unsigned short>>& w_up_ushort,
    const std::vector<ROOT::RVec<unsigned short>>& w_dn_ushort,
    int nbins, double xmin, double xmax,
    bool fold_overflow,
    bool use_up_variations) {
  const size_t n_evt = std::min(x.size(), w_nom.size());

  size_t n_knobs = 0;
  const size_t n_up = std::min(n_evt, w_up_ushort.size());
  const size_t n_dn = std::min(n_evt, w_dn_ushort.size());
  for (size_t i = 0; i < n_up; ++i) n_knobs = std::max(n_knobs, w_up_ushort[i].size());
  for (size_t i = 0; i < n_dn; ++i) n_knobs = std::max(n_knobs, w_dn_ushort[i].size());

  std::vector<std::vector<double>> out(n_knobs, std::vector<double>(static_cast<size_t>(nbins), 0.0));
  if (n_knobs == 0) return out;

  for (size_t i = 0; i < n_evt; ++i) {
    const int b = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (b < 0) continue;

    static const ROOT::RVec<unsigned short> empty;
    const auto& wup = (i < w_up_ushort.size()) ? w_up_ushort[i] : empty;
    const auto& wdn = (i < w_dn_ushort.size()) ? w_dn_ushort[i] : empty;

    for (size_t k = 0; k < n_knobs; ++k) {
      double ratio = 1.0;
      if (use_up_variations) {
        if (k < wup.size()) ratio = 1.0e-3 * static_cast<double>(wup[k]);
      } else {
        if (k < wdn.size()) ratio = 1.0e-3 * static_cast<double>(wdn[k]);
      }
      out[k][static_cast<size_t>(b)] += w_nom[i] * ratio;
    }
  }

  return out;
}

struct HistBand {
  std::vector<double> lo;
  std::vector<double> hi;
};

double quantile_from_sorted(const std::vector<double>& sorted, double q) {
  if (sorted.empty()) return 0.0;
  if (q <= 0.0) return sorted.front();
  if (q >= 1.0) return sorted.back();

  const double pos = q * static_cast<double>(sorted.size() - 1);
  const size_t i = static_cast<size_t>(std::floor(pos));
  const double frac = pos - static_cast<double>(i);
  if (i + 1 >= sorted.size()) return sorted.back();
  return (1.0 - frac) * sorted[i] + frac * sorted[i + 1];
}

HistBand build_quantile_band(const std::vector<std::vector<double>>& univ_hists,
                             int nbins, double q_lo, double q_hi) {
  HistBand band;
  band.lo.assign(static_cast<size_t>(nbins), 0.0);
  band.hi.assign(static_cast<size_t>(nbins), 0.0);

  std::vector<double> vals;
  vals.reserve(univ_hists.size());

  for (int b = 0; b < nbins; ++b) {
    vals.clear();
    for (const auto& hu : univ_hists) {
      if (static_cast<size_t>(b) >= hu.size()) continue;
      const double v = hu[static_cast<size_t>(b)];
      if (std::isfinite(v)) vals.push_back(v);
    }
    if (vals.empty()) continue;

    std::sort(vals.begin(), vals.end());
    band.lo[static_cast<size_t>(b)] = quantile_from_sorted(vals, q_lo);
    band.hi[static_cast<size_t>(b)] = quantile_from_sorted(vals, q_hi);
  }

  return band;
}

HistBand build_central_band(const std::vector<std::vector<double>>& univ_hists,
                            int nbins, double coverage) {
  const double cov = std::max(0.0, std::min(1.0, coverage));
  const double tail = 0.5 * (1.0 - cov);
  return build_quantile_band(univ_hists, nbins, tail, 1.0 - tail);
}

double max_band_value(const HistBand& band) {
  if (band.hi.empty()) return 0.0;
  return *std::max_element(band.hi.begin(), band.hi.end());
}

TH1D* make_hist1d(const std::vector<double>& bins, const std::string& name,
                  int nbins, double xmin, double xmax) {
  auto* h = new TH1D(name.c_str(), "", nbins, xmin, xmax);
  h->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) h->SetBinContent(i + 1, bins[static_cast<size_t>(i)]);
  return h;
}

TH1D* make_band_hist(const HistBand& band, const std::string& name,
                     int nbins, double xmin, double xmax) {
  auto* h = new TH1D(name.c_str(), "", nbins, xmin, xmax);
  h->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) {
    const double lo = band.lo[static_cast<size_t>(i)];
    const double hi = band.hi[static_cast<size_t>(i)];
    const double mid = 0.5 * (lo + hi);
    const double err = 0.5 * std::max(0.0, hi - lo);
    h->SetBinContent(i + 1, mid);
    h->SetBinError(i + 1, err);
  }
  return h;
}

void apply_plot_style() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetEndErrorSize(0);
  gStyle->SetErrorX(0.5);
  gStyle->SetTitleSize(0.060, "XY");
  gStyle->SetLabelSize(0.050, "XY");
  gStyle->SetTitleOffset(0.85, "X");
  gStyle->SetTitleOffset(0.72, "Y");
  gStyle->SetNdivisions(506, "XY");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void style_cv_hist(TH1D* h, const std::string& x_title, const std::string& y_title = "Events") {
  h->SetTitle((";" + x_title + ";Events").c_str());
  h->GetYaxis()->SetTitle(y_title.c_str());
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.0);
  h->SetFillStyle(0);
}

void style_band_hist(TH1D* h, int color, float alpha) {
  h->SetLineColor(color);
  h->SetLineWidth(1);
  h->SetFillColorAlpha(color, alpha);
  h->SetFillStyle(1001);
  h->SetMarkerSize(0.0);
}

void style_band_edge(TH1D* h, int color, int line_width = 1) {
  h->SetLineColor(color);
  h->SetLineWidth(line_width);
  h->SetLineStyle(1);
  h->SetMarkerSize(0.0);
  h->SetFillStyle(0);
}

TLegend* make_band_legend(TH1D* h68, TH1D* h95, TH1D* h100, TH1D* h_cv,
                          TH1D* h_data = nullptr,
                          double x1 = 0.58, double y1 = 0.63,
                          double x2 = 0.92, double y2 = 0.89) {
  auto* leg = new TLegend(x1, y1, x2, y2);
  leg->SetNColumns(1);
  leg->SetTextSize(0.040);
  leg->AddEntry(h68, "Universe (68%)", "lf");
  leg->AddEntry(h95, "Universe (95%)", "lf");
  leg->AddEntry(h100, "Universe (100%)", "lf");
  leg->AddEntry(h_cv, "Central Value", "l");
  if (h_data != nullptr) leg->AddEntry(h_data, "Data", "ep");
  return leg;
}

} // namespace

/*
  Draw the central-value inference-score distribution in black and summarize the
  enabled universe variations as bin-by-bin 68%, 95%, and 100% central
  quantile bands.

  Supported sources:
    source = "genie"       -> weightsGenie
    source = "reint"       -> weightsReint
    source = "ppfx"        -> weightsPPFX (optionally divided by ppfx_cv)
    source = "flux"        -> weightsFlux

  Important:
    For GENIE, w_nominal already contains weightTune. The GENIE multisim
    universes therefore need to be applied relative to weightTune, not on top
    of it. The same holds for the GENIE up/down knob vectors.

    source = "genie_up"    -> weightsGenieUp
    source = "genie_dn"    -> weightsGenieDn

  The histogrammed observable is the chosen entry of `inf_scores`.
*/
int plot_inference_score_universes(const std::string& samples_tsv = "",
                                   const char* extra_sel = "sel_muon",
                                   const std::string& nominal_weight = "w_nominal",
                                   const std::string& source = "genie",
                                   int score_index = 0,
                                   int nbins = 24,
                                   double xmin = -15.0,
                                   double xmax = 15.0,
                                   bool fold_overflow = true,
                                   bool include_ext = false,
                                   bool include_data = false,
                                   bool ppfx_is_already_in_nominal = true,
                                   const std::string& output_stem = "rareproc_score0_universes") {
  const std::string extra_sel_expr = (extra_sel != nullptr) ? extra_sel : "";
  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plot_inference_score_universes] input=" << list_path << "\n";

  if (!looks_like_event_list_root(list_path)) {
    std::cerr << "[plot_inference_score_universes] input is not an event list ROOT file: " << list_path << "\n";
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
      {"inf_scores"});

  ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                 .Filter(
                                     [mask_ext](int sid) {
                                       return !(sid >= 0 && sid < static_cast<int>(mask_ext->size()) &&
                                                (*mask_ext)[static_cast<std::size_t>(sid)]);
                                     },
                                     {"sample_id"})
                                 .Define("__w__", nominal_weight);
  ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext).Define("__w__", nominal_weight);
  ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

  node_mc = apply_optional_filter(node_mc, extra_sel_expr, "extra selection");
  node_ext = apply_optional_filter(node_ext, extra_sel_expr, "extra selection");
  node_data = apply_optional_filter(node_data, extra_sel_expr, "extra selection");

  auto x_mc_h = node_mc.Take<float>("score0");
  auto w_mc_h = node_mc.Take<double>("__w__");
  const std::vector<double> cv_mc = build_cv_hist(*x_mc_h, *w_mc_h, nbins, xmin, xmax, fold_overflow);

  std::vector<std::vector<double>> univ_hists;
  std::string source_label = source;

  if (source == "genie") {
    if (!has_column(node_mc, "weightsGenie")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsGenie\n";
      return 1;
    }

    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenie");
    if (has_column(node_mc, "weightTune")) {
      auto tune_h = node_mc.Take<float>("weightTune");
      univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, &(*tune_h),
                                                      nbins, xmin, xmax, fold_overflow);
    } else {
      univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, nullptr,
                                                      nbins, xmin, xmax, fold_overflow);
    }
    source_label = "GENIE multisim";
  } else if (source == "reint") {
    if (!has_column(node_mc, "weightsReint")) {
      return 1;
    }
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsReint");
    univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins, xmin, xmax, fold_overflow);
    source_label = "Reinteraction multisim";
  } else if (source == "ppfx") {
    if (!has_column(node_mc, "weightsPPFX")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsPPFX\n";
      return 1;
    }
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsPPFX");
    if (ppfx_is_already_in_nominal && has_column(node_mc, "ppfx_cv")) {
      auto cv_h = node_mc.Take<float>("ppfx_cv");
      univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, &(*cv_h), nbins, xmin, xmax, fold_overflow);
    } else {
      univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins, xmin, xmax, fold_overflow);
    }
    source_label = "PPFX multisim";
  } else if (source == "flux") {
    if (!has_column(node_mc, "weightsFlux")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsFlux\n";
      return 1;
    }
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsFlux");
    univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins, xmin, xmax, fold_overflow);
    source_label = "Flux multisim";
  } else if (source == "genie_up") {
    if (!has_column(node_mc, "weightsGenieUp")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsGenieUp\n";
      return 1;
    }
    auto up_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenieUp");
    std::vector<ROOT::RVec<unsigned short>> dummy;
    if (has_column(node_mc, "weightTune")) {
      auto tune_h = node_mc.Take<float>("weightTune");
      const auto w_notune = divide_weights_by_cv_source(*w_mc_h, *tune_h);
      univ_hists = build_universe_hists_from_unisim(*x_mc_h, w_notune, *up_h, dummy,
                                                    nbins, xmin, xmax, fold_overflow, true);
    } else {
      univ_hists = build_universe_hists_from_unisim(*x_mc_h, *w_mc_h, *up_h, dummy,
                                                    nbins, xmin, xmax, fold_overflow, true);
    }
    source_label = "GENIE up unisims";
  } else if (source == "genie_dn") {
    if (!has_column(node_mc, "weightsGenieDn")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsGenieDn\n";
      return 1;
    }
    auto dn_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenieDn");
    std::vector<ROOT::RVec<unsigned short>> dummy;
    if (has_column(node_mc, "weightTune")) {
      auto tune_h = node_mc.Take<float>("weightTune");
      const auto w_notune = divide_weights_by_cv_source(*w_mc_h, *tune_h);
      univ_hists = build_universe_hists_from_unisim(*x_mc_h, w_notune, dummy, *dn_h,
                                                    nbins, xmin, xmax, fold_overflow, false);
    } else {
      univ_hists = build_universe_hists_from_unisim(*x_mc_h, *w_mc_h, dummy, *dn_h,
                                                    nbins, xmin, xmax, fold_overflow, false);
    }
    source_label = "GENIE down unisims";
  } else {
    std::cerr << "[plot_inference_score_universes] unsupported source '" << source
              << "'. Use one of: genie, reint, ppfx, flux, genie_up, genie_dn\n";
    return 1;
  }

  if (univ_hists.empty()) {
    std::cerr << "[plot_inference_score_universes] source '" << source << "' has zero universes.\n";
    return 1;
  }

  std::vector<double> cv_total = cv_mc;

  TH1D* h_ext = nullptr;
  if (include_ext) {
    auto x_ext_h = node_ext.Take<float>("score0");
    auto w_ext_h = node_ext.Take<double>("__w__");
    const std::vector<double> ext = build_cv_hist(*x_ext_h, *w_ext_h, nbins, xmin, xmax, fold_overflow);
    h_ext = make_hist1d(ext, "h_ext_overlay", nbins, xmin, xmax);
    for (int i = 0; i < nbins; ++i) cv_total[static_cast<size_t>(i)] += ext[static_cast<size_t>(i)];
    for (auto& hu : univ_hists) {
      for (int i = 0; i < nbins; ++i) hu[static_cast<size_t>(i)] += ext[static_cast<size_t>(i)];
    }
  }

  const HistBand band68 = build_central_band(univ_hists, nbins, 0.68);
  const HistBand band95 = build_central_band(univ_hists, nbins, 0.95);
  const HistBand band100 = build_central_band(univ_hists, nbins, 1.00);

  TH1D* h_cv = make_hist1d(cv_total, "h_cv", nbins, xmin, xmax);
  style_cv_hist(h_cv, "Inference score[" + std::to_string(score_index) + "]");

  const int color68 = TColor::GetColor("#EAD94C");
  const int color95 = TColor::GetColor("#67B7B0");
  const int color100 = TColor::GetColor("#7A4F8A");

  TH1D* h_band68 = make_band_hist(band68, "h_band_68", nbins, xmin, xmax);
  TH1D* h_band95 = make_band_hist(band95, "h_band_95", nbins, xmin, xmax);
  TH1D* h_band100 = make_band_hist(band100, "h_band_100", nbins, xmin, xmax);
  style_band_hist(h_band68, color68, 0.50f);
  style_band_hist(h_band95, color95, 0.30f);
  style_band_hist(h_band100, color100, 0.18f);

  TH1D* h_band68_lo = make_hist1d(band68.lo, "h_band_68_lo", nbins, xmin, xmax);
  TH1D* h_band68_hi = make_hist1d(band68.hi, "h_band_68_hi", nbins, xmin, xmax);
  TH1D* h_band95_lo = make_hist1d(band95.lo, "h_band_95_lo", nbins, xmin, xmax);
  TH1D* h_band95_hi = make_hist1d(band95.hi, "h_band_95_hi", nbins, xmin, xmax);
  TH1D* h_band100_lo = make_hist1d(band100.lo, "h_band_100_lo", nbins, xmin, xmax);
  TH1D* h_band100_hi = make_hist1d(band100.hi, "h_band_100_hi", nbins, xmin, xmax);
  style_band_edge(h_band68_lo, color68);
  style_band_edge(h_band68_hi, color68);
  style_band_edge(h_band95_lo, color95);
  style_band_edge(h_band95_hi, color95);
  style_band_edge(h_band100_lo, color100);
  style_band_edge(h_band100_hi, color100);

  double ymax = std::max(h_cv->GetMaximum(), max_band_value(band100));

  TH1D* h_data = nullptr;
  if (include_data) {
    auto x_data_h = node_data.Take<float>("score0");
    std::vector<double> data_bins(static_cast<size_t>(nbins), 0.0);
    for (const float x : *x_data_h) {
      const int b = find_bin(x, nbins, xmin, xmax, fold_overflow);
      if (b >= 0) data_bins[static_cast<size_t>(b)] += 1.0;
    }
    h_data = make_hist1d(data_bins, "h_data", nbins, xmin, xmax);
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.1);
    h_data->SetLineWidth(2);
    ymax = std::max(ymax, h_data->GetMaximum());
  }

  if (ymax <= 0.0) ymax = 1.0;

  apply_plot_style();
  const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
  const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
  gSystem->mkdir(out_dir.c_str(), true);

  const std::string stem = output_stem.empty() ? "rareproc_score0_universes" : output_stem;
  const std::string out_path = out_dir + "/" + stem + "_" + source + "." + out_fmt;
  const std::string root_path = out_dir + "/" + stem + "_" + source + ".root";

  TCanvas c("c_score_universes", "c_score_universes", 980, 760);
  c.SetLeftMargin(0.11);
  c.SetRightMargin(0.03);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.08);

  style_cv_hist(h_cv, "Inference score[" + std::to_string(score_index) + "]");
  h_cv->GetYaxis()->SetRangeUser(0.0, 1.18 * ymax);
  h_cv->Draw("HIST");

  h_band100->Draw("E2 SAME");
  h_band95->Draw("E2 SAME");
  h_band68->Draw("E2 SAME");
  h_band100_lo->Draw("HIST SAME");
  h_band100_hi->Draw("HIST SAME");
  h_band95_lo->Draw("HIST SAME");
  h_band95_hi->Draw("HIST SAME");
  h_band68_lo->Draw("HIST SAME");
  h_band68_hi->Draw("HIST SAME");
  h_cv->Draw("HIST SAME");
  if (h_data != nullptr) h_data->Draw("E1 SAME");

  TLegend* leg = make_band_legend(h_band68, h_band95, h_band100, h_cv, h_data);
  leg->Draw();

  if (include_ext && h_ext != nullptr) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.14, 0.84, "EXT included in all templates");
  }

  c.RedrawAxis();
  c.SaveAs(out_path.c_str());

  TFile fout(root_path.c_str(), "RECREATE");
  h_cv->Write("h_cv");
  if (h_ext != nullptr) h_ext->Write("h_ext");
  h_band68->Write("h_band_68");
  h_band95->Write("h_band_95");
  h_band100->Write("h_band_100");
  h_band68_lo->Write("h_band_68_lo");
  h_band68_hi->Write("h_band_68_hi");
  h_band95_lo->Write("h_band_95_lo");
  h_band95_hi->Write("h_band_95_hi");
  h_band100_lo->Write("h_band_100_lo");
  h_band100_hi->Write("h_band_100_hi");
  if (h_data != nullptr) h_data->Write("h_data");
  for (size_t u = 0; u < univ_hists.size(); ++u) {
    std::unique_ptr<TH1D> h_universe(
        make_hist1d(univ_hists[u],
                    "h_universe_tmp_" + std::to_string(static_cast<unsigned long long>(u)),
                    nbins, xmin, xmax));
    fout.WriteObject(
        h_universe.get(),
        ("h_universe_" + std::to_string(static_cast<unsigned long long>(u))).c_str());
  }
  fout.Close();

  std::cout << "[plot_inference_score_universes] source = " << source_label << "\n";
  std::cout << "[plot_inference_score_universes] universes summarized = " << univ_hists.size() << "\n";
  std::cout << "[plot_inference_score_universes] wrote " << out_path << "\n";
  std::cout << "[plot_inference_score_universes] wrote " << root_path << "\n";

  delete h_cv;
  delete h_ext;
  delete h_data;
  delete h_band68;
  delete h_band95;
  delete h_band100;
  delete h_band68_lo;
  delete h_band68_hi;
  delete h_band95_lo;
  delete h_band95_hi;
  delete h_band100_lo;
  delete h_band100_hi;
  delete leg;

  return 0;
}

#endif
