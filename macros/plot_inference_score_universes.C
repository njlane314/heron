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
#include <TFile.h>
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

TH1D* make_hist1d(const std::vector<double>& bins, const std::string& name,
                  int nbins, double xmin, double xmax) {
  auto* h = new TH1D(name.c_str(), "", nbins, xmin, xmax);
  h->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) h->SetBinContent(i + 1, bins[static_cast<size_t>(i)]);
  return h;
}

void style_cv_hist(TH1D* h, const std::string& x_title) {
  h->SetTitle((";" + x_title + ";Events").c_str());
  h->SetLineColor(kBlack);
  h->SetLineWidth(4);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.0);
}

void style_universe_hist(TH1D* h) {
  h->SetLineColor(kGreen + 1);
  h->SetLineWidth(2);
  h->SetLineStyle(1);
  h->SetMarkerSize(0.0);
}

} // namespace

/*
  Draw the central-value inference-score distribution in black and every enabled
  universe variation in bright green.

  Supported sources:
    source = "genie"       -> weightsGenie
    source = "reint"       -> weightsReint
    source = "ppfx"        -> weightsPPFX (optionally divided by ppfx_cv)
    source = "flux"        -> weightsFlux
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
    univ_hists = build_universe_hists_from_multisim(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins, xmin, xmax, fold_overflow);
    source_label = "GENIE multisim";
  } else if (source == "reint") {
    if (!has_column(node_mc, "weightsReint")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsReint\n";
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
    univ_hists = build_universe_hists_from_unisim(*x_mc_h, *w_mc_h, *up_h, dummy, nbins, xmin, xmax, fold_overflow, true);
    source_label = "GENIE up unisims";
  } else if (source == "genie_dn") {
    if (!has_column(node_mc, "weightsGenieDn")) {
      std::cerr << "[plot_inference_score_universes] missing branch weightsGenieDn\n";
      return 1;
    }
    auto dn_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenieDn");
    std::vector<ROOT::RVec<unsigned short>> dummy;
    univ_hists = build_universe_hists_from_unisim(*x_mc_h, *w_mc_h, dummy, *dn_h, nbins, xmin, xmax, fold_overflow, false);
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

  TH1D* h_cv = make_hist1d(cv_total, "h_cv", nbins, xmin, xmax);
  style_cv_hist(h_cv, "Inference score[" + std::to_string(score_index) + "]");

  std::vector<TH1D*> h_univ;
  h_univ.reserve(univ_hists.size());
  double ymax = h_cv->GetMaximum();

  for (size_t u = 0; u < univ_hists.size(); ++u) {
    TH1D* h = make_hist1d(univ_hists[u], "h_univ_" + std::to_string(static_cast<unsigned long long>(u)), nbins, xmin, xmax);
    style_universe_hist(h);
    ymax = std::max(ymax, h->GetMaximum());
    h_univ.push_back(h);
  }

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
    h_data->SetLineWidth(2);
    ymax = std::max(ymax, h_data->GetMaximum());
  }

  if (ymax <= 0.0) ymax = 1.0;

  gStyle->SetOptStat(0);
  const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
  const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
  gSystem->mkdir(out_dir.c_str(), true);

  const std::string stem = output_stem.empty() ? "rareproc_score0_universes" : output_stem;
  const std::string out_path = out_dir + "/" + stem + "_" + source + "." + out_fmt;
  const std::string root_path = out_dir + "/" + stem + "_" + source + ".root";

  TCanvas c("c_score_universes", "c_score_universes", 980, 760);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetBottomMargin(0.12);

  h_cv->GetYaxis()->SetRangeUser(0.0, 1.30 * ymax);
  h_cv->Draw("HIST");

  for (TH1D* h : h_univ) h->Draw("HIST SAME");
  h_cv->Draw("HIST SAME");
  if (h_data != nullptr) h_data->Draw("E1 SAME");

  TLegend leg(0.54, 0.68, 0.90, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(h_cv, "CV", "l");
  leg.AddEntry(h_univ.front(), (source_label + " universes").c_str(), "l");
  if (include_ext && h_ext != nullptr) leg.AddEntry((TObject*)nullptr, "EXT included in all templates", "");
  if (h_data != nullptr) leg.AddEntry(h_data, "Data", "ep");
  leg.Draw();

  c.RedrawAxis();
  c.SaveAs(out_path.c_str());

  TFile fout(root_path.c_str(), "RECREATE");
  h_cv->Write("h_cv");
  if (h_ext != nullptr) h_ext->Write("h_ext");
  if (h_data != nullptr) h_data->Write("h_data");
  for (size_t u = 0; u < h_univ.size(); ++u) {
    fout.WriteObject(h_univ[u], ("h_universe_" + std::to_string(static_cast<unsigned long long>(u))).c_str());
  }
  fout.Close();

  std::cout << "[plot_inference_score_universes] source = " << source_label << "\n";
  std::cout << "[plot_inference_score_universes] universes drawn = " << h_univ.size() << "\n";
  std::cout << "[plot_inference_score_universes] wrote " << out_path << "\n";
  std::cout << "[plot_inference_score_universes] wrote " << root_path << "\n";

  delete h_cv;
  delete h_ext;
  delete h_data;
  for (TH1D* h : h_univ) delete h;

  return 0;
}

#endif
