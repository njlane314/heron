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
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
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

std::string trim_copy(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  return s.substr(b, e - b);
}

std::vector<int> parse_csv_ints(const std::string& csv) {
  std::vector<int> out;
  std::stringstream ss(csv);
  std::string item;
  while (std::getline(ss, item, ',')) {
    item = trim_copy(item);
    if (item.empty()) continue;
    try {
      out.push_back(std::stoi(item));
    } catch (...) {
      std::cerr << "[plot_systematic_percentage_stack] ignoring non-integer boundary '" << item << "'\n";
    }
  }
  return out;
}

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
    std::cerr << "[plot_systematic_percentage_stack] " << label << " column '" << expr
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

HistSummary make_empty_hist(int nbins) {
  HistSummary h;
  h.sumw.assign(static_cast<size_t>(nbins), 0.0);
  h.sumw2.assign(static_cast<size_t>(nbins), 0.0);
  return h;
}

HistSummary build_hist_summary(const std::vector<float>& x, const std::vector<double>& w, int nbins, double xmin,
                               double xmax, bool fold_overflow) {
  HistSummary out = make_empty_hist(nbins);
  const size_t n = std::min(x.size(), w.size());
  for (size_t i = 0; i < n; ++i) {
    const int b = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (b < 0) continue;
    out.sumw[static_cast<size_t>(b)] += w[i];
    out.sumw2[static_cast<size_t>(b)] += w[i] * w[i];
  }
  return out;
}

std::vector<double> make_zero_vector(int nbins) {
  return std::vector<double>(static_cast<size_t>(nbins), 0.0);
}

void add_in_place(std::vector<double>& dst, const std::vector<double>& src) {
  if (dst.size() < src.size()) dst.resize(src.size(), 0.0);
  for (size_t i = 0; i < src.size(); ++i) dst[i] += src[i];
}

std::vector<double> make_fullcorr_diag(const HistSummary& h, double frac) {
  std::vector<double> out(h.sumw.size(), 0.0);
  const double f2 = frac * frac;
  for (size_t i = 0; i < h.sumw.size(); ++i) out[i] = h.sumw[i] * h.sumw[i] * f2;
  return out;
}

std::vector<double> build_diag_multisim_variance_from_vectors(
    const std::vector<float>& x,
    const std::vector<double>& w_nom,
    const std::vector<ROOT::RVec<unsigned short>>& w_univ_ushort,
    const std::vector<float>* w_cv_src,
    int nbins, double xmin, double xmax,
    bool fold_overflow, bool average_universes) {
  const size_t n_evt = std::min(x.size(), w_nom.size());

  size_t n_universes = 0;
  const size_t n_vec = std::min(n_evt, w_univ_ushort.size());
  for (size_t i = 0; i < n_vec; ++i) n_universes = std::max(n_universes, w_univ_ushort[i].size());

  std::vector<double> out(static_cast<size_t>(nbins), 0.0);
  if (n_universes == 0) return out;

  std::vector<double> cv_bin(static_cast<size_t>(nbins), 0.0);
  std::vector<std::vector<double>> univ_bin(n_universes, std::vector<double>(static_cast<size_t>(nbins), 0.0));

  for (size_t i = 0; i < n_evt; ++i) {
    const int b = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (b < 0) continue;

    const double w0 = w_nom[i];
    cv_bin[static_cast<size_t>(b)] += w0;

    static const ROOT::RVec<unsigned short> empty;
    const auto& wu = (i < w_univ_ushort.size()) ? w_univ_ushort[i] : empty;
    const double cv_src = (w_cv_src != nullptr && i < w_cv_src->size()) ? static_cast<double>((*w_cv_src)[i]) : 1.0;

    for (size_t u = 0; u < n_universes; ++u) {
      double ratio = 1.0;
      if (u < wu.size()) ratio = 1.0e-3 * static_cast<double>(wu[u]);
      if (w_cv_src != nullptr && cv_src > 0.0) ratio /= cv_src;
      univ_bin[u][static_cast<size_t>(b)] += w0 * ratio;
    }
  }

  for (size_t u = 0; u < n_universes; ++u) {
    for (int b = 0; b < nbins; ++b) {
      const double d = cv_bin[static_cast<size_t>(b)] - univ_bin[u][static_cast<size_t>(b)];
      out[static_cast<size_t>(b)] += d * d;
    }
  }

  if (average_universes && n_universes > 0) {
    const double inv = 1.0 / static_cast<double>(n_universes);
    for (double& v : out) v *= inv;
  }

  return out;
}

std::vector<double> build_diag_unisim_variance_from_vectors(
    const std::vector<float>& x,
    const std::vector<double>& w_nom,
    const std::vector<ROOT::RVec<unsigned short>>& w_up_ushort,
    const std::vector<ROOT::RVec<unsigned short>>& w_dn_ushort,
    const std::vector<float>* w_cv_src,
    int nbins, double xmin, double xmax,
    bool fold_overflow) {
  const size_t n_evt = std::min(x.size(), w_nom.size());

  size_t n_knobs = 0;
  const size_t n_up = std::min(n_evt, w_up_ushort.size());
  const size_t n_dn = std::min(n_evt, w_dn_ushort.size());
  for (size_t i = 0; i < n_up; ++i) n_knobs = std::max(n_knobs, w_up_ushort[i].size());
  for (size_t i = 0; i < n_dn; ++i) n_knobs = std::max(n_knobs, w_dn_ushort[i].size());

  std::vector<double> out(static_cast<size_t>(nbins), 0.0);
  if (n_knobs == 0) return out;

  std::vector<std::vector<double>> delta(n_knobs, std::vector<double>(static_cast<size_t>(nbins), 0.0));

  for (size_t i = 0; i < n_evt; ++i) {
    const int b = find_bin(x[i], nbins, xmin, xmax, fold_overflow);
    if (b < 0) continue;

    static const ROOT::RVec<unsigned short> empty;
    const auto& wup = (i < w_up_ushort.size()) ? w_up_ushort[i] : empty;
    const auto& wdn = (i < w_dn_ushort.size()) ? w_dn_ushort[i] : empty;
    const double cv_src = (w_cv_src != nullptr && i < w_cv_src->size()) ? static_cast<double>((*w_cv_src)[i]) : 1.0;

    for (size_t k = 0; k < n_knobs; ++k) {
      double rup = 1.0;
      double rdn = 1.0;
      if (k < wup.size()) rup = 1.0e-3 * static_cast<double>(wup[k]);
      if (k < wdn.size()) rdn = 1.0e-3 * static_cast<double>(wdn[k]);
      if (w_cv_src != nullptr && cv_src > 0.0) {
        rup /= cv_src;
        rdn /= cv_src;
      }
      const double d = 0.5 * w_nom[i] * (rup - rdn);
      delta[k][static_cast<size_t>(b)] += d;
    }
  }

  for (size_t k = 0; k < n_knobs; ++k) {
    for (int b = 0; b < nbins; ++b) {
      const double d = delta[k][static_cast<size_t>(b)];
      out[static_cast<size_t>(b)] += d * d;
    }
  }

  return out;
}

void apply_plot_style() {
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.055, "XY");
  gStyle->SetLabelSize(0.045, "XY");
  gStyle->SetTitleOffset(1.05, "X");
  gStyle->SetTitleOffset(1.10, "Y");
  gStyle->SetNdivisions(510, "XY");
}

TH1D* make_percent_hist(const std::vector<double>& var, const std::vector<double>& denom, const std::string& name,
                        const std::string& title, int color, int nbins) {
  auto* h = new TH1D(name.c_str(), title.c_str(), nbins, 0.0, static_cast<double>(nbins));
  h->SetDirectory(nullptr);
  h->SetFillColor(color);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);

  for (int i = 0; i < nbins; ++i) {
    const double d = (i < static_cast<int>(denom.size())) ? denom[static_cast<size_t>(i)] : 0.0;
    const double v = (i < static_cast<int>(var.size())) ? var[static_cast<size_t>(i)] : 0.0;
    const double pct = (d > 0.0) ? (100.0 * std::max(0.0, v) / d) : 0.0;
    h->SetBinContent(i + 1, pct);
    h->SetBinError(i + 1, 0.0);
  }

  return h;
}

}  // namespace

/*
  Build a stacked percentage plot directly from the event list, without reading
  an intermediate systematics ROOT file.

  For each bin i of the chosen inference-score histogram, the macro computes the
  diagonal variance contribution V_i^(k) for each systematic category k and plots

      p_i^(k) = 100 * V_i^(k) / sum_j V_i^(j)

  Categories in this implementation:
    Flux     <- PPFX multisim + Flux multisim
    Xs       <- GENIE multisim + reinteraction multisim + GENIE unisims + MC full corr
    Detector <- category full corr
    MC stat  <- diagonal MC statistical variance
    Dirt     <- diagonal EXT statistical variance
*/
int plot_systematic_percentage_stack(const std::string& samples_tsv = "",
                                     const char* extra_sel = "sel_muon",
                                     const std::string& nominal_weight = "w_nominal",
                                     int score_index = 0,
                                     int nbins = 24,
                                     double xmin = -15.0,
                                     double xmax = 15.0,
                                     bool fold_overflow = true,
                                     bool average_universes = true,
                                     bool add_mcstat = true,
                                     bool add_extstat = true,
                                     double mc_fullcorr_frac = 0.0,
                                     double category_fullcorr_frac = 0.0,
                                     const std::string& category_sel = "is_signal",
                                     bool use_genie = true,
                                     bool use_reint = true,
                                     bool use_ppfx = true,
                                     bool ppfx_is_already_in_nominal = true,
                                     bool use_flux = false,
                                     bool use_genie_unisims = false,
                                     const std::string& boundary_csv = "",
                                     const std::string& output_stem = "rareproc_score0_syst_percentage_stack") {
  apply_plot_style();

  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plot_systematic_percentage_stack] input=" << list_path << "\n";

  if (!looks_like_event_list_root(list_path)) {
    std::cerr << "[plot_systematic_percentage_stack] input is not an event list ROOT file: " << list_path << "\n";
    return 2;
  }

  if (nbins < 1) nbins = 1;
  if (xmax < xmin) std::swap(xmin, xmax);
  if (xmax == xmin) xmax = xmin + 1.0;

  EventListIO el(list_path);
  ROOT::RDataFrame rdf0 = el.rdf();

  auto mask_ext = el.mask_for_ext();
  auto mask_mc = el.mask_for_mc_like();

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

  const std::string extra_sel_expr = (extra_sel != nullptr) ? extra_sel : "";
  node_mc = apply_optional_filter(node_mc, extra_sel_expr, "extra selection");
  node_ext = apply_optional_filter(node_ext, extra_sel_expr, "extra selection");

  auto x_mc_h = node_mc.Take<float>("score0");
  auto w_mc_h = node_mc.Take<double>("__w__");
  const HistSummary h_mc = build_hist_summary(*x_mc_h, *w_mc_h, nbins, xmin, xmax, fold_overflow);

  std::vector<double> var_flux = make_zero_vector(nbins);
  std::vector<double> var_xs = make_zero_vector(nbins);
  std::vector<double> var_detector = make_zero_vector(nbins);
  std::vector<double> var_mcstat = make_zero_vector(nbins);
  std::vector<double> var_dirt = make_zero_vector(nbins);

  if (use_genie && has_column(node_mc, "weightsGenie")) {
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenie");
    if (has_column(node_mc, "weightTune")) {
      auto tune_h = node_mc.Take<float>("weightTune");
      add_in_place(var_xs, build_diag_multisim_variance_from_vectors(*x_mc_h, *w_mc_h, *u_h, &(*tune_h), nbins,
                                                                     xmin, xmax, fold_overflow, average_universes));
    } else {
      add_in_place(var_xs, build_diag_multisim_variance_from_vectors(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins,
                                                                     xmin, xmax, fold_overflow, average_universes));
    }
  }

  if (use_reint && has_column(node_mc, "weightsReint")) {
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsReint");
    add_in_place(var_xs, build_diag_multisim_variance_from_vectors(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins,
                                                                   xmin, xmax, fold_overflow, average_universes));
  }

  if (use_ppfx && has_column(node_mc, "weightsPPFX")) {
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsPPFX");
    if (ppfx_is_already_in_nominal && has_column(node_mc, "ppfx_cv")) {
      auto cv_h = node_mc.Take<float>("ppfx_cv");
      add_in_place(var_flux, build_diag_multisim_variance_from_vectors(*x_mc_h, *w_mc_h, *u_h, &(*cv_h), nbins,
                                                                       xmin, xmax, fold_overflow, average_universes));
    } else {
      add_in_place(var_flux, build_diag_multisim_variance_from_vectors(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins,
                                                                       xmin, xmax, fold_overflow, average_universes));
    }
  }

  if (use_flux && has_column(node_mc, "weightsFlux")) {
    auto u_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsFlux");
    add_in_place(var_flux, build_diag_multisim_variance_from_vectors(*x_mc_h, *w_mc_h, *u_h, nullptr, nbins,
                                                                     xmin, xmax, fold_overflow, average_universes));
  }

  if (use_genie_unisims && has_column(node_mc, "weightsGenieUp") && has_column(node_mc, "weightsGenieDn")) {
    auto up_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenieUp");
    auto dn_h = node_mc.Take<ROOT::RVec<unsigned short>>("weightsGenieDn");
    if (has_column(node_mc, "weightTune")) {
      auto tune_h = node_mc.Take<float>("weightTune");
      add_in_place(var_xs, build_diag_unisim_variance_from_vectors(*x_mc_h, *w_mc_h, *up_h, *dn_h, &(*tune_h),
                                                                   nbins, xmin, xmax, fold_overflow));
    } else {
      add_in_place(var_xs, build_diag_unisim_variance_from_vectors(*x_mc_h, *w_mc_h, *up_h, *dn_h, nullptr,
                                                                   nbins, xmin, xmax, fold_overflow));
    }
  }

  if (add_mcstat) add_in_place(var_mcstat, h_mc.sumw2);
  if (mc_fullcorr_frac > 0.0) add_in_place(var_xs, make_fullcorr_diag(h_mc, mc_fullcorr_frac));

  if (category_fullcorr_frac > 0.0) {
    ROOT::RDF::RNode node_cat = apply_optional_filter(node_mc, category_sel, "category selection");
    auto x_cat_h = node_cat.Take<float>("score0");
    auto w_cat_h = node_cat.Take<double>("__w__");
    const HistSummary h_cat = build_hist_summary(*x_cat_h, *w_cat_h, nbins, xmin, xmax, fold_overflow);
    add_in_place(var_detector, make_fullcorr_diag(h_cat, category_fullcorr_frac));
  }

  if (add_extstat) {
    auto x_ext_h = node_ext.Take<float>("score0");
    auto w_ext_h = node_ext.Take<double>("__w__");
    const HistSummary h_ext = build_hist_summary(*x_ext_h, *w_ext_h, nbins, xmin, xmax, fold_overflow);
    add_in_place(var_dirt, h_ext.sumw2);
  }

  std::vector<double> total_var(static_cast<size_t>(nbins), 0.0);
  for (int i = 0; i < nbins; ++i) {
    total_var[static_cast<size_t>(i)] = var_flux[static_cast<size_t>(i)] +
                                        var_xs[static_cast<size_t>(i)] +
                                        var_detector[static_cast<size_t>(i)] +
                                        var_mcstat[static_cast<size_t>(i)] +
                                        var_dirt[static_cast<size_t>(i)];
  }

  double sum_total_var = 0.0;
  for (double v : total_var) sum_total_var += v;
  if (!(sum_total_var > 0.0)) {
    std::cerr << "[plot_systematic_percentage_stack] total variance is zero for the requested configuration.\n";
    return 1;
  }

  auto boundaries = parse_csv_ints(boundary_csv);

  auto* h_flux = make_percent_hist(var_flux, total_var, "h_pct_flux", "Flux", TColor::GetColor("#f0301a"), nbins);
  auto* h_xs = make_percent_hist(var_xs, total_var, "h_pct_xs", "Xs", TColor::GetColor("#1010f0"), nbins);
  auto* h_det = make_percent_hist(var_detector, total_var, "h_pct_detector", "Detector", TColor::GetColor("#d933e6"), nbins);
  auto* h_mcstat = make_percent_hist(var_mcstat, total_var, "h_pct_mcstat", "MC stat", TColor::GetColor("#61d23f"), nbins);
  auto* h_dirt = make_percent_hist(var_dirt, total_var, "h_pct_dirt", "Dirt", TColor::GetColor("#f0a33a"), nbins);

  THStack hs("hs_pct", ";Bin index;Syst. percentage");
  hs.Add(h_flux);
  hs.Add(h_xs);
  hs.Add(h_det);
  hs.Add(h_mcstat);
  hs.Add(h_dirt);

  const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
  const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
  gSystem->mkdir(out_dir.c_str(), true);

  TCanvas c("c_pct_stack", "c_pct_stack", 1320, 760);
  c.SetLeftMargin(0.10);
  c.SetRightMargin(0.24);
  c.SetBottomMargin(0.14);
  c.SetTopMargin(0.06);

  hs.Draw("HIST");
  hs.SetMinimum(0.0);
  hs.SetMaximum(105.0);
  hs.GetXaxis()->SetLimits(0.0, static_cast<double>(nbins));
  hs.GetXaxis()->SetTitle("Bin index");
  hs.GetYaxis()->SetTitle("Syst. percentage");
  hs.GetXaxis()->CenterTitle(false);
  hs.GetYaxis()->CenterTitle(false);

  for (int b : boundaries) {
    TLine line(static_cast<double>(b), 0.0, static_cast<double>(b), 105.0);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.SetLineColor(kGray + 2);
    line.Draw("SAME");
  }

  TLegend leg(0.76, 0.58, 0.97, 0.96);
  leg.SetBorderSize(1);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.045);
  leg.AddEntry(h_flux, "Flux", "f");
  leg.AddEntry(h_xs, "Xs", "f");
  leg.AddEntry(h_det, "Detector", "f");
  leg.AddEntry(h_mcstat, "MC stat", "f");
  leg.AddEntry(h_dirt, "Dirt", "f");
  leg.Draw();

  c.RedrawAxis();

  const std::string pdf_path = out_dir + "/" + output_stem + "." + out_fmt;
  const std::string root_path = out_dir + "/" + output_stem + ".root";
  c.SaveAs(pdf_path.c_str());

  TFile fout(root_path.c_str(), "RECREATE");
  h_flux->Write();
  h_xs->Write();
  h_det->Write();
  h_mcstat->Write();
  h_dirt->Write();
  fout.Close();

  std::cout << "[plot_systematic_percentage_stack] wrote " << pdf_path << "\n";
  std::cout << "[plot_systematic_percentage_stack] wrote " << root_path << "\n";

  delete h_flux;
  delete h_xs;
  delete h_det;
  delete h_mcstat;
  delete h_dirt;

  return 0;
}
