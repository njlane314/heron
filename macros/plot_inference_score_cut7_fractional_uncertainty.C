#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TBox.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
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

using Cut7UShortRVec = ROOT::RVec<unsigned short>;

struct Cut7Component {
  std::string label;
  double frac = 0.0;  // sigma_i / total_prediction
  int color = kBlack;
};

template <class T>
using Cut7Result = ROOT::RDF::RResultPtr<T>;

bool cut7_has_column(ROOT::RDF::RNode node, const std::string& name) {
  const auto cols = node.GetColumnNames();
  return std::find(cols.begin(), cols.end(), name) != cols.end();
}

double cut7_multisim_variance_1bin(const std::vector<double>& w_nom,
                                   const std::vector<Cut7UShortRVec>& w_univ_ushort,
                                   const std::vector<float>* cv_src,
                                   bool average_universes) {
  const std::size_t n_evt = w_nom.size();
  if (n_evt == 0) return 0.0;

  std::size_t n_universes = 0;
  for (const auto& wu : w_univ_ushort) n_universes = std::max(n_universes, wu.size());
  if (n_universes == 0) return 0.0;

  const double cv = std::accumulate(w_nom.begin(), w_nom.end(), 0.0);
  std::vector<double> univ_sum(n_universes, 0.0);

  static const Cut7UShortRVec empty;
  for (std::size_t i = 0; i < n_evt; ++i) {
    const double w0 = w_nom[i];
    const auto& wu = (i < w_univ_ushort.size()) ? w_univ_ushort[i] : empty;
    const double cvw =
        (cv_src != nullptr && i < cv_src->size()) ? static_cast<double>((*cv_src)[i]) : 1.0;

    for (std::size_t u = 0; u < n_universes; ++u) {
      double ratio = 1.0;
      if (u < wu.size()) ratio = 1.0e-3 * static_cast<double>(wu[u]);
      if (cv_src != nullptr && cvw > 0.0) ratio /= cvw;
      univ_sum[u] += w0 * ratio;
    }
  }

  double var = 0.0;
  for (const double yu : univ_sum) {
    const double d = cv - yu;
    var += d * d;
  }

  if (average_universes) var /= static_cast<double>(n_universes);
  return var;
}

void draw_cut7_single_bin_fractional_uncertainty(const std::vector<Cut7Component>& components,
                                                 double total_frac,
                                                 const std::string& out_path) {
  gStyle->SetOptStat(0);

  TCanvas c("c_score0_cut7_frac", "c_score0_cut7_frac", 800, 650);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetBottomMargin(0.16);

  TH1D frame("h_score0_cut7_frac", ";Selection;Fractional uncertainty", 1, 0.0, 1.0);
  frame.SetDirectory(nullptr);
  frame.GetXaxis()->SetBinLabel(1, "score[0] >= 7");
  frame.SetMinimum(0.0);
  frame.SetMaximum(std::max(0.05, 1.35 * total_frac));
  frame.SetLineColor(kBlack);
  frame.SetLineWidth(1);
  frame.Draw();

  const double xlo = 0.08;
  const double xhi = 0.92;

  std::vector<std::unique_ptr<TBox>> boxes;
  boxes.reserve(components.size());

  TLegend leg(0.58, 0.60, 0.90, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);

  double qsum = 0.0;
  for (const auto& comp : components) {
    if (comp.frac <= 0.0) continue;

    const double y0 = std::sqrt(qsum);
    qsum += comp.frac * comp.frac;
    const double y1 = std::sqrt(qsum);

    boxes.emplace_back(std::make_unique<TBox>(xlo, y0, xhi, y1));
    boxes.back()->SetFillColorAlpha(comp.color, 0.55);
    boxes.back()->SetLineColor(comp.color);
    boxes.back()->SetLineWidth(2);
    boxes.back()->Draw();

    leg.AddEntry(boxes.back().get(), comp.label.c_str(), "f");
  }

  TLine total_line(xlo, total_frac, xhi, total_frac);
  total_line.SetLineColor(kBlack);
  total_line.SetLineWidth(3);
  total_line.Draw();

  leg.AddEntry(&total_line, "Total", "l");
  leg.Draw();

  frame.Draw("AXIS SAME");
  c.RedrawAxis();
  c.SaveAs(out_path.c_str());
}

}  // namespace

/*
  Fixed choices:
    - input: default_event_list_root()
    - selection: sel_muon
    - weight: w_nominal
    - cut: inf_scores[0] >= 7
    - uncertainties: GENIE + Reinteraction + PPFX + MC stat
    - output: ./scratch/out/score0_geq7_fractional_uncertainty.pdf
*/
int plot_inference_score_cut7_fractional_uncertainty() {
  constexpr bool kAverageUniverses = true;
  constexpr float kScoreCut = 7.0f;

  const std::string list_path = default_event_list_root();
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] input = " << list_path << "\n";

  if (!looks_like_event_list_root(list_path)) {
    std::cerr << "[plot_inference_score_cut7_fractional_uncertainty] "
              << "input is not an event list ROOT file: " << list_path << "\n";
    return 2;
  }

  EventListIO el(list_path);
  ROOT::RDataFrame rdf0 = el.rdf();

  auto mask_ext = el.mask_for_ext();
  auto mask_mc = el.mask_for_mc_like();

  auto filter_by_mask = [](ROOT::RDF::RNode node,
                           std::shared_ptr<const std::vector<char>> mask) {
    return node.Filter(
        [mask](int sid) {
          return sid >= 0 &&
                 sid < static_cast<int>(mask->size()) &&
                 (*mask)[static_cast<std::size_t>(sid)];
        },
        {"sample_id"});
  };

  ROOT::RDF::RNode base = SelectionService::decorate(rdf0).Define(
      "score0",
      [](const ROOT::RVec<float>& scores) {
        return scores.empty() ? -9999.0f : scores[0];
      },
      {"inf_scores"});

  if (!cut7_has_column(base, "sample_id")) {
    std::cerr << "[plot_inference_score_cut7_fractional_uncertainty] missing column 'sample_id'.\n";
    return 1;
  }
  if (!cut7_has_column(base, "w_nominal")) {
    std::cerr << "[plot_inference_score_cut7_fractional_uncertainty] missing column 'w_nominal'.\n";
    return 1;
  }
  if (!cut7_has_column(base, "sel_muon")) {
    std::cerr << "[plot_inference_score_cut7_fractional_uncertainty] missing column 'sel_muon'.\n";
    return 1;
  }

  const auto pass_score_cut = [=](float s) {
    return s >= kScoreCut;  // change to (s > kScoreCut) if you want strict >
  };

  ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                 .Filter(
                                     [mask_ext](int sid) {
                                       return !(sid >= 0 &&
                                                sid < static_cast<int>(mask_ext->size()) &&
                                                (*mask_ext)[static_cast<std::size_t>(sid)]);
                                     },
                                     {"sample_id"})
                                 .Filter([](bool pass) { return pass; }, {"sel_muon"})
                                 .Filter(pass_score_cut, {"score0"})
                                 .Define("__w__", "1.0 * w_nominal");

  ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext)
                                  .Filter([](bool pass) { return pass; }, {"sel_muon"})
                                  .Filter(pass_score_cut, {"score0"})
                                  .Define("__w__", "1.0 * w_nominal");

  auto w_mc_h = node_mc.Take<double>("__w__");
  auto w_ext_h = node_ext.Take<double>("__w__");

  const bool have_genie = cut7_has_column(node_mc, "weightsGenie");
  const bool have_reint = cut7_has_column(node_mc, "weightsReint");
  const bool have_ppfx = cut7_has_column(node_mc, "weightsPPFX");
  const bool have_ppfx_cv = cut7_has_column(node_mc, "ppfx_cv");

  std::unique_ptr<Cut7Result<std::vector<Cut7UShortRVec>>> genie_h;
  std::unique_ptr<Cut7Result<std::vector<Cut7UShortRVec>>> reint_h;
  std::unique_ptr<Cut7Result<std::vector<Cut7UShortRVec>>> ppfx_h;
  std::unique_ptr<Cut7Result<std::vector<float>>> ppfx_cv_h;

  if (have_genie) {
    genie_h.reset(new Cut7Result<std::vector<Cut7UShortRVec>>(
        node_mc.Take<Cut7UShortRVec>("weightsGenie")));
  }
  if (have_reint) {
    reint_h.reset(new Cut7Result<std::vector<Cut7UShortRVec>>(
        node_mc.Take<Cut7UShortRVec>("weightsReint")));
  }
  if (have_ppfx) {
    ppfx_h.reset(new Cut7Result<std::vector<Cut7UShortRVec>>(
        node_mc.Take<Cut7UShortRVec>("weightsPPFX")));
  }
  if (have_ppfx && have_ppfx_cv) {
    ppfx_cv_h.reset(new Cut7Result<std::vector<float>>(
        node_mc.Take<float>("ppfx_cv")));
  }

  const auto& w_mc = *w_mc_h;
  const auto& w_ext = *w_ext_h;

  const double y_mc = std::accumulate(w_mc.begin(), w_mc.end(), 0.0);
  const double y_ext = std::accumulate(w_ext.begin(), w_ext.end(), 0.0);
  const double y_pred = y_mc + y_ext;

  double var_mcstat = 0.0;
  for (const double w : w_mc) var_mcstat += w * w;

  double var_genie = 0.0;
  double var_reint = 0.0;
  double var_ppfx = 0.0;

  if (genie_h) {
    var_genie = cut7_multisim_variance_1bin(w_mc, **genie_h, nullptr, kAverageUniverses);
  }
  if (reint_h) {
    var_reint = cut7_multisim_variance_1bin(w_mc, **reint_h, nullptr, kAverageUniverses);
  }
  if (ppfx_h) {
    const std::vector<float>* ppfx_cv = ppfx_cv_h ? &(**ppfx_cv_h) : nullptr;
    var_ppfx = cut7_multisim_variance_1bin(w_mc, **ppfx_h, ppfx_cv, kAverageUniverses);
  }

  const double var_total = var_mcstat + var_genie + var_reint + var_ppfx;
  const double sigma = std::sqrt(std::max(0.0, var_total));
  const double frac = (y_pred > 0.0) ? (sigma / y_pred) : 0.0;
  const double frac_mcstat = (y_pred > 0.0) ? (std::sqrt(std::max(0.0, var_mcstat)) / y_pred) : 0.0;
  const double frac_genie  = (y_pred > 0.0) ? (std::sqrt(std::max(0.0, var_genie )) / y_pred) : 0.0;
  const double frac_reint  = (y_pred > 0.0) ? (std::sqrt(std::max(0.0, var_reint )) / y_pred) : 0.0;
  const double frac_ppfx   = (y_pred > 0.0) ? (std::sqrt(std::max(0.0, var_ppfx  )) / y_pred) : 0.0;

  std::vector<Cut7Component> components;
  components.push_back({"MC stat",        frac_mcstat, kGray + 2});
  components.push_back({"GENIE",          frac_genie,  kRed + 1});
  components.push_back({"Reinteraction",  frac_reint,  kGreen + 2});
  components.push_back({"PPFX",           frac_ppfx,   kBlue + 1});

  if (!have_genie) {
    std::cout << "[plot_inference_score_cut7_fractional_uncertainty] "
              << "weightsGenie missing; skipping GENIE.\n";
    components.erase(
        std::remove_if(components.begin(), components.end(),
                       [](const Cut7Component& c) { return c.label == "GENIE"; }),
        components.end());
  }
  if (!have_reint) {
    std::cout << "[plot_inference_score_cut7_fractional_uncertainty] "
              << "weightsReint missing; skipping reinteraction.\n";
    components.erase(
        std::remove_if(components.begin(), components.end(),
                       [](const Cut7Component& c) { return c.label == "Reinteraction"; }),
        components.end());
  }
  if (!have_ppfx) {
    std::cout << "[plot_inference_score_cut7_fractional_uncertainty] "
              << "weightsPPFX missing; skipping PPFX.\n";
    components.erase(
        std::remove_if(components.begin(), components.end(),
                       [](const Cut7Component& c) { return c.label == "PPFX"; }),
        components.end());
  } else if (!have_ppfx_cv) {
    std::cout << "[plot_inference_score_cut7_fractional_uncertainty] "
              << "ppfx_cv missing; using raw weightsPPFX.\n";
  }

  const std::string out_dir = "./scratch/out";
  const std::string out_path = out_dir + "/score0_geq7_fractional_uncertainty.pdf";
  gSystem->mkdir(out_dir.c_str(), true);
  draw_cut7_single_bin_fractional_uncertainty(components, frac, out_path);

  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] selected MC yield    = " << y_mc << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] selected EXT yield   = " << y_ext << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] total prediction     = " << y_pred << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] absolute uncertainty = " << sigma << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] MC stat frac         = " << frac_mcstat << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] GENIE frac           = " << frac_genie << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] Reinteraction frac   = " << frac_reint << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] PPFX frac            = " << frac_ppfx << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] fractional unc.      = " << frac << "\n";
  std::cout << "[plot_inference_score_cut7_fractional_uncertainty] wrote " << out_path << "\n";

  return 0;
}
