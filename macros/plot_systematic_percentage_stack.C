#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "PlotEnv.hh"
#include "PlottingHelper.hh"

using namespace nu;

namespace {

std::string lower_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

std::string trim_copy(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  return s.substr(b, e - b);
}

std::vector<std::string> parse_csv_strings(const std::string& csv) {
  std::vector<std::string> out;
  std::stringstream ss(csv);
  std::string item;
  while (std::getline(ss, item, ',')) {
    item = trim_copy(lower_copy(item));
    if (!item.empty()) out.push_back(item);
  }
  return out;
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

bool contains_any_token(const std::string& text_lower, const std::vector<std::string>& patterns_lower) {
  for (const auto& p : patterns_lower) {
    if (!p.empty() && text_lower.find(p) != std::string::npos) return true;
  }
  return false;
}

std::string strip_total_suffix(std::string s) {
  const std::string suffix = " total";
  const std::string s_lower = lower_copy(s);
  if (s_lower.size() >= suffix.size() &&
      s_lower.substr(s_lower.size() - suffix.size(), suffix.size()) == suffix) {
    s.resize(s.size() - suffix.size());
  }
  return s;
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

std::vector<double> diagonal_from_matrix_hist(const TH2D& h) {
  const int nbins = std::min(h.GetNbinsX(), h.GetNbinsY());
  std::vector<double> out(static_cast<size_t>(nbins), 0.0);
  for (int i = 1; i <= nbins; ++i) {
    out[static_cast<size_t>(i - 1)] = h.GetBinContent(i, i);
  }
  return out;
}

}  // namespace

/*
  Build a stacked percentage plot from the component covariance matrices written by
  plot_inference_score_systematics[_split].

  For each displayed category k and bin i, the plotted value is

      p_i^(k) = 100 * V_i^(k) / sum_j V_i^(j)

  where V_i^(k) is the diagonal variance contribution in bin i from category k.

  Categories are assigned by substring matching against the titles of the input
  histograms h_cov_component_total_*. By default:
    Flux     <- "ppfx", "flux"
    Xs       <- "genie", "reinteraction", "xsec", "xs", "mc full corr"
    Detector <- "detector", "category full corr"
    MC stat  <- "mc stat"
    Dirt     <- "ext stat", "dirt"

  If your component names differ, override the pattern CSV strings.
*/
int plot_systematic_percentage_stack(const std::string& input_root = "",
                                     const std::string& output_stem = "rareproc_score0_syst_percentage_stack",
                                     const std::string& boundary_csv = "",
                                     const std::string& flux_patterns = "ppfx,flux",
                                     const std::string& xs_patterns = "genie,reinteraction,xsec,xs,mc full corr",
                                     const std::string& detector_patterns = "detector,category full corr",
                                     const std::string& mcstat_patterns = "mc stat",
                                     const std::string& dirt_patterns = "ext stat,dirt",
                                     bool normalise_to_shown_categories_only = true) {
  apply_plot_style();

  const std::string out_dir = env_or("HERON_PLOT_OUT_DIR", "./scratch/out");
  const std::string out_fmt = env_or("HERON_PLOT_OUT_FMT", "pdf");
  gSystem->mkdir(out_dir.c_str(), true);

  const std::string in_path = input_root.empty() ? (out_dir + "/rareproc_score0_systematics_split.root") : input_root;
  TFile fin(in_path.c_str(), "READ");
  if (fin.IsZombie()) {
    std::cerr << "[plot_systematic_percentage_stack] failed to open input file: " << in_path << "\n";
    return 1;
  }

  auto flux_tokens = parse_csv_strings(flux_patterns);
  auto xs_tokens = parse_csv_strings(xs_patterns);
  auto detector_tokens = parse_csv_strings(detector_patterns);
  auto mcstat_tokens = parse_csv_strings(mcstat_patterns);
  auto dirt_tokens = parse_csv_strings(dirt_patterns);
  auto boundaries = parse_csv_ints(boundary_csv);

  std::map<std::string, std::vector<double>> cat_var;
  std::vector<std::string> cat_order = {"Flux", "Xs", "Detector", "MC stat", "Dirt"};
  for (const auto& cat : cat_order) cat_var[cat] = {};

  std::vector<std::string> unmatched_titles;
  int nbins = -1;

  TIter next(fin.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    const std::string name = key->GetName();
    if (name.rfind("h_cov_component_total_", 0) != 0) continue;

    auto* h = dynamic_cast<TH2D*>(key->ReadObj());
    if (h == nullptr) continue;

    const std::string raw_title = h->GetTitle();
    const std::string base_title = strip_total_suffix(raw_title);
    const std::string title_lower = lower_copy(base_title);
    const std::vector<double> diag = diagonal_from_matrix_hist(*h);
    if (nbins < 0) nbins = static_cast<int>(diag.size());

    auto add_to = [&](const std::string& cat) {
      if (cat_var[cat].empty()) cat_var[cat].assign(diag.size(), 0.0);
      for (size_t i = 0; i < diag.size(); ++i) cat_var[cat][i] += std::max(0.0, diag[i]);
    };

    bool matched = false;
    if (contains_any_token(title_lower, flux_tokens)) {
      add_to("Flux");
      matched = true;
    }
    if (contains_any_token(title_lower, xs_tokens)) {
      add_to("Xs");
      matched = true;
    }
    if (contains_any_token(title_lower, detector_tokens)) {
      add_to("Detector");
      matched = true;
    }
    if (contains_any_token(title_lower, mcstat_tokens)) {
      add_to("MC stat");
      matched = true;
    }
    if (contains_any_token(title_lower, dirt_tokens)) {
      add_to("Dirt");
      matched = true;
    }

    if (!matched) unmatched_titles.push_back(base_title);
    delete h;
  }

  if (nbins <= 0) {
    std::cerr << "[plot_systematic_percentage_stack] no h_cov_component_total_* histograms found in " << in_path << "\n";
    return 1;
  }

  for (const auto& cat : cat_order) {
    if (cat_var[cat].empty()) cat_var[cat].assign(static_cast<size_t>(nbins), 0.0);
  }

  std::vector<double> denom(static_cast<size_t>(nbins), 0.0);
  for (int i = 0; i < nbins; ++i) {
    double shown_sum = 0.0;
    for (const auto& cat : cat_order) shown_sum += cat_var[cat][static_cast<size_t>(i)];
    denom[static_cast<size_t>(i)] = shown_sum;
  }

  if (!normalise_to_shown_categories_only) {
    auto* htot = dynamic_cast<TH2D*>(fin.Get("h_score_cov_total"));
    if (htot != nullptr) {
      const std::vector<double> total_diag = diagonal_from_matrix_hist(*htot);
      for (int i = 0; i < nbins && i < static_cast<int>(total_diag.size()); ++i) {
        denom[static_cast<size_t>(i)] = std::max(0.0, total_diag[static_cast<size_t>(i)]);
      }
    }
  }

  auto* h_flux = make_percent_hist(cat_var["Flux"], denom, "h_pct_flux", "Flux", TColor::GetColor("#f0301a"), nbins);
  auto* h_xs = make_percent_hist(cat_var["Xs"], denom, "h_pct_xs", "Xs", TColor::GetColor("#1010f0"), nbins);
  auto* h_det = make_percent_hist(cat_var["Detector"], denom, "h_pct_detector", "Detector", TColor::GetColor("#d933e6"), nbins);
  auto* h_mcstat = make_percent_hist(cat_var["MC stat"], denom, "h_pct_mcstat", "MC stat", TColor::GetColor("#61d23f"), nbins);
  auto* h_dirt = make_percent_hist(cat_var["Dirt"], denom, "h_pct_dirt", "Dirt", TColor::GetColor("#f0a33a"), nbins);

  THStack hs("hs_pct", ";Bin index;Syst. percentage");
  hs.Add(h_flux);
  hs.Add(h_xs);
  hs.Add(h_det);
  hs.Add(h_mcstat);
  hs.Add(h_dirt);

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

  if (!unmatched_titles.empty()) {
    std::cout << "[plot_systematic_percentage_stack] unmatched component titles:" << "\n";
    for (const auto& t : unmatched_titles) std::cout << "  - " << t << "\n";
    std::cout << "[plot_systematic_percentage_stack] adjust the pattern CSV arguments if needed.\n";
  }

  delete h_flux;
  delete h_xs;
  delete h_det;
  delete h_mcstat;
  delete h_dirt;

  return 0;
}
