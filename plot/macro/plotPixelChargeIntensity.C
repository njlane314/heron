// plot/macro/plotPixelChargeIntensity.C
//
// Plot per-event pixel "charge intensity" distributions, split by semantic class:
//   - cosmic: semantic == Cosmic (index 1)
//   - neutrino: semantic >= Muon (index >= 2)
// The plotted variables are per event (summing all planes):
//   - mean_adc_* : (sum of adc over pixels in category) / (number of pixels in category)
//   - sum_adc_*  : sum of adc over pixels in category
//
// Notes:
//   - Requires semantic_image_u/v/w to be present (MC-like).
//   - This scans per-pixel arrays once per event per plane to accumulate sums.
//   - Histograms are normalised to "Fraction of Dataset [%]" like the example.
//
// Run with:
//   ./heron macro plotPixelChargeIntensity.C
//   ./heron macro plotPixelChargeIntensity.C 'plotPixelChargeIntensity("/path/to/event_list.root","true")'
//
// Output:
//   $HERON_PLOT_DIR (default: ./scratch/plots)
//   $HERON_PLOT_FORMAT (default: pdf)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
#include "PlottingHelper.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &p)
{
  const auto n = p.size();
  if (n < 5 || p.substr(n - 5) != ".root")
    return false;

  std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
  if (!f || f->IsZombie())
    return false;

  const bool has_refs = (f->Get("sample_refs") != nullptr);
  const bool has_events_tree = (f->Get("events") != nullptr);
  const bool has_event_tree_key = (f->Get("event_tree") != nullptr);

  return has_refs && (has_events_tree || has_event_tree_key);
}

std::string plot_out_dir()
{
  const char *v = std::getenv("HERON_PLOT_DIR");
  return v ? std::string(v) : std::string("./scratch/plots");
}

std::string plot_out_fmt()
{
  const char *v = std::getenv("HERON_PLOT_FORMAT");
  return v ? std::string(v) : std::string("pdf");
}

std::vector<double> log_edges(int nbins, double xmin, double xmax)
{
  std::vector<double> edges;
  edges.reserve(static_cast<size_t>(nbins) + 1);

  xmin = std::max(xmin, 1e-12);
  xmax = std::max(xmax, xmin * 1.0001);

  const double logmin = std::log10(xmin);
  const double logmax = std::log10(xmax);
  for (int i = 0; i <= nbins; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(nbins);
    edges.push_back(std::pow(10.0, logmin + t * (logmax - logmin)));
  }
  return edges;
}

void normalise_to_percent(TH1 &h)
{
  const double tot = h.Integral(0, h.GetNbinsX() + 1);
  if (tot > 0.0)
    h.Scale(100.0 / tot);
}

void style_hist(TH1 &h, int color, double alpha)
{
  h.SetLineColor(color);
  h.SetLineWidth(2);
  h.SetFillColorAlpha(color, alpha);
}

int require_columns(const std::unordered_set<std::string> &columns,
                    const std::vector<std::string> &required,
                    const std::string &label)
{
  std::vector<std::string> missing;
  missing.reserve(required.size());
  for (const auto &name : required)
  {
    if (columns.find(name) == columns.end())
      missing.push_back(name);
  }

  if (missing.empty())
    return 0;

  std::cerr << "[plotPixelChargeIntensity] missing required columns for " << label << ":\n";
  for (const auto &name : missing)
    std::cerr << "  - " << name << "\n";
  return 1;
}

int at_or_zero(const ROOT::VecOps::RVec<int> &v, int idx)
{
  if (idx < 0)
    return 0;
  const size_t u = static_cast<size_t>(idx);
  if (u >= v.size())
    return 0;
  return v[u];
}

int sum_from_or_zero(const ROOT::VecOps::RVec<int> &v, int idx_from)
{
  if (idx_from < 0)
    idx_from = 0;
  const size_t start = static_cast<size_t>(idx_from);
  int sum = 0;
  for (size_t i = start; i < v.size(); ++i)
    sum += v[i];
  return sum;
}

// Return {sum_all_pos, sum_cosmic, sum_neutrino}
std::array<double, 3> sum_adc_components(const ROOT::VecOps::RVec<float> &adc,
                                         const ROOT::VecOps::RVec<int32_t> &sem)
{
  double sum_all = 0.0;
  double sum_cos = 0.0;
  double sum_nu = 0.0;

  const size_t n = std::min(adc.size(), sem.size());
  for (size_t i = 0; i < n; ++i)
  {
    const float a = adc[i];
    if (a <= 0.0f)
      continue;

    sum_all += static_cast<double>(a);
    const int32_t s = sem[i];
    if (s == 1)
      sum_cos += static_cast<double>(a);
    else if (s > 1)
      sum_nu += static_cast<double>(a);
  }

  return {sum_all, sum_cos, sum_nu};
}

} // namespace

int plotPixelChargeIntensity(const std::string &samples_tsv = "",
                             const std::string &extra_sel = "true",
                             int nbins = 60,
                             double mean_xmin = 1e-1,
                             double mean_xmax = 1e3,
                             double sum_xmin = 1.0,
                             double sum_xmax = 1e7)
{
  ROOT::EnableImplicitMT();
  std::cout << "[plotPixelChargeIntensity] implicit MT enabled\n";

  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plotPixelChargeIntensity] input=" << list_path << "\n";
  std::cout << "[plotPixelChargeIntensity] extra_sel=" << extra_sel << "\n";

  if (!looks_like_event_list_root(list_path))
  {
    std::cerr << "[plotPixelChargeIntensity] input is not an event list root file: " << list_path << "\n";
    return 1;
  }

  EventListIO el(list_path);
  ROOT::RDataFrame rdf = el.rdf();

  auto mask_ext = el.mask_for_ext();
  auto mask_mc = el.mask_for_mc_like();

  auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) -> ROOT::RDF::RNode {
    if (!mask)
      return n;
    return n.Filter(
        [mask](int sid) {
          return sid >= 0 && sid < static_cast<int>(mask->size()) &&
                 (*mask)[static_cast<size_t>(sid)];
        },
        {"sample_id"});
  };

  // MC-only: MC-like minus EXT.
  ROOT::RDF::RNode node = filter_by_mask(rdf, mask_mc);
  if (mask_ext)
  {
    node = node.Filter(
        [mask_ext](int sid) {
          return !(sid >= 0 && sid < static_cast<int>(mask_ext->size()) &&
                   (*mask_ext)[static_cast<size_t>(sid)]);
        },
        {"sample_id"});
  }

  if (!extra_sel.empty())
    node = node.Filter(extra_sel);

  const auto col_names = node.GetColumnNames();
  const std::unordered_set<std::string> columns(col_names.begin(), col_names.end());

  const std::vector<std::string> required = {
      "sample_id",
      "detector_image_u",
      "detector_image_v",
      "detector_image_w",
      "semantic_image_u",
      "semantic_image_v",
      "semantic_image_w",
      "active_pixels_u",
      "active_pixels_v",
      "active_pixels_w",
      "slice_semantic_active_pixels_u",
      "slice_semantic_active_pixels_v",
      "slice_semantic_active_pixels_w"};

  if (require_columns(columns, required, "pixel charge intensity") != 0)
    return 1;

  const int kCosmicIdx = 1;
  const int kFirstNuIdx = 2;

  auto n = node
               // per-plane sums: {all, cosmic, neutrino}
               .Define("sum3_u", sum_adc_components, {"detector_image_u", "semantic_image_u"})
               .Define("sum3_v", sum_adc_components, {"detector_image_v", "semantic_image_v"})
               .Define("sum3_w", sum_adc_components, {"detector_image_w", "semantic_image_w"})

               .Define("sum_adc_all", "sum3_u[0] + sum3_v[0] + sum3_w[0]")
               .Define("sum_adc_cosmic", "sum3_u[1] + sum3_v[1] + sum3_w[1]")
               .Define("sum_adc_neutrino", "sum3_u[2] + sum3_v[2] + sum3_w[2]")

               .Define("n_active_all", "active_pixels_u + active_pixels_v + active_pixels_w")

               .Define("n_cosmic",
                       [=](const ROOT::VecOps::RVec<int> &u,
                           const ROOT::VecOps::RVec<int> &v,
                           const ROOT::VecOps::RVec<int> &w) {
                         return at_or_zero(u, kCosmicIdx) + at_or_zero(v, kCosmicIdx) + at_or_zero(w, kCosmicIdx);
                       },
                       {"slice_semantic_active_pixels_u", "slice_semantic_active_pixels_v", "slice_semantic_active_pixels_w"})
               .Define("n_neutrino",
                       [=](const ROOT::VecOps::RVec<int> &u,
                           const ROOT::VecOps::RVec<int> &v,
                           const ROOT::VecOps::RVec<int> &w) {
                         return sum_from_or_zero(u, kFirstNuIdx) + sum_from_or_zero(v, kFirstNuIdx) + sum_from_or_zero(w, kFirstNuIdx);
                       },
                       {"slice_semantic_active_pixels_u", "slice_semantic_active_pixels_v", "slice_semantic_active_pixels_w"})

               .Define("mean_adc_all",
                       [](double sum, int npx) {
                         if (npx <= 0)
                           return 0.0;
                         return sum / static_cast<double>(npx);
                       },
                       {"sum_adc_all", "n_active_all"})
               .Define("mean_adc_cosmic",
                       [](double sum, int npx) {
                         if (npx <= 0)
                           return 0.0;
                         return sum / static_cast<double>(npx);
                       },
                       {"sum_adc_cosmic", "n_cosmic"})
               .Define("mean_adc_neutrino",
                       [](double sum, int npx) {
                         if (npx <= 0)
                           return 0.0;
                         return sum / static_cast<double>(npx);
                       },
                       {"sum_adc_neutrino", "n_neutrino"});

  const auto mean_edges = log_edges(nbins, mean_xmin, mean_xmax);
  const auto sum_edges = log_edges(nbins, sum_xmin, sum_xmax);

  auto mean_model = [&](const std::string &name) {
    return ROOT::RDF::TH1DModel(name.c_str(),
                               ";Mean ADC per (active) pixel / event;Fraction of Dataset [%]",
                               nbins, mean_edges.data());
  };
  auto sum_model = [&](const std::string &name) {
    return ROOT::RDF::TH1DModel(name.c_str(),
                               ";Sum ADC over (active) pixels / event;Fraction of Dataset [%]",
                               nbins, sum_edges.data());
  };

  auto h_mean_cos = n.Histo1D(mean_model("h_mean_cos"), "mean_adc_cosmic");
  auto h_mean_nu = n.Histo1D(mean_model("h_mean_nu"), "mean_adc_neutrino");
  auto h_sum_cos = n.Histo1D(sum_model("h_sum_cos"), "sum_adc_cosmic");
  auto h_sum_nu = n.Histo1D(sum_model("h_sum_nu"), "sum_adc_neutrino");

  // Also produce "all active pixels" references (often useful).
  auto h_mean_all = n.Histo1D(mean_model("h_mean_all"), "mean_adc_all");
  auto h_sum_all = n.Histo1D(sum_model("h_sum_all"), "sum_adc_all");

  ROOT::RDF::RunGraphs({h_mean_cos, h_mean_nu, h_sum_cos, h_sum_nu, h_mean_all, h_sum_all});

  gStyle->SetOptStat(0);

  const std::string out_dir = plot_out_dir();
  const std::string fmt = plot_out_fmt();
  gSystem->mkdir(out_dir.c_str(), true);

  auto draw_two = [&](TH1D *h1, const char *l1,
                      TH1D *h2, const char *l2,
                      const std::string &tag,
                      bool logx) {
    normalise_to_percent(*h1);
    normalise_to_percent(*h2);

    style_hist(*h1, kAzure + 1, 0.35);
    style_hist(*h2, kOrange + 1, 0.35);

    const double ymax = 1.15 * std::max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMaximum(ymax);
    h1->SetTitle("");

    TCanvas c(("c_" + tag).c_str(), ("c_" + tag).c_str(), 900, 550);
    if (logx)
      c.SetLogx();

    h1->Draw("hist");
    h2->Draw("hist same");

    TLegend leg(0.14, 0.84, 0.44, 0.96);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h1, l1, "f");
    leg.AddEntry(h2, l2, "f");
    leg.Draw();

    const std::string out = out_dir + "/" + tag + "." + fmt;
    c.SaveAs(out.c_str());
    std::cout << "[plotPixelChargeIntensity] wrote " << out << "\n";
  };

  auto draw_one = [&](TH1D *h, const char *label, const std::string &tag, bool logx) {
    normalise_to_percent(*h);
    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->SetTitle("");

    TCanvas c(("c_" + tag).c_str(), ("c_" + tag).c_str(), 900, 550);
    if (logx)
      c.SetLogx();

    h->Draw("hist");

    TLegend leg(0.14, 0.84, 0.44, 0.96);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h, label, "l");
    leg.Draw();

    const std::string out = out_dir + "/" + tag + "." + fmt;
    c.SaveAs(out.c_str());
    std::cout << "[plotPixelChargeIntensity] wrote " << out << "\n";
  };

  draw_two(&(*h_mean_cos), "Cosmic Pixels", &(*h_mean_nu), "Neutrino Pixels",
           "img_mean_adc_cosmic_vs_neutrino", true);

  draw_two(&(*h_sum_cos), "Cosmic Pixels", &(*h_sum_nu), "Neutrino Pixels",
           "img_sum_adc_cosmic_vs_neutrino", true);

  draw_one(&(*h_mean_all), "All active pixels", "img_mean_adc_all_active", true);
  draw_one(&(*h_sum_all), "All active pixels", "img_sum_adc_all_active", true);

  return 0;
}
