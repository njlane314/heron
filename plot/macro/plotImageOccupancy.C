// plot/macro/plotImageOccupancy.C
//
// Plot per-event semantic pixel occupancy (Cosmic vs Neutrino) as a percentage
// of all pixels in the slice images, in the style of the example plot.
//
// Definitions (per event):
//   - total_pixels = |detector_image_u| + |detector_image_v| + |detector_image_w|
//   - cosmic_pixels = sum_planes slice_semantic_active_pixels_[u/v/w][Cosmic]
//   - neutrino_pixels = sum_planes sum_{label>=Muon} slice_semantic_active_pixels_[u/v/w][label]
//   - percentages are 100 * (pixels / total_pixels)
//
// Notes:
//   - Requires MC-like samples with semantic labels filled.
//   - For data/EXT, slice_semantic_active_pixels_* is typically empty (or all zeros).
//   - Uses event_list_<analysis>.root format (EventListIO).
//
// Run with:
//   ./nuxsec macro plotImageOccupancy.C
//   ./nuxsec macro plotImageOccupancy.C 'plotImageOccupancy("/path/to/event_list.root","true")'
//   ./nuxsec macro plotImageOccupancy.C 'plotImageOccupancy("/path/to/event_list.root","sel_triggered_slice")'
//
// Output:
//   Saves plots to:
//     $NUXSEC_PLOT_DIR (default: ./scratch/plots)
//   with format:
//     $NUXSEC_PLOT_FORMAT (default: pdf)

#include <algorithm>
#include <cmath>
#include <cstdlib>
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
  const char *v = std::getenv("NUXSEC_PLOT_DIR");
  return v ? std::string(v) : std::string("./scratch/plots");
}

std::string plot_out_fmt()
{
  const char *v = std::getenv("NUXSEC_PLOT_FORMAT");
  return v ? std::string(v) : std::string("pdf");
}

std::vector<double> log_edges(int nbins, double xmin, double xmax)
{
  std::vector<double> edges;
  edges.reserve(static_cast<size_t>(nbins) + 1);

  // Guard.
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

  std::cerr << "[plotImageOccupancy] missing required columns for " << label << ":\n";
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

} // namespace

int plotImageOccupancy(const std::string &samples_tsv = "",
                       const std::string &extra_sel = "true",
                       int nbins = 60,
                       double xmin_pct = 1e-4,
                       double xmax_pct = 1e1,
                       bool draw_planes = true)
{
  ROOT::EnableImplicitMT();
  std::cout << "[plotImageOccupancy] implicit MT enabled\n";

  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plotImageOccupancy] input=" << list_path << "\n";
  std::cout << "[plotImageOccupancy] extra_sel=" << extra_sel << "\n";

  if (!looks_like_event_list_root(list_path))
  {
    std::cerr << "[plotImageOccupancy] input is not an event list root file: " << list_path << "\n";
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
      "slice_semantic_active_pixels_u",
      "slice_semantic_active_pixels_v",
      "slice_semantic_active_pixels_w"};

  if (require_columns(columns, required, "semantic occupancy") != 0)
    return 1;

  // SemanticClassifier enum order:
  //   0 Empty, 1 Cosmic, 2 Muon, ...
  const int kCosmicIdx = 1;
  const int kFirstNuIdx = 2;

  auto n = node
               .Define("n_pix_u", "(int)detector_image_u.size()")
               .Define("n_pix_v", "(int)detector_image_v.size()")
               .Define("n_pix_w", "(int)detector_image_w.size()")
               .Define("n_pix_tot", "n_pix_u + n_pix_v + n_pix_w")

               .Define("cosmic_u", [=](const ROOT::VecOps::RVec<int> &c) { return at_or_zero(c, kCosmicIdx); }, {"slice_semantic_active_pixels_u"})
               .Define("cosmic_v", [=](const ROOT::VecOps::RVec<int> &c) { return at_or_zero(c, kCosmicIdx); }, {"slice_semantic_active_pixels_v"})
               .Define("cosmic_w", [=](const ROOT::VecOps::RVec<int> &c) { return at_or_zero(c, kCosmicIdx); }, {"slice_semantic_active_pixels_w"})

               .Define("nu_u", [=](const ROOT::VecOps::RVec<int> &c) { return sum_from_or_zero(c, kFirstNuIdx); }, {"slice_semantic_active_pixels_u"})
               .Define("nu_v", [=](const ROOT::VecOps::RVec<int> &c) { return sum_from_or_zero(c, kFirstNuIdx); }, {"slice_semantic_active_pixels_v"})
               .Define("nu_w", [=](const ROOT::VecOps::RVec<int> &c) { return sum_from_or_zero(c, kFirstNuIdx); }, {"slice_semantic_active_pixels_w"})

               .Define("cosmic_tot", "cosmic_u + cosmic_v + cosmic_w")
               .Define("nu_tot", "nu_u + nu_v + nu_w")

               .Define("cosmic_pct_tot",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"cosmic_tot", "n_pix_tot"})
               .Define("nu_pct_tot",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"nu_tot", "n_pix_tot"})

               .Define("cosmic_pct_u",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"cosmic_u", "n_pix_u"})
               .Define("nu_pct_u",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"nu_u", "n_pix_u"})

               .Define("cosmic_pct_v",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"cosmic_v", "n_pix_v"})
               .Define("nu_pct_v",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"nu_v", "n_pix_v"})

               .Define("cosmic_pct_w",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"cosmic_w", "n_pix_w"})
               .Define("nu_pct_w",
                       [](int num, int den) {
                         if (den <= 0)
                           return 0.0;
                         return 100.0 * static_cast<double>(num) / static_cast<double>(den);
                       },
                       {"nu_w", "n_pix_w"});

  const auto edges = log_edges(nbins, xmin_pct, xmax_pct);

  auto make_model = [&](const std::string &name) {
    return ROOT::RDF::TH1DModel(name.c_str(),
                               ";Percentage of pixels / event;Fraction of Dataset [%]",
                               nbins, edges.data());
  };

  auto h_cos_all = n.Histo1D(make_model("h_cos_all"), "cosmic_pct_tot");
  auto h_nu_all = n.Histo1D(make_model("h_nu_all"), "nu_pct_tot");

  ROOT::RDF::RResultPtr<TH1D> h_cos_u, h_nu_u, h_cos_v, h_nu_v, h_cos_w, h_nu_w;
  if (draw_planes)
  {
    h_cos_u = n.Histo1D(make_model("h_cos_u"), "cosmic_pct_u");
    h_nu_u = n.Histo1D(make_model("h_nu_u"), "nu_pct_u");
    h_cos_v = n.Histo1D(make_model("h_cos_v"), "cosmic_pct_v");
    h_nu_v = n.Histo1D(make_model("h_nu_v"), "nu_pct_v");
    h_cos_w = n.Histo1D(make_model("h_cos_w"), "cosmic_pct_w");
    h_nu_w = n.Histo1D(make_model("h_nu_w"), "nu_pct_w");
  }

  if (draw_planes)
    ROOT::RDF::RunGraphs({h_cos_all, h_nu_all, h_cos_u, h_nu_u, h_cos_v, h_nu_v, h_cos_w, h_nu_w});
  else
    ROOT::RDF::RunGraphs({h_cos_all, h_nu_all});

  gStyle->SetOptStat(0);

  const std::string out_dir = plot_out_dir();
  const std::string fmt = plot_out_fmt();
  gSystem->mkdir(out_dir.c_str(), true);

  auto draw_one = [&](TH1D *hcos, TH1D *hnu, const std::string &tag) {
    normalise_to_percent(*hcos);
    normalise_to_percent(*hnu);

    style_hist(*hcos, kAzure + 1, 0.35);
    style_hist(*hnu, kOrange + 1, 0.35);

    const double ymax = 1.15 * std::max(hcos->GetMaximum(), hnu->GetMaximum());
    hcos->SetMaximum(ymax);
    hcos->SetTitle("");

    TCanvas c(("c_" + tag).c_str(), ("c_" + tag).c_str(), 900, 550);
    c.SetLogx();

    hcos->Draw("hist");
    hnu->Draw("hist same");

    TLegend leg(0.14, 0.84, 0.44, 0.96);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hcos, "Cosmic Pixels", "f");
    leg.AddEntry(hnu, "Neutrino Pixels", "f");
    leg.Draw();

    const std::string out = out_dir + "/" + tag + "." + fmt;
    c.SaveAs(out.c_str());
    std::cout << "[plotImageOccupancy] wrote " << out << "\n";
  };

  draw_one(&(*h_cos_all), &(*h_nu_all), "img_occ_cosmic_vs_neutrino_all");

  if (draw_planes)
  {
    draw_one(&(*h_cos_u), &(*h_nu_u), "img_occ_cosmic_vs_neutrino_u");
    draw_one(&(*h_cos_v), &(*h_nu_v), "img_occ_cosmic_vs_neutrino_v");
    draw_one(&(*h_cos_w), &(*h_nu_w), "img_occ_cosmic_vs_neutrino_w");
  }

  return 0;
}
