// plot/macro/plotSemanticPixelDensity.C
//
// Plot per-event pixel density (percentage of pixels / event) for each semantic label.
//
// Definition (per event, per label i):
//   pct_i = 100 * (count_i_u + count_i_v + count_i_w) / (|det_u| + |det_v| + |det_w|)
//
// Uses slice_semantic_active_pixels_{u,v,w} counts, so it does not need to scan
// per-pixel semantic_image arrays.
//
// Run with:
//   ./nuxsec macro plotSemanticPixelDensity.C
//   ./nuxsec macro plotSemanticPixelDensity.C 'plotSemanticPixelDensity("/path/to/event_list.root","true")'
//
// Output:
//   One plot per label, saved to $NUXSEC_PLOT_DIR in $NUXSEC_PLOT_FORMAT.

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

std::string sanitise(const std::string &s)
{
  std::string out;
  out.reserve(s.size());
  for (unsigned char c : s)
  {
    if ((c >= 'a' && c <= 'z') ||
        (c >= 'A' && c <= 'Z') ||
        (c >= '0' && c <= '9'))
      out.push_back(static_cast<char>(c));
    else
      out.push_back('_');
  }
  return out;
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

  std::cerr << "[plotSemanticPixelDensity] missing required columns for " << label << ":\n";
  for (const auto &name : missing)
    std::cerr << "  - " << name << "\n";
  return 1;
}

int at_or_zero(const std::vector<int> &v, int idx)
{
  if (idx < 0)
    return 0;
  const size_t u = static_cast<size_t>(idx);
  if (u >= v.size())
    return 0;
  return v[u];
}

} // namespace

int plotSemanticPixelDensity(const std::string &samples_tsv = "",
                             const std::string &extra_sel = "true",
                             int nbins = 60,
                             double xmin_pct = 1e-4,
                             double xmax_pct = 1e1)
{
  ROOT::EnableImplicitMT();
  std::cout << "[plotSemanticPixelDensity] implicit MT enabled\n";

  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plotSemanticPixelDensity] input=" << list_path << "\n";
  std::cout << "[plotSemanticPixelDensity] extra_sel=" << extra_sel << "\n";

  if (!looks_like_event_list_root(list_path))
  {
    std::cerr << "[plotSemanticPixelDensity] input is not an event list root file: " << list_path << "\n";
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

  if (require_columns(columns, required, "semantic density") != 0)
    return 1;

  // Default label names matching your SemanticClassifier.
  std::vector<std::string> labels = {
      "Empty",
      "Cosmic",
      "Muon",
      "Electron",
      "Photon",
      "ChargedPion",
      "NeutralPion",
      "Neutron",
      "Proton",
      "ChargedKaon",
      "NeutralKaon",
      "Lambda",
      "ChargedSigma",
      "NeutralSigma",
      "Other",
      "Ambiguous"};

  // If semantic_label_names exists, prefer its first-entry contents.
  if (columns.find("semantic_label_names") != columns.end())
  {
    auto first = node.Range(1).Take<std::vector<std::string>>("semantic_label_names").GetValue();
    if (!first.empty() && !first[0].empty())
    {
      labels = first[0];
      std::cout << "[plotSemanticPixelDensity] using semantic_label_names from tree (n=" << labels.size() << ")\n";
    }
  }

  const int nlabels = static_cast<int>(labels.size());
  if (nlabels <= 0)
  {
    std::cerr << "[plotSemanticPixelDensity] no semantic labels found\n";
    return 1;
  }

  auto n = node.Define("n_pix_tot",
                       [](const std::vector<float> &u, const std::vector<float> &v, const std::vector<float> &w) {
                         return static_cast<long long>(u.size()) + static_cast<long long>(v.size()) + static_cast<long long>(w.size());
                       },
                       {"detector_image_u", "detector_image_v", "detector_image_w"});

  // Define one percentage column per label and book all histograms.
  const auto edges = log_edges(nbins, xmin_pct, xmax_pct);
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hists;
  hists.reserve(static_cast<size_t>(nlabels));

  for (int i = 0; i < nlabels; ++i)
  {
    const std::string col = "sem_pct_" + std::to_string(i);
    n = n.Define(
        col,
        [i](const std::vector<int> &cu,
            const std::vector<int> &cv,
            const std::vector<int> &cw,
            long long n_pix_tot) {
          if (n_pix_tot <= 0)
            return 0.0;
          const long long cnt =
              static_cast<long long>(at_or_zero(cu, i)) +
              static_cast<long long>(at_or_zero(cv, i)) +
              static_cast<long long>(at_or_zero(cw, i));
          return 100.0 * static_cast<double>(cnt) / static_cast<double>(n_pix_tot);
        },
        {"slice_semantic_active_pixels_u", "slice_semantic_active_pixels_v", "slice_semantic_active_pixels_w", "n_pix_tot"});

    const std::string hname = "h_" + col;
    ROOT::RDF::TH1DModel model(hname.c_str(),
                               ";Percentage of pixels / event;Fraction of Dataset [%]",
                               nbins, edges.data());

    hists.push_back(n.Histo1D(model, col));
  }

  // Single pass over the dataset for all label hists.
  std::vector<ROOT::RDF::RResultHandle> handles;
  handles.reserve(hists.size());
  for (const auto &hist : hists)
    handles.emplace_back(hist);
  ROOT::RDF::RunGraphs(handles);

  gStyle->SetOptStat(0);

  const std::string out_dir = plot_out_dir();
  const std::string fmt = plot_out_fmt();
  gSystem->mkdir(out_dir.c_str(), true);

  for (int i = 0; i < nlabels; ++i)
  {
    TH1D *h = &(*hists[static_cast<size_t>(i)]);
    normalise_to_percent(*h);

    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->SetTitle("");

    TCanvas c(("c_sem_" + std::to_string(i)).c_str(),
              ("c_sem_" + std::to_string(i)).c_str(),
              900, 550);
    c.SetLogx();

    h->Draw("hist");

    TLegend leg(0.14, 0.84, 0.70, 0.96);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    const std::string entry = "Label: " + labels[static_cast<size_t>(i)];
    leg.AddEntry(h, entry.c_str(), "l");
    leg.Draw();

    const std::string tag = "img_sem_density_" + std::to_string(i) + "_" + sanitise(labels[static_cast<size_t>(i)]);
    const std::string out = out_dir + "/" + tag + "." + fmt;
    c.SaveAs(out.c_str());
    std::cout << "[plotSemanticPixelDensity] wrote " << out << "\n";
  }

  return 0;
}
