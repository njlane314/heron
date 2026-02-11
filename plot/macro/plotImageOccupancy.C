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
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TSystem.h>

#include "ColumnDerivationService.hh"

#include "EventListIO.hh"
#include "Plotter.hh"
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

  gStyle->SetOptStat(0);

  const std::string out_dir = plot_out_dir();
  const std::string fmt = plot_out_fmt();
  gSystem->mkdir(out_dir.c_str(), true);

  const double n_events = static_cast<double>(*n.Count());
  const double evt_weight = (n_events > 0.0) ? (100.0 / n_events) : 0.0;

  auto draw_one = [&](const std::string &cos_col,
                      const std::string &nu_col,
                      const std::string &tag) {
    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;

    std::vector<Entry> entries;
    entries.reserve(2);

    auto cosmic_node = n
                           .Define("plot_occ_value", cos_col)
                           .Define("plot_occ_weight", [evt_weight]() { return evt_weight; })
                           .Define("plot_occ_channel", "97");
    entries.emplace_back(make_entry(std::move(cosmic_node), rec_mc));

    auto neutrino_node = n
                             .Define("plot_occ_value", nu_col)
                             .Define("plot_occ_weight", [evt_weight]() { return evt_weight; })
                             .Define("plot_occ_channel", "98");
    entries.emplace_back(make_entry(std::move(neutrino_node), rec_mc));

    std::vector<const Entry *> mc;
    mc.reserve(entries.size());
    for (auto &e : entries)
      mc.push_back(&e);

    TH1DModel spec;
    spec.id = tag;
    spec.name = tag;
    spec.title = ";Percentage of pixels / event;Fraction of Dataset [%]";
    spec.expr = "plot_occ_value";
    spec.weight = "plot_occ_weight";
    spec.nbins = nbins;
    spec.xmin = xmin_pct;
    spec.xmax = xmax_pct;
    spec.sel = Preset::Empty;

    Plotter plotter;
    auto &opt = plotter.options();
    opt.out_dir = out_dir;
    opt.image_format = fmt;
    opt.show_ratio = false;
    opt.show_ratio_band = false;
    opt.annotate_numbers = false;
    opt.overlay_signal = false;
    opt.legend_on_top = false;
    opt.show_legend = true;
    opt.use_log_y = false;
    opt.channel_column = "plot_occ_channel";
    opt.unstack_channel_keys = {97, 98};
    opt.unstack_channel_labels = {
        {97, "Cosmic Pixels"},
        {98, "Neutrino Pixels"}
    };
    opt.unstack_channel_colours = {
        {97, kAzure + 1},
        {98, kOrange + 1}
    };

    plotter.draw_unstack(spec, mc);
  };

  draw_one("cosmic_pct_tot", "nu_pct_tot", "img_occ_cosmic_vs_neutrino_all");

  if (draw_planes)
  {
    draw_one("cosmic_pct_u", "nu_pct_u", "img_occ_cosmic_vs_neutrino_u");
    draw_one("cosmic_pct_v", "nu_pct_v", "img_occ_cosmic_vs_neutrino_v");
    draw_one("cosmic_pct_w", "nu_pct_w", "img_occ_cosmic_vs_neutrino_w");
  }

  return 0;
}
