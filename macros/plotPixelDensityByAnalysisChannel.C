// plot/macro/plotPixelDensityByAnalysisChannel.C
//
// Plot pixel occupancy (active pixels / total pixels) vs analysis_channels.
// Produces a TProfile of mean occupancy per channel.
//
// Notes:
//   - analysis_channels must already exist in the event list.
//   - Occupancy is computed from (active_pixels_u/v/w) and detector_image sizes.
//   - By default it overlays MC and EXT (if present) and optionally Data.
//
// Run with:
//   ./heron macro plotPixelDensityByAnalysisChannel.C
//   ./heron macro plotPixelDensityByAnalysisChannel.C 'plotPixelDensityByAnalysisChannel("/path/to/event_list.root","true","",false)'
//   ./heron macro plotPixelDensityByAnalysisChannel.C 'plotPixelDensityByAnalysisChannel("/path/to/event_list.root","sel","w_nominal",true)'
//
// Args:
//   samples_tsv  : event_list root file path (default: default_event_list_root())
//   extra_sel    : additional selection string (default: "true")
//   weight_expr  : optional weight expression (default: "", meaning unweighted)
//   include_data : overlay data profile if true (default: false)

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
#include <TLegend.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
#include "PlottingHelper.hh"
#include "include/MacroGuard.hh"
#include "include/MacroEnv.hh"
#include "include/MacroColumns.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{



int require_columns(const std::unordered_set<std::string> &columns,
                    const std::vector<std::string> &required,
                    const std::string &label)
{
    return heron::macro::require_columns(columns, required, "plotPixelDensityByAnalysisChannel", label);
}

} // namespace

int plotPixelDensityByAnalysisChannel(const std::string &samples_tsv = "",
                                      const std::string &extra_sel = "true",
                                      const std::string &weight_expr = "",
                                      bool include_data = false)
{
  return heron::macro::run_with_guard("plotPixelDensityByAnalysisChannel", [&]() -> int {
    ROOT::EnableImplicitMT();
    std::cout << "[plotPixelDensityByAnalysisChannel] implicit MT enabled\n";

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotPixelDensityByAnalysisChannel] input=" << list_path << "\n";
    std::cout << "[plotPixelDensityByAnalysisChannel] extra_sel=" << extra_sel << "\n";
    std::cout << "[plotPixelDensityByAnalysisChannel] weight_expr=" << (weight_expr.empty() ? "(unweighted)" : weight_expr) << "\n";
    std::cout << "[plotPixelDensityByAnalysisChannel] include_data=" << (include_data ? "true" : "false") << "\n";

    if (!heron::macro::looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotPixelDensityByAnalysisChannel] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

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

    ROOT::RDF::RNode base = rdf;
    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc);
    if (mask_ext)
    {
        node_mc = node_mc.Filter([mask_ext](int sid) {
            return !(sid >= 0 &&
                     sid < static_cast<int>(mask_ext->size()) &&
                     (*mask_ext)[static_cast<size_t>(sid)]);
        },
        {"sample_id"});
    }
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

    if (!extra_sel.empty())
    {
        node_mc = node_mc.Filter(extra_sel);
        node_ext = node_ext.Filter(extra_sel);
        if (include_data)
            node_data = node_data.Filter(extra_sel);
    }

    const auto col_names = base.GetColumnNames();
    const std::unordered_set<std::string> columns(col_names.begin(), col_names.end());

    const std::vector<std::string> required = {
        "analysis_channels",
        "active_pixels_u",
        "active_pixels_v",
        "active_pixels_w",
        "detector_image_u",
        "detector_image_v",
        "detector_image_w",
        "sample_id"};

    if (require_columns(columns, required, "occupancy vs channels") != 0)
        return 1;

    auto add_occ = [&](ROOT::RDF::RNode n) {
        n = n.Define("n_pix_tot",
                     [](const ROOT::RVec<float> &u, const ROOT::RVec<float> &v, const ROOT::RVec<float> &w) {
                         return static_cast<long long>(u.size()) + static_cast<long long>(v.size()) + static_cast<long long>(w.size());
                     },
                     {"detector_image_u", "detector_image_v", "detector_image_w"});
        n = n.Define("n_active_tot", "active_pixels_u + active_pixels_v + active_pixels_w");
        n = n.Define("occ_pct",
                     [](int n_active, long long n_pix) {
                         if (n_pix <= 0)
                             return 0.0;
                         return 100.0 * static_cast<double>(n_active) / static_cast<double>(n_pix);
                     },
                     {"n_active_tot", "n_pix_tot"});

        if (!weight_expr.empty())
            n = n.Define("__w__", weight_expr);

        return n;
    };

    node_mc = add_occ(node_mc);
    node_ext = add_occ(node_ext);
    if (include_data)
        node_data = add_occ(node_data);

    const int nbins = 100;
    const double xmin = -0.5;
    const double xmax = 99.5;

    ROOT::RDF::TProfile1DModel pm_mc("p_mc", ";analysis_channels;Mean pixel occupancy [%]",
                                     nbins, xmin, xmax);
    ROOT::RDF::TProfile1DModel pm_ext("p_ext", ";analysis_channels;Mean pixel occupancy [%]",
                                      nbins, xmin, xmax);
    ROOT::RDF::TProfile1DModel pm_data("p_data", ";analysis_channels;Mean pixel occupancy [%]",
                                       nbins, xmin, xmax);

    ROOT::RDF::RResultPtr<TProfile> p_mc;
    ROOT::RDF::RResultPtr<TProfile> p_ext;
    ROOT::RDF::RResultPtr<TProfile> p_data;

    if (!weight_expr.empty())
    {
        p_mc = node_mc.Profile1D(pm_mc, "analysis_channels", "occ_pct", "__w__");
        p_ext = node_ext.Profile1D(pm_ext, "analysis_channels", "occ_pct", "__w__");
        if (include_data)
            p_data = node_data.Profile1D(pm_data, "analysis_channels", "occ_pct");
    }
    else
    {
        p_mc = node_mc.Profile1D(pm_mc, "analysis_channels", "occ_pct");
        p_ext = node_ext.Profile1D(pm_ext, "analysis_channels", "occ_pct");
        if (include_data)
            p_data = node_data.Profile1D(pm_data, "analysis_channels", "occ_pct");
    }

    if (include_data)
        ROOT::RDF::RunGraphs({p_mc, p_ext, p_data});
    else
        ROOT::RDF::RunGraphs({p_mc, p_ext});

    gStyle->SetOptStat(0);

    TProfile *pmc = &(*p_mc);
    TProfile *pext = &(*p_ext);

    pmc->SetLineColor(kBlue + 1);
    pmc->SetMarkerColor(kBlue + 1);
    pmc->SetMarkerStyle(20);

    pext->SetLineColor(kRed + 1);
    pext->SetMarkerColor(kRed + 1);
    pext->SetMarkerStyle(21);

    TCanvas c("c_occ_chan", "c_occ_chan", 1100, 600);
    pmc->SetTitle("");
    pmc->Draw("E1");
    pext->Draw("E1 same");

    TLegend leg(0.14, 0.84, 0.44, 0.96);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(pmc, "MC", "lp");
    leg.AddEntry(pext, "EXT", "lp");

    if (include_data)
    {
        TProfile *pdata = &(*p_data);
        pdata->SetLineColor(kBlack);
        pdata->SetMarkerColor(kBlack);
        pdata->SetMarkerStyle(22);
        pdata->Draw("E1 same");
        leg.AddEntry(pdata, "Data", "lp");
    }

    leg.Draw();

    const std::string out_dir = heron::macro::plot_out_dir();
    const std::string fmt = heron::macro::plot_out_fmt();
    gSystem->mkdir(out_dir.c_str(), true);

    const std::string out = out_dir + "/img_occ_pct_vs_analysis_channels." + fmt;
    c.SaveAs(out.c_str());
    std::cout << "[plotPixelDensityByAnalysisChannel] wrote " << out << "\n";

    return 0;

  });
}
