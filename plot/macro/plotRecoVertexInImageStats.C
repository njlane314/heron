// plot/macro/plotRecoVertexInImageStats.C
//
// Summarise and plot how often the reconstructed neutrino vertex falls inside
// the stored slice images, using the boolean branches:
//   - is_vtx_in_image_u/v/w
//
// Produces a fraction-per-plane plot (with binomial error bars).
//
// Run with:
//   ./heron macro plotRecoVertexInImageStats.C
//   ./heron macro plotRecoVertexInImageStats.C 'plotRecoVertexInImageStats("/path/to/event_list.root","true",true)'
//
// Args:
//   samples_tsv   : event_list root file (default: default_event_list_root())
//   extra_sel     : additional selection (default: "true")
//   include_data  : include data in overlay (default: false)

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
#include "SampleCLI.hh"

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

    std::cerr << "[plotRecoVertexInImageStats] missing required columns for " << label << ":\n";
    for (const auto &name : missing)
        std::cerr << "  - " << name << "\n";
    return 1;
}

struct Counts
{
    long long n_tot = 0;
    long long n_u = 0;
    long long n_v = 0;
    long long n_w = 0;
    long long n_all3 = 0;
};

Counts compute_counts(ROOT::RDF::RNode node)
{
    auto a_tot = node.Count();
    auto a_u = node.Filter("is_vtx_in_image_u").Count();
    auto a_v = node.Filter("is_vtx_in_image_v").Count();
    auto a_w = node.Filter("is_vtx_in_image_w").Count();
    auto a_all3 = node.Filter("is_vtx_in_image_u && is_vtx_in_image_v && is_vtx_in_image_w").Count();

    ROOT::RDF::RunGraphs({a_tot, a_u, a_v, a_w, a_all3});

    Counts c;
    c.n_tot = a_tot.GetValue();
    c.n_u = a_u.GetValue();
    c.n_v = a_v.GetValue();
    c.n_w = a_w.GetValue();
    c.n_all3 = a_all3.GetValue();
    return c;
}

void fill_frac_hist(TH1D &h, const Counts &c)
{
    auto set = [&](int bin, long long n_pass) {
        if (c.n_tot <= 0)
        {
            h.SetBinContent(bin, 0.0);
            h.SetBinError(bin, 0.0);
            return;
        }
        const double p = static_cast<double>(n_pass) / static_cast<double>(c.n_tot);
        const double err = std::sqrt(std::max(0.0, p * (1.0 - p) / static_cast<double>(c.n_tot)));
        h.SetBinContent(bin, 100.0 * p);
        h.SetBinError(bin, 100.0 * err);
    };

    set(1, c.n_u);
    set(2, c.n_v);
    set(3, c.n_w);
    set(4, c.n_all3);
}

} // namespace

int plotRecoVertexInImageStats(const std::string &samples_tsv = "",
                               const std::string &extra_sel = "true",
                               bool include_data = false)
{
    ROOT::EnableImplicitMT();
    std::cout << "[plotRecoVertexInImageStats] implicit MT enabled\n";

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotRecoVertexInImageStats] input=" << list_path << "\n";
    std::cout << "[plotRecoVertexInImageStats] extra_sel=" << extra_sel << "\n";
    std::cout << "[plotRecoVertexInImageStats] include_data=" << (include_data ? "true" : "false") << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotRecoVertexInImageStats] input is not an event list root file: " << list_path << "\n";
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
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0 &&
                                                sid < static_cast<int>(mask_ext->size()) &&
                                                (*mask_ext)[static_cast<size_t>(sid)]);
                                   },
                                           {"sample_id"});
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
        "is_vtx_in_image_u",
        "is_vtx_in_image_v",
        "is_vtx_in_image_w",
        "sample_id"};

    if (require_columns(columns, required, "vertex-in-image") != 0)
        return 1;

    // Compute counts and print.
    const Counts c_mc = compute_counts(node_mc);
    const Counts c_ext = compute_counts(node_ext);
    Counts c_data;
    if (include_data)
        c_data = compute_counts(node_data);

    auto print_one = [&](const char *label, const Counts &c) {
        std::cout << "\n[plotRecoVertexInImageStats] " << label << ":\n";
        std::cout << "  total=" << c.n_tot << "\n";
        if (c.n_tot > 0)
        {
            std::cout << "  U    = " << c.n_u << " (" << (100.0 * c.n_u / c.n_tot) << "%)\n";
            std::cout << "  V    = " << c.n_v << " (" << (100.0 * c.n_v / c.n_tot) << "%)\n";
            std::cout << "  W    = " << c.n_w << " (" << (100.0 * c.n_w / c.n_tot) << "%)\n";
            std::cout << "  all3 = " << c.n_all3 << " (" << (100.0 * c.n_all3 / c.n_tot) << "%)\n";
        }
    };

    print_one("MC", c_mc);
    print_one("EXT", c_ext);
    if (include_data)
        print_one("Data", c_data);

    // Build fraction hists for plotting.
    TH1D h_mc("h_mc", ";Plane;Events with reco #nu vtx in image [%]", 4, 0.5, 4.5);
    TH1D h_ext("h_ext", ";Plane;Events with reco #nu vtx in image [%]", 4, 0.5, 4.5);
    TH1D h_data("h_data", ";Plane;Events with reco #nu vtx in image [%]", 4, 0.5, 4.5);

    h_mc.GetXaxis()->SetBinLabel(1, "U");
    h_mc.GetXaxis()->SetBinLabel(2, "V");
    h_mc.GetXaxis()->SetBinLabel(3, "W");
    h_mc.GetXaxis()->SetBinLabel(4, "U&V&W");

    h_ext.GetXaxis()->SetBinLabel(1, "U");
    h_ext.GetXaxis()->SetBinLabel(2, "V");
    h_ext.GetXaxis()->SetBinLabel(3, "W");
    h_ext.GetXaxis()->SetBinLabel(4, "U&V&W");

    h_data.GetXaxis()->SetBinLabel(1, "U");
    h_data.GetXaxis()->SetBinLabel(2, "V");
    h_data.GetXaxis()->SetBinLabel(3, "W");
    h_data.GetXaxis()->SetBinLabel(4, "U&V&W");

    fill_frac_hist(h_mc, c_mc);
    fill_frac_hist(h_ext, c_ext);
    if (include_data)
        fill_frac_hist(h_data, c_data);

    gStyle->SetOptStat(0);

    h_mc.SetLineColor(kBlue + 1);
    h_mc.SetMarkerColor(kBlue + 1);
    h_mc.SetMarkerStyle(20);

    h_ext.SetLineColor(kRed + 1);
    h_ext.SetMarkerColor(kRed + 1);
    h_ext.SetMarkerStyle(21);

    if (include_data)
    {
        h_data.SetLineColor(kBlack);
        h_data.SetMarkerColor(kBlack);
        h_data.SetMarkerStyle(22);
    }

    TCanvas c("c_vtx_in_img", "c_vtx_in_img", 900, 550);
    h_mc.SetTitle("");
    h_mc.SetMinimum(0.0);
    h_mc.Draw("E1");
    h_ext.Draw("E1 same");
    if (include_data)
        h_data.Draw("E1 same");

    TLegend leg(0.14, 0.84, 0.44, 0.96);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&h_mc, "MC", "lp");
    leg.AddEntry(&h_ext, "EXT", "lp");
    if (include_data)
        leg.AddEntry(&h_data, "Data", "lp");
    leg.Draw();

    const std::string out_dir = plot_out_dir();
    const std::string fmt = plot_out_fmt();
    gSystem->mkdir(out_dir.c_str(), true);

    const std::string out = out_dir + "/reco_vtx_in_image_fractions." + fmt;
    c.SaveAs(out.c_str());
    std::cout << "[plotRecoVertexInImageStats] wrote " << out << "\n";

    return 0;
}
