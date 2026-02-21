// plot/macro/plotInclusiveMuCCCutFlow.C
//
// Build a selection cut-flow plot for the inclusive νμ CC analysis.
//
// Default stages:
//   Stage 0: no cuts
//   Stage 1: sel_trigger
//   Stage 2: sel_slice
//   Stage 3: sel_fiducial
//   Stage 4: sel_muon
//
// The macro reports and plots:
//   - efficiency   = selected signal / total signal
//   - purity       = selected signal / selected (MC + EXT)
//   - MC purity    = selected signal / selected MC
//
// Run with:
//   ./heron macro plotInclusiveMuCCCutFlow.C
//   ./heron macro plotInclusiveMuCCCutFlow.C \
//     'plotInclusiveMuCCCutFlow("./scratch/out/event_list_myana.root")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TStyle.h>

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "SelectionService.hh"

using namespace nu;

namespace
{
struct CutFlowPoint
{
    std::string label;
    std::string expr;
    double efficiency = 0.0;
    double purity = 0.0;
    double mc_purity = 0.0;
};

std::shared_ptr<const std::vector<char>> build_truth_mc_mask(const EventListIO &el)
{
    const auto &refs = el.sample_refs();
    int max_sample_id = -1;
    for (const auto &kv : refs)
    {
        if (kv.first > max_sample_id)
            max_sample_id = kv.first;
    }

    auto mask = std::make_shared<std::vector<char>>(static_cast<size_t>(max_sample_id + 1), 0);

    const int overlay = static_cast<int>(SampleIO::SampleOrigin::kOverlay);
    const int dirt = static_cast<int>(SampleIO::SampleOrigin::kDirt);
    const int strangeness = static_cast<int>(SampleIO::SampleOrigin::kStrangeness);

    for (const auto &kv : refs)
    {
        const int sid = kv.first;
        const int origin = kv.second.sample_origin;

        if (sid < 0 || sid > max_sample_id)
            continue;

        (*mask)[static_cast<size_t>(sid)] =
            (origin == overlay || origin == dirt || origin == strangeness) ? 1 : 0;
    }

    return mask;
}

std::string cumulative_expr(const std::vector<std::string> &flags, std::size_t upto)
{
    if (upto == 0)
        return "true";

    std::string expr;
    for (std::size_t i = 0; i < upto; ++i)
    {
        if (!expr.empty())
            expr += " && ";
        expr += "(" + flags[i] + ")";
    }

    return expr.empty() ? std::string("true") : expr;
}

} // namespace

int plotInclusiveMuCCCutFlow(const std::string &event_list_path = "",
                             const std::string &signal_sel = "is_signal",
                             const std::string &mc_weight = "w_nominal",
                             const std::string &output_stem = "inclusive_mucc_cutflow")
{
    ROOT::EnableImplicitMT();

    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
    if (!looks_like_event_list_root(input_path))
    {
        std::cerr << "[plotInclusiveMuCCCutFlow] input is not an event-list root file: " << input_path << "\n";
        return 1;
    }

    const std::vector<std::string> cut_flags = {
        "sel_trigger",
        "sel_slice",
        "sel_fiducial",
        "sel_muon"};

    const std::vector<std::string> cut_labels = {
        "no cuts",
        "trigger",
        "slice",
        "fiducial",
        "inclusive #nu_{#mu} CC"};

    EventListIO el(input_path);
    ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf());

    auto mask_mc = build_truth_mc_mask(el);
    auto mask_ext = el.mask_for_ext();

    ROOT::RDF::RNode node_mc = filter_by_sample_mask(rdf, mask_mc).Define("__w__", mc_weight);
    ROOT::RDF::RNode node_ext = filter_by_sample_mask(rdf, mask_ext).Define("__w__", mc_weight);

    const double signal_total = *(node_mc.Filter(signal_sel).Sum<double>("__w__"));
    if (signal_total <= 0.0)
    {
        std::cerr << "[plotInclusiveMuCCCutFlow] signal denominator is <= 0 for signal_sel='"
                  << signal_sel << "'.\n";
        return 1;
    }

    std::vector<CutFlowPoint> points;
    points.reserve(cut_flags.size() + 1);

    std::cout << "\n[plotInclusiveMuCCCutFlow] Inclusive νμ CC selection cut-flow\n";
    std::cout << "stage\tlabel\tefficiency\tpurity\tmc_purity\n";

    for (std::size_t i = 0; i <= cut_flags.size(); ++i)
    {
        CutFlowPoint p;
        p.label = cut_labels[i];
        p.expr = cumulative_expr(cut_flags, i);

        const std::string stage_sel = "(" + p.expr + ")";
        const std::string signal_and_stage = "(" + signal_sel + ") && " + stage_sel;

        const double signal_pass = *(node_mc.Filter(signal_and_stage).Sum<double>("__w__"));
        const double selected_mc = *(node_mc.Filter(stage_sel).Sum<double>("__w__"));
        const double selected_ext = *(node_ext.Filter(stage_sel).Sum<double>("__w__"));
        const double selected_all = selected_mc + selected_ext;

        p.efficiency = signal_pass / signal_total;
        p.purity = (selected_all > 0.0) ? signal_pass / selected_all : 0.0;
        p.mc_purity = (selected_mc > 0.0) ? signal_pass / selected_mc : 0.0;

        points.push_back(p);

        std::cout << i << "\t" << p.label << "\t"
                  << p.efficiency << "\t"
                  << p.purity << "\t"
                  << p.mc_purity << "\n";
    }

    const int n = static_cast<int>(points.size());
    std::vector<double> x(n), y_eff(n), y_pur(n), y_mc_pur(n);
    for (int i = 0; i < n; ++i)
    {
        x[i] = static_cast<double>(i + 1);
        y_eff[i] = points[static_cast<std::size_t>(i)].efficiency;
        y_pur[i] = points[static_cast<std::size_t>(i)].purity;
        y_mc_pur[i] = points[static_cast<std::size_t>(i)].mc_purity;
    }

    TH1D h_axis("h_axis", ";selection stage;efficiency or purity", n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i)
        h_axis.GetXaxis()->SetBinLabel(i + 1, points[static_cast<std::size_t>(i)].label.c_str());

    TGraph g_eff(n, x.data(), y_eff.data());
    TGraph g_pur(n, x.data(), y_pur.data());
    TGraph g_mc_pur(n, x.data(), y_mc_pur.data());

    g_eff.SetLineColor(kBlue + 1);
    g_eff.SetMarkerColor(kBlue + 1);
    g_eff.SetLineWidth(2);
    g_eff.SetMarkerStyle(20);

    g_pur.SetLineColor(kRed + 1);
    g_pur.SetMarkerColor(kRed + 1);
    g_pur.SetLineWidth(2);
    g_pur.SetMarkerStyle(20);

    g_mc_pur.SetLineColor(kBlack);
    g_mc_pur.SetMarkerColor(kBlack);
    g_mc_pur.SetLineWidth(2);
    g_mc_pur.SetMarkerStyle(20);

    TCanvas c("c_inclusive_mucc_cutflow", "Inclusive #nu_{#mu} CC cut-flow", 1000, 700);
    gStyle->SetOptStat(0);
    c.SetBottomMargin(0.22);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.07);

    h_axis.SetMinimum(0.0);
    h_axis.SetMaximum(1.05);
    h_axis.GetXaxis()->LabelsOption("v");
    h_axis.Draw("AXIS");

    g_eff.Draw("LP SAME");
    g_pur.Draw("LP SAME");
    g_mc_pur.Draw("LP SAME");

    TLegend leg(0.62, 0.18, 0.88, 0.34);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&g_eff, "efficiency", "lp");
    leg.AddEntry(&g_pur, "purity", "lp");
    leg.AddEntry(&g_mc_pur, "MC purity", "lp");
    leg.Draw();

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());

    std::cout << "\n[plotInclusiveMuCCCutFlow] saved plot: " << out << "\n";

    return 0;
}
