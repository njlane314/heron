// plot/macro/plotSelectionEvolutionAndTable.C
//
// Build a stage-by-stage selection-evolution plot (efficiency, purity, MC purity)
// and print a table suitable for direct copy into notes/papers.
//
// Default stage sequence follows the inclusive νμ CC selection flags in this codebase:
//   Stage 0: no cuts
//   Stage 1: sel_trigger
//   Stage 2: ... && sel_slice
//   Stage 3: ... && sel_fiducial
//   Stage 4: ... && sel_muon
//
// Run with:
//   ./heron macro plotSelectionEvolutionAndTable.C
//   ./heron macro plotSelectionEvolutionAndTable.C \
//     'plotSelectionEvolutionAndTable("./scratch/out/event_list_myana.root", "is_signal")'
//
//   // Custom cut sequence (comma-separated cumulative flags and labels)
//   ./heron macro plotSelectionEvolutionAndTable.C \
//     'plotSelectionEvolutionAndTable("./scratch/out/event_list_myana.root", "is_signal", "sel_trigger,sel_slice,sel_fiducial,sel_muon", "trigger,slice,fiducial,muon")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TPad.h>
#include <TStyle.h>

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"
#include "SelectionService.hh"

using namespace nu;

namespace
{
struct StageStats
{
    std::string cut_name;
    std::string stage_flag;
    std::string expr;
    double efficiency = 0.0;
    double purity = 0.0;
    double mc_purity = 0.0;
    double rel_efficiency = 0.0;
};

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

std::vector<std::string> split_csv(const std::string &csv)
{
    std::vector<std::string> out;
    std::stringstream ss(csv);
    std::string token;
    while (std::getline(ss, token, ','))
    {
        const auto ltrim = token.find_first_not_of(" \t\n\r");
        if (ltrim == std::string::npos)
            continue;
        const auto rtrim = token.find_last_not_of(" \t\n\r");
        out.push_back(token.substr(ltrim, rtrim - ltrim + 1));
    }
    return out;
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

std::string tex_number(double x)
{
    if (!std::isfinite(x))
        return "\\infty";

    std::ostringstream os;
    os << std::fixed << std::setprecision(4) << x;
    return os.str();
}

} // namespace

int plotSelectionEvolutionAndTable(const std::string &event_list_path = "",
                                   const std::string &signal_sel = "is_signal",
                                   const std::string &cuts_csv = "sel_trigger,sel_slice,sel_fiducial,sel_muon",
                                   const std::string &labels_csv = "trigger,slice,fiducial,muon",
                                   const std::string &mc_weight = "w_nominal",
                                   const std::string &output_stem = "selection_evolution")
{
    ROOT::EnableImplicitMT();

    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;

    if (signal_sel.find("analysis_channels") != std::string::npos)
    {
        std::cout << "[plotSelectionEvolutionAndTable] note: signal_sel uses 'analysis_channels'. "
                  << "For truth-level signal efficiency prefer 'is_signal' (or an explicit truth selector).\n";
    }

    if (!looks_like_event_list_root(input_path))
    {
        std::cerr << "[plotSelectionEvolutionAndTable] input is not an event-list root file: " << input_path << "\n";
        return 1;
    }

    const std::vector<std::string> cut_flags = split_csv(cuts_csv);
    std::vector<std::string> cut_labels = split_csv(labels_csv);

    if (cut_flags.empty())
    {
        std::cerr << "[plotSelectionEvolutionAndTable] cuts_csv is empty.\n";
        return 1;
    }

    if (cut_labels.size() < cut_flags.size())
    {
        for (std::size_t i = cut_labels.size(); i < cut_flags.size(); ++i)
            cut_labels.push_back(cut_flags[i]);
    }

    EventListIO el(input_path);
    ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf());

    auto mask_mc = build_truth_mc_mask(el);
    auto mask_ext = el.mask_for_ext();

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) {
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode node_mc = filter_by_mask(rdf, mask_mc).Define("__w__", mc_weight);
    ROOT::RDF::RNode node_ext = filter_by_mask(rdf, mask_ext).Define("__w__", mc_weight);

    const double signal_total = *(node_mc.Filter(signal_sel).Sum<double>("__w__"));
    if (signal_total <= 0.0)
    {
        std::cerr << "[plotSelectionEvolutionAndTable] signal denominator is <= 0 for signal_sel='" << signal_sel << "'.\n";
        return 1;
    }

    std::vector<StageStats> rows;
    rows.reserve(cut_flags.size() + 1);

    for (std::size_t i = 0; i <= cut_flags.size(); ++i)
    {
        StageStats r;
        r.cut_name = (i == 0) ? "no cuts" : cut_labels[i - 1];
        r.stage_flag = (i == 0) ? "true" : cut_flags[i - 1];
        r.expr = cumulative_expr(cut_flags, i);

        const std::string stage_sel = "(" + r.expr + ")";
        const std::string signal_and_stage = "(" + signal_sel + ") && " + stage_sel;

        const double signal_pass = *(node_mc.Filter(signal_and_stage).Sum<double>("__w__"));
        const double selected_mc = *(node_mc.Filter(stage_sel).Sum<double>("__w__"));
        const double selected_ext = *(node_ext.Filter(stage_sel).Sum<double>("__w__"));

        const double selected_all = selected_mc + selected_ext;

        r.efficiency = signal_pass / signal_total;
        r.purity = (selected_all > 0.0) ? signal_pass / selected_all : 0.0;
        r.mc_purity = (selected_mc > 0.0) ? signal_pass / selected_mc : 0.0;

        if (i == 0)
            r.rel_efficiency = std::numeric_limits<double>::infinity();
        else
            r.rel_efficiency = (rows.back().efficiency > 0.0) ? r.efficiency / rows.back().efficiency : 0.0;

        rows.push_back(r);
    }

    // Print an easy-to-copy markdown table.
    std::cout << "\n[plotSelectionEvolutionAndTable] Selection evolution table\n";
    std::cout << "| Cut | Stage flag | Applied selection | Efficiency | Purity | MC purity | Rel. eff. |\n";
    std::cout << "|---|---|---|---:|---:|---:|---:|\n";
    for (const auto &r : rows)
    {
        std::cout << "| " << r.cut_name
                  << " | `" << r.stage_flag
                  << "` | `" << r.expr
                  << "` | " << tex_number(r.efficiency)
                  << " | " << tex_number(r.purity)
                  << " | " << tex_number(r.mc_purity)
                  << " | " << tex_number(r.rel_efficiency)
                  << " |\n";
    }

    // Plot evolution curves.
    const int n = static_cast<int>(rows.size());
    std::vector<double> x(n), y_eff(n), y_pur(n), y_mcp(n);

    for (int i = 0; i < n; ++i)
    {
        x[i] = static_cast<double>(i + 1);
        y_eff[i] = rows[static_cast<std::size_t>(i)].efficiency;
        y_pur[i] = rows[static_cast<std::size_t>(i)].purity;
        y_mcp[i] = rows[static_cast<std::size_t>(i)].mc_purity;
    }

    TH1D haxis("haxis", ";cut stage;efficiency", n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i)
        haxis.GetXaxis()->SetBinLabel(i + 1, rows[static_cast<std::size_t>(i)].cut_name.c_str());

    TGraph geff(n, x.data(), y_eff.data());
    TGraph gpur(n, x.data(), y_pur.data());
    TGraph gmcp(n, x.data(), y_mcp.data());

    geff.SetLineColor(kBlue + 1);
    geff.SetMarkerColor(kBlue + 1);
    geff.SetMarkerStyle(20);
    geff.SetLineWidth(2);

    gpur.SetLineColor(kRed + 1);
    gpur.SetMarkerColor(kRed + 1);
    gpur.SetMarkerStyle(20);
    gpur.SetLineWidth(2);

    gmcp.SetLineColor(kBlack);
    gmcp.SetMarkerColor(kBlack);
    gmcp.SetMarkerStyle(20);
    gmcp.SetLineWidth(2);

    TCanvas c("c_selection_evolution", "Selection evolution", 1000, 700);
    gStyle->SetOptStat(0);
    c.SetBottomMargin(0.18);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.10);
    c.SetTopMargin(0.07);

    haxis.SetMinimum(0.0);
    haxis.SetMaximum(1.02);
    haxis.GetXaxis()->LabelsOption("v");
    haxis.GetXaxis()->SetLabelSize(0.038);
    haxis.GetXaxis()->SetTitleOffset(1.20);
    haxis.GetYaxis()->SetTitleOffset(1.2);
    haxis.Draw("AXIS");

    geff.Draw("LP SAME");

    double purity_min = std::numeric_limits<double>::infinity();
    double purity_max = 0.0;
    for (int i = 0; i < n; ++i)
    {
        if (y_pur[static_cast<std::size_t>(i)] > 0.0)
        {
            purity_min = std::min(purity_min, y_pur[static_cast<std::size_t>(i)]);
            purity_max = std::max(purity_max, y_pur[static_cast<std::size_t>(i)]);
        }

        if (y_mcp[static_cast<std::size_t>(i)] > 0.0)
        {
            purity_min = std::min(purity_min, y_mcp[static_cast<std::size_t>(i)]);
            purity_max = std::max(purity_max, y_mcp[static_cast<std::size_t>(i)]);
        }
    }

    if (!std::isfinite(purity_min) || purity_min <= 0.0 || purity_max <= 0.0)
    {
        purity_min = 1e-6;
        purity_max = 1.0;
    }
    else
    {
        purity_min = std::pow(10.0, std::floor(std::log10(purity_min)));
        purity_max = std::pow(10.0, std::ceil(std::log10(purity_max)));
        if (!(purity_max > purity_min))
            purity_max = 10.0 * purity_min;
    }

    TPad p_purity("p_purity", "p_purity", 0.0, 0.0, 1.0, 1.0);
    p_purity.SetFillStyle(4000);
    p_purity.SetFrameFillStyle(0);
    p_purity.SetFrameLineColor(0);
    p_purity.SetBottomMargin(c.GetBottomMargin());
    p_purity.SetLeftMargin(c.GetLeftMargin());
    p_purity.SetRightMargin(c.GetRightMargin());
    p_purity.SetTopMargin(c.GetTopMargin());
    p_purity.SetLogy();
    p_purity.Draw();
    p_purity.cd();

    TH1D haxis_purity("haxis_purity", ";;purity", n, 0.5, n + 0.5);
    haxis_purity.SetMinimum(purity_min);
    haxis_purity.SetMaximum(purity_max);
    haxis_purity.GetXaxis()->SetLabelSize(0.0);
    haxis_purity.GetXaxis()->SetTickLength(0.0);
    haxis_purity.GetYaxis()->SetTitleOffset(1.2);
    haxis_purity.GetYaxis()->SetLabelSize(haxis.GetYaxis()->GetLabelSize());
    haxis_purity.Draw("AXIS Y+");

    gpur.Draw("LP SAME");
    gmcp.Draw("LP SAME");

    c.cd();

    TLegend leg(0.16, 0.92, 0.90, 0.985);
    leg.SetNColumns(3);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.042);
    leg.AddEntry(&geff, "efficiency", "lp");
    leg.AddEntry(&gpur, "purity", "lp");
    leg.AddEntry(&gmcp, "MC purity", "lp");
    leg.Draw();

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());

    std::cout << "\n[plotSelectionEvolutionAndTable] saved plot: " << out << "\n";
    std::cout << "[plotSelectionEvolutionAndTable] signal_sel='" << signal_sel << "'  mc_weight='" << mc_weight << "'\n";

    return 0;
}
