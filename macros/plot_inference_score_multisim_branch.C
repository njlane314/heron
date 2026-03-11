#ifndef PLOT_INFERENCE_SCORE_MULTISIM_BRANCH_CXX
#define PLOT_INFERENCE_SCORE_MULTISIM_BRANCH_CXX

// macros/plot_inference_score_multisim_branch.C
//
// Generic helper for multisim vector branches such as:
//   - weightsGenie
//   - weightsFlux
//   - weightsReint
//   - weightsPPFX
//
// If the direct vector branch is missing but the weights map is present,
// the macro can also fall back to named map keys such as:
//   - All_UBGenie
//   - reint_all
//   - ppfx_all
//   - flux_all
//
// It builds the total inference-score[0] histogram for:
//   - nominal total = MC(nominal) + EXT(nominal)        [if include_ext=true]
//   - multisim universes = MC(w_nominal * universe_w) + EXT(nominal)
//
// The plot shows:
//   - nominal total histogram
//   - mean universe histogram
//   - mean +/- RMS band across universes
//   - min/max envelope across universes
//   - a small overlay subset of individual universes
//
// The vector branches are stored as std::vector<unsigned short> with an
// internal scale factor of 1000, so each universe weight is decoded as
//   w = 0.001 * weightsBranch[i]
// and any non-positive entry is treated as 1.0, matching the sanitising
// convention used elsewhere.
//
// Typical use:
//   ./heron macro plot_inference_score_multisim_branch.C \
//     'plot_inference_score_multisim_branch("./scratch/out/event_list.root", "sel_muon", "w_nominal", "weightsGenie", "GENIE")'
//
// Wrapper entry points are provided at the bottom of this file:
//   plot_inference_score_genie_universes
//   plot_inference_score_flux_universes
//   plot_inference_score_reint_universes
//   plot_inference_score_ppfx_universes

#include <algorithm>
#include <cmath>
#include <cctype>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

ROOT::RDF::RNode apply_selection(ROOT::RDF::RNode node, const std::string &sel)
{
    if (sel.empty())
        return node;

    if (has_column(node, sel))
        return node.Filter([](bool pass) { return pass; }, {sel});

    return node.Filter(sel);
}

ROOT::RDF::RNode define_score0(ROOT::RDF::RNode node)
{
    if (has_column(node, "inf_score_0"))
    {
        return node.Define(
            "score0_for_plot",
            [](double x)
            {
                return std::isfinite(x) ? x : -1.0e9;
            },
            {"inf_score_0"});
    }

    if (has_column(node, "inf_scores"))
    {
        return node.Define(
            "score0_for_plot",
            [](const ROOT::RVec<float> &scores)
            {
                return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
            },
            {"inf_scores"});
    }

    if (has_column(node, "inference_score"))
    {
        return node.Define(
            "score0_for_plot",
            [](const ROOT::RVec<float> &scores)
            {
                return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
            },
            {"inference_score"});
    }

    if (has_column(node, "inference_scores"))
    {
        return node.Define(
            "score0_for_plot",
            [](const ROOT::RVec<float> &scores)
            {
                return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
            },
            {"inference_scores"});
    }

    return node;
}

std::string make_safe_name(std::string s)
{
    for (char &c : s)
    {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (!(std::isalnum(uc) || c == '_'))
            c = '_';
    }
    return s;
}

std::unique_ptr<TH1D> clone_hist(const TH1D &h, const std::string &name)
{
    auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(name.c_str())));
    out->SetDirectory(nullptr);
    return out;
}

void style_line_hist(TH1D &h, int color, int style = 1, int width = 2)
{
    h.SetDirectory(nullptr);
    h.SetLineColor(color);
    h.SetMarkerColor(color);
    h.SetLineStyle(style);
    h.SetLineWidth(width);
    h.SetMarkerSize(0.0);
    h.SetFillStyle(0);
}

void style_band_hist(TH1D &h, int line_color, int fill_color)
{
    h.SetDirectory(nullptr);
    h.SetLineColor(line_color);
    h.SetLineWidth(1);
    h.SetFillColor(fill_color);
    h.SetFillStyle(1001);
    h.SetMarkerSize(0.0);
}

double sanitise_multisim_weight(double w)
{
    if (!std::isfinite(w) || w <= 0.0)
        return 1.0;
    return w;
}

template <class VecT>
double universe_weight_from_vec(const VecT &ws, int idx)
{
    if (idx < 0)
        return 1.0;
    const auto n = static_cast<int>(ws.size());
    if (idx >= n)
        return 1.0;

    const double w = 0.001 * static_cast<double>(ws[static_cast<std::size_t>(idx)]);
    return sanitise_multisim_weight(w);
}

double universe_weight_from_map(const std::map<std::string, std::vector<double>> &wmap,
                                const std::string &key,
                                int idx)
{
    auto it = wmap.find(key);
    if (it == wmap.end() || idx < 0)
        return 1.0;

    const auto &ws = it->second;
    const auto n = static_cast<int>(ws.size());
    if (idx >= n)
        return 1.0;

    return sanitise_multisim_weight(ws[static_cast<std::size_t>(idx)]);
}

std::vector<int> choose_overlay_indices(int first, int last, int max_overlay)
{
    std::vector<int> out;
    if (max_overlay <= 0 || last < first)
        return out;

    const int nsel = last - first + 1;
    if (nsel <= max_overlay)
    {
        out.reserve(static_cast<std::size_t>(nsel));
        for (int i = first; i <= last; ++i)
            out.push_back(i);
        return out;
    }

    out.reserve(static_cast<std::size_t>(max_overlay));
    for (int k = 0; k < max_overlay; ++k)
    {
        const double frac = (max_overlay == 1) ? 0.0 : static_cast<double>(k) / static_cast<double>(max_overlay - 1);
        const int idx = first + static_cast<int>(std::llround(frac * static_cast<double>(nsel - 1)));
        if (out.empty() || idx != out.back())
            out.push_back(idx);
    }
    return out;
}

struct UniverseBandSummary
{
    std::unique_ptr<TH1D> h_nominal_total;
    std::unique_ptr<TH1D> h_mean_total;
    std::unique_ptr<TH1D> h_band_total;
    std::unique_ptr<TH1D> h_env_lo_total;
    std::unique_ptr<TH1D> h_env_hi_total;

    std::unique_ptr<TH1D> h_ratio_mean;
    std::unique_ptr<TH1D> h_ratio_band;
    std::unique_ptr<TH1D> h_ratio_env_lo;
    std::unique_ptr<TH1D> h_ratio_env_hi;

    std::vector<std::unique_ptr<TH1D>> overlay_total;
    std::vector<std::unique_ptr<TH1D>> overlay_ratio;
};

UniverseBandSummary build_universe_summary(const TH1D &h_nom_mc,
                                           const TH1D *h_nom_ext,
                                           const std::vector<int> &overlay_indices,
                                           const std::vector<ROOT::RDF::RResultPtr<TH1D>> &booked_universe_hists,
                                           int first_universe,
                                           int last_universe,
                                           const std::string &safe_label)
{
    UniverseBandSummary out;

    out.h_nominal_total = clone_hist(h_nom_mc, "h_nominal_total_" + safe_label);
    if (h_nom_ext != nullptr)
        out.h_nominal_total->Add(h_nom_ext);

    out.h_mean_total = clone_hist(*out.h_nominal_total, "h_mean_total_" + safe_label);
    out.h_band_total = clone_hist(*out.h_nominal_total, "h_band_total_" + safe_label);
    out.h_env_lo_total = clone_hist(*out.h_nominal_total, "h_env_lo_total_" + safe_label);
    out.h_env_hi_total = clone_hist(*out.h_nominal_total, "h_env_hi_total_" + safe_label);

    out.h_ratio_mean = clone_hist(*out.h_nominal_total, "h_ratio_mean_" + safe_label);
    out.h_ratio_band = clone_hist(*out.h_nominal_total, "h_ratio_band_" + safe_label);
    out.h_ratio_env_lo = clone_hist(*out.h_nominal_total, "h_ratio_env_lo_" + safe_label);
    out.h_ratio_env_hi = clone_hist(*out.h_nominal_total, "h_ratio_env_hi_" + safe_label);

    out.h_mean_total->Reset("ICES");
    out.h_band_total->Reset("ICES");
    out.h_env_lo_total->Reset("ICES");
    out.h_env_hi_total->Reset("ICES");
    out.h_ratio_mean->Reset("ICES");
    out.h_ratio_band->Reset("ICES");
    out.h_ratio_env_lo->Reset("ICES");
    out.h_ratio_env_hi->Reset("ICES");

    const int n_bins_all = out.h_nominal_total->GetNbinsX() + 2;
    const int n_universes = std::max(0, last_universe - first_universe + 1);

    for (int b = 0; b < n_bins_all; ++b)
    {
        if (n_universes == 0)
            continue;

        double sum = 0.0;
        double sumsq = 0.0;
        double vmin = std::numeric_limits<double>::infinity();
        double vmax = -std::numeric_limits<double>::infinity();

        for (int u = first_universe; u <= last_universe; ++u)
        {
            const TH1D &h_mc_u = booked_universe_hists[static_cast<std::size_t>(u)].GetValue();
            double val = h_mc_u.GetBinContent(b);
            if (h_nom_ext != nullptr)
                val += h_nom_ext->GetBinContent(b);

            sum += val;
            sumsq += val * val;
            vmin = std::min(vmin, val);
            vmax = std::max(vmax, val);
        }

        const double mean = sum / static_cast<double>(n_universes);
        const double var = std::max(0.0, sumsq / static_cast<double>(n_universes) - mean * mean);
        const double rms = std::sqrt(var);

        out.h_mean_total->SetBinContent(b, mean);
        out.h_mean_total->SetBinError(b, 0.0);

        out.h_band_total->SetBinContent(b, mean);
        out.h_band_total->SetBinError(b, rms);

        out.h_env_lo_total->SetBinContent(b, vmin);
        out.h_env_lo_total->SetBinError(b, 0.0);
        out.h_env_hi_total->SetBinContent(b, vmax);
        out.h_env_hi_total->SetBinError(b, 0.0);

        const double nom = out.h_nominal_total->GetBinContent(b);
        if (nom > 0.0)
        {
            out.h_ratio_mean->SetBinContent(b, mean / nom);
            out.h_ratio_band->SetBinContent(b, mean / nom);
            out.h_ratio_band->SetBinError(b, rms / nom);
            out.h_ratio_env_lo->SetBinContent(b, vmin / nom);
            out.h_ratio_env_hi->SetBinContent(b, vmax / nom);
        }
        else
        {
            out.h_ratio_mean->SetBinContent(b, 1.0);
            out.h_ratio_band->SetBinContent(b, 1.0);
            out.h_ratio_band->SetBinError(b, 0.0);
            out.h_ratio_env_lo->SetBinContent(b, 1.0);
            out.h_ratio_env_hi->SetBinContent(b, 1.0);
        }
    }

    out.overlay_total.reserve(overlay_indices.size());
    out.overlay_ratio.reserve(overlay_indices.size());
    for (int u : overlay_indices)
    {
        auto htot = clone_hist(booked_universe_hists[static_cast<std::size_t>(u)].GetValue(),
                               "h_overlay_total_" + safe_label + "_u" + std::to_string(u));
        if (h_nom_ext != nullptr)
            htot->Add(h_nom_ext);

        auto hrat = clone_hist(*htot,
                               "h_overlay_ratio_" + safe_label + "_u" + std::to_string(u));
        hrat->Divide(out.h_nominal_total.get());

        out.overlay_total.push_back(std::move(htot));
        out.overlay_ratio.push_back(std::move(hrat));
    }

    return out;
}

std::pair<double, double> find_ratio_range(const UniverseBandSummary &sum)
{
    double rmin = std::numeric_limits<double>::infinity();
    double rmax = -std::numeric_limits<double>::infinity();

    const auto scan_hist = [&](const TH1D &h, bool use_errors)
    {
        const int nb = h.GetNbinsX();
        for (int b = 1; b <= nb; ++b)
        {
            const double c = h.GetBinContent(b);
            const double e = use_errors ? h.GetBinError(b) : 0.0;
            rmin = std::min(rmin, c - e);
            rmax = std::max(rmax, c + e);
        }
    };

    scan_hist(*sum.h_ratio_band, true);
    scan_hist(*sum.h_ratio_env_lo, false);
    scan_hist(*sum.h_ratio_env_hi, false);
    for (const auto &h : sum.overlay_ratio)
        scan_hist(*h, false);

    if (!std::isfinite(rmin) || !std::isfinite(rmax))
        return {0.8, 1.2};

    const double pad = 0.10 * std::max(0.05, rmax - rmin);
    return {std::max(0.0, rmin - pad), rmax + pad};
}

void draw_multisim_summary(const UniverseBandSummary &sum,
                           const std::string &family_label,
                           const std::vector<int> &overlay_indices,
                           int first_universe,
                           int last_universe,
                           double xmin,
                           double xmax,
                           bool use_logy,
                           const std::string &output_stem)
{
    auto h_nom = clone_hist(*sum.h_nominal_total, "h_nom_draw_" + make_safe_name(family_label));
    auto h_mean = clone_hist(*sum.h_mean_total, "h_mean_draw_" + make_safe_name(family_label));
    auto h_band = clone_hist(*sum.h_band_total, "h_band_draw_" + make_safe_name(family_label));
    auto h_env_lo = clone_hist(*sum.h_env_lo_total, "h_env_lo_draw_" + make_safe_name(family_label));
    auto h_env_hi = clone_hist(*sum.h_env_hi_total, "h_env_hi_draw_" + make_safe_name(family_label));
    auto h_ratio_mean = clone_hist(*sum.h_ratio_mean, "h_ratio_mean_draw_" + make_safe_name(family_label));
    auto h_ratio_band = clone_hist(*sum.h_ratio_band, "h_ratio_band_draw_" + make_safe_name(family_label));
    auto h_ratio_env_lo = clone_hist(*sum.h_ratio_env_lo, "h_ratio_env_lo_draw_" + make_safe_name(family_label));
    auto h_ratio_env_hi = clone_hist(*sum.h_ratio_env_hi, "h_ratio_env_hi_draw_" + make_safe_name(family_label));

    style_line_hist(*h_nom, kBlack, 1, 3);
    style_line_hist(*h_mean, kBlue + 1, 1, 3);
    style_band_hist(*h_band, kBlue + 1, kAzure - 9);
    style_line_hist(*h_env_lo, kBlue + 2, 2, 2);
    style_line_hist(*h_env_hi, kBlue + 2, 2, 2);
    style_line_hist(*h_ratio_mean, kBlue + 1, 1, 3);
    style_band_hist(*h_ratio_band, kBlue + 1, kAzure - 9);
    style_line_hist(*h_ratio_env_lo, kBlue + 2, 2, 2);
    style_line_hist(*h_ratio_env_hi, kBlue + 2, 2, 2);

    std::vector<std::unique_ptr<TH1D>> overlay_total;
    std::vector<std::unique_ptr<TH1D>> overlay_ratio;
    overlay_total.reserve(sum.overlay_total.size());
    overlay_ratio.reserve(sum.overlay_ratio.size());
    for (std::size_t i = 0; i < sum.overlay_total.size(); ++i)
    {
        overlay_total.push_back(clone_hist(*sum.overlay_total[i],
                                           "h_overlay_total_draw_" + std::to_string(i) + "_" + make_safe_name(family_label)));
        overlay_ratio.push_back(clone_hist(*sum.overlay_ratio[i],
                                           "h_overlay_ratio_draw_" + std::to_string(i) + "_" + make_safe_name(family_label)));
        style_line_hist(*overlay_total.back(), kGray + 1, 1, 1);
        style_line_hist(*overlay_ratio.back(), kGray + 1, 1, 1);
    }

    double ymax = h_nom->GetMaximum();
    ymax = std::max(ymax, h_mean->GetMaximum() + h_band->GetBinError(h_band->GetMaximumBin()));
    ymax = std::max(ymax, h_env_hi->GetMaximum());
    for (const auto &h : overlay_total)
        ymax = std::max(ymax, h->GetMaximum());
    if (!(ymax > 0.0))
        ymax = 1.0;

    const auto rr = find_ratio_range(sum);
    const std::string safe = make_safe_name(family_label);

    TCanvas c(("c_multisim_" + safe).c_str(), family_label.c_str(), 1000, 850);
    c.SetMargin(0.0, 0.0, 0.0, 0.0);

    TPad p1(("p1_multisim_" + safe).c_str(), "", 0.0, 0.30, 1.0, 1.0);
    TPad p2(("p2_multisim_" + safe).c_str(), "", 0.0, 0.00, 1.0, 0.30);

    p1.SetLeftMargin(0.12);
    p1.SetRightMargin(0.05);
    p1.SetBottomMargin(0.02);
    p1.SetTopMargin(0.06);
    if (use_logy)
        p1.SetLogy();

    p2.SetLeftMargin(0.12);
    p2.SetRightMargin(0.05);
    p2.SetTopMargin(0.02);
    p2.SetBottomMargin(0.35);

    p1.Draw();
    p2.Draw();

    p1.cd();
    TH1D hframe(("hframe_multisim_" + safe).c_str(), ";Inference score [0];Events", 100, xmin, xmax);
    hframe.SetDirectory(nullptr);
    hframe.SetMinimum(use_logy ? 0.3 : 0.0);
    hframe.SetMaximum(use_logy ? 20.0 * ymax : 1.35 * ymax);
    hframe.GetYaxis()->SetTitleSize(0.050);
    hframe.GetYaxis()->SetLabelSize(0.045);
    hframe.GetYaxis()->SetTitleOffset(1.10);
    hframe.GetXaxis()->SetTitleSize(0.0);
    hframe.GetXaxis()->SetLabelSize(0.0);
    hframe.Draw("AXIS");

    for (const auto &h : overlay_total)
        h->Draw("HIST SAME");
    h_band->Draw("E2 SAME");
    h_env_lo->Draw("HIST SAME");
    h_env_hi->Draw("HIST SAME");
    h_mean->Draw("HIST SAME");
    h_nom->Draw("HIST SAME");

    TLegend leg(0.50, 0.62, 0.93, 0.90);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h_nom.get(), "nominal total", "l");
    leg.AddEntry(h_mean.get(), (family_label + " mean").c_str(), "l");
    leg.AddEntry(h_band.get(), (family_label + " mean #pm RMS").c_str(), "f");
    leg.AddEntry(h_env_hi.get(), (family_label + " min/max envelope").c_str(), "l");
    if (!overlay_total.empty())
    {
        std::string ov = family_label + " overlay universes";
        if (!overlay_indices.empty())
            ov += " (" + std::to_string(overlay_indices.front()) + "..." + std::to_string(overlay_indices.back()) + ")";
        leg.AddEntry(overlay_total.front().get(), ov.c_str(), "l");
    }
    leg.Draw();

    p2.cd();
    TH1D rframe(("rframe_multisim_" + safe).c_str(), ";Inference score [0];universe / nominal", 100, xmin, xmax);
    rframe.SetDirectory(nullptr);
    rframe.SetMinimum(rr.first);
    rframe.SetMaximum(rr.second);
    rframe.GetXaxis()->SetTitleSize(0.12);
    rframe.GetXaxis()->SetLabelSize(0.10);
    rframe.GetYaxis()->SetTitleSize(0.10);
    rframe.GetYaxis()->SetLabelSize(0.08);
    rframe.GetYaxis()->SetTitleOffset(0.55);
    rframe.GetYaxis()->SetNdivisions(505);
    rframe.Draw("AXIS");

    TLine unity(xmin, 1.0, xmax, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();

    for (const auto &h : overlay_ratio)
        h->Draw("HIST SAME");
    h_ratio_band->Draw("E2 SAME");
    h_ratio_env_lo->Draw("HIST SAME");
    h_ratio_env_hi->Draw("HIST SAME");
    h_ratio_mean->Draw("HIST SAME");

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());

    const int last_bin = h_nom->GetNbinsX() + 1;
    std::cout << "[plot_inference_score_multisim_branch] "
              << family_label
              << " universes=" << (last_universe - first_universe + 1)
              << " selected_range=[" << first_universe << "," << last_universe << "]"
              << " nominal_total=" << h_nom->Integral(0, last_bin)
              << " mean_total=" << h_mean->Integral(0, last_bin)
              << " saved=" << out << "\n";
}

} // namespace

int plot_inference_score_multisim_branch(const std::string &event_list_path = "",
                                         const std::string &base_sel = "sel_muon",
                                         const std::string &mc_weight = "w_nominal",
                                         const std::string &weight_branch = "weightsGenie",
                                         const std::string &map_key = "All_UBGenie",
                                         const std::string &family_label = "GENIE multisims",
                                         int nbins = 60,
                                         double xmin = -15.0,
                                         double xmax = 15.0,
                                         bool include_ext = true,
                                         bool use_logy = false,
                                         int max_overlay_universes = 10,
                                         int first_universe = 0,
                                         int last_universe = -1,
                                         const std::string &output_stem = "inference_score_multisim_branch")
{
    return heron::macro::run_with_guard("plot_inference_score_multisim_branch", [&]() -> int {
        ROOT::EnableImplicitMT();
        TH1::SetDefaultSumw2();
        gStyle->SetOptStat(0);

        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_inference_score_multisim_branch] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inference_score_multisim_branch] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode base = define_score0(SelectionService::decorate(el.rdf()));
        if (!has_column(base, "score0_for_plot"))
        {
            std::cerr << "[plot_inference_score_multisim_branch] could not resolve a score column. "
                      << "Expected one of: inf_score_0, inf_scores, inference_score, inference_scores.\n";
            return 1;
        }

        base = apply_selection(base, base_sel);

        auto mask_ext = el.mask_for_ext();
        auto mask_mc_like = el.mask_for_mc_like();

        ROOT::RDF::RNode node_all = filter_by_sample_mask(base, mask_mc_like, "sample_id");
        ROOT::RDF::RNode node_mc = filter_not_sample_mask(node_all, mask_ext, "sample_id");
        ROOT::RDF::RNode node_ext = node_all;

        const bool has_direct_branch = has_column(node_mc, weight_branch);
        const bool has_weights_map = has_column(node_mc, "weights");
        const bool use_map_fallback = (!has_direct_branch && has_weights_map && !map_key.empty());

        if (!has_direct_branch && !use_map_fallback)
        {
            std::cerr << "[plot_inference_score_multisim_branch] missing multisim branch '" << weight_branch << "'";
            if (!map_key.empty())
                std::cerr << " and no usable weights-map fallback for key '" << map_key << "'";
            std::cerr << ".\n";
            return 1;
        }

        node_mc = node_mc.Define("__w_nom__", mc_weight);
        if (has_direct_branch)
        {
            node_mc = node_mc.Define("__n_univ__", "static_cast<int>(" + weight_branch + ".size())");
        }
        else
        {
            node_mc = node_mc.Define(
                "__n_univ__",
                [map_key](const std::map<std::string, std::vector<double>> &wmap)
                {
                    auto it = wmap.find(map_key);
                    return (it == wmap.end()) ? 0 : static_cast<int>(it->second.size());
                },
                {"weights"});
        }

        if (include_ext)
            node_ext = filter_by_sample_mask(node_all, mask_ext, "sample_id").Define("__w_nom__", mc_weight);

        int n_universes = 0;
        {
            auto max_n = node_mc.Max<int>("__n_univ__");
            n_universes = *max_n;
        }

        if (n_universes <= 0)
        {
            std::cerr << "[plot_inference_score_multisim_branch] branch '" << weight_branch
                      << "' has no universes after selection.\n";
            return 1;
        }

        first_universe = std::max(0, first_universe);
        if (last_universe < 0 || last_universe >= n_universes)
            last_universe = n_universes - 1;

        if (first_universe > last_universe)
        {
            std::cerr << "[plot_inference_score_multisim_branch] invalid universe range: ["
                      << first_universe << ", " << last_universe << "] for n_universes=" << n_universes << "\n";
            return 1;
        }

        std::cout << "[plot_inference_score_multisim_branch] family='" << family_label
                  << "' source='" << (has_direct_branch ? weight_branch : (std::string("weights[") + map_key + "]"))
                  << "' n_universes=" << n_universes
                  << " using range [" << first_universe << ", " << last_universe << "]\n";

        ROOT::RDF::TH1DModel h_nom_mc_model(("h_nom_mc_" + make_safe_name(weight_branch)).c_str(), "", nbins, xmin, xmax);
        auto h_nom_mc = node_mc.Histo1D(h_nom_mc_model, "score0_for_plot", "__w_nom__");

        ROOT::RDF::RResultPtr<TH1D> h_nom_ext;
        if (include_ext)
        {
            ROOT::RDF::TH1DModel h_nom_ext_model(("h_nom_ext_" + make_safe_name(weight_branch)).c_str(), "", nbins, xmin, xmax);
            h_nom_ext = node_ext.Histo1D(h_nom_ext_model, "score0_for_plot", "__w_nom__");
        }

        std::vector<ROOT::RDF::RResultPtr<TH1D>> booked_universe_hists;
        booked_universe_hists.reserve(static_cast<std::size_t>(n_universes));

        for (int u = 0; u < n_universes; ++u)
        {
            const std::string safe_u = make_safe_name(weight_branch) + "_u" + std::to_string(u);
            const std::string w_name = "__w_" + safe_u;

            ROOT::RDF::RNode node_u = node_mc;
            if (has_direct_branch)
            {
                node_u = node_u.Define(
                    w_name,
                    "__w_nom__ * universe_weight_from_vec(" + weight_branch + ", " + std::to_string(u) + ")");
            }
            else
            {
                node_u = node_u.Define(
                    w_name,
                    "__w_nom__ * universe_weight_from_map(weights, \"" + map_key + "\", " + std::to_string(u) + ")");
            }

            ROOT::RDF::TH1DModel h_model(("h_mc_" + safe_u).c_str(), "", nbins, xmin, xmax);
            booked_universe_hists.push_back(node_u.Histo1D(h_model, "score0_for_plot", w_name));
        }

        const std::vector<int> overlay_indices = choose_overlay_indices(first_universe, last_universe, max_overlay_universes);

        const TH1D &nom_mc = *h_nom_mc;
        const TH1D *nom_ext = include_ext ? &(*h_nom_ext) : nullptr;

        UniverseBandSummary summary = build_universe_summary(nom_mc,
                                                             nom_ext,
                                                             overlay_indices,
                                                             booked_universe_hists,
                                                             first_universe,
                                                             last_universe,
                                                             make_safe_name(family_label));

        draw_multisim_summary(summary,
                              family_label,
                              overlay_indices,
                              first_universe,
                              last_universe,
                              xmin,
                              xmax,
                              use_logy,
                              output_stem);

        std::cout << "[plot_inference_score_multisim_branch] done\n";
        return 0;
    });
}

int plot_inference_score_genie_universes(const std::string &event_list_path = "",
                                         const std::string &base_sel = "sel_muon",
                                         const std::string &mc_weight = "w_nominal",
                                         int nbins = 60,
                                         double xmin = -15.0,
                                         double xmax = 15.0,
                                         bool include_ext = true,
                                         bool use_logy = false,
                                         int max_overlay_universes = 10,
                                         int first_universe = 0,
                                         int last_universe = -1,
                                         const std::string &output_stem = "inference_score_genie_universes")
{
    return plot_inference_score_multisim_branch(event_list_path,
                                                base_sel,
                                                mc_weight,
                                                "weightsGenie",
                                                "All_UBGenie",
                                                "GENIE universes",
                                                nbins,
                                                xmin,
                                                xmax,
                                                include_ext,
                                                use_logy,
                                                max_overlay_universes,
                                                first_universe,
                                                last_universe,
                                                output_stem);
}

int plot_inference_score_flux_universes(const std::string &event_list_path = "",
                                        const std::string &base_sel = "sel_muon",
                                        const std::string &mc_weight = "w_nominal",
                                        int nbins = 60,
                                        double xmin = -15.0,
                                        double xmax = 15.0,
                                        bool include_ext = true,
                                        bool use_logy = false,
                                        int max_overlay_universes = 10,
                                        int first_universe = 0,
                                        int last_universe = -1,
                                        const std::string &output_stem = "inference_score_flux_universes")
{
    return plot_inference_score_multisim_branch(event_list_path,
                                                base_sel,
                                                mc_weight,
                                                "weightsFlux",
                                                "flux_all",
                                                "Flux universes",
                                                nbins,
                                                xmin,
                                                xmax,
                                                include_ext,
                                                use_logy,
                                                max_overlay_universes,
                                                first_universe,
                                                last_universe,
                                                output_stem);
}

int plot_inference_score_reint_universes(const std::string &event_list_path = "",
                                         const std::string &base_sel = "sel_muon",
                                         const std::string &mc_weight = "w_nominal",
                                         int nbins = 60,
                                         double xmin = -15.0,
                                         double xmax = 15.0,
                                         bool include_ext = true,
                                         bool use_logy = false,
                                         int max_overlay_universes = 10,
                                         int first_universe = 0,
                                         int last_universe = -1,
                                         const std::string &output_stem = "inference_score_reint_universes")
{
    return plot_inference_score_multisim_branch(event_list_path,
                                                base_sel,
                                                mc_weight,
                                                "weightsReint",
                                                "reint_all",
                                                "Reinteraction universes",
                                                nbins,
                                                xmin,
                                                xmax,
                                                include_ext,
                                                use_logy,
                                                max_overlay_universes,
                                                first_universe,
                                                last_universe,
                                                output_stem);
}

int plot_inference_score_ppfx_universes(const std::string &event_list_path = "",
                                        const std::string &base_sel = "sel_muon",
                                        const std::string &mc_weight = "w_nominal",
                                        int nbins = 60,
                                        double xmin = -15.0,
                                        double xmax = 15.0,
                                        bool include_ext = true,
                                        bool use_logy = false,
                                        int max_overlay_universes = 10,
                                        int first_universe = 0,
                                        int last_universe = -1,
                                        const std::string &output_stem = "inference_score_ppfx_universes")
{
    return plot_inference_score_multisim_branch(event_list_path,
                                                base_sel,
                                                mc_weight,
                                                "weightsPPFX",
                                                "ppfx_all",
                                                "PPFX universes",
                                                nbins,
                                                xmin,
                                                xmax,
                                                include_ext,
                                                use_logy,
                                                max_overlay_universes,
                                                first_universe,
                                                last_universe,
                                                output_stem);
}

#endif // PLOT_INFERENCE_SCORE_MULTISIM_BRANCH_CXX
