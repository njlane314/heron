#ifndef PLOT_INF_SCORE0_MULTISIM_VARIATIONS_CXX
#define PLOT_INF_SCORE0_MULTISIM_VARIATIONS_CXX

// plot_inf_score0_multisim_variations.C
//
// Draw the nominal inf_scores[0] histogram together with many multisim
// variation histograms for a chosen packed-vector branch.
//
// Example:
//   ./heron macro macros/plot_inf_score0_multisim_variations.C \
//   'plot_inf_score0_multisim_variations("", "sel_muon", "true", "w_nominal", "weightsGenie", 12, -15.0, 15.0, false, 40, true, false, "inf_score0_genie_variations")'
//
// Useful branches:
//   weightsGenie
//   weightsPPFX
//   weightsFlux
//   weightsReint

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
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
#include "Plotter.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

using UShortVec = ROOT::VecOps::RVec<unsigned short>;
using DoubleTake = ROOT::RDF::RResultPtr<std::vector<double>>;
using PackedTake = ROOT::RDF::RResultPtr<std::vector<UShortVec>>;

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

ROOT::RDF::RNode apply_optional_filter(ROOT::RDF::RNode node,
                                       const std::string &expr)
{
    if (expr.empty() || expr == "true")
        return node;

    if (has_column(node, expr))
        return node.Filter([](bool pass) { return pass; }, {expr});

    return node.Filter(expr);
}

int score_bin(double x, int nbins, double xmin, double xmax)
{
    if (!(x >= xmin) || !(x < xmax))
        return -1;

    const double span = xmax - xmin;
    if (!(span > 0.0))
        return -1;

    int b = static_cast<int>(std::floor((x - xmin) / span * nbins));
    if (b < 0)
        return -1;
    if (b >= nbins)
        b = nbins - 1;
    return b;
}

double unpack_packed_weight(unsigned short packed)
{
    return 0.001 * static_cast<double>(packed);
}

std::vector<std::size_t> choose_overlay_ids(const std::vector<std::size_t> &active_ids,
                                            int max_overlay)
{
    if (max_overlay <= 0 || active_ids.empty())
        return {};

    if (static_cast<int>(active_ids.size()) <= max_overlay)
        return active_ids;

    std::vector<std::size_t> out;
    const double step = static_cast<double>(active_ids.size()) /
                        static_cast<double>(max_overlay);

    for (int i = 0; i < max_overlay; ++i)
    {
        const std::size_t idx =
            static_cast<std::size_t>(std::floor(i * step));
        out.push_back(active_ids[std::min(idx, active_ids.size() - 1)]);
    }

    return out;
}

void fill_hist_from_bins(TH1D &h, const std::vector<double> &bins)
{
    for (int b = 1; b <= h.GetNbinsX(); ++b)
        h.SetBinContent(b, bins[static_cast<std::size_t>(b - 1)]);
}

double max_content(const std::vector<std::vector<double>> &all_bins)
{
    double out = 0.0;
    for (const auto &bins : all_bins)
        for (double x : bins)
            out = std::max(out, x);
    return out;
}

} // namespace

int plot_inf_score0_multisim_variations(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &extra_sel = "true",
    const std::string &nominal_weight = "w_nominal",
    const std::string &multisim_col = "weightsGenie",
    int nbins = 12,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_nominal = false,
    int max_overlay = 40,
    bool draw_ratio = true,
    bool normalise = false,
    const std::string &output_stem = "inf_score0_multisim_variations")
{
    return heron::macro::run_with_guard("plot_inf_score0_multisim_variations", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        if (implicit_mt_enabled())
            ROOT::EnableImplicitMT();

        TH1::SetDefaultSumw2();

        const std::string input_path =
            event_list_path.empty() ? default_event_list_root() : event_list_path;

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inf_score0_multisim_variations] input is not an event-list ROOT file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode node = SelectionService::decorate(el.rdf());

        if (!has_column(node, "inf_scores"))
        {
            std::cerr << "[plot_inf_score0_multisim_variations] missing inf_scores column.\n";
            return 1;
        }

        if (!has_column(node, multisim_col))
        {
            std::cerr << "[plot_inf_score0_multisim_variations] missing multisim column: "
                      << multisim_col << "\n";
            return 1;
        }

        node = node.Define(
            "inf_score_0",
            [](const ROOT::VecOps::RVec<float> &scores) {
                return scores.empty()
                           ? std::numeric_limits<double>::quiet_NaN()
                           : static_cast<double>(scores[0]);
            },
            {"inf_scores"});

        node = node.Define("__nom_w__", nominal_weight);

        const auto mask_mc_like = el.mask_for_mc_like();
        const auto mask_ext = el.mask_for_ext();

        node = filter_by_sample_mask(node, mask_mc_like, "sample_id");
        if (!include_ext_in_nominal)
            node = filter_not_sample_mask(node, mask_ext, "sample_id");

        node = apply_optional_filter(node, base_sel);
        node = apply_optional_filter(node, extra_sel);

        node = node.Filter(
            [xmin, xmax](double s, double w) {
                return std::isfinite(s) && std::isfinite(w) && (s >= xmin) && (s < xmax);
            },
            {"inf_score_0", "__nom_w__"});

        DoubleTake scores_take = node.Take<double>("inf_score_0");
        DoubleTake nom_take = node.Take<double>("__nom_w__");
        PackedTake multisim_take = node.Take<UShortVec>(multisim_col);
        auto n_selected = node.Count();

        const ULong64_t n_evt = *n_selected;
        if (n_evt == 0)
        {
            std::cerr << "[plot_inf_score0_multisim_variations] no selected events.\n";
            return 1;
        }

        const auto &scores = *scores_take;
        const auto &nom_w = *nom_take;
        const auto &packed_vecs = *multisim_take;

        if (scores.size() != nom_w.size() || scores.size() != packed_vecs.size())
        {
            std::cerr << "[plot_inf_score0_multisim_variations] internal size mismatch.\n";
            return 1;
        }

        std::size_t n_total = 0;
        for (const auto &v : packed_vecs)
            n_total = std::max(n_total, static_cast<std::size_t>(v.size()));

        if (n_total == 0)
        {
            std::cerr << "[plot_inf_score0_multisim_variations] no universes found in "
                      << multisim_col << "\n";
            return 1;
        }

        std::vector<double> nom_bins(static_cast<std::size_t>(nbins), 0.0);
        std::vector<std::vector<double>> univ_bins(
            n_total, std::vector<double>(static_cast<std::size_t>(nbins), 0.0));
        std::vector<char> active(n_total, 0);

        for (std::size_t i = 0; i < scores.size(); ++i)
        {
            const int b = score_bin(scores[i], nbins, xmin, xmax);
            if (b < 0)
                continue;

            const std::size_t ib = static_cast<std::size_t>(b);
            const double w = nom_w[i];
            nom_bins[ib] += w;

            const auto &packed = packed_vecs[i];

            for (std::size_t u = 0; u < n_total; ++u)
            {
                const double wu =
                    (u < packed.size()) ? unpack_packed_weight(packed[u]) : 1.0;

                if (u < packed.size() && packed[u] != 1000)
                    active[u] = 1;

                univ_bins[u][ib] += w * wu;
            }
        }

        if (normalise)
        {
            double nom_int = 0.0;
            for (double x : nom_bins)
                nom_int += x;

            if (nom_int > 0.0)
                for (double &x : nom_bins)
                    x /= nom_int;

            for (auto &bins : univ_bins)
            {
                double s = 0.0;
                for (double x : bins)
                    s += x;
                if (s > 0.0)
                    for (double &x : bins)
                        x /= s;
            }
        }

        std::vector<std::size_t> active_ids;
        for (std::size_t u = 0; u < n_total; ++u)
            if (active[u])
                active_ids.push_back(u);

        if (active_ids.empty())
        {
            std::cerr << "[plot_inf_score0_multisim_variations] all universes are nominal in the selected sample.\n";
            return 1;
        }

        const auto draw_ids = choose_overlay_ids(active_ids, max_overlay);

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        TCanvas c("c_inf_score0_multisim_variations",
                  "Inference score multisim variations",
                  1100,
                  draw_ratio ? 900 : 700);

        TPad *pad_top = nullptr;
        TPad *pad_bot = nullptr;

        if (draw_ratio)
        {
            pad_top = new TPad("pad_top", "", 0.0, 0.30, 1.0, 1.0);
            pad_bot = new TPad("pad_bot", "", 0.0, 0.00, 1.0, 0.30);

            pad_top->SetBottomMargin(0.02);
            pad_top->SetLeftMargin(0.11);
            pad_top->SetRightMargin(0.04);

            pad_bot->SetTopMargin(0.03);
            pad_bot->SetBottomMargin(0.33);
            pad_bot->SetLeftMargin(0.11);
            pad_bot->SetRightMargin(0.04);

            pad_top->Draw();
            pad_bot->Draw();
            pad_top->cd();
        }
        else
        {
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);
        }

        TH1D h_nom("h_nom", "", nbins, xmin, xmax);
        fill_hist_from_bins(h_nom, nom_bins);
        h_nom.SetDirectory(nullptr);
        h_nom.SetLineColor(kBlack);
        h_nom.SetLineWidth(3);
        h_nom.SetFillStyle(0);

        std::vector<std::unique_ptr<TH1D>> h_vars;
        std::vector<std::vector<double>> drawn_bins;
        drawn_bins.push_back(nom_bins);

        for (std::size_t k = 0; k < draw_ids.size(); ++k)
        {
            const auto u = draw_ids[k];
            auto h = std::make_unique<TH1D>(
                ("h_var_" + std::to_string(u)).c_str(), "", nbins, xmin, xmax);
            fill_hist_from_bins(*h, univ_bins[u]);
            h->SetDirectory(nullptr);
            h->SetLineColor(kGreen + 1);
            h->SetLineWidth(1);
            h->SetFillStyle(0);
            h_vars.push_back(std::move(h));
            drawn_bins.push_back(univ_bins[u]);
        }

        const double ymax = 1.15 * max_content(drawn_bins);

        TH1D hframe("hframe", "", 100, xmin, xmax);
        hframe.SetMinimum(0.0);
        hframe.SetMaximum(ymax > 0.0 ? ymax : 1.0);
        hframe.GetXaxis()->SetTitle(draw_ratio ? "" : "Inference score [0]");
        hframe.GetYaxis()->SetTitle(normalise ? "Arbitrary units" : "Events");
        hframe.Draw("AXIS");

        for (const auto &h : h_vars)
            h->Draw("HIST SAME");

        h_nom.Draw("HIST SAME");

        TLegend leg(0.50, 0.77, 0.89, 0.89);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(&h_nom, "Central Value", "l");
        leg.AddEntry(h_vars.empty() ? &h_nom : h_vars.front().get(), "Variations", "l");
        leg.Draw();

        gPad->RedrawAxis();

        if (draw_ratio)
        {
            c.cd();
            pad_bot->cd();

            TH1D hframe_ratio("hframe_ratio", "", 100, xmin, xmax);
            hframe_ratio.SetMinimum(0.5);
            hframe_ratio.SetMaximum(1.5);
            hframe_ratio.GetXaxis()->SetTitle("Inference score [0]");
            hframe_ratio.GetYaxis()->SetTitle("Var./CV");
            hframe_ratio.GetXaxis()->SetTitleSize(0.11);
            hframe_ratio.GetXaxis()->SetLabelSize(0.10);
            hframe_ratio.GetYaxis()->SetTitleSize(0.10);
            hframe_ratio.GetYaxis()->SetLabelSize(0.08);
            hframe_ratio.GetYaxis()->SetTitleOffset(0.50);
            hframe_ratio.GetYaxis()->SetNdivisions(505);
            hframe_ratio.Draw("AXIS");

            for (std::size_t k = 0; k < draw_ids.size(); ++k)
            {
                const auto u = draw_ids[k];
                TH1D h_ratio(("h_ratio_" + std::to_string(u)).c_str(), "", nbins, xmin, xmax);
                h_ratio.SetDirectory(nullptr);
                h_ratio.SetLineColor(kGreen + 1);
                h_ratio.SetLineWidth(1);
                h_ratio.SetFillStyle(0);

                for (int b = 1; b <= nbins; ++b)
                {
                    const double nom = nom_bins[static_cast<std::size_t>(b - 1)];
                    const double var = univ_bins[u][static_cast<std::size_t>(b - 1)];
                    h_ratio.SetBinContent(b, nom > 0.0 ? var / nom : 0.0);
                }

                h_ratio.Draw("HIST SAME");
            }

            TLine line1(xmin, 1.0, xmax, 1.0);
            line1.SetLineColor(kBlack);
            line1.SetLineStyle(2);
            line1.Draw("SAME");

            gPad->RedrawAxis();
        }

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());

        std::cout << "[plot_inf_score0_multisim_variations] selected events = " << n_evt << "\n";
        std::cout << "[plot_inf_score0_multisim_variations] active universes = "
                  << active_ids.size() << "/" << n_total << "\n";
        std::cout << "[plot_inf_score0_multisim_variations] drawn universes  = "
                  << draw_ids.size() << "\n";
        std::cout << "[plot_inf_score0_multisim_variations] saved: " << out << "\n";

        return 0;
    });
}

#endif
