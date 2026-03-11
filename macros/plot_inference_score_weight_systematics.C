// macros/plot_inference_score_weight_systematics.C
//
// Overlay the total inference-score[0] distribution for nominal and
// nuisance-varied event weights.
//
// Conventions used here:
//   - total = MC + EXT by default (set include_ext=false for MC-only)
//   - nuisance variations are applied to MC only:
//         w_up = w_nominal * knob...up
//         w_dn = w_nominal * knob...dn
//   - EXT stays at nominal weight
//   - weightSpline, weightTune, ppfx_cv, and RootinoFix are NOT multiplied
//     again because they are already part of w_nominal
//
// Run with:
//   ./heron macro plot_inference_score_weight_systematics.C
//   ./heron macro plot_inference_score_weight_systematics.C \
//     'plot_inference_score_weight_systematics("./scratch/out/event_list_myana.root")'
//
// Single nuisance example:
//   ./heron macro plot_inference_score_weight_systematics.C \
//     'plot_inference_score_weight_systematics("./scratch/out/event_list_myana.root", "sel_muon", "w_nominal", 60, -15, 15, true, false, "RPA_CCQE")'

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
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
#include "PlotEnv.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

struct VariationSpec
{
    std::string tag;
    std::string up_col;
    std::string dn_col;
};

struct BookedVariation
{
    VariationSpec spec;
    ROOT::RDF::RResultPtr<TH1D> h_up_mc;
    ROOT::RDF::RResultPtr<TH1D> h_dn_mc;
};

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

std::string make_safe_name(std::string s)
{
    for (char &c : s)
    {
        if (!std::isalnum(static_cast<unsigned char>(c)))
            c = '_';
    }
    return s;
}

bool matches_filter(const VariationSpec &spec, const std::string &filter)
{
    if (filter.empty())
        return true;

    return spec.tag.find(filter) != std::string::npos ||
           spec.up_col.find(filter) != std::string::npos ||
           spec.dn_col.find(filter) != std::string::npos;
}

double sanitise_variation_weight(double w)
{
    if (!std::isfinite(w) || w <= 0.0)
        return 1.0;
    return w;
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

std::unique_ptr<TH1D> clone_hist(const TH1D &h, const std::string &name)
{
    auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(name.c_str())));
    out->SetDirectory(nullptr);
    return out;
}

void style_hist(TH1D &h, int color, int line_style = 1, int line_width = 3)
{
    h.SetDirectory(nullptr);
    h.SetLineColor(color);
    h.SetMarkerColor(color);
    h.SetLineStyle(line_style);
    h.SetLineWidth(line_width);
    h.SetMarkerSize(0.0);
}

std::pair<double, double> ratio_range(const TH1D &num1, const TH1D &num2, const TH1D &den)
{
    double rmin = 1.0;
    double rmax = 1.0;
    bool have = false;

    const int nb = den.GetNbinsX();
    for (int b = 1; b <= nb; ++b)
    {
        const double d = den.GetBinContent(b);
        if (!(d > 0.0))
            continue;

        const double r1 = num1.GetBinContent(b) / d;
        const double r2 = num2.GetBinContent(b) / d;

        if (!have)
        {
            rmin = std::min(r1, r2);
            rmax = std::max(r1, r2);
            have = true;
        }
        else
        {
            rmin = std::min(rmin, std::min(r1, r2));
            rmax = std::max(rmax, std::max(r1, r2));
        }
    }

    if (!have)
        return {0.8, 1.2};

    const double pad = 0.10 * std::max(0.02, rmax - rmin);
    const double lo = std::max(0.0, rmin - pad);
    const double hi = rmax + pad;
    return {lo, hi};
}

void draw_one_variation(const VariationSpec &spec,
                        const TH1D &h_nom_mc,
                        const TH1D *h_nom_ext,
                        const TH1D &h_up_mc,
                        const TH1D &h_dn_mc,
                        double xmin,
                        double xmax,
                        bool use_logy,
                        const std::string &output_stem)
{
    const std::string safe = make_safe_name(spec.tag);

    auto h_nom = clone_hist(h_nom_mc, "h_nom_total_" + safe);
    auto h_up = clone_hist(h_up_mc, "h_up_total_" + safe);
    auto h_dn = clone_hist(h_dn_mc, "h_dn_total_" + safe);

    if (h_nom_ext != nullptr)
    {
        h_nom->Add(h_nom_ext);
        h_up->Add(h_nom_ext);
        h_dn->Add(h_nom_ext);
    }

    style_hist(*h_nom, kBlack, 1, 3);
    style_hist(*h_up, kRed + 1, 1, 3);
    style_hist(*h_dn, kBlue + 1, 2, 3);

    auto h_ratio_up = clone_hist(*h_up, "h_ratio_up_" + safe);
    auto h_ratio_dn = clone_hist(*h_dn, "h_ratio_dn_" + safe);
    h_ratio_up->Divide(h_nom.get());
    h_ratio_dn->Divide(h_nom.get());
    style_hist(*h_ratio_up, kRed + 1, 1, 3);
    style_hist(*h_ratio_dn, kBlue + 1, 2, 3);

    const auto rr = ratio_range(*h_up, *h_dn, *h_nom);

    double ymax = std::max({h_nom->GetMaximum(), h_up->GetMaximum(), h_dn->GetMaximum()});
    if (!(ymax > 0.0))
        ymax = 1.0;

    TCanvas c(("c_" + safe).c_str(), spec.tag.c_str(), 1000, 850);
    c.SetMargin(0.0, 0.0, 0.0, 0.0);

    TPad p1(("p1_" + safe).c_str(), "", 0.0, 0.30, 1.0, 1.0);
    TPad p2(("p2_" + safe).c_str(), "", 0.0, 0.00, 1.0, 0.30);

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
    TH1D hframe(("hframe_" + safe).c_str(), ";Inference score [0];Events", 100, xmin, xmax);
    hframe.SetDirectory(nullptr);
    hframe.SetMinimum(use_logy ? 0.3 : 0.0);
    hframe.SetMaximum(use_logy ? 20.0 * ymax : 1.35 * ymax);
    hframe.GetYaxis()->SetTitleSize(0.050);
    hframe.GetYaxis()->SetLabelSize(0.045);
    hframe.GetYaxis()->SetTitleOffset(1.10);
    hframe.GetXaxis()->SetTitleSize(0.0);
    hframe.GetXaxis()->SetLabelSize(0.0);
    hframe.Draw("AXIS");

    h_nom->Draw("HIST SAME");
    h_up->Draw("HIST SAME");
    h_dn->Draw("HIST SAME");

    TLegend leg(0.56, 0.71, 0.93, 0.90);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h_nom.get(), "nominal", "l");
    leg.AddEntry(h_up.get(), (spec.tag + " up").c_str(), "l");
    leg.AddEntry(h_dn.get(), (spec.tag + " down").c_str(), "l");
    leg.Draw();

    p2.cd();
    TH1D rframe(("rframe_" + safe).c_str(), ";Inference score [0];variation / nominal", 100, xmin, xmax);
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

    h_ratio_up->Draw("HIST SAME");
    h_ratio_dn->Draw("HIST SAME");

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem + "_" + safe).string();
    c.SaveAs(out.c_str());

    const int last_bin = h_nom->GetNbinsX() + 1;
    const double nom_int = h_nom->Integral(0, last_bin);
    const double up_int = h_up->Integral(0, last_bin);
    const double dn_int = h_dn->Integral(0, last_bin);

    std::cout << "[plot_inference_score_weight_systematics] "
              << spec.tag
              << " nominal=" << nom_int
              << " up=" << up_int
              << " down=" << dn_int
              << " saved=" << out << "\n";
}

} // namespace

int plot_inference_score_weight_systematics(const std::string &event_list_path = "",
                                            const std::string &base_sel = "sel_muon",
                                            const std::string &mc_weight = "w_nominal",
                                            int nbins = 60,
                                            double xmin = -15.0,
                                            double xmax = 15.0,
                                            bool include_ext = true,
                                            bool use_logy = false,
                                            const std::string &only_nuisance = "",
                                            const std::string &output_stem = "inference_score_weight_systematics")
{
    return heron::macro::run_with_guard("plot_inference_score_weight_systematics", [&]() -> int {
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
        std::cout << "[plot_inference_score_weight_systematics] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inference_score_weight_systematics] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode base = define_score0(SelectionService::decorate(el.rdf()));
        if (!has_column(base, "score0_for_plot"))
        {
            std::cerr << "[plot_inference_score_weight_systematics] could not resolve a score column. "
                      << "Expected one of: inf_score_0, inf_scores, inference_score, inference_scores.\n";
            return 1;
        }

        base = apply_selection(base, base_sel);

        auto mask_ext = el.mask_for_ext();
        auto mask_mc_like = el.mask_for_mc_like();

        ROOT::RDF::RNode node_all = filter_by_sample_mask(base, mask_mc_like, "sample_id");
        ROOT::RDF::RNode node_mc = filter_not_sample_mask(node_all, mask_ext, "sample_id")
                                      .Define("__w_nom__", mc_weight);

        ROOT::RDF::RNode node_ext = node_all;
        if (include_ext)
            node_ext = filter_by_sample_mask(node_all, mask_ext, "sample_id").Define("__w_nom__", mc_weight);

        const std::vector<VariationSpec> all_specs = {
            {"RPA_CCQE", "knobRPAup", "knobRPAdn"},
            {"XSecShape_CCMEC", "knobCCMECup", "knobCCMECdn"},
            {"AxFFCCQEshape", "knobAxFFCCQEup", "knobAxFFCCQEdn"},
            {"VecFFCCQEshape", "knobVecFFCCQEup", "knobVecFFCCQEdn"},
            {"DecayAngMEC", "knobDecayAngMECup", "knobDecayAngMECdn"},
            {"Theta_Delta2Npi", "knobThetaDelta2Npiup", "knobThetaDelta2Npidn"},
            {"ThetaDelta2NRad", "knobThetaDelta2NRadup", "knobThetaDelta2NRaddn"},
            {"NormCCCOH", "knobNormCCCOHup", "knobNormCCCOHdn"},
            {"NormNCCOH", "knobNormNCCOHup", "knobNormNCCOHdn"},
            {"xsr_scc_Fv3", "knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn"},
            {"xsr_scc_Fa3", "knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn"},
        };

        ROOT::RDF::TH1DModel h_nom_mc_model("h_nom_mc_total", "", nbins, xmin, xmax);
        auto h_nom_mc = node_mc.Histo1D(h_nom_mc_model, "score0_for_plot", "__w_nom__");

        ROOT::RDF::RResultPtr<TH1D> h_nom_ext;
        if (include_ext)
        {
            ROOT::RDF::TH1DModel h_nom_ext_model("h_nom_ext_total", "", nbins, xmin, xmax);
            h_nom_ext = node_ext.Histo1D(h_nom_ext_model, "score0_for_plot", "__w_nom__");
        }

        std::vector<BookedVariation> booked;
        booked.reserve(all_specs.size());

        for (const auto &spec : all_specs)
        {
            if (!matches_filter(spec, only_nuisance))
                continue;

            if (!has_column(node_mc, spec.up_col) || !has_column(node_mc, spec.dn_col))
            {
                std::cout << "[plot_inference_score_weight_systematics] skipping " << spec.tag
                          << " because one or both columns are missing: "
                          << spec.up_col << ", " << spec.dn_col << "\n";
                continue;
            }

            const std::string safe = make_safe_name(spec.tag);
            const std::string w_up_name = "__w_up_" + safe;
            const std::string w_dn_name = "__w_dn_" + safe;

            ROOT::RDF::RNode node_var = node_mc
                                            .Define(w_up_name,
                                                    [](double w_nom, double w_var)
                                                    {
                                                        return w_nom * sanitise_variation_weight(w_var);
                                                    },
                                                    {"__w_nom__", spec.up_col})
                                            .Define(w_dn_name,
                                                    [](double w_nom, double w_var)
                                                    {
                                                        return w_nom * sanitise_variation_weight(w_var);
                                                    },
                                                    {"__w_nom__", spec.dn_col});

            ROOT::RDF::TH1DModel h_up_model(("h_up_mc_" + safe).c_str(), "", nbins, xmin, xmax);
            ROOT::RDF::TH1DModel h_dn_model(("h_dn_mc_" + safe).c_str(), "", nbins, xmin, xmax);

            booked.push_back(BookedVariation{
                spec,
                node_var.Histo1D(h_up_model, "score0_for_plot", w_up_name),
                node_var.Histo1D(h_dn_model, "score0_for_plot", w_dn_name)});
        }

        if (booked.empty())
        {
            std::cerr << "[plot_inference_score_weight_systematics] no nuisance variations were booked.\n";
            return 1;
        }

        const TH1D &nom_mc = *h_nom_mc;
        const TH1D *nom_ext = include_ext ? &(*h_nom_ext) : nullptr;

        for (const auto &entry : booked)
        {
            draw_one_variation(entry.spec,
                               nom_mc,
                               nom_ext,
                               *entry.h_up_mc,
                               *entry.h_dn_mc,
                               xmin,
                               xmax,
                               use_logy,
                               output_stem);
        }

        std::cout << "[plot_inference_score_weight_systematics] done\n";
        return 0;
    });
}
