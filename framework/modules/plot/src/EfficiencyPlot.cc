/* -- C++ -- */
/**
 *  @file  framework/plot/src/EfficiencyPlot.cc
 *
 *  @brief Efficiency plotting helper implementation.
 */

#include "EfficiencyPlot.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "PlotEnv.hh"


namespace nu
{
namespace
{

void apply_env_defaults(Options &opt)
{
    if (opt.out_dir.empty())
    {
        opt.out_dir = plot_output_dir();
    }
    if (opt.image_format.empty())
    {
        opt.image_format = plot_image_format();
    }
}

void set_global_style_like_plotter()
{
    const int font_style = 42;
    TStyle *style = gROOT->GetStyle("PlotterStyle");
    if (style == nullptr)
    {
        style = new TStyle("PlotterStyle", "Plotter Style");
    }

    style->SetTitleFont(font_style, "X");
    style->SetTitleFont(font_style, "Y");
    style->SetTitleFont(font_style, "Z");
    style->SetTitleSize(0.04, "X");
    style->SetTitleSize(0.04, "Y");
    style->SetTitleSize(0.05, "Z");
    style->SetLabelFont(font_style, "X");
    style->SetLabelFont(font_style, "Y");
    style->SetLabelFont(font_style, "Z");
    style->SetLabelSize(0.045, "X");
    style->SetLabelSize(0.045, "Y");
    style->SetLabelSize(0.045, "Z");
    style->SetLabelOffset(0.005, "X");
    style->SetLabelOffset(0.005, "Y");
    style->SetLabelOffset(0.005, "Z");
    style->SetTitleOffset(1.10, "X");
    style->SetTitleOffset(1.25, "Y");
    style->SetOptStat(0);
    style->SetOptTitle(0);
    style->SetPadTickX(1);
    style->SetPadTickY(1);
    TGaxis::SetMaxDigits(4);
    style->SetPadLeftMargin(0.15);
    style->SetPadRightMargin(0.06);
    style->SetPadTopMargin(0.07);
    style->SetPadBottomMargin(0.12);
    style->SetMarkerSize(1.0);
    style->SetCanvasColor(0);
    style->SetPadColor(0);
    style->SetFrameFillColor(0);
    style->SetCanvasBorderMode(0);
    style->SetPadBorderMode(0);
    style->SetStatColor(0);
    style->SetFrameBorderMode(0);
    style->SetTitleFillColor(0);
    style->SetTitleBorderSize(0);

    gROOT->SetStyle("PlotterStyle");
    gROOT->ForceStyle();
}

} // namespace

EfficiencyPlot::EfficiencyPlot(TH1DModel spec, Options opt)
    : EfficiencyPlot(std::move(spec), std::move(opt), Config())
{
}

EfficiencyPlot::EfficiencyPlot(TH1DModel spec, Options opt, Config cfg)
    : spec_(std::move(spec)), opt_(std::move(opt)), cfg_(std::move(cfg))
{
    apply_env_defaults(opt_);
}

std::string EfficiencyPlot::sanitise_(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s)
    {
        if (std::isalnum(c) || c == '_' || c == '-')
        {
            out.push_back(static_cast<char>(c));
        }
        else if (c == '.' || c == ' ')
        {
            out.push_back('_');
        }
        else
        {
            out.push_back('_');
        }
    }
    if (out.empty())
    {
        out = "plot";
    }
    return out;
}

int EfficiencyPlot::compute(ROOT::RDF::RNode base,
                            const std::string &denom_sel,
                            const std::string &pass_sel,
                            const std::string &extra_sel)
{
    ROOT::RDF::RNode denom = base;
    if (!extra_sel.empty())
    {
        denom = denom.Filter(extra_sel);
    }
    if (!denom_sel.empty())
    {
        denom = denom.Filter(denom_sel);
    }

    ROOT::RDF::RNode numer = denom;
    if (!pass_sel.empty())
    {
        numer = numer.Filter(pass_sel);
    }

    return compute(denom, numer);
}

int EfficiencyPlot::compute(ROOT::RDF::RNode denom_node, ROOT::RDF::RNode pass_node)
{
    ready_ = false;
    h_total_.reset();
    h_passed_.reset();
    g_eff_.reset();
    n_denom_ = 0;
    n_pass_ = 0;

    const std::string nan_guard = spec_.expr + " == " + spec_.expr;
    ROOT::RDF::RNode denom_finite = denom_node.Filter(nan_guard);
    ROOT::RDF::RNode pass_finite = pass_node.Filter(nan_guard);

    n_denom_ = denom_finite.Count().GetValue();
    if (n_denom_ == 0)
    {
        std::cout << "[EfficiencyPlot] skip " << spec_.expr
                  << " (no denom entries after selection)\n";
        return 0;
    }
    n_pass_ = pass_finite.Count().GetValue();

    int nbins = spec_.nbins;
    double xmin = spec_.xmin;
    double xmax = spec_.xmax;

    if (cfg_.auto_x_range)
    {
        const double vmin = static_cast<double>(denom_finite.Min(spec_.expr).GetValue());
        const double vmax = static_cast<double>(denom_finite.Max(spec_.expr).GetValue());
        if (std::isfinite(vmin) && std::isfinite(vmax) && vmax > vmin)
        {
            const double span = vmax - vmin;
            const double pad = (span > 0.0 ? cfg_.x_pad_fraction * span : 1.0);
            xmin = vmin - pad;
            xmax = vmax + pad;
        }
    }

    const std::string tag = sanitise_(spec_.expr);
    const std::string htot_name = "h_eff_total_" + tag;
    const std::string hpas_name = "h_eff_passed_" + tag;

    ROOT::RDF::RResultPtr<TH1D> htot_r;
    ROOT::RDF::RResultPtr<TH1D> hpas_r;

    if (!spec_.weight.empty())
    {
        htot_r = denom_finite.Histo1D({htot_name.c_str(), "", nbins, xmin, xmax},
                                      spec_.expr,
                                      spec_.weight);
        hpas_r = pass_finite.Histo1D({hpas_name.c_str(), "", nbins, xmin, xmax},
                                     spec_.expr,
                                     spec_.weight);
    }
    else
    {
        htot_r = denom_finite.Histo1D({htot_name.c_str(), "", nbins, xmin, xmax},
                                      spec_.expr);
        hpas_r = pass_finite.Histo1D({hpas_name.c_str(), "", nbins, xmin, xmax},
                                     spec_.expr);
    }

    TH1D *htot = htot_r.GetPtr();
    TH1D *hpas = hpas_r.GetPtr();
    if (htot == nullptr || hpas == nullptr)
    {
        std::cerr << "[EfficiencyPlot] null histogram pointer for " << spec_.expr
                  << " (htot=" << htot << ", hpas=" << hpas << ")\n";
        return 1;
    }

    h_total_.reset(static_cast<TH1D *>(htot->Clone((htot_name + "_clone").c_str())));
    h_passed_.reset(static_cast<TH1D *>(hpas->Clone((hpas_name + "_clone").c_str())));
    if (!h_total_ || !h_passed_)
    {
        std::cerr << "[EfficiencyPlot] failed to clone histograms for " << spec_.expr << "\n";
        return 1;
    }
    h_total_->SetDirectory(nullptr);
    h_passed_->SetDirectory(nullptr);

    if (!TEfficiency::CheckConsistency(*h_passed_, *h_total_))
    {
        std::cerr << "[EfficiencyPlot] TEfficiency consistency check failed for " << spec_.expr
                  << " (passed > total in at least one bin, or negative contents)\n";
        return 1;
    }

    TEfficiency eff(*h_passed_, *h_total_);
    eff.SetConfidenceLevel(cfg_.conf_level);
    eff.SetStatisticOption(cfg_.stat);
    if (cfg_.use_weighted_events || !spec_.weight.empty())
    {
        eff.SetUseWeightedEvents();
    }

    if (cfg_.stat == TEfficiency::kBJeffrey ||
        cfg_.stat == TEfficiency::kBUniform ||
        cfg_.stat == TEfficiency::kBBayesian)
    {
        eff.SetPosteriorMode();
    }

    g_eff_.reset(eff.CreateGraph());
    if (!g_eff_)
    {
        std::cerr << "[EfficiencyPlot] failed to create graph for " << spec_.expr << "\n";
        return 1;
    }

    ready_ = true;
    return 0;
}

int EfficiencyPlot::draw_and_save(const std::string &file_stem,
                                  const std::string &format) const
{
    if (!ready_)
    {
        return 0;
    }

    if (!h_total_ || !g_eff_)
    {
        std::cerr << "[EfficiencyPlot] draw_and_save called before compute()\n";
        return 1;
    }

    set_global_style_like_plotter();
    gROOT->SetBatch(true);

    const std::string out_dir = opt_.out_dir.empty() ? plot_output_dir() : opt_.out_dir;
    const std::string out_fmt = !format.empty() ? format : opt_.image_format;
    const std::string stem = !file_stem.empty() ? sanitise_(file_stem)
                                                : ("eff_" + sanitise_(spec_.expr));

    gSystem->mkdir(out_dir.c_str(), true);

    const std::string out_path = out_dir + "/" + stem + "." + out_fmt;

    const double eff_all =
        (n_denom_ > 0 ? static_cast<double>(n_pass_) / static_cast<double>(n_denom_) : 0.0);

    const int bright_green = TColor::GetColor("#00cc00");
    const int bright_red = TColor::GetColor("#ff0000");

    TCanvas c(("c_" + stem).c_str(), "", 990, 700);

    TPad p_plot("p_plot", "p_plot", 0.0, 0.0, 1.0, 0.88);
    TPad p_leg("p_leg", "p_leg", 0.0, 0.88, 1.0, 1.0);

    p_leg.SetBottomMargin(0.0);
    // Too much top margin makes the legend pad look "floaty".
    p_leg.SetTopMargin(0.05);

    // Keep the top gap small so the plot frame sits close to the legend strip.
    // (Large values create an obvious white band between legend and frame.)
    p_plot.SetTopMargin(0.01);
    p_plot.SetBottomMargin(0.12);
    // Reduce side padding so the x-axis frame uses more horizontal space,
    // while still keeping axis labels/titles readable.
    p_plot.SetLeftMargin(0.11);
    // Leave enough space for the right-hand (efficiency) axis title/labels.
    // (Too small and ROOT will clip them.)
    p_plot.SetRightMargin(0.10);

    p_leg.Draw();
    p_plot.Draw();

    p_leg.cd();
    TLegend leg(0.10, 0.00, 0.90, 1.00);
    leg.SetBorderSize(0);
    leg.SetNColumns(2);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);

    p_plot.cd();
    if (cfg_.logy)
    {
        p_plot.SetLogy(true);
    }

    // Keep left-axis exponent text clear of the legend pad boundary.
    // ROOT draws exponents near the frame corner by default, which can make
    // "#times10^{n}" appear clipped when the top plot margin is tight.
    TGaxis::SetExponentOffset(-0.055, 0.02, "y");

    // Disable mirrored right-hand ticks from the primary (black) y-axis.
    // The dedicated red efficiency axis is drawn explicitly on the right.
    p_plot.SetTicky(0);
    p_plot.SetGrid(0, 0);

    std::unique_ptr<TH1D> h_tot_draw(static_cast<TH1D *>(h_total_->Clone("h_tot_draw")));
    std::unique_ptr<TH1D> h_pas_draw(h_passed_ ? static_cast<TH1D *>(h_passed_->Clone("h_pas_draw")) : nullptr);
    h_tot_draw->SetDirectory(nullptr);
    h_tot_draw->SetStats(0);
    if (h_pas_draw)
    {
        h_pas_draw->SetDirectory(nullptr);
        h_pas_draw->SetStats(0);
    }

    const std::string x_title = cfg_.x_title.empty() ? spec_.expr : cfg_.x_title;

    THStack hs("hs", (";" + x_title + ";" + cfg_.y_counts_title).c_str());

    if (cfg_.draw_distributions)
    {
        if (cfg_.draw_total_hist && h_tot_draw)
        {
            h_tot_draw->SetLineColor(1);
            h_tot_draw->SetLineWidth(2);
            h_tot_draw->SetFillStyle(0);
            // "HIST" looks much cleaner than "E0" (which produces tick-mark errorbar styling).
            hs.Add(h_tot_draw.get(), "HIST");
            leg.AddEntry(h_tot_draw.get(), cfg_.legend_total.c_str(), "L");
        }

        if (cfg_.draw_passed_hist && h_pas_draw)
        {
            h_pas_draw->SetLineColor(bright_green);
            h_pas_draw->SetLineWidth(2);
            h_pas_draw->SetFillStyle(0);
            hs.Add(h_pas_draw.get(), "HIST");
            leg.AddEntry(h_pas_draw.get(), cfg_.legend_passed.c_str(), "L");
        }

        hs.Draw("nostack HIST");

        hs.GetYaxis()->SetNoExponent(cfg_.no_exponent_y);

        if (!cfg_.logy)
        {
            const double ymax = std::max(1.0, hs.GetMaximum("nostack"));
            hs.SetMaximum(1.20 * ymax);
            hs.SetMinimum(0.0);
        }
        else
        {
            // log-y requires strictly positive minimum. Use a small floor to avoid empty pads.
            // (We avoid scanning bins here; this is a simple, robust default.)
            hs.SetMinimum(0.5);
        }

        p_plot.Update();
    }
    else
    {
        g_eff_->SetTitle((";" + x_title + ";" + cfg_.y_eff_title).c_str());
        g_eff_->Draw("AP");

        g_eff_->GetYaxis()->SetRangeUser(cfg_.eff_ymin, cfg_.eff_ymax);
        p_plot.Update();
    }

    if (cfg_.draw_distributions)
    {
        // Anchor the right-hand axis to the histogram axis limits rather than pad
        // user coordinates so the red axis always spans the full frame height.
        // (Using GetUymin/GetUymax can yield a shorter axis after ROOT applies
        // internal padding/scaling.)
        const double left_min = hs.GetMinimum("nostack");
        const double left_max = hs.GetMaximum("nostack");

        double eff_min = cfg_.eff_ymin;
        double eff_max = cfg_.eff_ymax;
        if (!(eff_max > eff_min))
        {
            eff_max = eff_min + 1.0;
        }

        const double scale = (left_max - left_min) / (eff_max - eff_min);

        TGraphAsymmErrors g_scaled(*g_eff_);
        for (int i = 0; i < g_scaled.GetN(); ++i)
        {
            double x = 0.0;
            double y = 0.0;
            g_scaled.GetPoint(i, x, y);

            const double y_left = left_min + (y - eff_min) * scale;
            const double eyl = g_scaled.GetErrorYlow(i) * scale;
            const double eyh = g_scaled.GetErrorYhigh(i) * scale;

            g_scaled.SetPoint(i, x, y_left);
            g_scaled.SetPointError(i,
                                   g_scaled.GetErrorXlow(i),
                                   g_scaled.GetErrorXhigh(i),
                                   eyl,
                                   eyh);
        }

        g_scaled.SetLineColor(bright_red);
        g_scaled.SetMarkerColor(bright_red);
        g_scaled.SetMarkerStyle(5);
        g_scaled.SetMarkerSize(1.2);
        g_scaled.SetLineWidth(2);

        // IMPORTANT:
        // ROOT keeps pointers to drawn primitives in the pad. If we draw stack objects
        // (like TGaxis/TGraphAsymmErrors locals) and they go out of scope before SaveAs(),
        // the efficiency axis/graph can disappear (or crash). DrawClone() makes the pad
        // own a heap clone with a safe lifetime.
        TGaxis axis_tmp(p_plot.GetUxmax(), left_min,
                        p_plot.GetUxmax(), left_max,
                        eff_min, eff_max, 510, "+L");
        axis_tmp.SetTitle(cfg_.y_eff_title.c_str());
        axis_tmp.SetTitleColor(bright_red);
        axis_tmp.SetLabelColor(bright_red);
        axis_tmp.SetLineColor(bright_red);
        axis_tmp.SetTitleSize(0.04);
        axis_tmp.SetLabelSize(0.045);
        axis_tmp.SetTitleOffset(1.10);
        axis_tmp.DrawClone();

        auto *g_draw = static_cast<TGraphAsymmErrors *>(g_scaled.DrawClone("P SAME"));
        if (g_draw != nullptr)
        {
            leg.AddEntry(g_draw, cfg_.legend_eff.c_str(), "LP");
        }
        else
        {
            // Fallback (shouldn't happen, but keeps legend sane if DrawClone fails)
            leg.AddEntry(g_eff_.get(), cfg_.legend_eff.c_str(), "LP");
        }
    }
    else
    {
        g_eff_->SetLineColor(bright_red);
        g_eff_->SetMarkerColor(bright_red);
        g_eff_->SetMarkerStyle(5);
        g_eff_->SetMarkerSize(1.2);
        g_eff_->SetLineWidth(2);
        leg.AddEntry(g_eff_.get(), cfg_.legend_eff.c_str(), "LP");
    }

    // Remove "stats text" from the plot (cleaner) and print to stdout instead.
    if (cfg_.print_stats)
    {
        std::cout << "[EfficiencyPlot] " << stem
                  << "  denom=" << n_denom_
                  << "  numer=" << n_pass_
                  << "  overall_eff=" << eff_all
                  << "  -> " << out_path << "\n";
    }

    // Optional user-provided annotations (kept separate from stats).
    if (!cfg_.extra_text_lines.empty() || cfg_.draw_stats_text)
    {
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.035);

        double y = 0.92;
        if (cfg_.draw_stats_text)
        {
            std::ostringstream os1;
            os1 << "Denom: " << n_denom_ << "   Numer: " << n_pass_;
            lat.DrawLatex(0.16, y, os1.str().c_str());
            y -= 0.04;

            std::ostringstream os2;
            os2 << "Overall eff: " << eff_all;
            lat.DrawLatex(0.16, y, os2.str().c_str());
            y -= 0.04;
        }

        for (const auto &line : cfg_.extra_text_lines)
        {
            lat.DrawLatex(0.16, y, line.c_str());
            y -= 0.04;
        }
    }

    p_plot.RedrawAxis();

    p_leg.cd();
    leg.Draw();

    c.SaveAs(out_path.c_str());
    return 0;
}

} // namespace nu
