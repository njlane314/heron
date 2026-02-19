/* -- C++ -- */
/**
 *  @file  plot/src/UnstackedHist.cpp
 *
 *  @brief Unstacked (overlay) histogram plotting helper.
 */

#include "UnstackedHist.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RDFHelpers.hxx"
#include "TArrow.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TLatex.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixDSym.h"
#include "TPad.h"

#include "AdaptiveBinningService.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"

namespace nu
{

// NOTE:
// These utilities are already defined in plot/src/StackedHist.cpp in namespace nu.
// We forward-declare them here to reuse the same implementation and avoid duplicate
// symbol definitions.
void apply_total_errors(TH1D &h,
                        const TMatrixDSym *cov,
                        const std::vector<double> *syst_bin,
                        bool density_mode);

double integral_in_visible_range(const TH1D &h, double xmin, double xmax);

double maximum_in_visible_range(const TH1D &h, double xmin, double xmax, bool include_error);

namespace
{
bool unstack_debug_enabled()
{
    const char *env = std::getenv("HERON_DEBUG_PLOT_UNSTACK");
    return env != nullptr && std::string(env) != "0";
}

void unstack_debug_log(const std::string &msg)
{
    if (!unstack_debug_enabled())
    {
        return;
    }
    std::cout << "[UnstackedHist][debug] " << msg << "\n";
    std::cout.flush();
}

int line_style_for_index(std::size_t i)
{
    // ROOT common line styles: 1=solid, 2=dashed, 3=dotted, 4=dash-dot, ...
    static const int styles[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    return styles[i % (sizeof(styles) / sizeof(styles[0]))];
}
} // namespace

UnstackedHist::UnstackedHist(TH1DModel spec,
                             Options opt,
                             std::vector<const Entry *> mc,
                             std::vector<const Entry *> data)
    : spec_(std::move(spec)),
      opt_(std::move(opt)),
      mc_(std::move(mc)),
      data_(std::move(data)),
      plot_name_(Plotter::sanitise(spec_.id)),
      output_directory_(opt_.out_dir)
{
}

void UnstackedHist::setup_pads(TCanvas &c, TPad *&p_main, TPad *&p_ratio, TPad *&p_legend) const
{
    auto disable_primitive_ownership = [](TPad *pad) {
        if (!pad)
        {
            return;
        }
        auto *primitives = pad->GetListOfPrimitives();
        if (primitives)
        {
            primitives->SetOwner(kFALSE);
        }
    };

    c.cd();
    p_main = nullptr;
    p_ratio = nullptr;
    p_legend = nullptr;

    const double split = 0.85;

    if (opt_.legend_on_top)
    {
        if (want_ratio())
        {
            const std::string ratio_name = plot_name_ + "_pad_ratio";
            const std::string main_name = plot_name_ + "_pad_main";
            const std::string legend_name = plot_name_ + "_pad_legend";
            p_ratio = new TPad(ratio_name.c_str(), ratio_name.c_str(), 0., 0.00, 1., 0.30);
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.30, 1., split);
            p_legend = new TPad(legend_name.c_str(), legend_name.c_str(), 0., split, 1., 1.00);

            p_main->SetTopMargin(0.02);
            p_main->SetBottomMargin(0.02);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);

            p_ratio->SetTopMargin(0.05);
            p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12);
            p_ratio->SetRightMargin(0.05);

            p_legend->SetTopMargin(0.05);
            p_legend->SetBottomMargin(0.01);
            p_legend->SetLeftMargin(0.00);
            p_legend->SetRightMargin(0.00);
        }
        else
        {
            const std::string main_name = plot_name_ + "_pad_main";
            const std::string legend_name = plot_name_ + "_pad_legend";
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.00, 1., split);
            p_legend = new TPad(legend_name.c_str(), legend_name.c_str(), 0., split, 1., 1.00);

            p_main->SetTopMargin(0.01);
            p_main->SetBottomMargin(0.12);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);

            p_legend->SetTopMargin(0.05);
            p_legend->SetBottomMargin(0.01);
        }
        if (opt_.use_log_y && p_main)
        {
            p_main->SetLogy();
        }
        if (opt_.use_log_x && p_main)
        {
            p_main->SetLogx();
        }
        if (p_ratio)
        {
            p_ratio->Draw();
        }
        if (p_main)
        {
            p_main->Draw();
            disable_primitive_ownership(p_main);
        }
        if (p_legend)
        {
            p_legend->Draw();
        }
        disable_primitive_ownership(p_ratio);
        disable_primitive_ownership(p_main);
        disable_primitive_ownership(p_legend);
    }
    else
    {
        if (want_ratio())
        {
            const std::string main_name = plot_name_ + "_pad_main";
            const std::string ratio_name = plot_name_ + "_pad_ratio";
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.30, 1., 1.);
            p_ratio = new TPad(ratio_name.c_str(), ratio_name.c_str(), 0., 0., 1., 0.30);
            p_main->SetTopMargin(0.06);
            p_main->SetBottomMargin(0.02);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);
            p_ratio->SetTopMargin(0.05);
            p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12);
            p_ratio->SetRightMargin(0.05);
            if (opt_.use_log_y)
            {
                p_main->SetLogy();
            }
            if (opt_.use_log_x)
            {
                p_main->SetLogx();
            }
            p_ratio->Draw();
            p_main->Draw();
            disable_primitive_ownership(p_ratio);
            disable_primitive_ownership(p_main);
        }
        else
        {
            const std::string main_name = plot_name_ + "_pad_main";
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0., 1., 1.);
            p_main->SetTopMargin(0.06);
            p_main->SetBottomMargin(0.12);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);
            if (opt_.use_log_y)
            {
                p_main->SetLogy();
            }
            if (opt_.use_log_x)
            {
                p_main->SetLogx();
            }
            p_main->Draw();
            disable_primitive_ownership(p_main);
        }
    }
}

void UnstackedHist::build_histograms()
{
    unstack_debug_log("build_histograms: begin spec=" + spec_.id +
                      " expr=" + (spec_.expr.empty() ? std::string("<id>") : spec_.expr) +
                      " mc_entries=" + std::to_string(mc_.size()) +
                      " data_entries=" + std::to_string(data_.size()));

    const auto axes = spec_.axis_title();
    overlay_ = std::make_unique<THStack>((spec_.id + "_overlay").c_str(), axes.c_str());
    if (overlay_->GetHists())
    {
        overlay_->GetHists()->SetOwner(kFALSE);
    }

    mc_ch_hists_.clear();
    chan_order_.clear();
    chan_event_yields_.clear();
    mc_total_.reset();
    mc_unc_band_.reset();
    data_hist_.reset();
    sig_hist_.reset();
    ratio_hist_.reset();
    ratio_band_.reset();
    legend_.reset();
    legend_proxies_.clear();

    signal_events_ = 0.0;
    signal_scale_ = 1.0;
    total_mc_events_ = 0.0;
    density_mode_ = false;

    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> booked;
    const std::vector<int> channels = opt_.unstack_channel_keys.empty()
                                          ? Channels::mc_keys()
                                          : opt_.unstack_channel_keys;
    const auto cfg = AdaptiveBinningService::config_from(opt_);

    // If adaptive binning is enabled, fill a finer histogram first and then merge.
    TH1DModel fill_spec = spec_;
    if (cfg.enabled)
    {
        const int factor = std::max(1, opt_.adaptive_fine_bin_factor);
        long long nb = 1LL * fill_spec.nbins * factor;
        nb = std::max<long long>(1, std::min<long long>(nb, 5000));
        fill_spec.nbins = static_cast<int>(nb);
    }

    for (size_t ie = 0; ie < mc_.size(); ++ie)
    {
        unstack_debug_log("build_histograms: booking MC source index=" + std::to_string(ie));
        const Entry *e = mc_[ie];
        if (!e)
        {
            continue;
        }
        auto n0 = apply(e->rnode(), spec_.sel);
        auto n = (spec_.expr.empty() ? n0 : n0.Define("_nx_expr_", spec_.expr));
        const std::string var = spec_.expr.empty() ? spec_.id : "_nx_expr_";

        for (int ch : channels)
        {
            auto nf = n.Filter([ch](int c) { return c == ch; }, {opt_.channel_column});
            auto h = nf.Histo1D(fill_spec.model("_mc_ch" + std::to_string(ch) + "_src" + std::to_string(ie)),
                                var,
                                spec_.weight);
            booked[ch].push_back(h);
        }
    }

    unstack_debug_log("build_histograms: booked channels=" + std::to_string(booked.size()));

    std::map<int, std::unique_ptr<TH1D>> sum_by_channel;

    for (int ch : channels)
    {
        auto it = booked.find(ch);
        if (it == booked.end() || it->second.empty())
        {
            continue;
        }
        std::unique_ptr<TH1D> sum;
        for (auto &rr : it->second)
        {
            const TH1D &h = rr.GetValue();
            if (!sum)
            {
                sum.reset(static_cast<TH1D *>(h.Clone((spec_.id + "_mc_sum_ch" + std::to_string(ch)).c_str())));
                sum->SetDirectory(nullptr);
            }
            else
            {
                sum->Add(&h);
            }
        }
        if (sum)
        {
            sum_by_channel.emplace(ch, std::move(sum));
        }
    }

    // ---- Adaptive min-stat-per-bin rebin (derived from TOTAL MC) ----
    std::vector<double> adaptive_edges;
    if (cfg.enabled && !sum_by_channel.empty())
    {
        std::vector<const TH1D *> parts;
        parts.reserve(sum_by_channel.size());
        for (auto &kv : sum_by_channel)
        {
            if (kv.second)
            {
                parts.push_back(kv.second.get());
            }
        }

        auto mc_tot_fine = AdaptiveBinningService::instance().sum_hists(
            parts, spec_.id + "_mc_total_fine_for_binning", cfg);

        if (mc_tot_fine)
        {
            adaptive_edges = AdaptiveBinningService::instance().edges_min_stat(*mc_tot_fine, cfg);
            if (adaptive_edges.size() >= 2)
            {
                for (auto &kv : sum_by_channel)
                {
                    if (!kv.second)
                    {
                        continue;
                    }
                    kv.second = AdaptiveBinningService::instance().rebin_to_edges(
                        *kv.second,
                        adaptive_edges,
                        spec_.id + "_mc_sum_ch" + std::to_string(kv.first) + "_rebin",
                        cfg);
                }
            }
        }
    }

    // ---- Compute yields AFTER rebin so legend ordering matches what you plot ----
    std::vector<std::pair<int, double>> yields;
    yields.reserve(sum_by_channel.size());
    for (auto &kv : sum_by_channel)
    {
        if (!kv.second)
        {
            continue;
        }
        yields.emplace_back(kv.first, kv.second->Integral());
    }

    std::stable_sort(yields.begin(), yields.end(), [](const auto &a, const auto &b) {
        if (a.second == b.second)
        {
            return a.first < b.first;
        }
        return a.second > b.second;
    });

    // Style + add to overlay stack (nostack draw later).
    unstack_debug_log("build_histograms: sorted channel yields count=" + std::to_string(yields.size()));

    for (size_t i = 0; i < yields.size(); ++i)
    {
        const int ch = yields[i].first;
        const double yld = yields[i].second;

        auto it = sum_by_channel.find(ch);
        if (it == sum_by_channel.end() || !it->second)
        {
            continue;
        }

        auto &sum = it->second;

        // Unstacked visual: lines only to avoid filled-overlap hiding earlier curves.
        sum->SetFillStyle(0);
        int line_colour = Channels::colour(ch);
        auto colour_it = opt_.unstack_channel_colours.find(ch);
        if (colour_it != opt_.unstack_channel_colours.end())
        {
            line_colour = colour_it->second;
        }
        sum->SetLineColor(line_colour);
        sum->SetLineWidth(2);
        sum->SetLineStyle(line_style_for_index(i));

        overlay_->Add(sum.get(), "HIST");
        if (auto *hists = overlay_->GetHists())
        {
            hists->SetOwner(kFALSE);
        }

        mc_ch_hists_.push_back(std::move(sum));
        chan_order_.push_back(ch);
        chan_event_yields_.push_back(yld); // pre-density event yield
    }

    // Total MC
    for (auto &uptr : mc_ch_hists_)
    {
        if (!uptr)
        {
            continue;
        }
        if (!mc_total_)
        {
            mc_total_.reset(static_cast<TH1D *>(uptr->Clone((spec_.id + "_mc_total").c_str())));
            mc_total_->SetDirectory(nullptr);
        }
        else
        {
            mc_total_->Add(uptr.get());
        }
    }

    total_mc_events_ = mc_total_ ? mc_total_->Integral() : 0.0;
    unstack_debug_log("build_histograms: mc channels=" + std::to_string(mc_ch_hists_.size()) +
                      " total_mc_events=" + std::to_string(total_mc_events_));

    // Data
    if (!data_.empty())
    {
        std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
        for (size_t ie = 0; ie < data_.size(); ++ie)
        {
            const Entry *e = data_[ie];
            if (!e)
            {
                continue;
            }
            auto n0 = apply(e->rnode(), spec_.sel);
            auto n = (spec_.expr.empty() ? n0 : n0.Define("_nx_expr_", spec_.expr));
            const std::string var = spec_.expr.empty() ? spec_.id : "_nx_expr_";
            parts.push_back(n.Histo1D(fill_spec.model("_data_src" + std::to_string(ie)), var));
        }
        for (auto &rr : parts)
        {
            const TH1D &h = rr.GetValue();
            if (!data_hist_)
            {
                data_hist_.reset(static_cast<TH1D *>(h.Clone((spec_.id + "_data").c_str())));
                data_hist_->SetDirectory(nullptr);
            }
            else
            {
                data_hist_->Add(&h);
            }
        }
        if (data_hist_ && !adaptive_edges.empty())
        {
            data_hist_ = AdaptiveBinningService::instance().rebin_to_edges(
                *data_hist_, adaptive_edges, spec_.id + "_data_rebin", cfg);
        }
        if (data_hist_)
        {
            data_hist_->SetMarkerStyle(kFullCircle);
            data_hist_->SetMarkerSize(0.9);
            data_hist_->SetLineColor(kBlack);
            data_hist_->SetFillStyle(0);
        }
    }

    // Signal overlay (scaled to total in visible range), same logic as StackedHist.
    if (opt_.overlay_signal && !opt_.signal_channels.empty() && !mc_ch_hists_.empty() && mc_total_)
    {
        const double tot_sum = integral_in_visible_range(*mc_total_, spec_.xmin, spec_.xmax);

        auto sig = std::make_unique<TH1D>(*mc_ch_hists_.front());
        sig->Reset();

        for (size_t i = 0; i < mc_ch_hists_.size(); ++i)
        {
            const int ch = chan_order_.at(i);
            if (std::find(opt_.signal_channels.begin(), opt_.signal_channels.end(), ch) != opt_.signal_channels.end())
            {
                sig->Add(mc_ch_hists_[i].get());
            }
        }

        signal_events_ = integral_in_visible_range(*sig, spec_.xmin, spec_.xmax);
        const double sig_sum = signal_events_;

        if (sig_sum > 0.0 && tot_sum > 0.0)
        {
            signal_scale_ = tot_sum / sig_sum;
            sig->Scale(signal_scale_);
        }

        sig->SetLineColor(kGreen + 2);
        sig->SetLineStyle(kDashed);
        sig->SetLineWidth(3);
        sig->SetFillStyle(0);

        sig_hist_ = std::move(sig);
    }

    // If we used adaptive edges, bins are generally not uniform: draw densities.
    const bool do_density = (cfg.enabled && adaptive_edges.size() >= 2);
    if (do_density)
    {
        auto scale_width = [](TH1D &h) { h.Scale(1.0, "width"); };

        for (auto &h : mc_ch_hists_)
        {
            if (h)
            {
                scale_width(*h);
            }
        }
        if (mc_total_)
        {
            scale_width(*mc_total_);
        }
        if (data_hist_)
        {
            scale_width(*data_hist_);
        }
        if (sig_hist_)
        {
            scale_width(*sig_hist_);
        }

        density_mode_ = true;
    }
}

void UnstackedHist::draw_overlay_and_unc(TPad *p_main, double &max_y)
{
    if (!p_main)
    {
        unstack_debug_log("draw_overlay_and_unc: main pad is null");
        return;
    }
    p_main->cd();

    if (!overlay_)
    {
        unstack_debug_log("draw_overlay_and_unc: overlay is null");
        return;
    }

    unstack_debug_log("draw_overlay_and_unc: drawing overlay with " + std::to_string(mc_ch_hists_.size()) +
                      " channel histograms");

    // Titles
    if (!opt_.x_title.empty() || !opt_.y_title.empty())
    {
        const std::string default_x = !spec_.name.empty() ? spec_.name : spec_.id;
        const std::string x = opt_.x_title.empty() ? default_x : opt_.x_title;
        const std::string y = opt_.y_title.empty() ? "Events" : opt_.y_title;
        overlay_->SetTitle((";" + x + ";" + y).c_str());
    }

    // Apply total errors if requested (needed for band + ratio band).
    if (mc_total_ && (opt_.total_cov || !opt_.syst_bin.empty()))
    {
        apply_total_errors(*mc_total_, opt_.total_cov.get(),
                           opt_.syst_bin.empty() ? nullptr : &opt_.syst_bin,
                           density_mode_);
    }

    // Determine y-max in visible range.
    max_y = 0.0;
    for (const auto &h : mc_ch_hists_)
    {
        if (!h)
        {
            continue;
        }
        max_y = std::max(max_y, maximum_in_visible_range(*h, spec_.xmin, spec_.xmax, false));
    }
    if (mc_total_)
    {
        max_y = std::max(max_y, maximum_in_visible_range(*mc_total_, spec_.xmin, spec_.xmax, true));
    }
    if (sig_hist_)
    {
        max_y = std::max(max_y, maximum_in_visible_range(*sig_hist_, spec_.xmin, spec_.xmax, false));
    }
    if (has_data())
    {
        max_y = std::max(max_y, maximum_in_visible_range(*data_hist_, spec_.xmin, spec_.xmax, true));
    }
    if (opt_.y_max > 0.0)
    {
        max_y = opt_.y_max;
    }
    if (max_y <= 0.0)
    {
        max_y = 1.0;
    }

    overlay_->SetMaximum(max_y * (opt_.use_log_y ? 10. : 1.3));
    overlay_->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);

    // Draw the overlay (THStack in "nostack" mode).
    overlay_->Draw("NOSTACK HIST");

    TH1 *frame = overlay_->GetHistogram();
    if (frame)
    {
        frame->SetLineWidth(2);
        if (spec_.xmin < spec_.xmax)
        {
            frame->GetXaxis()->SetRangeUser(spec_.xmin, spec_.xmax);
        }
        frame->GetXaxis()->SetNdivisions(510);
        frame->GetXaxis()->SetTickLength(0.02);
        frame->GetXaxis()->CenterTitle(false);

        if (!opt_.x_title.empty())
        {
            frame->GetXaxis()->SetTitle(opt_.x_title.c_str());
        }
        if (!opt_.y_title.empty())
        {
            frame->GetYaxis()->SetTitle(opt_.y_title.c_str());
        }
    }

    // Uncertainty band for total MC (draw first, then redraw lines on top).
    mc_unc_band_.reset();
    if (mc_total_)
    {
        mc_unc_band_.reset(static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_mc_totband").c_str())));
        mc_unc_band_->SetDirectory(nullptr);
        mc_unc_band_->SetFillColor(kBlack);
        mc_unc_band_->SetFillStyle(3004);
        mc_unc_band_->SetMarkerSize(0);
        mc_unc_band_->SetLineColor(kBlack);
        mc_unc_band_->SetLineWidth(1);
        mc_unc_band_->Draw("E2 SAME");
    }

    // Redraw channel curves so they stay on top of the band.
    for (auto &h : mc_ch_hists_)
    {
        if (h)
        {
            h->Draw("HIST SAME");
        }
    }

    // Total MC central value (black line)
    if (mc_total_)
    {
        mc_total_->SetFillStyle(0);
        mc_total_->SetLineColor(kBlack);
        mc_total_->SetLineWidth(2);
        mc_total_->SetLineStyle(1);
        mc_total_->Draw("HIST SAME");
    }

    // Signal overlay
    if (sig_hist_)
    {
        sig_hist_->Draw("HIST SAME");
    }

    // Data
    if (has_data())
    {
        data_hist_->Draw("E1 SAME");
    }

    // THStack can refresh internal axis state; re-apply explicit axis titles if set.
    if (!opt_.x_title.empty() && overlay_->GetXaxis())
    {
        overlay_->GetXaxis()->SetTitle(opt_.x_title.c_str());
        overlay_->GetXaxis()->CenterTitle(false);
    }
    if (!opt_.y_title.empty() && overlay_->GetYaxis())
    {
        overlay_->GetYaxis()->SetTitle(opt_.y_title.c_str());
    }
}

void UnstackedHist::draw_ratio(TPad *p_ratio)
{
    if (!p_ratio || !has_data() || !mc_total_)
    {
        return;
    }
    p_ratio->cd();

    ratio_hist_.reset(static_cast<TH1D *>(data_hist_->Clone((spec_.id + "_ratio").c_str())));
    ratio_hist_->SetDirectory(nullptr);
    ratio_hist_->Divide(mc_total_.get());
    ratio_hist_->SetTitle("; ;Data / MC");
    ratio_hist_->SetMaximum(1.99);
    ratio_hist_->SetMinimum(0.01);
    ratio_hist_->GetYaxis()->SetNdivisions(505);
    ratio_hist_->GetXaxis()->SetLabelSize(0.10);
    ratio_hist_->GetYaxis()->SetLabelSize(0.10);
    ratio_hist_->GetYaxis()->SetTitleSize(0.10);
    ratio_hist_->GetYaxis()->SetTitleOffset(0.4);
    ratio_hist_->GetYaxis()->SetTitle("Data / MC");
    ratio_hist_->GetXaxis()->CenterTitle(false);
    ratio_hist_->GetXaxis()->SetTitle(opt_.x_title.empty() ? data_hist_->GetXaxis()->GetTitle() : opt_.x_title.c_str());

    ratio_hist_->Draw("E1");

    if (opt_.show_ratio_band)
    {
        ratio_band_.reset(static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_ratio_band").c_str())));
        ratio_band_->SetDirectory(nullptr);
        const int nb = ratio_band_->GetNbinsX();
        for (int i = 1; i <= nb; ++i)
        {
            const double m = mc_total_->GetBinContent(i);
            const double em = mc_total_->GetBinError(i);
            ratio_band_->SetBinContent(i, 1.0);
            ratio_band_->SetBinError(i, (m > 0 ? em / m : 0.0));
        }
        ratio_band_->SetFillColor(kBlack);
        ratio_band_->SetFillStyle(3004);
        ratio_band_->SetLineColor(kBlack);
        ratio_band_->SetMarkerSize(0);
        ratio_band_->Draw("E2 SAME");
    }
    else
    {
        ratio_band_.reset();
    }

    ratio_hist_->Draw("E1 SAME");
}

void UnstackedHist::draw_legend(TPad *p)
{
    if (!p)
    {
        unstack_debug_log("draw_legend: legend pad is null");
        return;
    }
    unstack_debug_log("draw_legend: begin mc_ch_hists=" + std::to_string(mc_ch_hists_.size()) +
                      " chan_order=" + std::to_string(chan_order_.size()) +
                      " chan_yields=" + std::to_string(chan_event_yields_.size()));
    p->cd();

    legend_ = std::make_unique<TLegend>(0.12, 0.0, 0.95, 0.75);
    auto *leg = legend_.get();
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

    int n_entries = static_cast<int>(mc_ch_hists_.size());
    if (mc_total_)
    {
        n_entries += 2; // total line + unc band
    }
    if (sig_hist_)
    {
        ++n_entries;
    }
    if (has_data())
    {
        ++n_entries;
    }
    if (n_entries > 0)
    {
        const int n_cols = (n_entries > 4) ? 3 : 2;
        leg->SetNColumns(n_cols);
    }

    legend_proxies_.clear();

    for (size_t i = 0; i < mc_ch_hists_.size(); ++i)
    {
        if (i >= chan_order_.size())
        {
            unstack_debug_log("draw_legend: skipping histogram index " + std::to_string(i) +
                              " because chan_order size is " + std::to_string(chan_order_.size()));
            continue;
        }

        if (!mc_ch_hists_[i])
        {
            unstack_debug_log("draw_legend: skipping null histogram at index " + std::to_string(i));
            continue;
        }

        const int ch = chan_order_.at(i);

        double sum = 0.0;
        if (i < chan_event_yields_.size())
        {
            sum = chan_event_yields_[i];
        }
        else if (mc_ch_hists_[i])
        {
            sum = density_mode_ ? mc_ch_hists_[i]->Integral("width") : mc_ch_hists_[i]->Integral();
        }

        std::string label = Channels::label(ch);
        auto label_it = opt_.unstack_channel_labels.find(ch);
        if (label_it != opt_.unstack_channel_labels.end())
        {
            label = label_it->second;
        }
        if (label == "#emptyset")
        {
            label = "\xE2\x88\x85";
        }
        if (opt_.annotate_numbers)
        {
            label += " : " + Plotter::fmt_commas(sum, 2);
        }

        auto proxy = std::unique_ptr<TH1D>(
            static_cast<TH1D *>(mc_ch_hists_[i]->Clone((spec_.id + "_leg_ch" + std::to_string(ch)).c_str())));
        proxy->SetDirectory(nullptr);
        proxy->Reset("ICES");
        proxy->SetFillStyle(0);
        proxy->SetLineColor(mc_ch_hists_[i]->GetLineColor());
        proxy->SetLineWidth(mc_ch_hists_[i]->GetLineWidth());
        proxy->SetLineStyle(mc_ch_hists_[i]->GetLineStyle());

        leg->AddEntry(proxy.get(), label.c_str(), "l");
        leg->SetEntrySeparation(0.01);
        legend_proxies_.push_back(std::move(proxy));
    }

    // Total MC line + unc band
    if (mc_total_)
    {
        auto proxy_line = std::unique_ptr<TH1D>(
            static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_leg_total").c_str())));
        proxy_line->SetDirectory(nullptr);
        proxy_line->Reset("ICES");
        proxy_line->SetFillStyle(0);
        proxy_line->SetLineColor(kBlack);
        proxy_line->SetLineWidth(2);
        proxy_line->SetLineStyle(1);
        leg->AddEntry(proxy_line.get(), "Total MC", "l");
        legend_proxies_.push_back(std::move(proxy_line));

        auto proxy_band = std::unique_ptr<TH1D>(
            static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_leg_unc").c_str())));
        proxy_band->SetDirectory(nullptr);
        proxy_band->Reset("ICES");
        proxy_band->SetFillColor(kBlack);
        proxy_band->SetFillStyle(3004);
        proxy_band->SetLineColor(kBlack);
        proxy_band->SetLineWidth(1);
        leg->AddEntry(proxy_band.get(), "Stat. #oplus Syst. Unc.", "f");
        legend_proxies_.push_back(std::move(proxy_band));
    }

    if (sig_hist_)
    {
        std::ostringstream signal_label;
        signal_label << "Signal"
                     << " : " << Plotter::fmt_commas(signal_events_, 2)
                     << " (x" << std::fixed << std::setprecision(2) << signal_scale_ << ")";
        leg->AddEntry(sig_hist_.get(), signal_label.str().c_str(), "l");
    }

    if (has_data())
    {
        leg->AddEntry(data_hist_.get(), "Data", "lep");
    }

    unstack_debug_log("draw_legend: final legend entries=" + std::to_string(leg->GetNRows()));
    leg->Draw();
}

void UnstackedHist::draw_cuts(TPad *p, double max_y)
{
    if (!opt_.show_cuts || opt_.cuts.empty())
    {
        return;
    }
    if (!p)
    {
        return;
    }
    p->cd();

    TH1 *frame = overlay_ ? overlay_->GetHistogram() : nullptr;
    if (!frame)
    {
        return;
    }

    const double y = max_y * 0.80;
    const double xmin = frame->GetXaxis()->GetXmin();
    const double xmax = frame->GetXaxis()->GetXmax();
    const double xr = xmax - xmin;
    const double alen = xr * 0.04;

    for (const auto &c : opt_.cuts)
    {
        auto *line = new TLine(c.x, 0., c.x, max_y * 1.3);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(kDashed);
        line->Draw("same");

        const double x1 = c.x;
        const double x2 = (c.dir == CutDir::GreaterThan) ? c.x + alen : c.x - alen;
        auto *arr = new TArrow(x1, y, x2, y, 0.025, ">");
        arr->SetLineColor(kRed);
        arr->SetFillColor(kRed);
        arr->SetLineWidth(2);
        arr->Draw("same");
    }
}

void UnstackedHist::draw_watermark(TPad *p, double total_mc) const
{
    if (!p)
    {
        return;
    }
    p->cd();

    const std::string line1 = "#bf{#muBooNE Simulation, Preliminary}";

    const auto sum_pot = [](const std::vector<const Entry *> &entries) {
        constexpr double rel_tol = 1e-6;
        double total = 0.0;
        std::map<std::pair<std::string, std::string>, std::vector<double>> seen;
        for (const auto *e : entries)
        {
            if (!e || e->pot_nom <= 0.0)
            {
                continue;
            }
            const auto key = std::make_pair(e->beamline, e->period);
            auto &values = seen[key];
            const double pot = e->pot_nom;
            const bool already_accounted = std::any_of(values.begin(), values.end(), [&](double v) {
                const double scale = std::max(1.0, std::max(std::abs(v), std::abs(pot)));
                return std::abs(v - pot) <= rel_tol * scale;
            });
            if (!already_accounted)
            {
                values.push_back(pot);
                total += pot;
            }
        }
        return total;
    };

    double pot_value = opt_.total_protons_on_target;
    if (pot_value <= 0.0)
    {
        pot_value = sum_pot(data_);
    }
    if (pot_value <= 0.0)
    {
        pot_value = sum_pot(mc_);
    }

    auto format_pot = [](double value) {
        std::ostringstream ss;
        ss << std::scientific << std::setprecision(2) << value;
        auto text = ss.str();
        const auto pos = text.find('e');
        if (pos != std::string::npos)
        {
            const int exponent = std::stoi(text.substr(pos + 1));
            text = text.substr(0, pos) + " #times 10^{" + std::to_string(exponent) + "}";
        }
        return text;
    };

    const std::string pot_str = pot_value > 0.0 ? format_pot(pot_value) : "N/A";

    const auto first_non_empty = [](const std::vector<const Entry *> &entries, auto getter) -> std::string {
        for (const auto *e : entries)
        {
            if (!e)
            {
                continue;
            }
            auto value = getter(*e);
            if (!value.empty())
            {
                return value;
            }
        }
        return {};
    };

    auto beam_name = opt_.beamline;
    if (beam_name.empty())
    {
        beam_name = first_non_empty(data_, [](const auto &e) { return e.beamline; });
    }
    if (beam_name.empty())
    {
        beam_name = first_non_empty(mc_, [](const auto &e) { return e.beamline; });
    }
    std::string beam_key = beam_name;
    std::transform(beam_key.begin(), beam_key.end(), beam_key.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (beam_key == "numi")
    {
        beam_name = "NuMI";
    }
    else if (beam_key == "numi_fhc")
    {
        beam_name = "NuMI FHC";
    }
    else if (beam_key == "numi_rhc")
    {
        beam_name = "NuMI RHC";
    }
    if (beam_name.empty())
    {
        beam_name = "N/A";
    }

    std::vector<std::string> runs = opt_.run_numbers;
    if (runs.empty())
    {
        runs = opt_.periods;
    }
    if (runs.empty())
    {
        auto run = first_non_empty(data_, [](const auto &e) { return e.period; });
        if (run.empty())
        {
            run = first_non_empty(mc_, [](const auto &e) { return e.period; });
        }
        if (!run.empty())
        {
            runs.push_back(std::move(run));
        }
    }

    auto format_run = [](std::string label) {
        if (label.rfind("run", 0) == 0)
        {
            label.erase(0, 3);
        }
        try
        {
            label = Plotter::fmt_commas(std::stod(label), 0);
        }
        catch (...)
        {
        }
        return label;
    };

    std::string runs_str = "N/A";
    if (!runs.empty())
    {
        std::ostringstream ss;
        for (size_t i = 0; i < runs.size(); ++i)
        {
            if (i)
            {
                ss << ", ";
            }
            ss << format_run(runs[i]);
        }
        runs_str = ss.str();
    }

    std::string region_label = opt_.analysis_region_label;
    if (region_label.empty())
    {
        region_label = SelectionService::selection_label(spec_.sel);
    }

    const std::string line2 = "Beam(s), Run(s): " + beam_name + ", " + runs_str +
                              " (" + pot_str + " POT)";
    const std::string line3 = "Analysis Region: " + region_label + " (" +
                              Plotter::fmt_commas(total_mc, 2) + " events)";

    TLatex watermark;
    watermark.SetNDC();
    watermark.SetTextAlign(33);
    watermark.SetTextFont(62);
    watermark.SetTextSize(0.05);
    const double x = 1 - p->GetRightMargin() - 0.03;
    const double top = 1 - p->GetTopMargin();
    double y = top - 0.03;
    watermark.DrawLatex(x, y, line1.c_str());
    y -= 0.06;
    watermark.SetTextFont(42);
    watermark.SetTextSize(0.05 * 0.8);
    watermark.DrawLatex(x, y, line2.c_str());
    watermark.DrawLatex(x, y - 0.06, line3.c_str());
}

void UnstackedHist::draw(TCanvas &canvas)
{
    build_histograms();
    unstack_debug_log("draw: build_histograms complete");

    TPad *p_main = nullptr;
    TPad *p_ratio = nullptr;
    TPad *p_legend = nullptr;

    setup_pads(canvas, p_main, p_ratio, p_legend);

    double max_y = 1.0;
    draw_overlay_and_unc(p_main, max_y);
    unstack_debug_log("draw: draw_overlay_and_unc complete");
    draw_cuts(p_main, max_y);
    if (opt_.show_watermark)
    {
        draw_watermark(p_main, total_mc_events_);
    }

    if (opt_.show_legend)
    {
        draw_legend(p_legend ? p_legend : p_main);
    }
    if (want_ratio())
    {
        draw_ratio(p_ratio);
    }

    if (p_main)
    {
        p_main->RedrawAxis();
    }
    canvas.Update();
}

void UnstackedHist::draw_and_save(const std::string &image_format)
{
    std::filesystem::create_directories(output_directory_);
    unstack_debug_log("draw_and_save enter: plot='" + plot_name_ +
                      "', out_dir='" + output_directory_ + "'");

    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
    unstack_debug_log("canvas constructed: plot='" + plot_name_ + "'");

    draw(canvas);
    unstack_debug_log("draw finished: plot='" + plot_name_ + "'");

    const std::string fmt = image_format.empty() ? "png" : image_format;
    const std::string out = output_directory_ + "/" + plot_name_ + "." + fmt;

    unstack_debug_log("SaveAs start: plot='" + plot_name_ + "', file='" + out + "'");
    canvas.SaveAs(out.c_str());
    unstack_debug_log("SaveAs done: plot='" + plot_name_ + "'");
}

} // namespace nu
