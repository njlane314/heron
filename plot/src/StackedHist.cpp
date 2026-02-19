/* -- C++ -- */
/**
 *  @file  plot/src/StackedHist.cpp
 *
 *  @brief Stacked histogram plotting helper.
 */

#include "StackedHist.hh"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"
#include "TArrow.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TImage.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixDSym.h"

#include "AdaptiveBinningService.hh"
#include "PlotChannels.hh"
#include "ParticleChannels.hh"
#include "Plotter.hh"

namespace nu
{

namespace
{
constexpr int k_panel_fill_colour = 17;
constexpr int k_uncertainty_fill_style = 3002;

bool stack_debug_enabled()
{
    const char *env = std::getenv("HERON_DEBUG_PLOT_STACK");
    return env != nullptr && std::string(env) != "0";
}

void stack_debug_log(const std::string &msg)
{
    if (!stack_debug_enabled())
    {
        return;
    }
    std::cout << "[StackedHist][debug] " << msg << "\n";
    std::cout.flush();
}
ROOT::VecOps::RVec<double> select_particle_values(const ROOT::VecOps::RVec<double> &values,
                                                  const ROOT::VecOps::RVec<int> &pdg_codes,
                                                  int key,
                                                  bool drop_nan)
{
    ROOT::VecOps::RVec<double> out;
    out.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        const double v = values[i];
        if (drop_nan && !std::isfinite(v))
        {
            continue;
        }
        // If the PDG vector is shorter than the value vector (possible if a producer
        // did not push placeholders), treat missing PDGs as "unmatched" (0).
        const int pdg = (i < pdg_codes.size()) ? pdg_codes[i] : 0;
        if (ParticleChannels::matches(key, pdg))
        {
            out.push_back(v);
        }
    }
    return out;
}

} // namespace

void apply_total_errors(TH1D &h,
                        const TMatrixDSym *cov,
                        const std::vector<double> *syst_bin,
                        bool density_mode)
{
    const int nb = h.GetNbinsX();
    for (int i = 1; i <= nb; ++i)
    {
        const double stat = h.GetBinError(i);
        double syst = 0.0;
        if (cov && i - 1 < cov->GetNrows())
        {
            syst = std::sqrt((*cov)(i - 1, i - 1));
        }
        else if (syst_bin && i - 1 < static_cast<int>(syst_bin->size()))
        {
            syst = std::max(0.0, (*syst_bin)[i - 1]);
        }

        // If the histogram content has been scaled to a density (events / unit-x)
        // via h.Scale(1, "width"), then the systematic uncertainties must be put
        // in the same units. For diagonal-only usage, that is a 1/width factor.
        // (Full covariance needs V'_{ij} = V_{ij} / (w_i w_j).)
        if (density_mode)
        {
            const double w = h.GetXaxis() ? h.GetXaxis()->GetBinWidth(i) : 0.0;
            if (w > 0.0)
            {
                syst /= w;
            }
        }

        const double tot = std::sqrt(stat * stat + syst * syst);
        h.SetBinError(i, tot);
    }
}

double integral_in_visible_range(const TH1D &h, double xmin, double xmax)
{
    const TAxis *axis = h.GetXaxis();
    if (!axis)
    {
        return h.Integral();
    }

    int first_bin = 1;
    int last_bin = h.GetNbinsX();
    if (xmin < xmax)
    {
        first_bin = std::max(1, axis->FindFixBin(xmin));
        last_bin = std::min(h.GetNbinsX(), axis->FindFixBin(xmax));
        if (xmax <= axis->GetBinLowEdge(last_bin))
        {
            last_bin = std::max(first_bin, last_bin - 1);
        }
    }

    if (first_bin > last_bin)
    {
        return 0.0;
    }
    return h.Integral(first_bin, last_bin);
}

double maximum_in_visible_range(const TH1D &h, double xmin, double xmax, bool include_error)
{
    const TAxis *axis = h.GetXaxis();
    if (!axis)
    {
        return include_error ? h.GetMaximum() + h.GetBinError(h.GetMaximumBin()) : h.GetMaximum();
    }

    int first_bin = 1;
    int last_bin = h.GetNbinsX();
    if (xmin < xmax)
    {
        first_bin = std::max(1, axis->FindFixBin(xmin));
        last_bin = std::min(h.GetNbinsX(), axis->FindFixBin(xmax));
        if (xmax <= axis->GetBinLowEdge(last_bin))
        {
            last_bin = std::max(first_bin, last_bin - 1);
        }
    }

    if (first_bin > last_bin)
    {
        return 0.0;
    }

    double max_y = 0.0;
    for (int i = first_bin; i <= last_bin; ++i)
    {
        const double y = h.GetBinContent(i) + (include_error ? h.GetBinError(i) : 0.0);
        max_y = std::max(max_y, y);
    }
    return max_y;
}


StackedHist::StackedHist(TH1DModel spec, Options opt, std::vector<const Entry *> mc, std::vector<const Entry *> data)
    : spec_(std::move(spec)),
      opt_(std::move(opt)),
      mc_(std::move(mc)),
      data_(std::move(data)),
      plot_name_(Plotter::sanitise(spec_.id)),
      output_directory_(opt_.out_dir)
{
}

void StackedHist::setup_pads(TCanvas &c, TPad *&p_main, TPad *&p_ratio, TPad *&p_legend) const
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
    c.SetFillColor(k_panel_fill_colour);
    c.SetFrameFillColor(k_panel_fill_colour);
    p_main = nullptr;
    p_ratio = nullptr;
    p_legend = nullptr;

    const double split_top = 0.85;
    const bool separate_legend = opt_.show_legend;

    if (opt_.legend_on_top)
    {
        if (want_ratio())
        {
            const std::string ratio_name = plot_name_ + "_pad_ratio";
            const std::string main_name = plot_name_ + "_pad_main";
            const std::string legend_name = plot_name_ + "_pad_legend";
            p_ratio = new TPad(ratio_name.c_str(), ratio_name.c_str(), 0., 0.00, 1., 0.30);
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.30, 1., split_top);
            p_legend = new TPad(legend_name.c_str(), legend_name.c_str(), 0., split_top, 1., 1.00);

            p_main->SetTopMargin(0.02);
            p_main->SetBottomMargin(0.02);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);
            p_main->SetFillColor(k_panel_fill_colour);
            p_main->SetFrameFillColor(k_panel_fill_colour);

            p_ratio->SetTopMargin(0.05);
            p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12);
            p_ratio->SetRightMargin(0.05);
            p_ratio->SetFillColor(k_panel_fill_colour);
            p_ratio->SetFrameFillColor(k_panel_fill_colour);

            p_legend->SetTopMargin(0.05);
            p_legend->SetBottomMargin(0.01);
            p_legend->SetLeftMargin(0.00);
            p_legend->SetRightMargin(0.00);
            p_legend->SetFillColor(k_panel_fill_colour);
            p_legend->SetFrameFillColor(k_panel_fill_colour);
        }
        else
        {
            const std::string main_name = plot_name_ + "_pad_main";
            const std::string legend_name = plot_name_ + "_pad_legend";
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.00, 1., split_top);
            p_legend = new TPad(legend_name.c_str(), legend_name.c_str(), 0., split_top, 1., 1.00);

            p_main->SetTopMargin(0.01);
            p_main->SetBottomMargin(0.12);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);
            p_main->SetFillColor(k_panel_fill_colour);
            p_main->SetFrameFillColor(k_panel_fill_colour);

            p_legend->SetTopMargin(0.05);
            p_legend->SetBottomMargin(0.01);
            p_legend->SetFillColor(k_panel_fill_colour);
            p_legend->SetFrameFillColor(k_panel_fill_colour);
        }
        if (opt_.use_log_y && p_main)
        {
            p_main->SetLogy();
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
        return;
    }

    const double legend_split = separate_legend ? 0.78 : 1.00;

    if (want_ratio())
    {
        const std::string main_name = plot_name_ + "_pad_main";
        const std::string ratio_name = plot_name_ + "_pad_ratio";
        p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.30, legend_split, 1.);
        p_ratio = new TPad(ratio_name.c_str(), ratio_name.c_str(), 0., 0., legend_split, 0.30);

        p_main->SetTopMargin(0.06);
        p_main->SetBottomMargin(0.02);
        p_main->SetLeftMargin(0.12);
        p_main->SetRightMargin(separate_legend ? 0.03 : 0.05);
        p_main->SetFillColor(k_panel_fill_colour);
        p_main->SetFrameFillColor(k_panel_fill_colour);

        p_ratio->SetTopMargin(0.05);
        p_ratio->SetBottomMargin(0.35);
        p_ratio->SetLeftMargin(0.12);
        p_ratio->SetRightMargin(separate_legend ? 0.03 : 0.05);
        p_ratio->SetFillColor(k_panel_fill_colour);
        p_ratio->SetFrameFillColor(k_panel_fill_colour);

        if (opt_.use_log_y)
        {
            p_main->SetLogy();
        }

        if (separate_legend)
        {
            const std::string legend_name = plot_name_ + "_pad_legend";
            p_legend = new TPad(legend_name.c_str(), legend_name.c_str(), legend_split, 0., 1., 1.);
            p_legend->SetTopMargin(0.06);
            p_legend->SetBottomMargin(0.08);
            p_legend->SetLeftMargin(0.02);
            p_legend->SetRightMargin(0.02);
            p_legend->SetFillColor(k_panel_fill_colour);
            p_legend->SetFrameFillColor(k_panel_fill_colour);
        }

        p_ratio->Draw();
        p_main->Draw();
        if (p_legend)
        {
            p_legend->Draw();
        }
        disable_primitive_ownership(p_ratio);
        disable_primitive_ownership(p_main);
        disable_primitive_ownership(p_legend);
        return;
    }

    const std::string main_name = plot_name_ + "_pad_main";
    p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0., legend_split, 1.);
    p_main->SetTopMargin(0.06);
    p_main->SetBottomMargin(0.12);
    p_main->SetLeftMargin(0.12);
    p_main->SetRightMargin(separate_legend ? 0.03 : 0.05);
    p_main->SetFillColor(k_panel_fill_colour);
    p_main->SetFrameFillColor(k_panel_fill_colour);
    if (opt_.use_log_y)
    {
        p_main->SetLogy();
    }

    if (separate_legend)
    {
        const std::string legend_name = plot_name_ + "_pad_legend";
        p_legend = new TPad(legend_name.c_str(), legend_name.c_str(), legend_split, 0., 1., 1.);
        p_legend->SetTopMargin(0.06);
        p_legend->SetBottomMargin(0.08);
        p_legend->SetLeftMargin(0.02);
        p_legend->SetRightMargin(0.02);
        p_legend->SetFillColor(k_panel_fill_colour);
        p_legend->SetFrameFillColor(k_panel_fill_colour);
    }

    p_main->Draw();
    if (p_legend)
    {
        p_legend->Draw();
    }
    disable_primitive_ownership(p_main);
    disable_primitive_ownership(p_legend);
}


void StackedHist::build_histograms()
{
    const auto axes = spec_.axis_title();
    stack_ = std::make_unique<THStack>((spec_.id + "_stack").c_str(), axes.c_str());
    if (stack_->GetHists())
    {
        stack_->GetHists()->SetOwner(kFALSE);
    }
    mc_ch_hists_.clear();
    mc_total_.reset();
    data_hist_.reset();
    sig_hist_.reset();
    ratio_hist_.reset();
    ratio_band_.reset();
    signal_events_ = 0.0;
    signal_scale_ = 1.0;
    chan_event_yields_.clear();
    total_mc_events_ = 0.0;
    density_mode_ = false;
    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> booked;
    const bool particle_level = opt_.particle_level;
    const auto &channels = particle_level ? ParticleChannels::keys() : Channels::mc_keys();
    const std::string pdg_branch = particle_level
                                       ? (opt_.particle_pdg_branch.empty() ? "backtracked_pdg_codes" : opt_.particle_pdg_branch)
                                       : std::string{};
    const auto cfg = AdaptiveBinningService::config_from(opt_);

    // If adaptive binning is enabled, fill a finer histogram first and then merge.
    // Otherwise the merge algorithm has nothing to do and bins will look uniform.
    TH1DModel fill_spec = spec_;
    if (cfg.enabled)
    {
        const int factor = std::max(1, opt_.adaptive_fine_bin_factor);
        long long nb = 1LL * fill_spec.nbins * factor;
        // Clamp to something sane.
        nb = std::max<long long>(1, std::min<long long>(nb, 5000));
        fill_spec.nbins = static_cast<int>(nb);
    }

    for (size_t ie = 0; ie < mc_.size(); ++ie)
    {
        const Entry *e = mc_[ie];
        if (!e)
        {
            continue;
        }
        auto n0 = apply(e->rnode(), spec_.sel, e->selection);
        auto n = (spec_.expr.empty() ? n0 : n0.Define("_nx_expr_", spec_.expr));
        const std::string var = spec_.expr.empty() ? spec_.id : "_nx_expr_";

        if (!particle_level)
        {
            for (int ch : channels)
            {
                auto nf = n.Filter([ch](int c) { return c == ch; }, {"analysis_channels"});
                auto h = nf.Histo1D(fill_spec.model("_mc_ch" + std::to_string(ch) + "_src" + std::to_string(ie)), var, spec_.weight);
                booked[ch].push_back(h);
            }
        }
        else
        {
            // Convert any vector-like numeric column to RVec<double> once, then
            // per-category select based on the truth-matched PDG vector.
            const std::string val_d = "_nx_val_d_";
            const std::string val_expr = "ROOT::VecOps::RVec<double>(" + var + ".begin(), " + var + ".end())";
            auto nvec = n.Define(val_d, val_expr);

            for (int ch : channels)
            {
                const std::string sel = "_nx_part_sel_" + std::to_string(ch) + "_src" + std::to_string(ie);
                auto nsel = nvec.Define(
                    sel,
                    [ch, drop_nan = opt_.particle_drop_nan](const ROOT::VecOps::RVec<double> &values,
                                                            const ROOT::VecOps::RVec<int> &pdg_codes) {
                        return select_particle_values(values, pdg_codes, ch, drop_nan);
                    },
                    std::vector<std::string>{val_d, pdg_branch});

                auto h = nsel.Histo1D(fill_spec.model("_mc_pdg" + std::to_string(ch) + "_src" + std::to_string(ie)),
                                      sel, spec_.weight);
                booked[ch].push_back(h);
            }
        }
    }

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

    // ---- Compute yields AFTER rebin so ordering matches what you plot ----
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

    chan_order_.clear();
    chan_event_yields_.clear();
    for (auto &cy : yields)
    {
        const int ch = cy.first;
        auto it = sum_by_channel.find(ch);
        if (it == sum_by_channel.end() || !it->second)
        {
            continue;
        }
        auto &sum = it->second;
        if (particle_level)
        {
            sum->SetFillColor(ParticleChannels::colour(ch));
            sum->SetFillStyle(ParticleChannels::fill_style(ch));
        }
        else
        {
            sum->SetFillColor(Channels::colour(ch));
            sum->SetFillStyle(Channels::fill_style(ch));
        }
        sum->SetLineColor(kBlack);
        sum->SetLineWidth(1);
        stack_->Add(sum.get(), "HIST");
        if (auto *hists = stack_->GetHists())
        {
            // THStack owns its histogram list by default. We also own these
            // histograms via std::unique_ptr in mc_ch_hists_, so disable list
            // ownership to avoid double-deletes when StackedHist tears down.
            hists->SetOwner(kFALSE);
        }
        mc_ch_hists_.push_back(std::move(sum));
        chan_order_.push_back(ch);
        chan_event_yields_.push_back(cy.second); // pre-density event yield
    }

    for (auto &uptr : mc_ch_hists_)
    {
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

    // Cache the event-count total before any density scaling.
    total_mc_events_ = mc_total_ ? mc_total_->Integral() : 0.0;

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
            auto n0 = apply(e->rnode(), spec_.sel, e->selection);
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

    if (opt_.overlay_signal && !opt_.signal_channels.empty() && !mc_ch_hists_.empty())
    {
        double tot_sum = mc_total_ ? integral_in_visible_range(*mc_total_, spec_.xmin, spec_.xmax) : 0.0;
        auto sig = std::make_unique<TH1D>(*mc_ch_hists_.front());
        sig->Reset();
        for (size_t i = 0; i < mc_ch_hists_.size(); ++i)
        {
            int ch = chan_order_.at(i);
            if (std::find(opt_.signal_channels.begin(), opt_.signal_channels.end(), ch) != opt_.signal_channels.end())
            {
                sig->Add(mc_ch_hists_[i].get());
            }
        }
        signal_events_ = integral_in_visible_range(*sig, spec_.xmin, spec_.xmax);
        double sig_sum = signal_events_;
        if (sig_sum > 0.0 && tot_sum > 0.0)
        {
            signal_scale_ = tot_sum / sig_sum;
            sig->Scale(signal_scale_);
        }
        sig->SetLineColor(kGreen + 2);
        sig->SetLineStyle(kDashed);
        sig->SetLineWidth(2);
        sig->SetFillStyle(0);
        sig_hist_ = std::move(sig);
    }

    // ---- If we used adaptive edges, the bin widths are generally not uniform.
    // Plot densities (events / unit-x) so bar heights are comparable across bins.
    //
    // NOTE: This changes histogram units; any later error model must match those units.
    const bool do_density = (cfg.enabled && adaptive_edges.size() >= 2);
    if (do_density)
    {
        auto scale_width = [](TH1D &h) { h.Scale(1.0, "width"); };

        for (auto &h : mc_ch_hists_)
        {
            if (h)
                scale_width(*h);
        }
        if (mc_total_)
            scale_width(*mc_total_);
        if (data_hist_)
            scale_width(*data_hist_);
        if (sig_hist_)
            scale_width(*sig_hist_);

        density_mode_ = true;
    }
}

void StackedHist::draw_stack_and_unc(TPad *p_main, double &max_y)
{
    if (!p_main)
    {
        return;
    }
    p_main->cd();

    if (auto *hists = stack_->GetHists())
    {
        for (TObject *obj = hists->First(); obj != nullptr; obj = hists->After(obj))
        {
            if (auto *hist = dynamic_cast<TH1 *>(obj))
            {
                hist->SetLineColor(kBlack);
            }
        }
    }

    if (!opt_.x_title.empty() || !opt_.y_title.empty())
    {
        const std::string default_x = !spec_.name.empty() ? spec_.name : spec_.id;
        const std::string x = opt_.x_title.empty() ? default_x : opt_.x_title;
        const std::string default_y = opt_.particle_level ? "Particles" : "Events";
        const std::string y = opt_.y_title.empty() ? default_y : opt_.y_title;
        stack_->SetTitle((";" + x + ";" + y).c_str());
    }

    stack_->Draw("HIST");
    TH1 *frame = stack_->GetHistogram();
    if (frame)
    {
        frame->SetLineWidth(2);
        frame->SetFillColor(k_panel_fill_colour);
        frame->SetMarkerSize(0.0);
    }
    if (frame && spec_.xmin < spec_.xmax)
    {
        frame->GetXaxis()->SetRangeUser(spec_.xmin, spec_.xmax);
    }
    if (frame)
    {
        frame->GetXaxis()->SetNdivisions(510);
        frame->GetXaxis()->SetTickLength(0.02);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->CenterTitle(false);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleOffset(1.2);

        // If we didn't override titles explicitly, THStack's default comes from
        // TH1DModel::axis_title() (which defaults to "Events"). For particle-level
        // plots, update the default Y label to avoid misleading units.
        if (opt_.particle_level && opt_.y_title.empty())
        {
            const std::string cur = frame->GetYaxis() ? frame->GetYaxis()->GetTitle() : "";
            if (cur.empty() || cur == "Events")
            {
                frame->GetYaxis()->SetTitle("Particles");
            }
        }
        if (!opt_.x_title.empty())
        {
            frame->GetXaxis()->SetTitle(opt_.x_title.c_str());
        }
        if (!opt_.y_title.empty())
        {
            frame->GetYaxis()->SetTitle(opt_.y_title.c_str());
        }
    }
    if (mc_total_)
    {
        if (opt_.total_cov || !opt_.syst_bin.empty())
        {
            apply_total_errors(*mc_total_, opt_.total_cov.get(),
                               opt_.syst_bin.empty() ? nullptr : &opt_.syst_bin, density_mode_);
        }

        max_y = maximum_in_visible_range(*mc_total_, spec_.xmin, spec_.xmax, true);
        if (sig_hist_)
        {
            max_y = std::max(max_y, maximum_in_visible_range(*sig_hist_, spec_.xmin, spec_.xmax, false));
        }
        if (opt_.y_max > 0)
        {
            max_y = opt_.y_max;
        }

        stack_->SetMaximum(max_y * (opt_.use_log_y ? 10. : 1.3));
        stack_->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);

        auto *h = static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_mc_totband").c_str()));
        h->SetDirectory(nullptr);
        h->SetFillColor(kBlack);
        h->SetFillStyle(k_uncertainty_fill_style);
        h->SetMarkerSize(0);
        h->SetLineColor(kBlack);
        h->SetLineStyle(kDashed);
        h->SetLineWidth(2);
        h->Draw("E2 SAME");
    }

    // THStack can refresh its internal frame after range/max updates.
    // Apply axis titles and alignment after those updates so labels persist.
    if (!opt_.x_title.empty() && stack_->GetXaxis())
    {
        stack_->GetXaxis()->SetTitle(opt_.x_title.c_str());
        stack_->GetXaxis()->CenterTitle(false);
    }
    if (!opt_.y_title.empty() && stack_->GetYaxis())
    {
        stack_->GetYaxis()->SetTitle(opt_.y_title.c_str());
    }

    if (sig_hist_)
    {
        sig_hist_->Draw("HIST SAME");
    }
    if (has_data())
    {
        data_hist_->Draw("E1 SAME");
    }
}

void StackedHist::draw_ratio(TPad *p_ratio)
{
    if (!p_ratio || !has_data() || !mc_total_)
    {
        return;
    }
    p_ratio->cd();

    ratio_hist_.reset(static_cast<TH1D *>(
        data_hist_->Clone((spec_.id + "_ratio").c_str())));
    ratio_hist_->SetDirectory(nullptr);
    ratio_hist_->Divide(mc_total_.get());
    ratio_hist_->SetTitle("; ;Data / MC");
    ratio_hist_->SetMaximum(1.20);
    ratio_hist_->SetMinimum(0.80);
    ratio_hist_->GetYaxis()->SetNdivisions(505);
    ratio_hist_->GetXaxis()->SetLabelSize(0.10);
    ratio_hist_->GetXaxis()->SetTitleSize(0.12);
    ratio_hist_->GetXaxis()->SetTitleOffset(1.05);
    ratio_hist_->GetYaxis()->SetLabelSize(0.10);
    ratio_hist_->GetYaxis()->SetTitleSize(0.10);
    ratio_hist_->GetYaxis()->SetTitleOffset(0.55);
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
        ratio_band_->SetFillStyle(k_uncertainty_fill_style);
        ratio_band_->SetLineColor(kBlack);
        ratio_band_->SetLineStyle(kDashed);
        ratio_band_->SetLineWidth(2);
        ratio_band_->SetMarkerSize(0);
        ratio_band_->Draw("E2 SAME");
    }
    else
    {
        ratio_band_.reset();
    }

    auto *unity = new TLine(ratio_hist_->GetXaxis()->GetXmin(), 1.0,
                            ratio_hist_->GetXaxis()->GetXmax(), 1.0);
    unity->SetLineColor(kBlack);
    unity->SetLineStyle(kDashed);
    unity->SetLineWidth(1);
    unity->Draw("SAME");

    ratio_hist_->Draw("E1 SAME");
}

void StackedHist::draw_legend(TPad *p)
{
    if (!p)
    {
        return;
    }
    p->cd();
    const bool compact_legend = p->GetWNDC() < 0.30;
    const double x1 = compact_legend ? 0.04 : 0.12;
    const double y1 = compact_legend ? 0.03 : 0.00;
    const double x2 = 0.98;
    const double y2 = compact_legend ? 0.97 : 0.75;
    legend_ = std::make_unique<TLegend>(x1, y1, x2, y2);
    auto *leg = legend_.get();
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetMargin(0.25);

    int n_entries = static_cast<int>(mc_ch_hists_.size());
    if (mc_total_)
    {
        ++n_entries;
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
        int n_cols = (n_entries > 4) ? 3 : 2;
        if (compact_legend)
        {
            n_cols = 1;
        }
        leg->SetNColumns(n_cols);
    }

    legend_proxies_.clear();

    for (size_t i = 0; i < mc_ch_hists_.size(); ++i)
    {
        int ch = chan_order_.at(i);
        double sum = 0.0;
        if (i < chan_event_yields_.size())
        {
            sum = chan_event_yields_[i];
        }
        else if (mc_ch_hists_[i])
        {
            // Fallback if yields were not cached for some reason.
            // If density_mode_ is true, Integral("width") restores event counts.
            sum = density_mode_ ? mc_ch_hists_[i]->Integral("width") : mc_ch_hists_[i]->Integral();
        }
        std::string label = opt_.particle_level ? ParticleChannels::label(ch) : Channels::label(ch);
        if (label == "#emptyset")
        {
            label = "\xE2\x88\x85";
        }
        if (opt_.annotate_numbers)
        {
            label += " : " + Plotter::fmt_commas(sum, 2);
        }
        auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D *>(
            mc_ch_hists_[i]->Clone((spec_.id + "_leg_ch" + std::to_string(ch)).c_str())));
        proxy->SetDirectory(nullptr);
        proxy->Reset("ICES");

        auto *entry = leg->AddEntry(proxy.get(), label.c_str(), "f");
        leg->SetEntrySeparation(0.01);
        legend_proxies_.push_back(std::move(proxy));
        (void)entry;
    }

    if (mc_total_)
    {
        auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D *>(
            mc_total_->Clone((spec_.id + "_leg_unc").c_str())));
        proxy->SetDirectory(nullptr);
        proxy->Reset("ICES");
        proxy->SetFillColor(kBlack);
        proxy->SetFillStyle(k_uncertainty_fill_style);
        proxy->SetLineColor(kBlack);
        proxy->SetLineStyle(kDashed);
        proxy->SetLineWidth(2);
        leg->AddEntry(proxy.get(), "Stat. #oplus Syst. Unc.", "f");
        legend_proxies_.push_back(std::move(proxy));
    }

    if (sig_hist_)
    {
        std::ostringstream signal_label;
        signal_label << "Signal";
        signal_label << " : "
                     << Plotter::fmt_commas(signal_events_, 2)
                     << " (x" << std::fixed << std::setprecision(2) << signal_scale_ << ")";
        leg->AddEntry(sig_hist_.get(), signal_label.str().c_str(), "l");
    }

    if (has_data())
    {
        leg->AddEntry(data_hist_.get(), "Data", "lep");
    }

    leg->Draw();
}

void StackedHist::draw_cuts(TPad *p, double max_y)
{
    if (!opt_.show_cuts || opt_.cuts.empty())
    {
        return;
    }
    p->cd();
    TH1 *frame = stack_->GetHistogram();
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

void StackedHist::draw(TCanvas &canvas)
{
    build_histograms();
    TPad *p_main = nullptr, *p_ratio = nullptr, *p_legend = nullptr;
    setup_pads(canvas, p_main, p_ratio, p_legend);
    double max_y = 1.;
    draw_stack_and_unc(p_main, max_y);
    draw_cuts(p_main, max_y);
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

void StackedHist::draw_and_save(const std::string &image_format)
{
    std::filesystem::create_directories(output_directory_);
    stack_debug_log("draw_and_save enter: plot='" + plot_name_ +
                    "', out_dir='" + output_directory_ + "'");
    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
    stack_debug_log("canvas constructed: plot='" + plot_name_ + "'");
    draw(canvas);
    stack_debug_log("draw finished: plot='" + plot_name_ + "'");
    const std::string fmt = image_format.empty() ? "png" : image_format;
    const std::string out = output_directory_ + "/" + plot_name_ + "." + fmt;

    stack_debug_log("SaveAs start: plot='" + plot_name_ + "', file='" + out + "'");
    canvas.SaveAs(out.c_str());
    stack_debug_log("SaveAs done: plot='" + plot_name_ + "'");
}

} // namespace nu
