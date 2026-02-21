/* -- C++ -- */
/**
 *  @file  framework/plot/src/StackedHist.cc
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
#include "TDecompChol.h"
#include "TImage.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixDSym.h"
#include "TPaveText.h"
#include "TVectorD.h"

#include "PlotChannels.hh"
#include "ParticleChannels.hh"
#include "Plotter.hh"

namespace nu
{

namespace
{
constexpr int k_panel_fill_colour = kWhite;
// Default (non-STV) uncertainty appearance.
constexpr int k_uncertainty_fill_colour = kGray + 2;
constexpr int k_uncertainty_fill_style = 3345;
// STV-style uncertainty outline (matches the reference plots: dotted line, no fill).
constexpr int k_uncertainty_outline_style = kDotted;
constexpr int k_uncertainty_line_colour = kGray + 2;
constexpr int k_uncertainty_line_style = k_uncertainty_outline_style;
constexpr int k_uncertainty_line_width = 2;

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

std::pair<int, int> visible_bin_range(const TH1D &h, double xmin, double xmax)
{
    const TAxis *axis = h.GetXaxis();
    int first_bin = 1;
    int last_bin = h.GetNbinsX();
    if (axis && xmin < xmax)
    {
        first_bin = std::max(1, axis->FindFixBin(xmin));
        last_bin = std::min(h.GetNbinsX(), axis->FindFixBin(xmax));
        if (xmax <= axis->GetBinLowEdge(last_bin))
        {
            last_bin = std::max(first_bin, last_bin - 1);
        }
    }
    return {first_bin, last_bin};
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

    auto use_single_sided_ticks = [](TPad *pad) {
        if (!pad)
        {
            return;
        }
        pad->SetTickx(0);
        pad->SetTicky(0);
    };

    c.cd();
    c.SetFillColor(k_panel_fill_colour);
    c.SetFrameFillColor(k_panel_fill_colour);
    c.SetTickx(0);
    c.SetTicky(0);

    p_main = nullptr;
    p_ratio = nullptr;
    p_legend = nullptr;

    if (opt_.stv_style)
    {
        const std::string main_name = plot_name_ + "_pad_main";
        if (want_ratio())
        {
            const std::string ratio_name = plot_name_ + "_pad_ratio";
            p_ratio = new TPad(ratio_name.c_str(), ratio_name.c_str(), 0.0, 0.01, 1.0, 0.23);
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0.0, 0.23, 1.0, 1.0);

            p_ratio->SetTopMargin(0.00);
            p_ratio->SetBottomMargin(0.38);
            p_ratio->SetLeftMargin(0.13);
            p_ratio->SetRightMargin(0.06);
            p_ratio->SetFillColor(k_panel_fill_colour);
            p_ratio->SetFrameFillColor(k_panel_fill_colour);
            p_ratio->SetFrameFillStyle(4000);
            use_single_sided_ticks(p_ratio);
            p_ratio->Draw();
            disable_primitive_ownership(p_ratio);
        }
        else
        {
            p_main = new TPad(main_name.c_str(), main_name.c_str(), 0.0, 0.0, 1.0, 1.0);
        }

        p_main->SetTopMargin(0.02);
        p_main->SetBottomMargin(want_ratio() ? 0.00 : 0.12);
        p_main->SetLeftMargin(0.13);
        p_main->SetRightMargin(0.06);
        p_main->SetFillColor(k_panel_fill_colour);
        p_main->SetFrameFillColor(k_panel_fill_colour);
        use_single_sided_ticks(p_main);
        if (opt_.use_log_y)
        {
            p_main->SetLogy();
        }

        p_main->Draw();
        disable_primitive_ownership(p_main);
        return;
    }

    if (want_ratio())
    {
        const std::string ratio_name = plot_name_ + "_pad_ratio";
        const std::string main_name = plot_name_ + "_pad_main";
        p_ratio = new TPad(ratio_name.c_str(), ratio_name.c_str(), 0., 0.00, 1., 0.30);
        p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.30, 1., 1.00);

        p_main->SetTopMargin(0.06);
        p_main->SetBottomMargin(0.02);
        p_main->SetLeftMargin(0.12);
        p_main->SetRightMargin(0.05);
        p_main->SetFillColor(k_panel_fill_colour);
        p_main->SetFrameFillColor(k_panel_fill_colour);
        use_single_sided_ticks(p_main);

        p_ratio->SetTopMargin(0.04);
        p_ratio->SetBottomMargin(0.35);
        p_ratio->SetLeftMargin(0.12);
        p_ratio->SetRightMargin(0.05);
        p_ratio->SetFillColor(k_panel_fill_colour);
        p_ratio->SetFrameFillColor(k_panel_fill_colour);
        use_single_sided_ticks(p_ratio);
    }
    else
    {
        const std::string main_name = plot_name_ + "_pad_main";
        p_main = new TPad(main_name.c_str(), main_name.c_str(), 0., 0.00, 1., 1.00);
        p_main->SetTopMargin(0.06);
        p_main->SetBottomMargin(0.12);
        p_main->SetLeftMargin(0.12);
        p_main->SetRightMargin(0.05);
        p_main->SetFillColor(k_panel_fill_colour);
        p_main->SetFrameFillColor(k_panel_fill_colour);
        use_single_sided_ticks(p_main);
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
    }
    disable_primitive_ownership(p_ratio);
    disable_primitive_ownership(p_main);
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
    mc_total_stat_.reset();
    data_hist_.reset();
    mc_unc_hist_.reset();
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
    const std::string channel_column = opt_.channel_column.empty() ? "analysis_channels" : opt_.channel_column;

    for (size_t ie = 0; ie < mc_.size(); ++ie)
    {
        const Entry *e = mc_[ie];
        if (!e)
        {
            continue;
        }

        auto n0 = apply(e->rnode(), spec_.sel);
        auto n = (spec_.expr.empty() ? n0 : n0.Define("_nx_expr_", spec_.expr));
        const std::string var = spec_.expr.empty() ? spec_.id : "_nx_expr_";

        if (!particle_level)
        {
            for (int ch : channels)
            {
                auto nf = n.Filter([ch](int c) { return c == ch; }, {channel_column});
                auto h = nf.Histo1D(spec_.model("_mc_ch" + std::to_string(ch) + "_src" + std::to_string(ie)), var, spec_.weight);
                booked[ch].push_back(h);
            }
            continue;
        }

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

            auto h = nsel.Histo1D(spec_.model("_mc_pdg" + std::to_string(ch) + "_src" + std::to_string(ie)),
                                  sel, spec_.weight);
            booked[ch].push_back(h);
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
            hists->SetOwner(kFALSE);
        }
        mc_ch_hists_.push_back(std::move(sum));
        chan_order_.push_back(ch);
        chan_event_yields_.push_back(cy.second);
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
            auto n0 = apply(e->rnode(), spec_.sel);
            auto n = (spec_.expr.empty() ? n0 : n0.Define("_nx_expr_", spec_.expr));
            const std::string var = spec_.expr.empty() ? spec_.id : "_nx_expr_";
            parts.push_back(n.Histo1D(spec_.model("_data_src" + std::to_string(ie)), var));
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
        const double total_sum = mc_total_ ? integral_in_visible_range(*mc_total_, spec_.xmin, spec_.xmax) : 0.0;
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
        if (signal_events_ > 0.0 && total_sum > 0.0)
        {
            signal_scale_ = total_sum / signal_events_;
            sig->Scale(signal_scale_);
        }
        sig->SetLineColor(kGreen + 2);
        sig->SetLineStyle(kDashed);
        sig->SetLineWidth(2);
        sig->SetFillStyle(0);
        sig_hist_ = std::move(sig);
    }

    if (opt_.normalise_by_bin_width)
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

void StackedHist::draw_stack_and_unc(TPad *p_main, double &max_y)
{
    if (!p_main)
    {
        return;
    }
    p_main->cd();
    p_main->SetTickx(0);
    p_main->SetTicky(0);

    if (auto *hists = stack_->GetHists())
    {
        for (TObject *obj = hists->First(); obj != nullptr; obj = hists->After(obj))
        {
            if (auto *hist = dynamic_cast<TH1 *>(obj))
            {
                hist->SetLineColor(hist->GetFillColor());
                hist->SetLineWidth(1);
            }
        }
    }

    if (!opt_.x_title.empty() || !opt_.y_title.empty())
    {
        const std::string default_x = !spec_.name.empty() ? spec_.name : spec_.id;
        const std::string x = opt_.x_title.empty() ? default_x : opt_.x_title;
        const std::string default_y = opt_.particle_level
                                          ? (opt_.stv_style ? "Particles/bin" : (density_mode_ ? "Particles / bin width" : "Particles"))
                                          : (opt_.stv_style ? "Events/bin" : (density_mode_ ? "Events / bin width" : "Events"));
        const std::string y = opt_.y_title.empty() ? default_y : opt_.y_title;
        stack_->SetTitle((";" + x + ";" + y).c_str());
    }

    stack_->Draw("HIST");
    TH1 *frame = stack_->GetHistogram();
    if (frame)
    {
        frame->SetLineWidth(1);
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
        frame->GetXaxis()->SetTickLength(opt_.stv_style ? 0.03 : 0.02);
        // When a ratio pad exists, the x-axis labels/titles live on the ratio pad.
        if (want_ratio())
        {
            frame->GetXaxis()->SetLabelSize(0.0);
            frame->GetXaxis()->SetTitleSize(0.0);
        }
        else
        {
            frame->GetXaxis()->SetLabelSize(opt_.stv_style ? 0.035 : 0.035);
            frame->GetXaxis()->SetTitleSize(opt_.stv_style ? 0.045 : 0.045);
        }
        frame->GetXaxis()->CenterTitle(opt_.stv_style);
        frame->GetYaxis()->SetLabelSize(opt_.stv_style ? 0.035 : 0.035);
        frame->GetYaxis()->SetTitleSize(opt_.stv_style ? 0.045 : 0.045);
        frame->GetYaxis()->SetTitleOffset(opt_.stv_style ? 0.95 : 1.2);
        frame->GetYaxis()->CenterTitle(opt_.stv_style);

        // If we didn't override titles explicitly, THStack's default comes from
        // TH1DModel::axis_title() (which defaults to "Events"). For particle-level
        // plots, update the default Y label to avoid misleading units.
        if (opt_.particle_level && opt_.y_title.empty())
        {
            const std::string cur = frame->GetYaxis() ? frame->GetYaxis()->GetTitle() : "";
            if (cur.empty() || cur == "Events")
            {
                frame->GetYaxis()->SetTitle(density_mode_ ? "Particles / bin width" : "Particles");
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

        const double headroom = opt_.use_log_y ? 10.0 : (opt_.stv_style ? 1.05 : 1.3);
        stack_->SetMaximum(max_y * headroom);
        stack_->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);

        // Uncertainty visual: draw a line-only stat+syst template.
        mc_unc_hist_.reset(static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_mc_unc").c_str())));
        mc_unc_hist_->SetDirectory(nullptr);
        mc_unc_hist_->SetMarkerSize(0);
        mc_unc_hist_->SetFillStyle(0);
        mc_unc_hist_->SetLineColor(k_uncertainty_line_colour);
        mc_unc_hist_->SetLineStyle(k_uncertainty_line_style);
        mc_unc_hist_->SetLineWidth(k_uncertainty_line_width);
        mc_unc_hist_->Draw("E2 SAME");
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
    p_ratio->SetTickx(0);
    p_ratio->SetTicky(0);

    ratio_hist_.reset(static_cast<TH1D *>(
        data_hist_->Clone((spec_.id + "_ratio").c_str())));
    ratio_hist_->SetDirectory(nullptr);
    ratio_hist_->Divide(mc_total_.get());
    ratio_hist_->SetTitle("");
    ratio_hist_->SetLineWidth(opt_.stv_style ? 2 : 1);
    ratio_hist_->SetLineColor(kBlack);
    ratio_hist_->SetMarkerStyle(kFullCircle);
    ratio_hist_->SetMarkerSize(0.8);

    if (opt_.stv_style)
    {
        ratio_hist_->GetXaxis()->SetTitle(opt_.x_title.empty() ? data_hist_->GetXaxis()->GetTitle() : opt_.x_title.c_str());
        ratio_hist_->GetXaxis()->CenterTitle(true);
        ratio_hist_->GetXaxis()->SetLabelSize(0.12);
        ratio_hist_->GetXaxis()->SetTitleSize(0.18);
        ratio_hist_->GetXaxis()->SetTickLength(0.05);
        ratio_hist_->GetXaxis()->SetTitleOffset(0.9);

        ratio_hist_->GetYaxis()->SetTitle("ratio");
        ratio_hist_->GetYaxis()->CenterTitle(true);
        ratio_hist_->GetYaxis()->SetLabelSize(0.08);
        ratio_hist_->GetYaxis()->SetTitleSize(0.15);
        ratio_hist_->GetYaxis()->SetTitleOffset(0.35);
        ratio_hist_->GetYaxis()->SetNdivisions(505);

        double rmax = 1.0;
        double rmin = 1.0;
        bool seeded = false;
        const int nb = ratio_hist_->GetNbinsX();
        for (int i = 1; i <= nb; ++i)
        {
            const double r = ratio_hist_->GetBinContent(i);
            if (!std::isfinite(r) || ratio_hist_->GetBinError(i) <= 0.0)
            {
                continue;
            }
            if (!seeded)
            {
                rmax = rmin = r;
                seeded = true;
            }
            rmax = std::max(rmax, r);
            rmin = std::min(rmin, r);
        }
        ratio_hist_->SetMaximum(rmax + 0.15 * std::abs(rmax));
        ratio_hist_->SetMinimum(rmin - 0.20 * std::abs(rmin));
    }
    else
    {
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
    }

    ratio_hist_->Draw("E1");

    if (opt_.show_ratio_band)
    {
        ratio_band_.reset(static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_ratio_band").c_str())));
        ratio_band_->SetDirectory(nullptr);
        const int nb2 = ratio_band_->GetNbinsX();
        for (int i = 1; i <= nb2; ++i)
        {
            const double m = mc_total_->GetBinContent(i);
            const double em = mc_total_->GetBinError(i);
            ratio_band_->SetBinContent(i, 1.0);
            ratio_band_->SetBinError(i, (m > 0 ? em / m : 0.0));
        }
        ratio_band_->SetFillColor(k_uncertainty_fill_colour);
        ratio_band_->SetFillStyle(k_uncertainty_fill_style);
        ratio_band_->SetLineColor(kGray + 2);
        ratio_band_->SetLineStyle(kSolid);
        ratio_band_->SetLineWidth(1);
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
    unity->SetLineStyle(opt_.stv_style ? 9 : kDashed);
    unity->SetLineWidth(1);
    unity->Draw("SAME");

    ratio_hist_->Draw("E1 SAME");
}

bool StackedHist::compute_chi2(double &chi2_out, int &ndf_out) const
{
    chi2_out = 0.0;
    ndf_out = 0;
    if (!has_data() || !mc_total_ || !data_hist_)
    {
        return false;
    }

    const TH1D &d = *data_hist_;
    const TH1D &m = *mc_total_;
    const TH1D &mstat = (mc_total_stat_ ? *mc_total_stat_ : *mc_total_);

    const auto [first_bin, last_bin] = visible_bin_range(m, spec_.xmin, spec_.xmax);
    std::vector<int> bins;
    bins.reserve(std::max(0, last_bin - first_bin + 1));

    for (int i = first_bin; i <= last_bin; ++i)
    {
        const double v_data = std::pow(d.GetBinError(i), 2);
        const double v_mcstat = std::pow(mstat.GetBinError(i), 2);
        double v_syst = 0.0;
        if (opt_.total_cov && (i - 1) < opt_.total_cov->GetNrows())
        {
            v_syst = (*opt_.total_cov)(i - 1, i - 1);
            if (density_mode_)
            {
                const double w = m.GetXaxis() ? m.GetXaxis()->GetBinWidth(i) : 0.0;
                if (w > 0.0)
                {
                    v_syst /= (w * w);
                }
            }
        }
        else if (!opt_.syst_bin.empty() && (i - 1) < static_cast<int>(opt_.syst_bin.size()))
        {
            double s = std::max(0.0, opt_.syst_bin[i - 1]);
            if (density_mode_)
            {
                const double w = m.GetXaxis() ? m.GetXaxis()->GetBinWidth(i) : 0.0;
                if (w > 0.0)
                {
                    s /= w;
                }
            }
            v_syst = s * s;
        }
        const double v_tot = v_data + v_mcstat + v_syst;
        if (v_tot > 0.0)
        {
            bins.push_back(i);
        }
    }

    if (bins.empty())
    {
        return false;
    }
    ndf_out = static_cast<int>(bins.size());

    if (opt_.total_cov)
    {
        const int n = static_cast<int>(bins.size());
        TMatrixDSym V(n);
        for (int a = 0; a < n; ++a)
        {
            const int i = bins[a];
            const double wi = density_mode_ && m.GetXaxis() ? m.GetXaxis()->GetBinWidth(i) : 1.0;
            for (int b = 0; b < n; ++b)
            {
                const int j = bins[b];
                const double wj = density_mode_ && m.GetXaxis() ? m.GetXaxis()->GetBinWidth(j) : 1.0;

                double v = 0.0;
                if (a == b)
                {
                    v += std::pow(d.GetBinError(i), 2);
                    v += std::pow(mstat.GetBinError(i), 2);
                }

                if ((i - 1) < opt_.total_cov->GetNrows() && (j - 1) < opt_.total_cov->GetNrows())
                {
                    double c = (*opt_.total_cov)(i - 1, j - 1);
                    if (density_mode_ && wi > 0.0 && wj > 0.0)
                    {
                        c /= (wi * wj);
                    }
                    v += c;
                }
                V(a, b) = v;
            }
        }

        TVectorD r(n);
        for (int a = 0; a < n; ++a)
        {
            const int i = bins[a];
            r[a] = d.GetBinContent(i) - m.GetBinContent(i);
        }

        TDecompChol chol(V);
        if (chol.Decompose())
        {
            TVectorD x(r);
            chol.Solve(x);
            chi2_out = r * x;
            return true;
        }
    }

    double chi2 = 0.0;
    for (int i : bins)
    {
        const double diff = d.GetBinContent(i) - m.GetBinContent(i);
        const double v_data = std::pow(d.GetBinError(i), 2);
        const double v_mcstat = std::pow(mstat.GetBinError(i), 2);
        double v_syst = 0.0;
        if (!opt_.syst_bin.empty() && (i - 1) < static_cast<int>(opt_.syst_bin.size()))
        {
            double s = std::max(0.0, opt_.syst_bin[i - 1]);
            if (density_mode_)
            {
                const double w = m.GetXaxis() ? m.GetXaxis()->GetBinWidth(i) : 0.0;
                if (w > 0.0)
                {
                    s /= w;
                }
            }
            v_syst = s * s;
        }
        else if (opt_.total_cov && (i - 1) < opt_.total_cov->GetNrows())
        {
            v_syst = (*opt_.total_cov)(i - 1, i - 1);
            if (density_mode_)
            {
                const double w = m.GetXaxis() ? m.GetXaxis()->GetBinWidth(i) : 0.0;
                if (w > 0.0)
                {
                    v_syst /= (w * w);
                }
            }
        }
        const double v = v_data + v_mcstat + v_syst;
        if (v > 0.0)
        {
            chi2 += (diff * diff) / v;
        }
    }
    chi2_out = chi2;
    return true;
}

void StackedHist::draw_chi2_box(TPad *p_main)
{
    if (!p_main || !opt_.show_chi2)
    {
        return;
    }
    double chi2 = 0.0;
    int ndf = 0;
    if (!compute_chi2(chi2, ndf) || ndf <= 0)
    {
        return;
    }
    p_main->cd();
    chi2_box_ = std::make_unique<TPaveText>(0.16, 0.83, 0.44, 0.91, "NDC");
    chi2_box_->SetFillColor(kWhite);
    chi2_box_->SetFillStyle(1001);
    chi2_box_->SetLineColor(kBlack);
    chi2_box_->SetBorderSize(1);
    chi2_box_->SetTextFont(42);
    chi2_box_->SetTextAlign(12);
    chi2_box_->SetTextSize(0.04);
    std::ostringstream ss;
    ss << "#chi^{2} = " << std::fixed << std::setprecision(2) << chi2 << "/" << ndf << " bins";
    chi2_box_->AddText(ss.str().c_str());
    chi2_box_->Draw("SAME");
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
    draw_chi2_box(p_main);
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
    // Match the reference stacked-plot canvas shape (ROOT default canvas size).
    const int cw = opt_.stv_style ? 700 : 800;
    const int ch = opt_.stv_style ? 500 : 600;
    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), cw, ch);
    stack_debug_log("canvas constructed: plot='" + plot_name_ + "'");
    draw(canvas);
    stack_debug_log("draw finished: plot='" + plot_name_ + "'");
    const std::string fmt = image_format.empty() ? "png" : image_format;
    const std::string out = output_directory_ + "/" + plot_name_ + "." + fmt;

    stack_debug_log("SaveAs start: plot='" + plot_name_ + "', file='" + out + "'");
    canvas.SaveAs(out.c_str());
    stack_debug_log("SaveAs done: plot='" + plot_name_ + "'");

    // The stacked histograms are owned by this StackedHist instance (unique_ptr).
    // Prevent TCanvas from attempting to delete drawn pad primitives during
    // teardown, which can otherwise double-delete histogram objects after SaveAs.
    auto *canvas_primitives = canvas.GetListOfPrimitives();
    if (canvas_primitives)
    {
        canvas_primitives->SetOwner(kFALSE);
    }
}

} // namespace nu
