/* -- C++ -- */
/**
 *  @file  plot/src/StackedHist.cc
 *
 *  @brief Stacked histogram plotting helper.
 */

#include "StackedHist.hh"

#include <algorithm>
#include <cctype>
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
#include "TArrow.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TLatex.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixDSym.h"

#include "PlotChannels.hh"
#include "Plotter.hh"

namespace
{

constexpr double kEdgeEps = 1e-12;

void fold_overflow_into_edges(TH1D &h)
{
    const int nb = h.GetNbinsX();

    // Underflow -> first bin
    {
        const double c0 = h.GetBinContent(0);
        const double e0 = h.GetBinError(0);
        const double c1 = h.GetBinContent(1);
        const double e1 = h.GetBinError(1);
        h.SetBinContent(1, c1 + c0);
        h.SetBinError(1, std::hypot(e1, e0));
        h.SetBinContent(0, 0.0);
        h.SetBinError(0, 0.0);
    }

    // Overflow -> last bin
    {
        const double co = h.GetBinContent(nb + 1);
        const double eo = h.GetBinError(nb + 1);
        const double cn = h.GetBinContent(nb);
        const double en = h.GetBinError(nb);
        h.SetBinContent(nb, cn + co);
        h.SetBinError(nb, std::hypot(en, eo));
        h.SetBinContent(nb + 1, 0.0);
        h.SetBinError(nb + 1, 0.0);
    }
}

std::vector<double> make_minstat_edges(const TH1D &h,
                                       double wmin,
                                       double relerrMax,
                                       bool keep_edge_bins)
{
    std::vector<double> edges;
    const int nb = h.GetNbinsX();
    const auto *ax = h.GetXaxis();

    if (nb <= 0)
    {
        edges.push_back(ax->GetXmin());
        edges.push_back(ax->GetXmax());
        return edges;
    }

    edges.reserve(nb + 1);
    edges.push_back(ax->GetBinLowEdge(1));

    int first_merge_bin = 1;
    int last_merge_bin = nb;

    if (keep_edge_bins && nb >= 3)
    {
        edges.push_back(ax->GetBinUpEdge(1));
        first_merge_bin = 2;
        last_merge_bin = nb - 1;
    }

    double sumw = 0.0;
    double sumw2 = 0.0;

    const bool use_wmin = (wmin > 0.0);
    const bool use_relerr = (relerrMax > 0.0);

    for (int i = first_merge_bin; i <= last_merge_bin; ++i)
    {
        const double w = h.GetBinContent(i);
        const double e = h.GetBinError(i);
        sumw += w;
        sumw2 += e * e;

        const double denom = std::abs(sumw);

        const bool pass_w = (!use_wmin) ? true : (denom >= wmin);
        const bool pass_re = (!use_relerr) ? true
                                           : (denom > 0.0 && std::sqrt(sumw2) / denom <= relerrMax);

        if (pass_w && pass_re)
        {
            const double up = ax->GetBinUpEdge(i);
            if (edges.empty() || up > edges.back() + kEdgeEps)
            {
                edges.push_back(up);
            }
            sumw = 0.0;
            sumw2 = 0.0;
        }
    }

    if (keep_edge_bins && nb >= 3)
    {
        const double last_low = ax->GetBinLowEdge(nb);
        if (edges.empty() || last_low > edges.back() + kEdgeEps)
        {
            edges.push_back(last_low);
        }
    }

    // Ensure we end exactly at xmax
    const double xmax = ax->GetBinUpEdge(nb);
    if (edges.empty() || xmax > edges.back() + kEdgeEps)
    {
        edges.push_back(xmax);
    }

    // De-dup in case of numerical weirdness
    edges.erase(std::unique(edges.begin(), edges.end(), [](double a, double b) { return std::abs(a - b) <= kEdgeEps; }),
                edges.end());

    if (edges.size() < 2)
    {
        edges.clear();
        edges.push_back(ax->GetXmin());
        edges.push_back(ax->GetXmax());
    }

    return edges;
}

std::unique_ptr<TH1D> rebin_to_edges(const TH1D &h, const std::vector<double> &edges, const std::string &new_name)
{
    // If no edges provided, just clone.
    if (edges.size() < 2)
    {
        auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone(new_name.c_str())));
        out->SetDirectory(nullptr);
        return out;
    }

    // Rebin is non-const, so operate on a detached clone.
    auto tmp = std::unique_ptr<TH1D>(static_cast<TH1D *>(h.Clone((new_name + "_tmp").c_str())));
    tmp->SetDirectory(nullptr);
    if (tmp->GetSumw2N() == 0)
    {
        tmp->Sumw2(kTRUE);
    }

    const int nnew = static_cast<int>(edges.size()) - 1;
    auto *reb = tmp->Rebin(nnew, new_name.c_str(), edges.data());
    auto out = std::unique_ptr<TH1D>(static_cast<TH1D *>(reb));
    out->SetDirectory(nullptr);
    return out;
}

void log_adaptive_bin_widths(const std::string &hist_id,
                             const std::vector<double> &edges)
{
    if (edges.size() < 2)
    {
        return;
    }

    std::ostringstream msg;
    msg << "[StackedHist] Adaptive bins settled for '" << hist_id << "': "
        << (edges.size() - 1) << " bins; widths [";

    for (std::size_t i = 1; i < edges.size(); ++i)
    {
        if (i > 1)
        {
            msg << ", ";
        }
        msg << (edges[i] - edges[i - 1]);
    }

    msg << "]";
    std::clog << msg.str() << '\n';
}

} // namespace

namespace nu
{

void apply_total_errors(TH1D &h, const TMatrixDSym *cov, const std::vector<double> *syst_bin)
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
            p_main->Draw();
            disable_primitive_ownership(p_main);
        }
    }
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
    signal_events_ = 0.0;
    signal_scale_ = 1.0;
    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> booked;
    const auto &channels = Channels::mc_keys();

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
        for (int ch : channels)
        {
            auto nf = n.Filter([ch](int c) { return c == ch; }, {"analysis_channels"});
            auto h = nf.Histo1D(spec_.model("_mc_ch" + std::to_string(ch) + "_src" + std::to_string(ie)), var, spec_.weight);
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

    // ---- Adaptive min-stat-per-bin rebin (derived from TOTAL MC) ----
    std::vector<double> adaptive_edges;
    if (opt_.adaptive_binning && !sum_by_channel.empty())
    {
        std::unique_ptr<TH1D> mc_tot;
        for (auto &kv : sum_by_channel)
        {
            if (!kv.second)
            {
                continue;
            }

            if (!mc_tot)
            {
                mc_tot.reset(static_cast<TH1D *>(kv.second->Clone((spec_.id + "_mc_total_fine").c_str())));
                mc_tot->SetDirectory(nullptr);
            }
            else
            {
                mc_tot->Add(kv.second.get());
            }
        }

        if (mc_tot)
        {
            if (opt_.adaptive_fold_overflow)
            {
                fold_overflow_into_edges(*mc_tot);
            }

            adaptive_edges = make_minstat_edges(*mc_tot,
                                                opt_.adaptive_min_sumw,
                                                opt_.adaptive_max_relerr,
                                                opt_.adaptive_keep_edge_bins);

            log_adaptive_bin_widths(spec_.id, adaptive_edges);

            if (adaptive_edges.size() >= 2)
            {
                for (auto &kv : sum_by_channel)
                {
                    if (!kv.second)
                    {
                        continue;
                    }

                    if (opt_.adaptive_fold_overflow)
                    {
                        fold_overflow_into_edges(*kv.second);
                    }

                    kv.second = rebin_to_edges(*kv.second,
                                               adaptive_edges,
                                               spec_.id + "_mc_sum_ch" + std::to_string(kv.first) + "_rebin");
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
    for (auto &cy : yields)
    {
        const int ch = cy.first;
        auto it = sum_by_channel.find(ch);
        if (it == sum_by_channel.end() || !it->second)
        {
            continue;
        }
        auto &sum = it->second;
        sum->SetFillColor(Channels::colour(ch));
        sum->SetFillStyle(Channels::fill_style(ch));
        sum->SetLineColor(kBlack);
        sum->SetLineWidth(1);
        stack_->Add(sum.get(), "HIST");
        mc_ch_hists_.push_back(std::move(sum));
        chan_order_.push_back(ch);
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
        if (data_hist_ && !adaptive_edges.empty())
        {
            if (opt_.adaptive_fold_overflow)
            {
                fold_overflow_into_edges(*data_hist_);
            }

            data_hist_ = rebin_to_edges(*data_hist_, adaptive_edges, (spec_.id + "_data_rebin"));
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
        const std::string y = opt_.y_title.empty() ? "Events" : opt_.y_title;
        stack_->SetTitle((";" + x + ";" + y).c_str());
    }

    stack_->Draw("HIST");
    TH1 *frame = stack_->GetHistogram();
    if (frame)
    {
        frame->SetLineWidth(2);
    }
    if (frame && spec_.xmin < spec_.xmax)
    {
        frame->GetXaxis()->SetRangeUser(spec_.xmin, spec_.xmax);
    }
    if (frame)
    {
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
    if (mc_total_)
    {
        if (opt_.total_cov || !opt_.syst_bin.empty())
        {
            apply_total_errors(*mc_total_, opt_.total_cov.get(), opt_.syst_bin.empty() ? nullptr : &opt_.syst_bin);
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
        h->SetFillStyle(3004);
        h->SetMarkerSize(0);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
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

    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D *>(
        data_hist_->Clone((spec_.id + "_ratio").c_str())));
    ratio->SetDirectory(nullptr);
    ratio->Divide(mc_total_.get());
    ratio->SetTitle("; ;Data / MC");
    ratio->SetMaximum(1.99);
    ratio->SetMinimum(0.01);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleSize(0.10);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetYaxis()->SetTitle("Data / MC");
    ratio->GetXaxis()->CenterTitle(false);
    ratio->GetXaxis()->SetTitle(opt_.x_title.empty() ? data_hist_->GetXaxis()->GetTitle() : opt_.x_title.c_str());

    ratio->Draw("E1");

    std::unique_ptr<TH1D> band;
    if (opt_.show_ratio_band)
    {
        band.reset(static_cast<TH1D *>(mc_total_->Clone((spec_.id + "_ratio_band").c_str())));
        band->SetDirectory(nullptr);
        const int nb = band->GetNbinsX();
        for (int i = 1; i <= nb; ++i)
        {
            const double m = mc_total_->GetBinContent(i);
            const double em = mc_total_->GetBinError(i);
            band->SetBinContent(i, 1.0);
            band->SetBinError(i, (m > 0 ? em / m : 0.0));
        }
        band->SetFillColor(kBlack);
        band->SetFillStyle(3004);
        band->SetLineColor(kBlack);
        band->SetMarkerSize(0);
        band->Draw("E2 SAME");
    }

    ratio->Draw("E1 SAME");
}

void StackedHist::draw_legend(TPad *p)
{
    if (!p)
    {
        return;
    }
    p->cd();
    legend_ = std::make_unique<TLegend>(0.12, 0.0, 0.95, 0.75);
    auto *leg = legend_.get();
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

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
        const int n_cols = (n_entries > 4) ? 3 : 2;
        leg->SetNColumns(n_cols);
    }

    legend_proxies_.clear();

    for (size_t i = 0; i < mc_ch_hists_.size(); ++i)
    {
        int ch = chan_order_.at(i);
        double sum = mc_ch_hists_[i]->Integral();
        std::string label = Channels::label(ch);
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
        proxy->SetFillStyle(3004);
        proxy->SetLineColor(kBlack);
        proxy->SetLineWidth(1);
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

void StackedHist::draw_watermark(TPad *p, double total_mc) const
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
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });

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

void StackedHist::draw(TCanvas &canvas)
{
    build_histograms();
    TPad *p_main = nullptr, *p_ratio = nullptr, *p_legend = nullptr;
    setup_pads(canvas, p_main, p_ratio, p_legend);
    double max_y = 1.;
    draw_stack_and_unc(p_main, max_y);
    draw_cuts(p_main, max_y);
    draw_watermark(p_main, mc_total_ ? mc_total_->Integral() : 0.0);
    draw_legend(p_legend ? p_legend : p_main);
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
    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
    draw(canvas);
    const std::string fmt = image_format.empty() ? "png" : image_format;
    const std::string out = output_directory_ + "/" + plot_name_ + "." + fmt;

    canvas.SaveAs(out.c_str());
}

} // namespace nu
