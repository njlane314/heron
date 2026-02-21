/* -- C++ -- */
/**
 *  @file  framework/plot/include/UnstackedHist.hh
 *
 *  @brief Unstacked (overlay) histogram plotting helper.
 */

#ifndef HERON_PLOT_UNSTACKED_H
#define HERON_PLOT_UNSTACKED_H

#include <memory>
#include <string>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>

#include "PlotDescriptors.hh"

class TCanvas;
class TPad;

namespace nu
{

class UnstackedHist
{
  public:
    UnstackedHist(TH1DModel spec,
                  Options opt,
                  std::vector<const Entry *> mc,
                  std::vector<const Entry *> data);

    void draw(TCanvas &canvas);
    void draw_and_save(const std::string &image_format);

  private:
    bool has_data() const noexcept { return static_cast<bool>(data_hist_); }
    bool want_ratio() const noexcept { return opt_.show_ratio && has_data(); }

    void setup_pads(TCanvas &c, TPad *&p_main, TPad *&p_ratio, TPad *&p_legend) const;

    void build_histograms();
    void draw_overlay_and_unc(TPad *p_main, double &max_y);
    void draw_ratio(TPad *p_ratio);
    void draw_legend(TPad *p);
    void draw_cuts(TPad *p, double max_y);
    void draw_watermark(TPad *p, double total_mc) const;

  private:
    TH1DModel spec_;
    Options opt_;
    std::vector<const Entry *> mc_;
    std::vector<const Entry *> data_;

    std::string plot_name_;
    std::string output_directory_;

    std::unique_ptr<THStack> overlay_;
    std::vector<std::unique_ptr<TH1D>> mc_ch_hists_;
    std::vector<int> chan_order_;
    std::vector<double> chan_event_yields_;

    std::unique_ptr<TH1D> mc_total_;
    std::unique_ptr<TH1D> mc_unc_band_;
    std::unique_ptr<TH1D> data_hist_;
    std::unique_ptr<TH1D> sig_hist_;

    std::unique_ptr<TH1D> ratio_hist_;
    std::unique_ptr<TH1D> ratio_band_;

    std::unique_ptr<TLegend> legend_;
    std::vector<std::unique_ptr<TH1D>> legend_proxies_;

    double total_mc_events_ = 0.0;
    double signal_events_ = 0.0;
    double signal_scale_ = 1.0;
    bool density_mode_ = false;
};

} // namespace nu

#endif // HERON_PLOT_UNSTACKED_H
