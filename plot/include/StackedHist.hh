/* -- C++ -- */
/**
 *  @file  plot/include/StackedHist.hh
 *
 *  @brief Stacked histogram plotting helper that manages histogram grouping,
 *         styling, and render orchestration.
 */

#ifndef HERON_PLOT_STACKED_HIST_H
#define HERON_PLOT_STACKED_HIST_H

#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TImage.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPad.h"

#include "PlotDescriptors.hh"


namespace nu
{

class StackedHist
{
  public:
    StackedHist(TH1DModel spec, Options opt, std::vector<const Entry *> mc, std::vector<const Entry *> data);
    ~StackedHist() = default;

    void draw_and_save(const std::string &image_format);

  protected:
    void draw(TCanvas &canvas);

  private:
    bool has_data() const { return data_hist_ && data_hist_->GetEntries() > 0.0; }
    bool want_ratio() const { return opt_.show_ratio && has_data() && mc_total_; }
    void build_histograms();
    void setup_pads(TCanvas &c, TPad *&p_main, TPad *&p_ratio, TPad *&p_legend) const;
    void draw_stack_and_unc(TPad *p_main, double &max_y);
    void draw_ratio(TPad *p_ratio);
    void draw_legend(TPad *p);
    void draw_cuts(TPad *p, double max_y);

    TH1DModel spec_;
    Options opt_;
    std::vector<const Entry *> mc_;
    std::vector<const Entry *> data_;
    std::string plot_name_;
    std::string output_directory_;
    std::unique_ptr<THStack> stack_;
    std::vector<std::unique_ptr<TH1D>> mc_ch_hists_;
    std::unique_ptr<TH1D> mc_total_;
    std::unique_ptr<TH1D> data_hist_;
    std::unique_ptr<TH1D> sig_hist_;
    std::unique_ptr<TH1D> ratio_hist_;
    std::unique_ptr<TH1D> ratio_band_;
    std::vector<int> chan_order_;
    // When we draw variable-width bins as a density (Scale("width")), TH1::Integral()
    // no longer returns an event count. Keep the pre-density yields around for legend
    // and watermark text.
    std::vector<double> chan_event_yields_; // aligned with mc_ch_hists_/chan_order_
    double total_mc_events_ = 0.0;          // event count (pre-density scaling)
    bool density_mode_ = false;             // true if we applied Scale("width")
    double signal_events_ = 0.0;
    double signal_scale_ = 1.0;
    mutable std::vector<std::unique_ptr<TH1D>> legend_proxies_;
    mutable std::unique_ptr<TLegend> legend_;
};

} // namespace nu


#endif // HERON_PLOT_STACKED_HIST_H
