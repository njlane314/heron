/* -- C++ -- */
/**
 *  @file  plot/include/StackedHist.hh
 *
 *  @brief Stacked histogram plotting helper.
 */

#ifndef NUXSEC_PLOT_STACKED_HIST_H
#define NUXSEC_PLOT_STACKED_HIST_H

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

namespace nuxsec
{
namespace plot
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
    bool want_ratio() const { return opt_.show_ratio && data_hist_ && mc_total_; }
    void build_histograms();
    void setup_pads(TCanvas &c, TPad *&p_main, TPad *&p_ratio, TPad *&p_legend) const;
    void draw_stack_and_unc(TPad *p_main, double &max_y);
    void draw_ratio(TPad *p_ratio);
    void draw_legend(TPad *p);
    void draw_cuts(TPad *p, double max_y);
    void draw_watermark(TPad *p, double total_mc) const;

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
    std::vector<int> chan_order_;
    double signal_scale_ = 1.0;
    mutable std::unique_ptr<TLegend> legend_;
    mutable std::vector<std::unique_ptr<TH1D>> legend_proxies_;
};

} // namespace plot
} // namespace nuxsec

#endif // NUXSEC_PLOT_STACKED_HIST_H
