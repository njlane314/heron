/* -- C++ -- */
/**
 *  @file  plot/include/Plotter.hh
 *
 *  @brief Plot orchestration helpers.
 */

#ifndef Nuxsec_PLOT_PLOTTER_H_INCLUDED
#define Nuxsec_PLOT_PLOTTER_H_INCLUDED

#include <string>
#include <vector>

#include "PlotDescriptors.hh"

namespace nuxsec
{
namespace plot
{

class StackedHist;

class Plotter
{
  public:
    Plotter();
    explicit Plotter(Options opt);

    const Options &options() const noexcept;
    Options &options() noexcept;
    void set_options(Options opt);

    void draw_stack_by_channel(const TH1DModel &spec, const std::vector<const Entry *> &mc) const;
    void draw_stack_by_channel(const TH1DModel &spec,
                               const std::vector<const Entry *> &mc,
                               const std::vector<const Entry *> &data) const;
    void draw_stack_by_channel_with_cov(const TH1DModel &spec,
                                        const std::vector<const Entry *> &mc,
                                        const std::vector<const Entry *> &data,
                                        const TMatrixDSym &total_cov) const;

    static std::string sanitise(const std::string &name);
    static std::string fmt_commas(double value, int precision = -1);

  private:
    void set_global_style() const;

    Options opt_;
};

} // namespace plot
} // namespace nuxsec

#endif // Nuxsec_PLOT_PLOTTER_H_INCLUDED
