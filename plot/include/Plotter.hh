/* -- C++ -- */
/**
 *  @file  plot/include/Plotter.hh
 *
 *  @brief Plot orchestration helpers.
 */

#ifndef NUXSEC_PLOT_PLOTTER_H
#define NUXSEC_PLOT_PLOTTER_H

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

    void draw_stack(const TH1DModel &spec, const std::vector<const Entry *> &mc) const;
    void draw_stack(const TH1DModel &spec,
                    const std::vector<const Entry *> &mc,
                    const std::vector<const Entry *> &data) const;
    void draw_stack_cov(const TH1DModel &spec,
                        const std::vector<const Entry *> &mc,
                        const std::vector<const Entry *> &data,
                        const TMatrixDSym &total_cov) const;

    void set_global_style() const;

    static std::string sanitise(const std::string &name);
    static std::string fmt_commas(double value, int precision = -1);

  private:
    Options opt_;
};

} // namespace plot
} // namespace nuxsec

#endif // NUXSEC_PLOT_PLOTTER_H
