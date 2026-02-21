/* -- C++ -- */
/**
 *  @file  framework/plot/src/Plotter.cc
 *
 *  @brief Plot orchestration helpers.
 */

#include "Plotter.hh"

#include <cstdlib>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

#include <TGaxis.h>
#include <TMatrixDSym.h>
#include <TROOT.h>
#include <TStyle.h>

#include "PlotEnv.hh"
#include "StackedHist.hh"
#include "UnstackedHist.hh"

namespace nu
{

namespace
{
bool debug_enabled(const char *env_name)
{
    const char *env = std::getenv(env_name);
    return env != nullptr && std::string(env) != "0";
}

void debug_log(bool enabled, const char *prefix, const std::string &msg)
{
    if (!enabled)
    {
        return;
    }
    std::cout << prefix << msg << "\n";
    std::cout.flush();
}

template <typename PlotType>
void draw_plot(const TH1DModel &spec,
               const Options &opt,
               const std::vector<const Entry *> &mc,
               const std::vector<const Entry *> &data,
               bool debug,
               const char *debug_prefix,
               const char *plot_label,
               const char *constructed_label)
{
    debug_log(debug,
              debug_prefix,
              std::string("draw_") + plot_label + " enter: hist='" + spec.name +
                  "', expr='" + spec.expr +
                  "', mc_entries=" + std::to_string(mc.size()) +
                  ", data_entries=" + std::to_string(data.size()));
    PlotType plot(spec, opt, mc, data);
    debug_log(debug, debug_prefix, std::string(constructed_label) + " constructed for hist='" + spec.name + "'");
    plot.draw_and_save(opt.image_format);
    debug_log(debug,
              debug_prefix,
              std::string("draw_") + plot_label + " exit: hist='" + spec.name + "'");
}

template <typename PlotType>
void draw_plot_cov(const TH1DModel &spec,
                   const Options &opt,
                   const std::vector<const Entry *> &mc,
                   const std::vector<const Entry *> &data,
                   const TMatrixDSym &total_cov)
{
    Plotter plotter(opt);
    plotter.set_global_style();
    auto cov_opt = opt;
    cov_opt.total_cov = std::make_shared<TMatrixDSym>(total_cov);
    PlotType plot(spec, std::move(cov_opt), mc, data);
    plot.draw_and_save(cov_opt.image_format);
}
} // namespace

void apply_env_defaults(Options &opt)
{
    // Keep explicit caller choices; only fill empty fields.
    if (opt.out_dir.empty())
    {
        opt.out_dir = plot_output_dir();
    }
    if (opt.image_format.empty())
    {
        opt.image_format = plot_image_format();
    }
}


Plotter::Plotter()
{
    apply_env_defaults(opt_);
}

Plotter::Plotter(Options opt)
    : opt_(std::move(opt))
{
    apply_env_defaults(opt_);
}

const Options &Plotter::options() const noexcept { return opt_; }

Options &Plotter::options() noexcept { return opt_; }

void Plotter::set_options(Options opt)
{
    opt_ = std::move(opt);
    apply_env_defaults(opt_);
}

void Plotter::draw_stack(const TH1DModel &spec, const std::vector<const Entry *> &mc) const
{
    static const std::vector<const Entry *> empty_data{};
    draw_stack(spec, mc, empty_data);
}

void Plotter::draw_stack(const TH1DModel &spec,
                         const std::vector<const Entry *> &mc,
                         const std::vector<const Entry *> &data) const
{
    set_global_style();
    draw_plot<StackedHist>(spec,
                           opt_,
                           mc,
                           data,
                           debug_enabled("HERON_DEBUG_PLOT_STACK"),
                           "[Plotter][debug] ",
                           "stack",
                           "StackedHist");
}

void Plotter::draw_stack_cov(const TH1DModel &spec,
                             const std::vector<const Entry *> &mc,
                             const std::vector<const Entry *> &data,
                             const TMatrixDSym &total_cov) const
{
    draw_plot_cov<StackedHist>(spec, opt_, mc, data, total_cov);
}

void Plotter::draw_unstack(const TH1DModel &spec, const std::vector<const Entry *> &mc) const
{
    static const std::vector<const Entry *> empty_data{};
    draw_unstack(spec, mc, empty_data);
}

void Plotter::draw_unstack(const TH1DModel &spec,
                           const std::vector<const Entry *> &mc,
                           const std::vector<const Entry *> &data) const
{
    set_global_style();
    draw_plot<UnstackedHist>(spec,
                             opt_,
                             mc,
                             data,
                             debug_enabled("HERON_DEBUG_PLOT_UNSTACK"),
                             "[Plotter][unstack-debug] ",
                             "unstack",
                             "UnstackedHist");
}

void Plotter::draw_unstack_cov(const TH1DModel &spec,
                               const std::vector<const Entry *> &mc,
                               const std::vector<const Entry *> &data,
                               const TMatrixDSym &total_cov) const
{
    draw_plot_cov<UnstackedHist>(spec, opt_, mc, data, total_cov);
}

std::string Plotter::sanitise(const std::string &name)
{
    std::string out;
    out.reserve(name.size());
    for (unsigned char c : name)
    {
        if (std::isalnum(c) || c == '_' || c == '-')
        {
            out.push_back(static_cast<char>(c));
        }
        else
        {
            out.push_back('_');
        }
    }
    if (out.empty())
    {
        return "plot";
    }
    return out;
}

std::string Plotter::fmt_commas(double value, int precision)
{
    std::ostringstream ss;
    if (precision >= 0)
    {
        ss << std::fixed << std::setprecision(precision);
    }
    ss << value;
    std::string text = ss.str();
    const auto pos = text.find('.');
    std::string integer = pos == std::string::npos ? text : text.substr(0, pos);
    std::string fraction = pos == std::string::npos ? std::string{} : text.substr(pos);
    bool negative = false;
    if (!integer.empty() && integer.front() == '-')
    {
        negative = true;
        integer.erase(integer.begin());
    }
    std::string with_commas;
    for (std::size_t i = 0; i < integer.size(); ++i)
    {
        if (i != 0 && (integer.size() - i) % 3 == 0)
        {
            with_commas.push_back(',');
        }
        with_commas.push_back(integer[i]);
    }
    if (negative)
    {
        with_commas.insert(with_commas.begin(), '-');
    }
    return with_commas + fraction;
}

void Plotter::set_global_style() const
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
    // Slightly lighter defaults; ratio pad overrides these explicitly.
    style->SetTitleSize(0.055, "X");
    style->SetTitleSize(0.055, "Y");
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
    style->SetTitleOffset(1.00, "X");
    style->SetTitleOffset(1.05, "Y");
    style->SetOptStat(0);
    style->SetOptTitle(0);
    style->SetPadTickX(1);
    style->SetPadTickY(1);
    // No horizontal x-error bars on binned points (matches reference).
    style->SetErrorX(0.0);
    // Use lighter axis and frame outlines.
    style->SetLineWidth(1);
    style->SetFrameLineWidth(1);
    style->SetHistLineWidth(2);
    // Grid appearance used by the reference code.
    style->SetGridColor(17);
    TGaxis::SetMaxDigits(4);
    style->SetPadLeftMargin(0.15);
    style->SetPadRightMargin(0.05);
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

} // namespace nu
