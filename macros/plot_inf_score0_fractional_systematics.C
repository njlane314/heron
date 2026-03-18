// plot_inf_score0_fractional_systematics.C
//
// Fractional systematic-uncertainty curves for the selected distribution of
// inf_scores[0].
//
// Definitions:
//   Nominal:
//      N_nom(bin) = Sum_events w_nominal
//
//   Scalar +/-1 sigma sources:
//      f(bin) = 0.5 * |N_up(bin) - N_dn(bin)| / N_nom(bin)
//
//   Multisim sources:
//      f(bin) = RMS_u[ N_u(bin) - N_nom(bin) ] / N_nom(bin)
//   where the RMS is taken only over "active" universes, i.e. universes that
//   differ from packed unity (1000) for at least one selected event.
//
//   GENIE total:
//      quadrature sum of GENIE multisim + all available scalar GENIE sources.
//
// Notes:
//   * weightsGenie / weightsPPFX are packed unsigned short weights, with
//     unpacked value = packed / 1000.
//   * With the current production XML, only fcl_evtw_00 is enabled, so the
//     multisim branches will likely report ~100 active universes rather than
//     the full 500 (GENIE) / 600 (PPFX).
//   * This macro plots fractional uncertainty on the predicted bin yield, not
//     shape-only uncertainty. If you want shape-only uncertainty, normalise
//     each varied histogram to its total before forming the ratios.
//
// Run with:
//   ./heron macro macros/plot_inf_score0_fractional_systematics.C
//
//   ./heron macro macros/plot_inf_score0_fractional_systematics.C \
//     'plot_inf_score0_fractional_systematics("", "sel_muon", "true", "w_nominal", 40, -15.0, 15.0, false, false, false, false, "inf_score0_frac_syst")'
//
// Signal-only example:
//   ./heron macro macros/plot_inf_score0_fractional_systematics.C \
//     'plot_inf_score0_fractional_systematics("", "sel_muon", "is_signal", "w_nominal", 40, -15.0, 15.0, false, false, false, false, "inf_score0_frac_syst_signal")'
//
// Include flux / reinteraction / PPFX lines:
//   ./heron macro macros/plot_inf_score0_fractional_systematics.C \
//     'plot_inf_score0_fractional_systematics("", "sel_muon", "true", "w_nominal", 40, -15.0, 15.0, false, true, true, true, "inf_score0_frac_syst_all")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

using UShortVec = ROOT::VecOps::RVec<unsigned short>;
using DoubleTake = ROOT::RDF::RResultPtr<std::vector<double>>;
using PackedTake = ROOT::RDF::RResultPtr<std::vector<UShortVec>>;

struct ScalarSpec
{
    std::string label;
    std::string up_col;
    std::string dn_col;
    int color = kBlack;
    int line_style = 1;
    bool show = true;
    std::vector<double> up_bins;
    std::vector<double> dn_bins;
    std::vector<double> frac_bins;
};

struct MultiSpec
{
    std::string label;
    std::string col;
    int color = kBlack;
    int line_style = 1;
    bool show = true;
    std::size_t n_total = 0;
    std::size_t n_active = 0;
    std::vector<std::vector<double>> univ_bins;
    std::vector<char> active;
    std::vector<double> frac_bins;
};

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

ROOT::RDF::RNode apply_optional_filter(ROOT::RDF::RNode node,
                                       const std::string &expr)
{
    if (expr.empty() || expr == "true")
        return node;

    if (has_column(node, expr))
        return node.Filter([](bool pass) { return pass; }, {expr});

    return node.Filter(expr);
}

int score_bin(double x, int nbins, double xmin, double xmax)
{
    if (!(x >= xmin) || !(x < xmax))
        return -1;

    const double span = xmax - xmin;
    if (!(span > 0.0))
        return -1;

    const double z = (x - xmin) / span;
    int b = static_cast<int>(std::floor(z * static_cast<double>(nbins)));

    if (b < 0)
        return -1;
    if (b >= nbins)
        b = nbins - 1;
    return b;
}

double unpack_packed_weight(unsigned short packed)
{
    return 0.001 * static_cast<double>(packed);
}

double max_in_bins(const std::vector<double> &bins)
{
    double out = 0.0;
    for (double x : bins)
    {
        if (std::isfinite(x))
            out = std::max(out, x);
    }
    return out;
}

void apply_curve_style(TH1D &h, int color, int line_style, int line_width = 3)
{
    h.SetLineColor(color);
    h.SetLineStyle(line_style);
    h.SetLineWidth(line_width);
    h.SetFillStyle(0);
    h.SetDirectory(nullptr);
}

std::unique_ptr<TH1D> make_curve_hist(const std::string &name,
                                      const std::vector<double> &bins,
                                      int nbins,
                                      double xmin,
                                      double xmax,
                                      int color,
                                      int line_style,
                                      int line_width = 3)
{
    auto h = std::make_unique<TH1D>(name.c_str(), "", nbins, xmin, xmax);
    for (int b = 1; b <= nbins && static_cast<std::size_t>(b) <= bins.size(); ++b)
    {
        h->SetBinContent(b, bins[static_cast<std::size_t>(b - 1)]);
    }
    apply_curve_style(*h, color, line_style, line_width);
    return h;
}

} // namespace

int plot_inf_score0_fractional_systematics(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &extra_sel = "true",
    const std::string &nominal_weight = "w_nominal",
    int nbins = 40,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_nominal = false,
    bool include_flux = false,
    bool include_reint = false,
    bool include_ppfx = false,
    const std::string &output_stem = "inf_score0_fractional_systematics")
{
    return heron::macro::run_with_guard("plot_inf_score0_fractional_systematics", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        if (implicit_mt_enabled())
            ROOT::EnableImplicitMT();

        TH1::SetDefaultSumw2();

        const std::string input_path =
            event_list_path.empty() ? default_event_list_root() : event_list_path;

        std::cout << "[plot_inf_score0_fractional_systematics] input="
                  << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inf_score0_fractional_systematics] input is not an event-list ROOT file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode node = SelectionService::decorate(el.rdf());
        if (!has_column(node, "inf_scores"))
        {
            std::cerr << "[plot_inf_score0_fractional_systematics] missing inf_scores column.\n";
            return 1;
        }

        node = node.Define(
            "inf_score_0",
            [](const ROOT::RVec<float> &scores) {
                return scores.empty()
                           ? std::numeric_limits<double>::quiet_NaN()
                           : static_cast<double>(scores[0]);
            },
            {"inf_scores"});

        node = node.Define("__nom_w__", nominal_weight);

        auto mask_mc_like = el.mask_for_mc_like();
        auto mask_ext = el.mask_for_ext();

        node = filter_by_sample_mask(node, mask_mc_like, "sample_id");
        if (!include_ext_in_nominal)
            node = filter_not_sample_mask(node, mask_ext, "sample_id");

        node = apply_optional_filter(node, base_sel);
        node = apply_optional_filter(node, extra_sel);

        node = node.Filter(
            [xmin, xmax](double s, double w) {
                return std::isfinite(s) && std::isfinite(w) &&
                       (s >= xmin) && (s < xmax);
            },
            {"inf_score_0", "__nom_w__"});

        auto n_selected = node.Count();
        auto scores_take = node.Take<double>("inf_score_0");
        auto nom_take = node.Take<double>("__nom_w__");

        std::vector<ScalarSpec> scalar_specs;
        std::vector<DoubleTake> scalar_up_takes;
        std::vector<DoubleTake> scalar_dn_takes;

        std::vector<MultiSpec> multi_specs;
        std::vector<PackedTake> multi_takes;

        auto add_scalar = [&](const std::string &label,
                              const std::string &up_col,
                              const std::string &dn_col,
                              int color,
                              int line_style,
                              bool show) {
            if (!(has_column(node, up_col) && has_column(node, dn_col)))
            {
                std::cout << "[plot_inf_score0_fractional_systematics] skip scalar source '"
                          << label << "' (missing " << up_col << "/" << dn_col << ")\n";
                return;
            }

            ScalarSpec spec;
            spec.label = label;
            spec.up_col = up_col;
            spec.dn_col = dn_col;
            spec.color = color;
            spec.line_style = line_style;
            spec.show = show;

            scalar_specs.push_back(std::move(spec));
            scalar_up_takes.emplace_back(node.Take<double>(up_col));
            scalar_dn_takes.emplace_back(node.Take<double>(dn_col));
        };

        auto add_multi = [&](const std::string &label,
                             const std::string &col,
                             int color,
                             int line_style,
                             bool show) -> int {
            if (!has_column(node, col))
            {
                std::cout << "[plot_inf_score0_fractional_systematics] skip multisim source '"
                          << label << "' (missing " << col << ")\n";
                return -1;
            }

            MultiSpec spec;
            spec.label = label;
            spec.col = col;
            spec.color = color;
            spec.line_style = line_style;
            spec.show = show;

            const int idx = static_cast<int>(multi_specs.size());
            multi_specs.push_back(std::move(spec));
            multi_takes.emplace_back(node.Take<UShortVec>(col));
            return idx;
        };

        // GENIE model sources
        add_scalar("NormCCCOH", "knobNormCCCOHup", "knobNormCCCOHdn", kBlue + 1, 2, true);
        add_scalar("RPA_CCQE", "knobRPAup", "knobRPAdn", kGreen + 2, 2, true);
        add_scalar("VecFFCCQEshape", "knobVecFFCCQEup", "knobVecFFCCQEdn", kOrange + 7, 1, true);
        add_scalar("DecayAngMEC", "knobDecayAngMECup", "knobDecayAngMECdn", kOrange + 1, 1, true);
        add_scalar("ThetaDelta2NRad", "knobThetaDelta2NRadup", "knobThetaDelta2NRaddn", kGray + 2, 3, true);
        add_scalar("XSecShape_CCMEC", "knobCCMECup", "knobCCMECdn", kMagenta + 2, 2, true);
        add_scalar("ThetaDelta2Npi", "knobThetaDelta2Npiup", "knobThetaDelta2Npidn", kRed + 1, 1, true);

        // Extra available branches, computed but not drawn by default.
        add_scalar("AxFFCCQEshape", "knobAxFFCCQEup", "knobAxFFCCQEdn", kOrange + 10, 7, false);
        add_scalar("NormNCCOH", "knobNormNCCOHup", "knobNormNCCOHdn", kBlue + 2, 7, false);
        add_scalar("xsr_scc_Fv3", "knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn", kCyan + 2, 7, false);
        add_scalar("xsr_scc_Fa3", "knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn", kViolet + 2, 7, false);
        add_scalar("RootinoFix", "RootinoFix", "RootinoFix", kPink + 6, 7, false);

        const int genie_multi_idx = add_multi("Multisim", "weightsGenie", kGray + 1, 3, true);

        int flux_multi_idx = -1;
        int reint_multi_idx = -1;
        int ppfx_multi_idx = -1;

        if (include_flux)
            flux_multi_idx = add_multi("Flux multisim", "weightsFlux", kAzure + 2, 7, true);
        if (include_reint)
            reint_multi_idx = add_multi("Reint multisim", "weightsReint", kTeal + 2, 7, true);
        if (include_ppfx)
            ppfx_multi_idx = add_multi("PPFX multisim", "weightsPPFX", kViolet + 1, 7, true);

        const ULong64_t n_evt = *n_selected;
        if (n_evt == 0)
        {
            std::cerr << "[plot_inf_score0_fractional_systematics] no selected events.\n";
            return 1;
        }

        const auto &scores = *scores_take;
        const auto &nom_w = *nom_take;

        if (scores.size() != nom_w.size() || scores.size() != static_cast<std::size_t>(n_evt))
        {
            std::cerr << "[plot_inf_score0_fractional_systematics] internal size mismatch after RDF actions.\n";
            return 1;
        }

        std::vector<double> nom_bins(static_cast<std::size_t>(nbins), 0.0);

        for (auto &spec : scalar_specs)
        {
            spec.up_bins.assign(static_cast<std::size_t>(nbins), 0.0);
            spec.dn_bins.assign(static_cast<std::size_t>(nbins), 0.0);
            spec.frac_bins.assign(static_cast<std::size_t>(nbins), 0.0);
        }

        for (std::size_t i = 0; i < multi_specs.size(); ++i)
        {
            const auto &per_event_vecs = *multi_takes[i];
            if (per_event_vecs.size() != scores.size())
            {
                std::cerr << "[plot_inf_score0_fractional_systematics] multisim column size mismatch for "
                          << multi_specs[i].label << ".\n";
                return 1;
            }

            std::size_t n_univ = 0;
            for (const auto &v : per_event_vecs)
                n_univ = std::max(n_univ, static_cast<std::size_t>(v.size()));

            multi_specs[i].n_total = n_univ;
            multi_specs[i].n_active = 0;
            multi_specs[i].univ_bins.assign(n_univ, std::vector<double>(static_cast<std::size_t>(nbins), 0.0));
            multi_specs[i].active.assign(n_univ, 0);
            multi_specs[i].frac_bins.assign(static_cast<std::size_t>(nbins), 0.0);
        }

        for (std::size_t i = 0; i < scores.size(); ++i)
        {
            const int b = score_bin(scores[i], nbins, xmin, xmax);
            if (b < 0)
                continue;

            const std::size_t ib = static_cast<std::size_t>(b);
            const double w = nom_w[i];
            nom_bins[ib] += w;

            for (std::size_t k = 0; k < scalar_specs.size(); ++k)
            {
                const auto &up_vals = *scalar_up_takes[k];
                const auto &dn_vals = *scalar_dn_takes[k];

                const double up_w = std::isfinite(up_vals[i]) ? up_vals[i] : 1.0;
                const double dn_w = std::isfinite(dn_vals[i]) ? dn_vals[i] : 1.0;

                scalar_specs[k].up_bins[ib] += w * up_w;
                scalar_specs[k].dn_bins[ib] += w * dn_w;
            }

            for (std::size_t m = 0; m < multi_specs.size(); ++m)
            {
                const auto &per_event_vecs = *multi_takes[m];
                const auto &packed = per_event_vecs[i];
                const std::size_t nu = std::min(static_cast<std::size_t>(packed.size()), multi_specs[m].n_total);

                for (std::size_t u = 0; u < nu; ++u)
                {
                    if (packed[u] != 1000)
                        multi_specs[m].active[u] = 1;

                    multi_specs[m].univ_bins[u][ib] += w * unpack_packed_weight(packed[u]);
                }
            }
        }

        for (auto &spec : scalar_specs)
        {
            for (int b = 0; b < nbins; ++b)
            {
                const std::size_t ib = static_cast<std::size_t>(b);
                const double nom = nom_bins[ib];
                if (nom > 0.0)
                {
                    spec.frac_bins[ib] =
                        0.5 * std::abs(spec.up_bins[ib] - spec.dn_bins[ib]) / nom;
                }
            }
        }

        for (auto &spec : multi_specs)
        {
            spec.n_active = static_cast<std::size_t>(
                std::count(spec.active.begin(), spec.active.end(), 1));

            if (spec.n_active == 0)
                continue;

            for (int b = 0; b < nbins; ++b)
            {
                const std::size_t ib = static_cast<std::size_t>(b);
                const double nom = nom_bins[ib];
                if (!(nom > 0.0))
                    continue;

                double sumsq = 0.0;
                for (std::size_t u = 0; u < spec.n_total; ++u)
                {
                    if (!spec.active[u])
                        continue;

                    const double d = spec.univ_bins[u][ib] - nom;
                    sumsq += d * d;
                }

                spec.frac_bins[ib] = std::sqrt(sumsq / static_cast<double>(spec.n_active)) / nom;
            }
        }

        std::vector<double> genie_total_bins(static_cast<std::size_t>(nbins), 0.0);
        bool have_genie_total = false;

        for (int b = 0; b < nbins; ++b)
        {
            const std::size_t ib = static_cast<std::size_t>(b);
            double sum2 = 0.0;

            if (genie_multi_idx >= 0)
            {
                const double x = multi_specs[static_cast<std::size_t>(genie_multi_idx)].frac_bins[ib];
                sum2 += x * x;
            }

            for (const auto &spec : scalar_specs)
            {
                const double x = spec.frac_bins[ib];
                sum2 += x * x;
            }

            genie_total_bins[ib] = std::sqrt(sum2);
            if (genie_total_bins[ib] > 0.0)
                have_genie_total = true;
        }

        std::cout << "[plot_inf_score0_fractional_systematics] selected events = "
                  << n_evt << "\n";
        std::cout << "[plot_inf_score0_fractional_systematics] nominal expression = "
                  << nominal_weight << "\n";

        for (const auto &spec : scalar_specs)
        {
            std::cout << "  scalar: " << std::setw(18) << std::left << spec.label
                      << " max_frac=" << std::setw(10) << std::right << max_in_bins(spec.frac_bins)
                      << (spec.show ? "" : "  (computed, not drawn by default)")
                      << "\n";
        }

        for (const auto &spec : multi_specs)
        {
            std::cout << "  multisim: " << std::setw(18) << std::left << spec.label
                      << " active=" << spec.n_active << "/" << spec.n_total
                      << "  max_frac=" << std::setw(10) << std::right << max_in_bins(spec.frac_bins)
                      << "\n";
        }

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        std::vector<std::unique_ptr<TH1D>> curves;
        std::vector<std::pair<TH1D *, std::string>> legend_entries;

        if (have_genie_total)
        {
            curves.push_back(make_curve_hist("h_genie_total",
                                             genie_total_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             kBlack,
                                             1,
                                             3));
            legend_entries.emplace_back(curves.back().get(), "GENIE total");
        }

        if (genie_multi_idx >= 0 &&
            max_in_bins(multi_specs[static_cast<std::size_t>(genie_multi_idx)].frac_bins) > 0.0)
        {
            const auto &spec = multi_specs[static_cast<std::size_t>(genie_multi_idx)];
            curves.push_back(make_curve_hist("h_genie_multisim",
                                             spec.frac_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             spec.color,
                                             spec.line_style,
                                             3));
            legend_entries.emplace_back(curves.back().get(), spec.label);
        }

        for (std::size_t i = 0; i < scalar_specs.size(); ++i)
        {
            const auto &spec = scalar_specs[i];
            if (!spec.show || max_in_bins(spec.frac_bins) <= 0.0)
                continue;

            std::ostringstream name;
            name << "h_scalar_" << i;
            curves.push_back(make_curve_hist(name.str(),
                                             spec.frac_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             spec.color,
                                             spec.line_style,
                                             3));
            legend_entries.emplace_back(curves.back().get(), spec.label);
        }

        for (std::size_t i = 0; i < multi_specs.size(); ++i)
        {
            if (static_cast<int>(i) == genie_multi_idx)
                continue;

            const auto &spec = multi_specs[i];
            if (!spec.show || max_in_bins(spec.frac_bins) <= 0.0)
                continue;

            std::ostringstream name;
            name << "h_multi_" << i;
            curves.push_back(make_curve_hist(name.str(),
                                             spec.frac_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             spec.color,
                                             spec.line_style,
                                             3));
            legend_entries.emplace_back(curves.back().get(), spec.label);
        }

        double ymax = 0.0;
        for (const auto &curve : curves)
            ymax = std::max(ymax, curve->GetMaximum());

        if (!(ymax > 0.0))
        {
            std::cerr << "[plot_inf_score0_fractional_systematics] all uncertainty curves are zero.\n";
            return 1;
        }

        TCanvas c("c_inf_score0_fractional_systematics",
                  "inf_scores[0] fractional systematics",
                  1150,
                  800);
        c.SetLeftMargin(0.11);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.12);

        TH1D hframe("hframe_inf_score0_frac_syst",
                    ";Inference score [0];Fractional uncertainty",
                    100,
                    xmin,
                    xmax);
        hframe.SetMinimum(0.0);
        hframe.SetMaximum(1.15 * ymax);
        hframe.Draw("AXIS");

        for (const auto &curve : curves)
            curve->Draw("HIST SAME");

        TLegend leg(0.43, 0.60, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        if (legend_entries.size() > 5)
            leg.SetNColumns(2);

        for (const auto &entry : legend_entries)
            leg.AddEntry(entry.first, entry.second.c_str(), "l");
        leg.Draw();

        c.RedrawAxis();

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());

        std::cout << "[plot_inf_score0_fractional_systematics] saved: "
                  << out << "\n";

        return 0;
    });
}

int plot_inference_score_fractional_systematics(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &extra_sel = "true",
    const std::string &nominal_weight = "w_nominal",
    int nbins = 40,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_nominal = false,
    bool include_flux = false,
    bool include_reint = false,
    bool include_ppfx = false,
    const std::string &output_stem = "inf_score0_fractional_systematics")
{
    return plot_inf_score0_fractional_systematics(event_list_path,
                                                  base_sel,
                                                  extra_sel,
                                                  nominal_weight,
                                                  nbins,
                                                  xmin,
                                                  xmax,
                                                  include_ext_in_nominal,
                                                  include_flux,
                                                  include_reint,
                                                  include_ppfx,
                                                  output_stem);
}
