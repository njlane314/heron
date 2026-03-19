// plot_inf_score0_uncertainty_summary.C
//
// Paper-style summary plot of per-bin fractional uncertainties for the
// histogram of inf_scores[0].
//
// Categories:
//   * Total            : quadrature sum of the systematic categories below
//   * MC Stat          : sqrt(sum w_nominal^2) / sum w_nominal in each bin
//   * X-Sec Total      : GENIE multisim + available dedicated GENIE unisims,
//                        combined in quadrature
//   * Flux Hadronic    : PPFX multisim from weightsPPFX
//   * Flux Beamline    : beamline/flux multisim from weightsFlux (if populated)
//   * Reinteractions   : Geant4 reinteraction multisim from weightsReint
//   * POT / NumTargets / Dirt Norm : optional flat fractional lines
//
// Definitions:
//   For scalar up/down sources:
//      frac(bin) = 0.5 * |N_up(bin) - N_dn(bin)| / N_nom(bin)
//
//   For multisim packed-vector sources:
//      frac(bin) = RMS_u[ N_u(bin) - N_nom(bin) ] / N_nom(bin)
//
// where N_nom(bin) is filled with w_nominal, and N_u(bin) is filled with
// w_nominal * w_u.
//
// Notes:
//   * weightsGenie / weightsPPFX / weightsFlux / weightsReint are packed
//     unsigned short vectors, with unpacked value = packed / 1000.
//   * For your current NuMI production, weightsFlux may be empty because the
//     beamline flux unisims are not flattened into a dedicated vector branch in
//     the makeNuMINtuple path.
//   * The current XML you showed only enables fcl_evtw_00 for the main beam
//     stages, so the active-universe counts for weightsGenie and weightsPPFX
//     may be partial rather than the full configured 500 / 600.
//
// Run with:
//   ./heron macro macros/plot_inf_score0_uncertainty_summary.C
//
//   ./heron macro macros/plot_inf_score0_uncertainty_summary.C \
//     'plot_inf_score0_uncertainty_summary("", "sel_muon", "true", "w_nominal", 12, -15.0, 15.0, false, 0.0, 0.0, 0.0, "inf_score0_uncertainty_summary")'
//
// Signal-only example:
//   ./heron macro macros/plot_inf_score0_uncertainty_summary.C \
//     'plot_inf_score0_uncertainty_summary("", "sel_muon", "is_signal", "w_nominal", 12, -15.0, 15.0, false, 0.0, 0.0, 0.0, "inf_score0_uncertainty_summary_signal")'

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
    std::vector<double> frac_bins;
};

struct MultiSpec
{
    std::string label;
    std::string col;
    int color = kBlack;
    int line_style = 1;
    std::size_t n_total = 0;
    std::size_t n_active = 0;
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

std::vector<double> constant_bins(int nbins, double value)
{
    return std::vector<double>(static_cast<std::size_t>(nbins), value);
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

std::vector<double> compute_mc_stat_bins(const std::vector<double> &nom_bins,
                                         const std::vector<double> &nom_sumw2_bins)
{
    std::vector<double> out(nom_bins.size(), 0.0);
    for (std::size_t i = 0; i < nom_bins.size(); ++i)
    {
        if (nom_bins[i] > 0.0 && nom_sumw2_bins[i] >= 0.0)
            out[i] = std::sqrt(nom_sumw2_bins[i]) / nom_bins[i];
    }
    return out;
}

std::vector<double> compute_scalar_frac_bins(const std::vector<double> &scores,
                                             const std::vector<double> &nom_w,
                                             const std::vector<double> &up_vals,
                                             const std::vector<double> &dn_vals,
                                             const std::vector<double> &nom_bins,
                                             int nbins,
                                             double xmin,
                                             double xmax)
{
    std::vector<double> up_bins(static_cast<std::size_t>(nbins), 0.0);
    std::vector<double> dn_bins(static_cast<std::size_t>(nbins), 0.0);
    std::vector<double> frac_bins(static_cast<std::size_t>(nbins), 0.0);

    for (std::size_t i = 0; i < scores.size(); ++i)
    {
        const int b = score_bin(scores[i], nbins, xmin, xmax);
        if (b < 0)
            continue;

        const std::size_t ib = static_cast<std::size_t>(b);
        const double w = nom_w[i];
        const double up = std::isfinite(up_vals[i]) ? up_vals[i] : 1.0;
        const double dn = std::isfinite(dn_vals[i]) ? dn_vals[i] : 1.0;

        up_bins[ib] += w * up;
        dn_bins[ib] += w * dn;
    }

    for (int b = 0; b < nbins; ++b)
    {
        const std::size_t ib = static_cast<std::size_t>(b);
        if (nom_bins[ib] > 0.0)
            frac_bins[ib] = 0.5 * std::abs(up_bins[ib] - dn_bins[ib]) / nom_bins[ib];
    }

    return frac_bins;
}

std::vector<double> compute_multisim_frac_bins(const std::vector<double> &scores,
                                               const std::vector<double> &nom_w,
                                               const std::vector<UShortVec> &packed_vecs,
                                               const std::vector<double> &nom_bins,
                                               int nbins,
                                               double xmin,
                                               double xmax,
                                               std::size_t &n_total,
                                               std::size_t &n_active)
{
    n_total = 0;
    n_active = 0;

    for (const auto &v : packed_vecs)
        n_total = std::max(n_total, static_cast<std::size_t>(v.size()));

    std::vector<double> frac_bins(static_cast<std::size_t>(nbins), 0.0);
    if (n_total == 0)
        return frac_bins;

    std::vector<std::vector<double>> univ_bins(
        n_total,
        std::vector<double>(static_cast<std::size_t>(nbins), 0.0));
    std::vector<char> active(n_total, 0);

    for (std::size_t i = 0; i < scores.size(); ++i)
    {
        const int b = score_bin(scores[i], nbins, xmin, xmax);
        if (b < 0)
            continue;

        const std::size_t ib = static_cast<std::size_t>(b);
        const double w = nom_w[i];
        const auto &packed = packed_vecs[i];
        const std::size_t nu = std::min(n_total, static_cast<std::size_t>(packed.size()));

        for (std::size_t u = 0; u < nu; ++u)
        {
            if (packed[u] != 1000)
                active[u] = 1;
            univ_bins[u][ib] += w * unpack_packed_weight(packed[u]);
        }
    }

    n_active = static_cast<std::size_t>(std::count(active.begin(), active.end(), 1));
    if (n_active == 0)
        return frac_bins;

    for (int b = 0; b < nbins; ++b)
    {
        const std::size_t ib = static_cast<std::size_t>(b);
        const double nom = nom_bins[ib];
        if (!(nom > 0.0))
            continue;

        double sumsq = 0.0;
        for (std::size_t u = 0; u < n_total; ++u)
        {
            if (!active[u])
                continue;
            const double d = univ_bins[u][ib] - nom;
            sumsq += d * d;
        }

        frac_bins[ib] = std::sqrt(sumsq / static_cast<double>(n_active)) / nom;
    }

    return frac_bins;
}

std::vector<double> quadrature_sum(const std::vector<std::vector<double>> &all_bins)
{
    if (all_bins.empty())
        return {};

    const std::size_t n = all_bins.front().size();
    std::vector<double> out(n, 0.0);

    for (const auto &bins : all_bins)
    {
        if (bins.size() != n)
            continue;
        for (std::size_t i = 0; i < n; ++i)
            out[i] += bins[i] * bins[i];
    }

    for (double &x : out)
        x = std::sqrt(x);

    return out;
}

} // namespace

int plot_inf_score0_uncertainty_summary(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &extra_sel = "true",
    const std::string &nominal_weight = "w_nominal",
    int nbins = 12,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_nominal = false,
    double pot_frac = 0.0,
    double num_targets_frac = 0.0,
    double dirt_frac = 0.0,
    const std::string &output_stem = "inf_score0_uncertainty_summary")
{
    return heron::macro::run_with_guard("plot_inf_score0_uncertainty_summary", [&]() -> int {
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
        std::cout << "[plot_inf_score0_uncertainty_summary] input="
                  << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_inf_score0_uncertainty_summary] input is not an event-list ROOT file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode node = SelectionService::decorate(el.rdf());
        if (!has_column(node, "inf_scores"))
        {
            std::cerr << "[plot_inf_score0_uncertainty_summary] missing inf_scores column.\n";
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
                return std::isfinite(s) && std::isfinite(w) && (s >= xmin) && (s < xmax);
            },
            {"inf_score_0", "__nom_w__"});

        auto n_selected = node.Count();
        auto scores_take = node.Take<double>("inf_score_0");
        auto nom_take = node.Take<double>("__nom_w__");

        std::vector<ScalarSpec> scalar_specs;
        std::vector<DoubleTake> scalar_up_takes;
        std::vector<DoubleTake> scalar_dn_takes;

        auto add_scalar = [&](const std::string &label,
                              const std::string &up_col,
                              const std::string &dn_col,
                              int color,
                              int line_style) {
            if (!(has_column(node, up_col) && has_column(node, dn_col)))
            {
                std::cout << "[plot_inf_score0_uncertainty_summary] skip scalar source '"
                          << label << "' (missing " << up_col << "/" << dn_col << ")\n";
                return;
            }

            ScalarSpec spec;
            spec.label = label;
            spec.up_col = up_col;
            spec.dn_col = dn_col;
            spec.color = color;
            spec.line_style = line_style;
            scalar_specs.push_back(std::move(spec));
            scalar_up_takes.emplace_back(node.Take<double>(up_col));
            scalar_dn_takes.emplace_back(node.Take<double>(dn_col));
        };

        add_scalar("RPA_CCQE", "knobRPAup", "knobRPAdn", kGreen + 2, 2);
        add_scalar("XSecShape_CCMEC", "knobCCMECup", "knobCCMECdn", kBlue + 1, 1);
        add_scalar("AxFFCCQEshape", "knobAxFFCCQEup", "knobAxFFCCQEdn", kOrange + 10, 7);
        add_scalar("VecFFCCQEshape", "knobVecFFCCQEup", "knobVecFFCCQEdn", kCyan + 2, 1);
        add_scalar("DecayAngMEC", "knobDecayAngMECup", "knobDecayAngMECdn", kOrange + 7, 1);
        add_scalar("ThetaDelta2Npi", "knobThetaDelta2Npiup", "knobThetaDelta2Npidn", kRed + 1, 1);
        add_scalar("ThetaDelta2NRad", "knobThetaDelta2NRadup", "knobThetaDelta2NRaddn", kGray + 2, 3);
        add_scalar("NormCCCOH", "knobNormCCCOHup", "knobNormCCCOHdn", kMagenta + 2, 2);
        add_scalar("NormNCCOH", "knobNormNCCOHup", "knobNormNCCOHdn", kViolet + 2, 7);
        add_scalar("xsr_scc_Fv3", "knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn", kAzure + 7, 7);
        add_scalar("xsr_scc_Fa3", "knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn", kPink + 7, 7);

        PackedTake genie_take;
        PackedTake ppfx_take;
        PackedTake flux_take;
        PackedTake reint_take;

        const bool have_genie = has_column(node, "weightsGenie");
        const bool have_ppfx = has_column(node, "weightsPPFX");
        const bool have_flux = has_column(node, "weightsFlux");
        const bool have_reint = has_column(node, "weightsReint");

        if (have_genie)
            genie_take = node.Take<UShortVec>("weightsGenie");
        if (have_ppfx)
            ppfx_take = node.Take<UShortVec>("weightsPPFX");
        if (have_flux)
            flux_take = node.Take<UShortVec>("weightsFlux");
        if (have_reint)
            reint_take = node.Take<UShortVec>("weightsReint");

        const ULong64_t n_evt = *n_selected;
        if (n_evt == 0)
        {
            std::cerr << "[plot_inf_score0_uncertainty_summary] no selected MC-like events.\n";
            return 1;
        }

        const auto &scores = *scores_take;
        const auto &nom_w = *nom_take;

        if (scores.size() != nom_w.size() || scores.size() != static_cast<std::size_t>(n_evt))
        {
            std::cerr << "[plot_inf_score0_uncertainty_summary] internal size mismatch after RDF actions.\n";
            return 1;
        }

        std::vector<double> nom_bins(static_cast<std::size_t>(nbins), 0.0);
        std::vector<double> nom_sumw2_bins(static_cast<std::size_t>(nbins), 0.0);

        for (std::size_t i = 0; i < scores.size(); ++i)
        {
            const int b = score_bin(scores[i], nbins, xmin, xmax);
            if (b < 0)
                continue;
            const std::size_t ib = static_cast<std::size_t>(b);
            const double w = nom_w[i];
            nom_bins[ib] += w;
            nom_sumw2_bins[ib] += w * w;
        }

        const std::vector<double> mc_stat_bins = compute_mc_stat_bins(nom_bins, nom_sumw2_bins);

        std::vector<std::vector<double>> xsec_parts;

        for (std::size_t i = 0; i < scalar_specs.size(); ++i)
        {
            const auto &up = *scalar_up_takes[i];
            const auto &dn = *scalar_dn_takes[i];
            if (up.size() != scores.size() || dn.size() != scores.size())
            {
                std::cerr << "[plot_inf_score0_uncertainty_summary] scalar size mismatch for "
                          << scalar_specs[i].label << "\n";
                return 1;
            }
            scalar_specs[i].frac_bins = compute_scalar_frac_bins(scores,
                                                                 nom_w,
                                                                 up,
                                                                 dn,
                                                                 nom_bins,
                                                                 nbins,
                                                                 xmin,
                                                                 xmax);
            xsec_parts.push_back(scalar_specs[i].frac_bins);
        }

        MultiSpec genie_multisim;
        genie_multisim.label = "GENIE Multisim";
        genie_multisim.col = "weightsGenie";
        genie_multisim.color = kGray + 1;
        genie_multisim.line_style = 3;
        if (have_genie)
        {
            const auto &vecs = *genie_take;
            if (vecs.size() != scores.size())
            {
                std::cerr << "[plot_inf_score0_uncertainty_summary] multisim size mismatch for weightsGenie\n";
                return 1;
            }
            genie_multisim.frac_bins = compute_multisim_frac_bins(scores,
                                                                  nom_w,
                                                                  vecs,
                                                                  nom_bins,
                                                                  nbins,
                                                                  xmin,
                                                                  xmax,
                                                                  genie_multisim.n_total,
                                                                  genie_multisim.n_active);
            xsec_parts.push_back(genie_multisim.frac_bins);
        }
        else
        {
            genie_multisim.frac_bins = constant_bins(nbins, 0.0);
        }

        std::vector<double> xsec_total_bins = quadrature_sum(xsec_parts);
        if (xsec_total_bins.empty())
            xsec_total_bins = constant_bins(nbins, 0.0);

        MultiSpec ppfx_multisim;
        ppfx_multisim.label = "Flux Hadronic";
        ppfx_multisim.col = "weightsPPFX";
        ppfx_multisim.color = kMagenta + 1;
        ppfx_multisim.line_style = 1;
        if (have_ppfx)
        {
            const auto &vecs = *ppfx_take;
            if (vecs.size() != scores.size())
            {
                std::cerr << "[plot_inf_score0_uncertainty_summary] multisim size mismatch for weightsPPFX\n";
                return 1;
            }
            ppfx_multisim.frac_bins = compute_multisim_frac_bins(scores,
                                                                 nom_w,
                                                                 vecs,
                                                                 nom_bins,
                                                                 nbins,
                                                                 xmin,
                                                                 xmax,
                                                                 ppfx_multisim.n_total,
                                                                 ppfx_multisim.n_active);
        }
        else
        {
            ppfx_multisim.frac_bins = constant_bins(nbins, 0.0);
        }

        MultiSpec flux_multisim;
        flux_multisim.label = "Flux Beamline";
        flux_multisim.col = "weightsFlux";
        flux_multisim.color = kCyan + 1;
        flux_multisim.line_style = 1;
        if (have_flux)
        {
            const auto &vecs = *flux_take;
            if (vecs.size() != scores.size())
            {
                std::cerr << "[plot_inf_score0_uncertainty_summary] multisim size mismatch for weightsFlux\n";
                return 1;
            }
            flux_multisim.frac_bins = compute_multisim_frac_bins(scores,
                                                                 nom_w,
                                                                 vecs,
                                                                 nom_bins,
                                                                 nbins,
                                                                 xmin,
                                                                 xmax,
                                                                 flux_multisim.n_total,
                                                                 flux_multisim.n_active);
        }
        else
        {
            flux_multisim.frac_bins = constant_bins(nbins, 0.0);
        }

        MultiSpec reint_multisim;
        reint_multisim.label = "Reinteractions";
        reint_multisim.col = "weightsReint";
        reint_multisim.color = kBlue + 2;
        reint_multisim.line_style = 1;
        if (have_reint)
        {
            const auto &vecs = *reint_take;
            if (vecs.size() != scores.size())
            {
                std::cerr << "[plot_inf_score0_uncertainty_summary] multisim size mismatch for weightsReint\n";
                return 1;
            }
            reint_multisim.frac_bins = compute_multisim_frac_bins(scores,
                                                                  nom_w,
                                                                  vecs,
                                                                  nom_bins,
                                                                  nbins,
                                                                  xmin,
                                                                  xmax,
                                                                  reint_multisim.n_total,
                                                                  reint_multisim.n_active);
        }
        else
        {
            reint_multisim.frac_bins = constant_bins(nbins, 0.0);
        }

        const std::vector<double> pot_bins = constant_bins(nbins, std::max(0.0, pot_frac));
        const std::vector<double> num_targets_bins = constant_bins(nbins, std::max(0.0, num_targets_frac));
        const std::vector<double> dirt_bins = constant_bins(nbins, std::max(0.0, dirt_frac));

        std::vector<std::vector<double>> sys_parts;
        sys_parts.push_back(xsec_total_bins);
        sys_parts.push_back(ppfx_multisim.frac_bins);
        sys_parts.push_back(flux_multisim.frac_bins);
        sys_parts.push_back(reint_multisim.frac_bins);
        if (pot_frac > 0.0)
            sys_parts.push_back(pot_bins);
        if (num_targets_frac > 0.0)
            sys_parts.push_back(num_targets_bins);
        if (dirt_frac > 0.0)
            sys_parts.push_back(dirt_bins);

        std::vector<double> sys_total_bins = quadrature_sum(sys_parts);
        if (sys_total_bins.empty())
            sys_total_bins = constant_bins(nbins, 0.0);

        std::vector<std::vector<double>> total_parts = sys_parts;
        total_parts.push_back(mc_stat_bins);

        std::vector<double> total_bins = quadrature_sum(total_parts);
        if (total_bins.empty())
            total_bins = constant_bins(nbins, 0.0);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[plot_inf_score0_uncertainty_summary] selected events = " << n_evt << "\n";
        std::cout << "  MC Stat        max = " << max_in_bins(mc_stat_bins) << "\n";
        std::cout << "  X-Sec Total    max = " << max_in_bins(xsec_total_bins) << "\n";
        std::cout << "  Sys Total      max = " << max_in_bins(sys_total_bins) << "\n";
        std::cout << "  Flux Hadronic  max = " << max_in_bins(ppfx_multisim.frac_bins)
                  << "  active = " << ppfx_multisim.n_active << "/" << ppfx_multisim.n_total << "\n";
        std::cout << "  Flux Beamline  max = " << max_in_bins(flux_multisim.frac_bins)
                  << "  active = " << flux_multisim.n_active << "/" << flux_multisim.n_total << "\n";
        std::cout << "  Reinteractions max = " << max_in_bins(reint_multisim.frac_bins)
                  << "  active = " << reint_multisim.n_active << "/" << reint_multisim.n_total << "\n";
        std::cout << "  Total          max = " << max_in_bins(total_bins) << "\n";
        if (genie_multisim.n_total > 0)
        {
            std::cout << "  GENIE Multisim max = " << max_in_bins(genie_multisim.frac_bins)
                      << "  active = " << genie_multisim.n_active << "/" << genie_multisim.n_total << "\n";
        }

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        std::vector<std::unique_ptr<TH1D>> curves;
        std::vector<std::pair<TH1D *, std::string>> legend_entries;

        curves.push_back(make_curve_hist("h_total",
                                         total_bins,
                                         nbins,
                                         xmin,
                                         xmax,
                                         kBlack,
                                         1,
                                         3));
        legend_entries.emplace_back(curves.back().get(), "Total");

        curves.push_back(make_curve_hist("h_mc_stat",
                                         mc_stat_bins,
                                         nbins,
                                         xmin,
                                         xmax,
                                         kBlack,
                                         2,
                                         3));
        legend_entries.emplace_back(curves.back().get(), "MC Stat");

        curves.push_back(make_curve_hist("h_sys_total",
                                         sys_total_bins,
                                         nbins,
                                         xmin,
                                         xmax,
                                         kGray + 2,
                                         7,
                                         3));
        legend_entries.emplace_back(curves.back().get(), "Sys Total");

        curves.push_back(make_curve_hist("h_xsec_total",
                                         xsec_total_bins,
                                         nbins,
                                         xmin,
                                         xmax,
                                         kBlue + 1,
                                         1,
                                         3));
        legend_entries.emplace_back(curves.back().get(), "X-Sec Total");

        if (max_in_bins(ppfx_multisim.frac_bins) > 0.0)
        {
            curves.push_back(make_curve_hist("h_flux_hadronic",
                                             ppfx_multisim.frac_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             ppfx_multisim.color,
                                             ppfx_multisim.line_style,
                                             3));
            legend_entries.emplace_back(curves.back().get(), "Flux Hadronic");
        }

        if (max_in_bins(flux_multisim.frac_bins) > 0.0)
        {
            curves.push_back(make_curve_hist("h_flux_beamline",
                                             flux_multisim.frac_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             flux_multisim.color,
                                             flux_multisim.line_style,
                                             3));
            legend_entries.emplace_back(curves.back().get(), "Flux Beamline");
        }

        if (max_in_bins(reint_multisim.frac_bins) > 0.0)
        {
            curves.push_back(make_curve_hist("h_reint",
                                             reint_multisim.frac_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             reint_multisim.color,
                                             reint_multisim.line_style,
                                             3));
            legend_entries.emplace_back(curves.back().get(), "Reinteractions");
        }

        if (pot_frac > 0.0)
        {
            curves.push_back(make_curve_hist("h_pot",
                                             pot_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             kRed + 1,
                                             1,
                                             2));
            legend_entries.emplace_back(curves.back().get(), "POT");
        }

        if (num_targets_frac > 0.0)
        {
            curves.push_back(make_curve_hist("h_num_targets",
                                             num_targets_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             kGreen + 2,
                                             1,
                                             2));
            legend_entries.emplace_back(curves.back().get(), "Num Targets");
        }

        if (dirt_frac > 0.0)
        {
            curves.push_back(make_curve_hist("h_dirt_norm",
                                             dirt_bins,
                                             nbins,
                                             xmin,
                                             xmax,
                                             kOrange + 7,
                                             1,
                                             2));
            legend_entries.emplace_back(curves.back().get(), "Dirt Norm");
        }

        double ymax = 0.0;
        for (const auto &curve : curves)
            ymax = std::max(ymax, curve->GetMaximum());

        if (!(ymax > 0.0))
        {
            std::cerr << "[plot_inf_score0_uncertainty_summary] all curves are zero.\n";
            return 1;
        }

        TCanvas c("c_inf_score0_uncertainty_summary",
                  "Inference score [0] uncertainty summary",
                  1150,
                  800);
        c.SetLeftMargin(0.11);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.12);

        TH1D hframe("hframe_inf_score0_uncertainty_summary",
                    ";Inference score [0];Fractional uncertainty",
                    100,
                    xmin,
                    xmax);
        hframe.SetMinimum(0.0);
        hframe.SetMaximum(1.15 * ymax);
        hframe.Draw("AXIS");

        for (const auto &curve : curves)
            curve->Draw("HIST SAME");

        TLegend leg(0.42, 0.59, 0.89, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        if (legend_entries.size() > 6)
            leg.SetNColumns(2);

        for (const auto &entry : legend_entries)
            leg.AddEntry(entry.first, entry.second.c_str(), "l");
        leg.Draw();

        c.RedrawAxis();

        const auto out = plot_output_file(output_stem).string();
        c.SaveAs(out.c_str());

        std::cout << "[plot_inf_score0_uncertainty_summary] saved: " << out << "\n";
        return 0;
    });
}

int plot_inference_score_uncertainty_summary(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &extra_sel = "true",
    const std::string &nominal_weight = "w_nominal",
    int nbins = 12,
    double xmin = -15.0,
    double xmax = 15.0,
    bool include_ext_in_nominal = false,
    double pot_frac = 0.0,
    double num_targets_frac = 0.0,
    double dirt_frac = 0.0,
    const std::string &output_stem = "inf_score0_uncertainty_summary")
{
    return plot_inf_score0_uncertainty_summary(event_list_path,
                                               base_sel,
                                               extra_sel,
                                               nominal_weight,
                                               nbins,
                                               xmin,
                                               xmax,
                                               include_ext_in_nominal,
                                               pot_frac,
                                               num_targets_frac,
                                               dirt_frac,
                                               output_stem);
}
