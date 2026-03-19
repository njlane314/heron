// plot_score_systematics.C
//
// Diagnostic plots for score-axis systematics built from event-weight branches.
//
// What it makes:
//   1) Nominal total prediction + data on inf_scores[0], with an absolute
//      systematic band from the selected weight variations.
//   2) Fractional systematic band sigma_i / N_i^0 on the total prediction.
//   3) The same fractional band for signal only.
//   4) The same fractional band for a user-chosen dominant background.
//   5) Cumulative tail-shift plots
//         R_k(t) = N(score > t)^k / N(score > t)^0
//      for total, signal, and dominant background.
//
// Default nuisance basis:
//   - named scalar GENIE up/down knobs
//   - flux multiverse (weightsPPFX if present, else weightsFlux)
//   - reinteraction multiverse (weightsReint)
//   - full GENIE multiverse (weightsGenie) is OFF by default to avoid
//     double counting with the named GENIE knobs
//
// Important assumptions:
//   - scalar knob branches are relative multipliers around the nominal weight
//   - multiverse branches store unsigned-short weights with scale factor 1000
//   - flux universes are treated as multiplicative relative variations on top
//     of w_nominal by default. If instead your PPFX universe weights replace
//     ppfx_cv, set flux_universes_replace_cv = true.
//
// Example:
//   ./heron macro plot_score_systematics.C
//
//   ./heron macro plot_score_systematics.C \
//     'plot_score_systematics("./scratch/out/event_list_myana.root",
//                             "sel_muon",
//                             "is_signal",
//                             "is_pi0_like",
//                             40, -15, 15,
//                             true,   // use named GENIE knobs
//                             false,  // use full GENIE multiverse
//                             true,   // use flux multiverse
//                             true,   // use reint multiverse
//                             false,  // flux universes replace ppfx_cv?
//                             6,      // max named tail curves to overlay
//                             "score_systematics")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
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

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

double sanitise_multiplier(double w)
{
    if (!std::isfinite(w) || w <= 0.0)
        return 1.0;
    return w;
}

double decode_u16_weight(const ROOT::RVec<unsigned short> &v, std::size_t i)
{
    if (i >= v.size())
        return 1.0;
    return sanitise_multiplier(0.001 * static_cast<double>(v[i]));
}

int score_bin(double x, int nbins, double xmin, double xmax)
{
    if (!(x >= xmin) || !(x < xmax) || !(xmax > xmin) || nbins <= 0)
        return -1;

    const double bw = (xmax - xmin) / static_cast<double>(nbins);
    const int b = static_cast<int>((x - xmin) / bw);
    if (b < 0 || b >= nbins)
        return -1;
    return b;
}

void fill_hist(std::vector<double> &h, int b, double w)
{
    if (b < 0 || b >= static_cast<int>(h.size()))
        return;
    h[static_cast<std::size_t>(b)] += w;
}

std::vector<double> cumulative_from_high(const std::vector<double> &h)
{
    std::vector<double> out(h.size(), 0.0);
    double run = 0.0;
    for (int b = static_cast<int>(h.size()) - 1; b >= 0; --b)
    {
        run += h[static_cast<std::size_t>(b)];
        out[static_cast<std::size_t>(b)] = run;
    }
    return out;
}

std::vector<double> add_hists(const std::vector<double> &a, const std::vector<double> &b)
{
    const std::size_t n = std::max(a.size(), b.size());
    std::vector<double> out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double va = (i < a.size()) ? a[i] : 0.0;
        const double vb = (i < b.size()) ? b[i] : 0.0;
        out[i] = va + vb;
    }
    return out;
}

double max_content(const std::vector<double> &y, const std::vector<double> *e = nullptr)
{
    double m = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        const double ey = (e && i < e->size()) ? (*e)[i] : 0.0;
        m = std::max(m, y[i] + ey);
    }
    return m;
}

void add_pair_variance(const std::vector<double> &up,
                       const std::vector<double> &dn,
                       std::vector<double> &diag2)
{
    const std::size_t n = std::min({up.size(), dn.size(), diag2.size()});
    for (std::size_t i = 0; i < n; ++i)
    {
        const double delta = 0.5 * (up[i] - dn[i]);
        diag2[i] += delta * delta;
    }
}

void add_multiverse_variance(const std::vector<std::vector<double>> &hu,
                             std::vector<double> &diag2)
{
    const std::size_t nu = hu.size();
    if (nu < 2 || diag2.empty())
        return;

    const std::size_t nbins = diag2.size();
    for (std::size_t b = 0; b < nbins; ++b)
    {
        double mean = 0.0;
        for (const auto &u : hu)
            mean += (b < u.size()) ? u[b] : 0.0;
        mean /= static_cast<double>(nu);

        double var = 0.0;
        for (const auto &u : hu)
        {
            const double x = (b < u.size()) ? u[b] : 0.0;
            const double d = x - mean;
            var += d * d;
        }
        var /= static_cast<double>(nu - 1);
        diag2[b] += var;
    }
}

std::vector<double> fractional_band(const std::vector<double> &nom,
                                    const std::vector<double> &diag2)
{
    const std::size_t n = std::min(nom.size(), diag2.size());
    std::vector<double> out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (nom[i] > 0.0 && diag2[i] >= 0.0)
            out[i] = std::sqrt(diag2[i]) / nom[i];
    }
    return out;
}

std::vector<double> absolute_sigma(const std::vector<double> &diag2)
{
    std::vector<double> out(diag2.size(), 0.0);
    for (std::size_t i = 0; i < diag2.size(); ++i)
        out[i] = (diag2[i] >= 0.0) ? std::sqrt(diag2[i]) : 0.0;
    return out;
}

std::unique_ptr<TH1D> make_hist(const std::string &name,
                                const std::string &title,
                                int nbins,
                                double xmin,
                                double xmax,
                                const std::vector<double> &y,
                                const std::vector<double> *e = nullptr)
{
    auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), nbins, xmin, xmax);
    h->SetDirectory(nullptr);
    for (int b = 1; b <= nbins; ++b)
    {
        const std::size_t i = static_cast<std::size_t>(b - 1);
        const double yi = (i < y.size()) ? y[i] : 0.0;
        const double ei = (e && i < e->size()) ? (*e)[i] : 0.0;
        h->SetBinContent(b, yi);
        h->SetBinError(b, ei);
    }
    return h;
}

struct NamedKnob
{
    std::string label;
    std::string up_col;
    std::string dn_col;

    std::vector<double> up_total;
    std::vector<double> dn_total;
    std::vector<double> up_sig;
    std::vector<double> dn_sig;
    std::vector<double> up_dom;
    std::vector<double> dn_dom;
};

struct TailCurve
{
    std::string label;
    std::vector<double> ratio;
    double impact = 0.0;
};

struct SampleBundle
{
    std::vector<double> mc_total;
    std::vector<double> ext_total;
    std::vector<double> signal;
    std::vector<double> dominant_bkg;
    std::vector<double> data;

    std::vector<double> diag2_total;
    std::vector<double> diag2_signal;
    std::vector<double> diag2_dombkg;

    std::vector<double> tail_diag2_total;
    std::vector<double> tail_diag2_signal;
    std::vector<double> tail_diag2_dombkg;
};

template <typename T, typename = void>
struct has_mask_for_data_like : std::false_type
{
};

template <typename T>
struct has_mask_for_data_like<T, std::void_t<decltype(std::declval<T &>().mask_for_data_like())>> : std::true_type
{
};

template <typename T, typename = void>
struct has_mask_for_data : std::false_type
{
};

template <typename T>
struct has_mask_for_data<T, std::void_t<decltype(std::declval<T &>().mask_for_data())>> : std::true_type
{
};

template <typename T, typename = void>
struct has_mask_for_ext_like : std::false_type
{
};

template <typename T>
struct has_mask_for_ext_like<T, std::void_t<decltype(std::declval<T &>().mask_for_ext_like())>> : std::true_type
{
};

template <typename T, typename = void>
struct has_mask_for_ext : std::false_type
{
};

template <typename T>
struct has_mask_for_ext<T, std::void_t<decltype(std::declval<T &>().mask_for_ext())>> : std::true_type
{
};

ROOT::RDF::RNode empty_node(ROOT::RDF::RNode node)
{
    return node.Filter([]() { return false; });
}

template <typename EL>
ROOT::RDF::RNode get_data_like_node(EL &el, ROOT::RDF::RNode node, bool &found)
{
    if constexpr (has_mask_for_data_like<EL>::value)
    {
        found = true;
        return filter_by_sample_mask(node, el.mask_for_data_like(), "sample_id");
    }
    else if constexpr (has_mask_for_data<EL>::value)
    {
        found = true;
        return filter_by_sample_mask(node, el.mask_for_data(), "sample_id");
    }
    else
    {
        found = false;
        return empty_node(node);
    }
}

template <typename EL>
ROOT::RDF::RNode get_ext_like_node(EL &el, ROOT::RDF::RNode node, bool &found)
{
    if constexpr (has_mask_for_ext_like<EL>::value)
    {
        found = true;
        return filter_by_sample_mask(node, el.mask_for_ext_like(), "sample_id");
    }
    else if constexpr (has_mask_for_ext<EL>::value)
    {
        found = true;
        return filter_by_sample_mask(node, el.mask_for_ext(), "sample_id");
    }
    else
    {
        found = false;
        return empty_node(node);
    }
}

std::vector<NamedKnob> make_named_knobs(int nbins,
                                        bool use_genie_named)
{
    std::vector<NamedKnob> knobs;
    if (!use_genie_named)
        return knobs;

    auto add = [&](const std::string &label,
                   const std::string &up_col,
                   const std::string &dn_col) {
        NamedKnob k;
        k.label = label;
        k.up_col = up_col;
        k.dn_col = dn_col;
        k.up_total.assign(static_cast<std::size_t>(nbins), 0.0);
        k.dn_total.assign(static_cast<std::size_t>(nbins), 0.0);
        k.up_sig.assign(static_cast<std::size_t>(nbins), 0.0);
        k.dn_sig.assign(static_cast<std::size_t>(nbins), 0.0);
        k.up_dom.assign(static_cast<std::size_t>(nbins), 0.0);
        k.dn_dom.assign(static_cast<std::size_t>(nbins), 0.0);
        knobs.push_back(std::move(k));
    };

    add("RPA_CCQE", "knobRPAup", "knobRPAdn");
    add("XSecShape_CCMEC", "knobCCMECup", "knobCCMECdn");
    add("AxFFCCQEshape", "knobAxFFCCQEup", "knobAxFFCCQEdn");
    add("VecFFCCQEshape", "knobVecFFCCQEup", "knobVecFFCCQEdn");
    add("DecayAngMEC", "knobDecayAngMECup", "knobDecayAngMECdn");
    add("Theta_Delta2Npi", "knobThetaDelta2Npiup", "knobThetaDelta2Npidn");
    add("ThetaDelta2NRad", "knobThetaDelta2NRadup", "knobThetaDelta2NRaddn");
    add("NormCCCOH", "knobNormCCCOHup", "knobNormCCCOHdn");
    add("NormNCCOH", "knobNormNCCOHup", "knobNormNCCOHdn");
    add("xsr_scc_Fv3", "knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn");
    add("xsr_scc_Fa3", "knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn");
    return knobs;
}

void ensure_scalar_column(ROOT::RDF::RNode &node, const std::string &name)
{
    if (!has_column(node, name))
        node = node.Define(name, [] { return 1.0; });
}

void ensure_vector_u16_column(ROOT::RDF::RNode &node, const std::string &name)
{
    if (!has_column(node, name))
        node = node.Define(name, [] { return ROOT::RVec<unsigned short>{}; });
}

int max_size_of_vector_column(ROOT::RDF::RNode node, const std::string &name)
{
    if (!has_column(node, name))
        return 0;

    auto r = node.Define("__tmp_n__", [](const ROOT::RVec<unsigned short> &v) {
                    return static_cast<int>(v.size());
                },
                {name})
                 .Max("__tmp_n__");
    return static_cast<int>(*r);
}

std::vector<TailCurve> top_tail_curves(const std::vector<NamedKnob> &knobs,
                                       const std::vector<double> &tail_nom_target,
                                       const std::vector<double> &tail_ext_total,
                                       int which_target, // 0 total, 1 signal, 2 dombkg
                                       std::size_t keep_n)
{
    std::vector<TailCurve> curves;
    for (const auto &k : knobs)
    {
        const auto &up_hist = (which_target == 0) ? k.up_total : (which_target == 1) ? k.up_sig
                                                                                      : k.up_dom;
        const auto &dn_hist = (which_target == 0) ? k.dn_total : (which_target == 1) ? k.dn_sig
                                                                                      : k.dn_dom;

        std::vector<double> up_tail = cumulative_from_high(up_hist);
        std::vector<double> dn_tail = cumulative_from_high(dn_hist);

        if (which_target == 0)
        {
            for (std::size_t i = 0; i < up_tail.size() && i < tail_ext_total.size(); ++i)
            {
                up_tail[i] += tail_ext_total[i];
                dn_tail[i] += tail_ext_total[i];
            }
        }

        TailCurve cu;
        cu.label = k.label + " up";
        cu.ratio.assign(tail_nom_target.size(), 1.0);

        TailCurve cd;
        cd.label = k.label + " dn";
        cd.ratio.assign(tail_nom_target.size(), 1.0);

        for (std::size_t i = 0; i < tail_nom_target.size(); ++i)
        {
            const double nom = tail_nom_target[i];
            if (nom > 0.0)
            {
                cu.ratio[i] = (i < up_tail.size()) ? up_tail[i] / nom : 1.0;
                cd.ratio[i] = (i < dn_tail.size()) ? dn_tail[i] / nom : 1.0;
                cu.impact = std::max(cu.impact, std::abs(cu.ratio[i] - 1.0));
                cd.impact = std::max(cd.impact, std::abs(cd.ratio[i] - 1.0));
            }
        }

        curves.push_back(std::move(cu));
        curves.push_back(std::move(cd));
    }

    std::sort(curves.begin(), curves.end(), [](const TailCurve &a, const TailCurve &b) {
        return a.impact > b.impact;
    });

    if (curves.size() > keep_n)
        curves.resize(keep_n);

    return curves;
}

void style_band_hist(TH1D &h)
{
    h.SetFillColorAlpha(kAzure - 9, 0.45);
    h.SetLineColor(kAzure - 3);
    h.SetLineWidth(1);
    h.SetMarkerSize(0.0);
}

void style_nominal_hist(TH1D &h)
{
    h.SetLineColor(kBlack);
    h.SetLineWidth(3);
    h.SetFillStyle(0);
}

void style_signal_hist(TH1D &h)
{
    h.SetLineColor(kRed + 1);
    h.SetLineStyle(2);
    h.SetLineWidth(3);
    h.SetFillStyle(0);
}

void style_dom_hist(TH1D &h)
{
    h.SetLineColor(kBlue + 1);
    h.SetLineStyle(3);
    h.SetLineWidth(3);
    h.SetFillStyle(0);
}

void style_data_hist(TH1D &h)
{
    h.SetLineColor(kBlack);
    h.SetMarkerColor(kBlack);
    h.SetMarkerStyle(20);
    h.SetMarkerSize(1.0);
    h.SetLineWidth(2);
}

void style_ratio_curve(TH1D &h, int idx)
{
    static const int cols[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2};
    static const int lsty[] = {1, 2, 3, 4, 5, 6, 7, 9};
    const int ncol = static_cast<int>(sizeof(cols) / sizeof(cols[0]));
    const int nls = static_cast<int>(sizeof(lsty) / sizeof(lsty[0]));
    h.SetLineColor(cols[idx % ncol]);
    h.SetLineStyle(lsty[idx % nls]);
    h.SetLineWidth(2);
    h.SetFillStyle(0);
}

double y_max_for_ratio_band(const std::vector<double> &frac, double floor = 0.05)
{
    double m = floor;
    for (double x : frac)
        m = std::max(m, x);
    return 1.35 * m;
}

void fill_poisson_errors(std::vector<double> &e, const std::vector<double> &n)
{
    e.assign(n.size(), 0.0);
    for (std::size_t i = 0; i < n.size(); ++i)
        e[i] = (n[i] >= 0.0) ? std::sqrt(n[i]) : 0.0;
}

} // namespace

int plot_score_systematics(const std::string &event_list_path = "",
                           const std::string &base_sel = "sel_muon",
                           const std::string &signal_sel = "is_signal",
                           const std::string &dominant_bkg_sel = "",
                           int nbins = 40,
                           double xmin = -15.0,
                           double xmax = 15.0,
                           bool use_genie_named = true,
                           bool use_genie_multiverse = false,
                           bool use_flux_multiverse = true,
                           bool use_reint_multiverse = true,
                           bool flux_universes_replace_cv = false,
                           int max_tail_curves = 6,
                           const std::string &output_stem = "score_systematics")
{
    return heron::macro::run_with_guard("plot_score_systematics", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        ROOT::DisableImplicitMT();
        TH1::SetDefaultSumw2();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[plot_score_systematics] input = " << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[plot_score_systematics] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);

        ROOT::RDF::RNode node = SelectionService::decorate(el.rdf())
                                    .Define("inf_score_0",
                                            [](const ROOT::RVec<float> &scores) {
                                                return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                            },
                                            {"inf_scores"});

        if (!base_sel.empty())
            node = node.Filter(base_sel);

        node = node.Filter(
            [xmin, xmax](double s) {
                return std::isfinite(s) && (s >= xmin) && (s < xmax);
            },
            {"inf_score_0"});

        if (!has_column(node, "w_nominal"))
        {
            std::cerr << "[plot_score_systematics] missing required column: w_nominal\n";
            return 1;
        }

        node = node.Define("__is_signal__", signal_sel);
        if (!dominant_bkg_sel.empty())
            node = node.Define("__is_dom__", dominant_bkg_sel);
        else
            node = node.Define("__is_dom__", "!__is_signal__");

        ensure_scalar_column(node, "ppfx_cv");
        node = node.Define(
            "__w_nominal_nofluxcv__",
            [](double w_nom, auto ppfx_cv) {
                const double cv = sanitise_multiplier(ppfx_cv);
                return (cv > 0.0) ? (w_nom / cv) : w_nom;
            },
            {"w_nominal", "ppfx_cv"});

        ROOT::RDF::RNode node_mc = filter_by_sample_mask(node, el.mask_for_mc_like(), "sample_id");

        bool have_data_node = false;
        bool have_ext_node = false;
        ROOT::RDF::RNode node_data = get_data_like_node(el, node, have_data_node);
        ROOT::RDF::RNode node_ext = get_ext_like_node(el, node, have_ext_node);

        const bool has_ppfx = has_column(node_mc, "weightsPPFX");
        const bool has_flux = has_column(node_mc, "weightsFlux");
        const bool has_reint = has_column(node_mc, "weightsReint");
        const bool has_genie = has_column(node_mc, "weightsGenie");

        std::string flux_branch = "";
        if (use_flux_multiverse)
        {
            if (has_ppfx)
                flux_branch = "weightsPPFX";
            else if (has_flux)
                flux_branch = "weightsFlux";
        }

        if (!flux_branch.empty())
            node_mc = node_mc.Define("__weights_flux__", [=](const ROOT::RVec<unsigned short> &v) { return v; }, {flux_branch});
        else
            node_mc = node_mc.Define("__weights_flux__", [] { return ROOT::RVec<unsigned short>{}; });

        if (use_reint_multiverse && has_reint)
            node_mc = node_mc.Define("__weights_reint__", [](const ROOT::RVec<unsigned short> &v) { return v; }, {"weightsReint"});
        else
            node_mc = node_mc.Define("__weights_reint__", [] { return ROOT::RVec<unsigned short>{}; });

        if (use_genie_multiverse && has_genie)
            node_mc = node_mc.Define("__weights_genie__", [](const ROOT::RVec<unsigned short> &v) { return v; }, {"weightsGenie"});
        else
            node_mc = node_mc.Define("__weights_genie__", [] { return ROOT::RVec<unsigned short>{}; });

        std::vector<NamedKnob> knobs = make_named_knobs(nbins, use_genie_named);

        // These columns appear in the Foreach signature below, so make sure they
        // exist even when the named-GENIE basis is disabled.
        ensure_scalar_column(node_mc, "knobRPAup");
        ensure_scalar_column(node_mc, "knobRPAdn");
        ensure_scalar_column(node_mc, "knobCCMECup");
        ensure_scalar_column(node_mc, "knobCCMECdn");
        ensure_scalar_column(node_mc, "knobAxFFCCQEup");
        ensure_scalar_column(node_mc, "knobAxFFCCQEdn");
        ensure_scalar_column(node_mc, "knobVecFFCCQEup");
        ensure_scalar_column(node_mc, "knobVecFFCCQEdn");
        ensure_scalar_column(node_mc, "knobDecayAngMECup");
        ensure_scalar_column(node_mc, "knobDecayAngMECdn");
        ensure_scalar_column(node_mc, "knobThetaDelta2Npiup");
        ensure_scalar_column(node_mc, "knobThetaDelta2Npidn");
        ensure_scalar_column(node_mc, "knobThetaDelta2NRadup");
        ensure_scalar_column(node_mc, "knobThetaDelta2NRaddn");
        ensure_scalar_column(node_mc, "knobNormCCCOHup");
        ensure_scalar_column(node_mc, "knobNormCCCOHdn");
        ensure_scalar_column(node_mc, "knobNormNCCOHup");
        ensure_scalar_column(node_mc, "knobNormNCCOHdn");
        ensure_scalar_column(node_mc, "knobxsr_scc_Fv3up");
        ensure_scalar_column(node_mc, "knobxsr_scc_Fv3dn");
        ensure_scalar_column(node_mc, "knobxsr_scc_Fa3up");
        ensure_scalar_column(node_mc, "knobxsr_scc_Fa3dn");

        const int n_flux_u = max_size_of_vector_column(node_mc, "__weights_flux__");
        const int n_reint_u = max_size_of_vector_column(node_mc, "__weights_reint__");
        const int n_genie_u = max_size_of_vector_column(node_mc, "__weights_genie__");

        std::cout << "[plot_score_systematics] nuisance basis summary\n";
        std::cout << "  named GENIE knobs       : " << (use_genie_named ? knobs.size() : 0) << "\n";
        std::cout << "  flux universes          : " << n_flux_u << " (" << (flux_branch.empty() ? "disabled" : flux_branch) << ")\n";
        std::cout << "  reinteraction universes : " << n_reint_u << "\n";
        std::cout << "  GENIE universes         : " << n_genie_u << " (enabled=" << std::boolalpha << use_genie_multiverse << std::noboolalpha << ")\n";
        std::cout << "  data node found         : " << std::boolalpha << have_data_node << std::noboolalpha << "\n";
        std::cout << "  ext node found          : " << std::boolalpha << have_ext_node << std::noboolalpha << "\n";

        if (use_genie_named && use_genie_multiverse)
            std::cout << "[plot_score_systematics] warning: both named GENIE knobs and weightsGenie multiverse are enabled; this usually double counts GENIE model uncertainty.\n";

        SampleBundle acc;
        const std::size_t nbinssz = static_cast<std::size_t>(nbins);
        acc.mc_total.assign(nbinssz, 0.0);
        acc.ext_total.assign(nbinssz, 0.0);
        acc.signal.assign(nbinssz, 0.0);
        acc.dominant_bkg.assign(nbinssz, 0.0);
        acc.data.assign(nbinssz, 0.0);

        acc.diag2_total.assign(nbinssz, 0.0);
        acc.diag2_signal.assign(nbinssz, 0.0);
        acc.diag2_dombkg.assign(nbinssz, 0.0);

        acc.tail_diag2_total.assign(nbinssz, 0.0);
        acc.tail_diag2_signal.assign(nbinssz, 0.0);
        acc.tail_diag2_dombkg.assign(nbinssz, 0.0);

        std::vector<std::vector<double>> h_flux_total(static_cast<std::size_t>(n_flux_u), std::vector<double>(nbinssz, 0.0));
        std::vector<std::vector<double>> h_flux_sig(static_cast<std::size_t>(n_flux_u), std::vector<double>(nbinssz, 0.0));
        std::vector<std::vector<double>> h_flux_dom(static_cast<std::size_t>(n_flux_u), std::vector<double>(nbinssz, 0.0));

        std::vector<std::vector<double>> h_reint_total(static_cast<std::size_t>(n_reint_u), std::vector<double>(nbinssz, 0.0));
        std::vector<std::vector<double>> h_reint_sig(static_cast<std::size_t>(n_reint_u), std::vector<double>(nbinssz, 0.0));
        std::vector<std::vector<double>> h_reint_dom(static_cast<std::size_t>(n_reint_u), std::vector<double>(nbinssz, 0.0));

        std::vector<std::vector<double>> h_genie_total(static_cast<std::size_t>(n_genie_u), std::vector<double>(nbinssz, 0.0));
        std::vector<std::vector<double>> h_genie_sig(static_cast<std::size_t>(n_genie_u), std::vector<double>(nbinssz, 0.0));
        std::vector<std::vector<double>> h_genie_dom(static_cast<std::size_t>(n_genie_u), std::vector<double>(nbinssz, 0.0));

        node_mc.Foreach(
            [&](double s,
                double w_nominal,
                double w_nominal_nofluxcv,
                bool is_signal,
                bool is_dom,
                auto knobRPAup,
                auto knobRPAdn,
                auto knobCCMECup,
                auto knobCCMECdn,
                auto knobAxFFCCQEup,
                auto knobAxFFCCQEdn,
                auto knobVecFFCCQEup,
                auto knobVecFFCCQEdn,
                auto knobDecayAngMECup,
                auto knobDecayAngMECdn,
                auto knobThetaDelta2Npiup,
                auto knobThetaDelta2Npidn,
                auto knobThetaDelta2NRadup,
                auto knobThetaDelta2NRaddn,
                auto knobNormCCCOHup,
                auto knobNormCCCOHdn,
                auto knobNormNCCOHup,
                auto knobNormNCCOHdn,
                auto knobxsr_scc_Fv3up,
                auto knobxsr_scc_Fv3dn,
                auto knobxsr_scc_Fa3up,
                auto knobxsr_scc_Fa3dn,
                const ROOT::RVec<unsigned short> &weights_flux,
                const ROOT::RVec<unsigned short> &weights_reint,
                const ROOT::RVec<unsigned short> &weights_genie) {
                const int b = score_bin(s, nbins, xmin, xmax);
                if (b < 0)
                    return;

                fill_hist(acc.mc_total, b, w_nominal);
                if (is_signal)
                    fill_hist(acc.signal, b, w_nominal);
                if (is_dom)
                    fill_hist(acc.dominant_bkg, b, w_nominal);

                if (!knobs.empty())
                {
                    const double upvals[] = {
                        knobRPAup, knobCCMECup, knobAxFFCCQEup, knobVecFFCCQEup,
                        knobDecayAngMECup, knobThetaDelta2Npiup, knobThetaDelta2NRadup,
                        knobNormCCCOHup, knobNormNCCOHup, knobxsr_scc_Fv3up, knobxsr_scc_Fa3up};
                    const double dnvals[] = {
                        knobRPAdn, knobCCMECdn, knobAxFFCCQEdn, knobVecFFCCQEdn,
                        knobDecayAngMECdn, knobThetaDelta2Npidn, knobThetaDelta2NRaddn,
                        knobNormCCCOHdn, knobNormNCCOHdn, knobxsr_scc_Fv3dn, knobxsr_scc_Fa3dn};

                    const std::size_t nfill = std::min(knobs.size(), static_cast<std::size_t>(11));
                    for (std::size_t i = 0; i < nfill; ++i)
                    {
                        const double w_up = w_nominal * sanitise_multiplier(upvals[i]);
                        const double w_dn = w_nominal * sanitise_multiplier(dnvals[i]);

                        fill_hist(knobs[i].up_total, b, w_up);
                        fill_hist(knobs[i].dn_total, b, w_dn);

                        if (is_signal)
                        {
                            fill_hist(knobs[i].up_sig, b, w_up);
                            fill_hist(knobs[i].dn_sig, b, w_dn);
                        }
                        if (is_dom)
                        {
                            fill_hist(knobs[i].up_dom, b, w_up);
                            fill_hist(knobs[i].dn_dom, b, w_dn);
                        }
                    }
                }

                for (int u = 0; u < n_flux_u; ++u)
                {
                    const double rel = decode_u16_weight(weights_flux, static_cast<std::size_t>(u));
                    const double w = (flux_universes_replace_cv ? w_nominal_nofluxcv : w_nominal) * rel;
                    fill_hist(h_flux_total[static_cast<std::size_t>(u)], b, w);
                    if (is_signal)
                        fill_hist(h_flux_sig[static_cast<std::size_t>(u)], b, w);
                    if (is_dom)
                        fill_hist(h_flux_dom[static_cast<std::size_t>(u)], b, w);
                }

                for (int u = 0; u < n_reint_u; ++u)
                {
                    const double w = w_nominal * decode_u16_weight(weights_reint, static_cast<std::size_t>(u));
                    fill_hist(h_reint_total[static_cast<std::size_t>(u)], b, w);
                    if (is_signal)
                        fill_hist(h_reint_sig[static_cast<std::size_t>(u)], b, w);
                    if (is_dom)
                        fill_hist(h_reint_dom[static_cast<std::size_t>(u)], b, w);
                }

                for (int u = 0; u < n_genie_u; ++u)
                {
                    const double w = w_nominal * decode_u16_weight(weights_genie, static_cast<std::size_t>(u));
                    fill_hist(h_genie_total[static_cast<std::size_t>(u)], b, w);
                    if (is_signal)
                        fill_hist(h_genie_sig[static_cast<std::size_t>(u)], b, w);
                    if (is_dom)
                        fill_hist(h_genie_dom[static_cast<std::size_t>(u)], b, w);
                }
            },
            {"inf_score_0",
             "w_nominal",
             "__w_nominal_nofluxcv__",
             "__is_signal__",
             "__is_dom__",
             "knobRPAup",
             "knobRPAdn",
             "knobCCMECup",
             "knobCCMECdn",
             "knobAxFFCCQEup",
             "knobAxFFCCQEdn",
             "knobVecFFCCQEup",
             "knobVecFFCCQEdn",
             "knobDecayAngMECup",
             "knobDecayAngMECdn",
             "knobThetaDelta2Npiup",
             "knobThetaDelta2Npidn",
             "knobThetaDelta2NRadup",
             "knobThetaDelta2NRaddn",
             "knobNormCCCOHup",
             "knobNormCCCOHdn",
             "knobNormNCCOHup",
             "knobNormNCCOHdn",
             "knobxsr_scc_Fv3up",
             "knobxsr_scc_Fv3dn",
             "knobxsr_scc_Fa3up",
             "knobxsr_scc_Fa3dn",
             "__weights_flux__",
             "__weights_reint__",
             "__weights_genie__"});

        if (have_ext_node)
        {
            node_ext.Foreach(
                [&](double s, double w_nominal) {
                    const int b = score_bin(s, nbins, xmin, xmax);
                    if (b < 0)
                        return;
                    fill_hist(acc.ext_total, b, w_nominal);
                },
                {"inf_score_0", "w_nominal"});
        }

        if (have_data_node)
        {
            node_data.Foreach(
                [&](double s) {
                    const int b = score_bin(s, nbins, xmin, xmax);
                    if (b < 0)
                        return;
                    fill_hist(acc.data, b, 1.0);
                },
                {"inf_score_0"});
        }

        for (const auto &k : knobs)
        {
            add_pair_variance(k.up_total, k.dn_total, acc.diag2_total);
            add_pair_variance(k.up_sig, k.dn_sig, acc.diag2_signal);
            add_pair_variance(k.up_dom, k.dn_dom, acc.diag2_dombkg);

            const auto up_tail_total = add_hists(cumulative_from_high(k.up_total), cumulative_from_high(acc.ext_total));
            const auto dn_tail_total = add_hists(cumulative_from_high(k.dn_total), cumulative_from_high(acc.ext_total));
            const auto up_tail_sig = cumulative_from_high(k.up_sig);
            const auto dn_tail_sig = cumulative_from_high(k.dn_sig);
            const auto up_tail_dom = cumulative_from_high(k.up_dom);
            const auto dn_tail_dom = cumulative_from_high(k.dn_dom);

            add_pair_variance(up_tail_total, dn_tail_total, acc.tail_diag2_total);
            add_pair_variance(up_tail_sig, dn_tail_sig, acc.tail_diag2_signal);
            add_pair_variance(up_tail_dom, dn_tail_dom, acc.tail_diag2_dombkg);
        }

        add_multiverse_variance(h_flux_total, acc.diag2_total);
        add_multiverse_variance(h_flux_sig, acc.diag2_signal);
        add_multiverse_variance(h_flux_dom, acc.diag2_dombkg);

        add_multiverse_variance(h_reint_total, acc.diag2_total);
        add_multiverse_variance(h_reint_sig, acc.diag2_signal);
        add_multiverse_variance(h_reint_dom, acc.diag2_dombkg);

        add_multiverse_variance(h_genie_total, acc.diag2_total);
        add_multiverse_variance(h_genie_sig, acc.diag2_signal);
        add_multiverse_variance(h_genie_dom, acc.diag2_dombkg);

        {
            std::vector<std::vector<double>> tail_flux_total, tail_flux_sig, tail_flux_dom;
            tail_flux_total.reserve(h_flux_total.size());
            tail_flux_sig.reserve(h_flux_sig.size());
            tail_flux_dom.reserve(h_flux_dom.size());

            for (const auto &h : h_flux_total)
                tail_flux_total.push_back(add_hists(cumulative_from_high(h), cumulative_from_high(acc.ext_total)));
            for (const auto &h : h_flux_sig)
                tail_flux_sig.push_back(cumulative_from_high(h));
            for (const auto &h : h_flux_dom)
                tail_flux_dom.push_back(cumulative_from_high(h));

            add_multiverse_variance(tail_flux_total, acc.tail_diag2_total);
            add_multiverse_variance(tail_flux_sig, acc.tail_diag2_signal);
            add_multiverse_variance(tail_flux_dom, acc.tail_diag2_dombkg);
        }

        {
            std::vector<std::vector<double>> tail_reint_total, tail_reint_sig, tail_reint_dom;
            tail_reint_total.reserve(h_reint_total.size());
            tail_reint_sig.reserve(h_reint_sig.size());
            tail_reint_dom.reserve(h_reint_dom.size());

            for (const auto &h : h_reint_total)
                tail_reint_total.push_back(add_hists(cumulative_from_high(h), cumulative_from_high(acc.ext_total)));
            for (const auto &h : h_reint_sig)
                tail_reint_sig.push_back(cumulative_from_high(h));
            for (const auto &h : h_reint_dom)
                tail_reint_dom.push_back(cumulative_from_high(h));

            add_multiverse_variance(tail_reint_total, acc.tail_diag2_total);
            add_multiverse_variance(tail_reint_sig, acc.tail_diag2_signal);
            add_multiverse_variance(tail_reint_dom, acc.tail_diag2_dombkg);
        }

        {
            std::vector<std::vector<double>> tail_genie_total, tail_genie_sig, tail_genie_dom;
            tail_genie_total.reserve(h_genie_total.size());
            tail_genie_sig.reserve(h_genie_sig.size());
            tail_genie_dom.reserve(h_genie_dom.size());

            for (const auto &h : h_genie_total)
                tail_genie_total.push_back(add_hists(cumulative_from_high(h), cumulative_from_high(acc.ext_total)));
            for (const auto &h : h_genie_sig)
                tail_genie_sig.push_back(cumulative_from_high(h));
            for (const auto &h : h_genie_dom)
                tail_genie_dom.push_back(cumulative_from_high(h));

            add_multiverse_variance(tail_genie_total, acc.tail_diag2_total);
            add_multiverse_variance(tail_genie_sig, acc.tail_diag2_signal);
            add_multiverse_variance(tail_genie_dom, acc.tail_diag2_dombkg);
        }

        const std::vector<double> nom_total = add_hists(acc.mc_total, acc.ext_total);
        const std::vector<double> tail_ext_total = cumulative_from_high(acc.ext_total);
        const std::vector<double> tail_nom_total = cumulative_from_high(nom_total);
        const std::vector<double> tail_nom_signal = cumulative_from_high(acc.signal);
        const std::vector<double> tail_nom_dombkg = cumulative_from_high(acc.dominant_bkg);

        const std::vector<double> abs_sigma_total = absolute_sigma(acc.diag2_total);
        const std::vector<double> frac_total = fractional_band(nom_total, acc.diag2_total);
        const std::vector<double> frac_signal = fractional_band(acc.signal, acc.diag2_signal);
        const std::vector<double> frac_dombkg = fractional_band(acc.dominant_bkg, acc.diag2_dombkg);

        const std::vector<double> tail_frac_total = fractional_band(tail_nom_total, acc.tail_diag2_total);
        const std::vector<double> tail_frac_signal = fractional_band(tail_nom_signal, acc.tail_diag2_signal);
        const std::vector<double> tail_frac_dombkg = fractional_band(tail_nom_dombkg, acc.tail_diag2_dombkg);

        std::vector<double> data_err;
        fill_poisson_errors(data_err, acc.data);

        auto tail_curves_total = top_tail_curves(knobs, tail_nom_total, tail_ext_total, 0, static_cast<std::size_t>(std::max(0, max_tail_curves)));
        auto tail_curves_signal = top_tail_curves(knobs, tail_nom_signal, {}, 1, static_cast<std::size_t>(std::max(0, max_tail_curves)));
        auto tail_curves_dombkg = top_tail_curves(knobs, tail_nom_dombkg, {}, 2, static_cast<std::size_t>(std::max(0, max_tail_curves)));

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        // 1) Nominal total prediction + data
        {
            auto h_nom_total = make_hist("h_nom_total", ";Inference score [0];Events / bin", nbins, xmin, xmax, nom_total);
            auto h_nom_band = make_hist("h_nom_band", ";Inference score [0];Events / bin", nbins, xmin, xmax, nom_total, &abs_sigma_total);
            auto h_sig = make_hist("h_sig_only", ";Inference score [0];Events / bin", nbins, xmin, xmax, acc.signal);
            auto h_dom = make_hist("h_dom_only", ";Inference score [0];Events / bin", nbins, xmin, xmax, acc.dominant_bkg);
            auto h_data = make_hist("h_data", ";Inference score [0];Events / bin", nbins, xmin, xmax, acc.data, &data_err);

            style_band_hist(*h_nom_band);
            style_nominal_hist(*h_nom_total);
            style_signal_hist(*h_sig);
            style_dom_hist(*h_dom);
            style_data_hist(*h_data);

            TCanvas c("c_score_nominal", "score nominal", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);
            c.SetLogy();

            double ymax = std::max(max_content(nom_total, &abs_sigma_total), max_content(acc.data, &data_err));
            ymax = std::max(ymax, max_content(acc.signal, nullptr));
            ymax = std::max(ymax, max_content(acc.dominant_bkg, nullptr));
            if (!(ymax > 0.0))
                ymax = 1.0;

            TH1D hframe("hframe_nom", ";Inference score [0];Events / bin", nbins, xmin, xmax);
            hframe.SetMinimum(0.25);
            hframe.SetMaximum(5.0 * ymax);
            hframe.Draw("AXIS");

            h_nom_band->Draw("E2 SAME");
            h_nom_total->Draw("HIST SAME");
            h_sig->Draw("HIST SAME");
            h_dom->Draw("HIST SAME");
            if (have_data_node)
                h_data->Draw("E1 SAME");

            TLegend leg(0.13, 0.68, 0.72, 0.89);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            if (have_data_node)
                leg.AddEntry(h_data.get(), "Data", "lep");
            leg.AddEntry(h_nom_total.get(), "Total prediction", "l");
            leg.AddEntry(h_nom_band.get(), "Selected syst. band", "f");
            leg.AddEntry(h_sig.get(), "Signal only", "l");
            leg.AddEntry(h_dom.get(), dominant_bkg_sel.empty() ? "Background only (current default)" : "Dominant background only", "l");
            leg.Draw();

            c.RedrawAxis();
            const auto out = plot_output_file(output_stem + "_nominal_total").string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_score_systematics] saved: " << out << "\n";
        }

        // 2) Fractional total band
        {
            auto h_frac = make_hist("h_frac_total", ";Inference score [0];#sigma_{syst} / N^{0}", nbins, xmin, xmax, frac_total);
            h_frac->SetLineWidth(3);
            h_frac->SetLineColor(kBlack);

            TCanvas c("c_frac_total", "fractional total band", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            TH1D hframe("hframe_frac_total", ";Inference score [0];#sigma_{syst} / N^{0}", nbins, xmin, xmax);
            hframe.SetMinimum(0.0);
            hframe.SetMaximum(y_max_for_ratio_band(frac_total));
            hframe.Draw("AXIS");
            h_frac->Draw("HIST SAME");

            c.RedrawAxis();
            const auto out = plot_output_file(output_stem + "_fractional_total").string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_score_systematics] saved: " << out << "\n";
        }

        // 3) Fractional signal band
        {
            auto h_frac = make_hist("h_frac_signal", ";Inference score [0];#sigma_{syst} / N^{0}_{signal}", nbins, xmin, xmax, frac_signal);
            h_frac->SetLineWidth(3);
            h_frac->SetLineColor(kRed + 1);

            TCanvas c("c_frac_signal", "fractional signal band", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            TH1D hframe("hframe_frac_signal", ";Inference score [0];#sigma_{syst} / N^{0}_{signal}", nbins, xmin, xmax);
            hframe.SetMinimum(0.0);
            hframe.SetMaximum(y_max_for_ratio_band(frac_signal));
            hframe.Draw("AXIS");
            h_frac->Draw("HIST SAME");

            c.RedrawAxis();
            const auto out = plot_output_file(output_stem + "_fractional_signal").string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_score_systematics] saved: " << out << "\n";
        }

        // 4) Fractional dominant-background band
        {
            auto h_frac = make_hist("h_frac_dombkg", ";Inference score [0];#sigma_{syst} / N^{0}_{dom. bkg}", nbins, xmin, xmax, frac_dombkg);
            h_frac->SetLineWidth(3);
            h_frac->SetLineColor(kBlue + 1);

            TCanvas c("c_frac_dombkg", "fractional dominant background band", 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            TH1D hframe("hframe_frac_dombkg", ";Inference score [0];#sigma_{syst} / N^{0}_{dom. bkg}", nbins, xmin, xmax);
            hframe.SetMinimum(0.0);
            hframe.SetMaximum(y_max_for_ratio_band(frac_dombkg));
            hframe.Draw("AXIS");
            h_frac->Draw("HIST SAME");

            c.RedrawAxis();
            const auto out = plot_output_file(output_stem + "_fractional_dombkg").string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_score_systematics] saved: " << out << "\n";
        }

        auto draw_tail_plot = [&](const std::string &canvas_name,
                                  const std::string &y_title,
                                  const std::vector<double> &tail_frac,
                                  const std::vector<TailCurve> &curves,
                                  const std::string &stem) {
            std::vector<double> ones(static_cast<std::size_t>(nbins), 1.0);
            auto h_band = make_hist("h_tail_band_" + stem,
                                    (";Threshold t on inference score [0];" + y_title).c_str(),
                                    nbins,
                                    xmin,
                                    xmax,
                                    ones,
                                    &tail_frac);
            style_band_hist(*h_band);

            TCanvas c(canvas_name.c_str(), canvas_name.c_str(), 1100, 800);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            double ymax = 1.25 + y_max_for_ratio_band(tail_frac, 0.05);
            double ymin = std::max(0.0, 1.0 - 1.25 * y_max_for_ratio_band(tail_frac, 0.05));

            for (const auto &tc : curves)
            {
                for (double r : tc.ratio)
                {
                    ymax = std::max(ymax, 1.10 * r);
                    ymin = std::min(ymin, 0.90 * r);
                }
            }

            TH1D hframe(("hframe_" + stem).c_str(),
                        (";Threshold t on inference score [0];" + y_title).c_str(),
                        nbins,
                        xmin,
                        xmax);
            hframe.SetMinimum(ymin);
            hframe.SetMaximum(ymax);
            hframe.Draw("AXIS");

            h_band->Draw("E2 SAME");

            TLine line(xmin, 1.0, xmax, 1.0);
            line.SetLineStyle(2);
            line.SetLineWidth(2);
            line.Draw("SAME");

            TLegend leg(0.13, 0.64, 0.78, 0.89);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(h_band.get(), "Combined 1#sigma tail band", "f");

            std::vector<std::unique_ptr<TH1D>> h_curves;
            h_curves.reserve(curves.size());

            for (std::size_t i = 0; i < curves.size(); ++i)
            {
                auto h = make_hist("hcurve_" + stem + "_" + std::to_string(i),
                                   "",
                                   nbins,
                                   xmin,
                                   xmax,
                                   curves[i].ratio);
                style_ratio_curve(*h, static_cast<int>(i));
                h->Draw("HIST SAME");
                leg.AddEntry(h.get(), curves[i].label.c_str(), "l");
                h_curves.push_back(std::move(h));
            }

            leg.Draw();
            c.RedrawAxis();

            const auto out = plot_output_file(output_stem + "_" + stem).string();
            c.SaveAs(out.c_str());
            std::cout << "[plot_score_systematics] saved: " << out << "\n";
        };

        draw_tail_plot("c_tail_total",
                       "R_{k}(t) = N(score > t)^{k} / N(score > t)^{0}",
                       tail_frac_total,
                       tail_curves_total,
                       "tail_total");

        draw_tail_plot("c_tail_signal",
                       "R_{k}(t) = N_{signal}(score > t)^{k} / N_{signal}(score > t)^{0}",
                       tail_frac_signal,
                       tail_curves_signal,
                       "tail_signal");

        draw_tail_plot("c_tail_dombkg",
                       "R_{k}(t) = N_{dom. bkg}(score > t)^{k} / N_{dom. bkg}(score > t)^{0}",
                       tail_frac_dombkg,
                       tail_curves_dombkg,
                       "tail_dombkg");

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[plot_score_systematics] integrals in plotted range\n";
        std::cout << "  total prediction = " << std::accumulate(nom_total.begin(), nom_total.end(), 0.0) << "\n";
        std::cout << "  signal           = " << std::accumulate(acc.signal.begin(), acc.signal.end(), 0.0) << "\n";
        std::cout << "  dominant bkg     = " << std::accumulate(acc.dominant_bkg.begin(), acc.dominant_bkg.end(), 0.0) << "\n";
        std::cout << "  data             = " << std::accumulate(acc.data.begin(), acc.data.end(), 0.0) << "\n";

        std::cout << "[plot_score_systematics] done\n";
        return 0;
    });
}
