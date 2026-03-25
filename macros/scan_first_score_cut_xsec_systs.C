// scan_first_score_cut_xsec_systs.C
//
// Cut-based expected systematic uncertainty study for a cross-section measurement
// using inf_scores[0] as the selection statistic.
//
// The macro treats the measured cross section as
//
//   sigma_hat(c) = (N_asimov(c) - B_pred(c)) / ( eps_pred(c) * Phi * N_targets )
//
// where c is a hard cut on inf_scores[0].
//
// Two Asimov modes are supported:
//   - fixed_cv_asimov = true:
//       N_asimov(c) = S_CV(c) + B_CV(c) is held fixed while the correction model
//       (background, efficiency, exposure factors) is varied. This is the better
//       default for cut optimisation.
//   - fixed_cv_asimov = false:
//       N_asimov(c) co-varies with the systematic variation. This mirrors the
//       "observable in universe space" covariance construction more directly.
//
// Systematic variances are built source-by-source as scalar analogues of the
// covariance construction in SystematicsCalculator:
//
//   V_k(c) = sum_u [ sigma_CV(c) - sigma_u(c) ]^2
//
// with an optional 1/N_u average when requested by the config.
//
// Important implementation choices:
//   - flux universes and detector variations use the CV truth-signal denominator
//     in the efficiency, matching the usual "vary reco response, keep truth denom
//     fixed" convention.
//   - xsec/reint/multisim universes use the varied truth-signal denominator.
//   - POT and numTargets are implemented as fully-correlated fractional errors.
//   - MCstat is propagated analytically from weighted sums in the nominal CV sample.
//   - BNB/EXT data statistical terms are intentionally not included.
//
// IMPORTANT TO CHECK IN YOUR OWN SETUP:
//   - For each weight systematic, set cv_component_branch correctly. If __w_cv__
//     already contains a CV factor that must be replaced (not multiplied) by the
//     universe weight, put that branch name there. For flux this is often
//     ppfx_cv_weight. For GENIE/reint knobs the correct choice depends on how
//     w_nominal is assembled in your event list.
//   - Fill in the detector-variation and alt-generator file paths.
//
// Example:
//   ./heron macro scan_first_score_cut_xsec_systs.C
//   ./heron macro scan_first_score_cut_xsec_systs.C \
//     'scan_first_score_cut_xsec_systs("./scratch/out/event_list_myana.root", "sel_muon", "is_signal", "is_signal", "w_nominal", 60, -15, 15, true, true, 1.0, 1.0, "score0_cut_xsec")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TGraph.h>
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

std::string unique_name(const std::string &stem)
{
    static unsigned long long counter = 0ULL;
    std::ostringstream ss;
    ss << stem << "_" << counter++;
    return ss.str();
}

double safe_div(double num, double den, double fallback = 0.0)
{
    if (!(den != 0.0) || !std::isfinite(num) || !std::isfinite(den))
        return fallback;
    return num / den;
}

enum class TruthDenomMode
{
    kUseCVTruthDenom,
    kUseVariedTruthDenom,
};

struct WeightSystematic
{
    std::string label;
    std::string universe_branch;
    std::string cv_component_branch; // divide this out of __w_cv__ before applying universe weight; empty => multiply directly
    bool average_over_universes = true;
    TruthDenomMode truth_denom_mode = TruthDenomMode::kUseVariedTruthDenom;
};

struct DetVarSystematic
{
    std::string label;
    std::string alt_event_list_path;
    std::string cv_reference_event_list_path; // optional override; leave empty to use the main CV file
    TruthDenomMode truth_denom_mode = TruthDenomMode::kUseCVTruthDenom;
};

struct AltFileUniverseSystematic
{
    std::string label;
    std::vector<std::string> universe_event_list_paths;
    bool average_over_universes = true;
    TruthDenomMode truth_denom_mode = TruthDenomMode::kUseVariedTruthDenom;
};

struct FullCorrSystematic
{
    std::string label;
    double frac = 0.0;
};

struct BookedYields
{
    ROOT::RDF::RResultPtr<TH1D> h_sig;
    ROOT::RDF::RResultPtr<TH1D> h_bkg;
    ROOT::RDF::RResultPtr<double> sum_truth_sig;
    ROOT::RDF::RResultPtr<double> sumw2_truth_sig;
};

struct EvaluatedYields
{
    double s_sel = 0.0;
    double b_sel = 0.0;
    double t_sig = 0.0;
    double var_s_sel = 0.0;
    double var_b_sel = 0.0;
    double var_t_sig = 0.0;
};

struct BookedWeightSystematic
{
    WeightSystematic cfg;
    std::vector<BookedYields> universes;
};

struct BookedDetVarSystematic
{
    DetVarSystematic cfg;
    std::unique_ptr<EventListIO> alt_io;
    BookedYields alt;
    std::unique_ptr<EventListIO> cv_ref_io;
    BookedYields cv_ref;
    bool has_cv_ref_override = false;
};

struct BookedAltFileSystematic
{
    AltFileUniverseSystematic cfg;
    std::vector<std::unique_ptr<EventListIO>> ios;
    std::vector<BookedYields> universes;
};

struct ScanRow
{
    double cut = 0.0;
    double sigma_cv = std::numeric_limits<double>::quiet_NaN();
    double eff_cv = std::numeric_limits<double>::quiet_NaN();
    double purity_cv = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, double> abs_var;
    std::map<std::string, double> rel_unc;
};

ROOT::RDF::RNode prepare_mc_node(EventListIO &el, const std::string &mc_weight)
{
    ROOT::RDF::RNode node = SelectionService::decorate(el.rdf())
                                .Define("inf_score_0",
                                        [](const ROOT::RVec<float> &scores) {
                                            return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                        },
                                        {"inf_scores"});

    auto mask_mc_like = el.mask_for_mc_like();
    node = filter_by_sample_mask(node, mask_mc_like, "sample_id")
               .Define("__w_cv__", mc_weight)
               .Filter("std::isfinite((double)inf_score_0)")
               .Filter("std::isfinite((double)__w_cv__)");

    return node;
}

BookedYields book_yields(ROOT::RDF::RNode node,
                         const std::string &base_sel,
                         const std::string &signal_sel,
                         const std::string &truth_denom_sel,
                         const std::string &weight_col,
                         const ROOT::RDF::TH1DModel &hmodel_sig,
                         const ROOT::RDF::TH1DModel &hmodel_bkg)
{
    ROOT::RDF::RNode node_sel = base_sel.empty() ? node : node.Filter(base_sel);
    ROOT::RDF::RNode node_sig = node_sel.Filter(signal_sel);
    ROOT::RDF::RNode node_bkg = node_sel.Filter("!(" + signal_sel + ")");
    ROOT::RDF::RNode node_truth_sig = node.Filter(truth_denom_sel);

    const std::string w2_col = unique_name("__w2_truth_sig");

    BookedYields out;
    out.h_sig = node_sig.Histo1D(hmodel_sig, "inf_score_0", weight_col);
    out.h_bkg = node_bkg.Histo1D(hmodel_bkg, "inf_score_0", weight_col);
    out.sum_truth_sig = node_truth_sig.Sum<double>(weight_col);
    out.sumw2_truth_sig = node_truth_sig.Define(w2_col, weight_col + " * " + weight_col).Sum<double>(w2_col);
    return out;
}

std::string make_universe_weight_expr(const std::string &base_weight_col,
                                      const std::string &universe_branch,
                                      int idx,
                                      const std::string &cv_component_branch)
{
    std::ostringstream ss;
    if (cv_component_branch.empty())
    {
        ss << "((" << universe_branch << ".size() > " << idx
           << ") && std::isfinite((double)" << universe_branch << "[" << idx << "])"
           << " && std::isfinite((double)" << base_weight_col << "))"
           << " ? ((double)" << base_weight_col << " * (double)" << universe_branch << "[" << idx << "])"
           << " : 0.0";
    }
    else
    {
        ss << "((" << universe_branch << ".size() > " << idx
           << ") && std::isfinite((double)" << universe_branch << "[" << idx << "])"
           << " && std::isfinite((double)" << cv_component_branch << ")"
           << " && std::abs((double)" << cv_component_branch << ") > 0.0"
           << " && std::isfinite((double)" << base_weight_col << "))"
           << " ? ((double)" << base_weight_col << " * (double)" << universe_branch << "[" << idx << "] / (double)" << cv_component_branch << ")"
           << " : 0.0";
    }
    return ss.str();
}

int get_num_universes(ROOT::RDF::RNode node, const std::string &branch)
{
    const std::string nuni_col = unique_name("__nuni");
    auto nuni = node.Define(nuni_col, "static_cast<int>(" + branch + ".size())").Max<int>(nuni_col);
    return *nuni;
}

int cut_to_first_selected_bin(const TH1D &h, double cut)
{
    const auto *ax = h.GetXaxis();
    if (cut <= ax->GetXmin())
        return 1;
    if (cut > ax->GetXmax())
        return h.GetNbinsX() + 1;

    int b = ax->FindFixBin(cut);
    if (b < 1)
        b = 1;
    if (b > h.GetNbinsX() + 1)
        b = h.GetNbinsX() + 1;
    return b;
}

double integral_hard_cut(const TH1D &h, double cut, bool keep_greater_than)
{
    const int nb = h.GetNbinsX();
    const int b0 = cut_to_first_selected_bin(h, cut);
    if (keep_greater_than)
        return h.Integral(b0, nb + 1);
    return h.Integral(0, b0 - 1);
}

double integral_err2_hard_cut(const TH1D &h, double cut, bool keep_greater_than)
{
    const int nb = h.GetNbinsX();
    const int b0 = cut_to_first_selected_bin(h, cut);
    double out = 0.0;
    if (keep_greater_than)
    {
        for (int b = b0; b <= nb + 1; ++b)
        {
            const double e = h.GetBinError(b);
            out += e * e;
        }
        return out;
    }

    for (int b = 0; b <= b0 - 1; ++b)
    {
        const double e = h.GetBinError(b);
        out += e * e;
    }
    return out;
}

EvaluatedYields evaluate_yields(const BookedYields &booked, double cut, bool keep_greater_than)
{
    const TH1D &hs = *booked.h_sig;
    const TH1D &hb = *booked.h_bkg;

    EvaluatedYields out;
    out.s_sel = integral_hard_cut(hs, cut, keep_greater_than);
    out.b_sel = integral_hard_cut(hb, cut, keep_greater_than);
    out.var_s_sel = integral_err2_hard_cut(hs, cut, keep_greater_than);
    out.var_b_sel = integral_err2_hard_cut(hb, cut, keep_greater_than);
    out.t_sig = *booked.sum_truth_sig;
    out.var_t_sig = *booked.sumw2_truth_sig;
    return out;
}

double compute_sigma_hat(const EvaluatedYields &cv_ref,
                         const EvaluatedYields &varied,
                         bool fixed_cv_asimov,
                         TruthDenomMode truth_mode,
                         double exposure_scale = 1.0)
{
    const double n_asimov = fixed_cv_asimov ? (cv_ref.s_sel + cv_ref.b_sel) : (varied.s_sel + varied.b_sel);
    const double truth_denom = (truth_mode == TruthDenomMode::kUseCVTruthDenom) ? cv_ref.t_sig : varied.t_sig;
    const double eps = safe_div(varied.s_sel, truth_denom, std::numeric_limits<double>::quiet_NaN());

    if (!(eps > 0.0) || !(exposure_scale > 0.0))
        return std::numeric_limits<double>::quiet_NaN();

    return (n_asimov - varied.b_sel) / (eps * exposure_scale);
}

double compute_mcstat_variance_nominal(const EvaluatedYields &cv,
                                       bool fixed_cv_asimov,
                                       double exposure_scale = 1.0)
{
    if (!(cv.t_sig > 0.0) || !(exposure_scale > 0.0))
        return 0.0;

    if (!fixed_cv_asimov)
        return cv.var_t_sig / (exposure_scale * exposure_scale);

    if (!(cv.s_sel > 0.0))
        return 0.0;

    const double r = cv.t_sig / cv.s_sel; // = 1/epsilon
    const double var = cv.var_t_sig + r * r * (cv.var_s_sel + cv.var_b_sel) - 2.0 * r * cv.var_s_sel;
    return std::max(0.0, var) / (exposure_scale * exposure_scale);
}

std::vector<double> sum_curves(const std::vector<std::vector<double>> &curves)
{
    std::vector<double> out;
    if (curves.empty())
        return out;

    out.assign(curves.front().size(), 0.0);
    for (const auto &c : curves)
    {
        if (c.size() != out.size())
            continue;
        for (std::size_t i = 0; i < c.size(); ++i)
            out[i] += c[i];
    }
    return out;
}

std::vector<double> rel_from_abs_var(const std::vector<double> &abs_var, const std::vector<double> &sigma_cv)
{
    std::vector<double> out(abs_var.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < abs_var.size() && i < sigma_cv.size(); ++i)
    {
        if (!(abs_var[i] >= 0.0) || !std::isfinite(abs_var[i]) || !(std::abs(sigma_cv[i]) > 0.0))
            continue;
        out[i] = std::sqrt(abs_var[i]) / std::abs(sigma_cv[i]);
    }
    return out;
}

std::unique_ptr<TGraph> make_graph(const std::vector<double> &x,
                                   const std::vector<double> &y,
                                   const std::string &name)
{
    std::vector<double> xx;
    std::vector<double> yy;
    xx.reserve(x.size());
    yy.reserve(y.size());

    for (std::size_t i = 0; i < x.size() && i < y.size(); ++i)
    {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
            continue;
        xx.push_back(x[i]);
        yy.push_back(y[i]);
    }

    auto g = std::make_unique<TGraph>(static_cast<int>(xx.size()), xx.data(), yy.data());
    g->SetName(name.c_str());
    g->SetLineWidth(2);
    return g;
}

} // namespace

int scan_first_score_cut_xsec_systs(const std::string &event_list_path = "",
                                    const std::string &base_sel = "sel_muon",
                                    const std::string &signal_sel = "is_signal",
                                    const std::string &truth_denom_sel = "is_signal",
                                    const std::string &mc_weight = "w_nominal",
                                    int nbins = 60,
                                    double xmin = -15.0,
                                    double xmax = 15.0,
                                    bool keep_greater_than = true,
                                    bool fixed_cv_asimov = true,
                                    double flux_times_targets_cv = 1.0,
                                    double flux_times_targets_var_default = 1.0,
                                    const std::string &output_stem = "scan_first_score_cut_xsec_systs")
{
    return heron::macro::run_with_guard("scan_first_score_cut_xsec_systs", [&]() -> int {
        if (nbins < 1)
            nbins = 1;
        if (xmax < xmin)
            std::swap(xmin, xmax);
        if (xmax == xmin)
            xmax = xmin + 1.0;

        if (implicit_mt_enabled())
            ROOT::EnableImplicitMT();

        TH1::SetDefaultSumw2();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[scan_first_score_cut_xsec_systs] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[scan_first_score_cut_xsec_systs] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        // ------------------------------------------------------------------
        // Edit this configuration block to match your MCC9 event-list setup.
        // ------------------------------------------------------------------
        const std::string kFluxCvComponent = "ppfx_cv_weight";
        const std::string kGenieCvComponent = ""; // e.g. "tuned_cv_weight" if you need replacement rather than multiplication

        const std::vector<WeightSystematic> weight_systs = {
            {"flux", "weight_flux_all", kFluxCvComponent, true, TruthDenomMode::kUseCVTruthDenom},
            {"reint", "weight_reint_all", kGenieCvComponent, true, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_multi", "weight_All_UBGenie", kGenieCvComponent, true, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_AxFFCCQEshape", "weight_AxFFCCQEshape", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_DecayAngMEC", "weight_DecayAngMEC", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_NormCCCOH", "weight_NormCCCOH", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_NormNCCOH", "weight_NormNCCOH", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_RPA_CCQE", "weight_RPA_CCQE", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_ThetaDelta2NRad", "weight_ThetaDelta2NRad", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_Theta_Delta2Npi", "weight_Theta_Delta2Npi", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_VecFFCCQEshape", "weight_VecFFCCQEshape", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_XSecShape_CCMEC", "weight_XSecShape_CCMEC", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_xsr_scc_Fa3_SCC", "weight_xsr_scc_Fa3_SCC", kGenieCvComponent, true, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_xsr_scc_Fv3_SCC", "weight_xsr_scc_Fv3_SCC", kGenieCvComponent, true, TruthDenomMode::kUseVariedTruthDenom},
        };

        const std::vector<DetVarSystematic> detvar_systs = {
            {"detVarLYatten", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarLYdown", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarLYrayl", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarRecomb2", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarSCE", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarWMAngleXZ", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarWMAngleYZ", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarWMdEdx", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarWMX", "", "", TruthDenomMode::kUseCVTruthDenom},
            {"detVarWMYZ", "", "", TruthDenomMode::kUseCVTruthDenom},
        };

        const std::vector<AltFileUniverseSystematic> altfile_systs = {
            {"NuWroGenie", {}, true, TruthDenomMode::kUseVariedTruthDenom},
        };

        const std::vector<FullCorrSystematic> fullcorr_systs = {
            {"POT", 0.02},
            {"numTargets", 0.01},
        };
        // ------------------------------------------------------------------

        EventListIO el_cv(input_path);
        ROOT::RDF::RNode node_cv = prepare_mc_node(el_cv, mc_weight);

        const ROOT::RDF::TH1DModel hmodel_sig("h_sig_cv", ";inf_scores[0];events", nbins, xmin, xmax);
        const ROOT::RDF::TH1DModel hmodel_bkg("h_bkg_cv", ";inf_scores[0];events", nbins, xmin, xmax);

        if (!has_column(node_cv, "inf_scores"))
        {
            std::cerr << "[scan_first_score_cut_xsec_systs] missing inf_scores branch.\n";
            return 1;
        }

        BookedYields cv_booked = book_yields(node_cv, base_sel, signal_sel, truth_denom_sel, "__w_cv__", hmodel_sig, hmodel_bkg);

        std::vector<BookedWeightSystematic> booked_weight_systs;
        booked_weight_systs.reserve(weight_systs.size());

        for (const auto &cfg : weight_systs)
        {
            if (!has_column(node_cv, cfg.universe_branch))
            {
                std::cout << "[scan_first_score_cut_xsec_systs] skip " << cfg.label
                          << ": missing branch '" << cfg.universe_branch << "'.\n";
                continue;
            }
            if (!cfg.cv_component_branch.empty() && !has_column(node_cv, cfg.cv_component_branch))
            {
                std::cout << "[scan_first_score_cut_xsec_systs] skip " << cfg.label
                          << ": missing cv component branch '" << cfg.cv_component_branch << "'.\n";
                continue;
            }

            const int nuni = get_num_universes(node_cv, cfg.universe_branch);
            std::cout << "[scan_first_score_cut_xsec_systs] booking " << cfg.label << " with " << nuni << " universes\n";

            BookedWeightSystematic booked;
            booked.cfg = cfg;
            booked.universes.reserve(static_cast<std::size_t>(std::max(0, nuni)));

            for (int iu = 0; iu < nuni; ++iu)
            {
                const std::string wcol = unique_name("__w_" + cfg.label);
                const std::string expr = make_universe_weight_expr("__w_cv__", cfg.universe_branch, iu, cfg.cv_component_branch);
                ROOT::RDF::RNode node_u = node_cv.Define(wcol, expr).Filter("std::isfinite((double)" + wcol + ")");
                booked.universes.push_back(book_yields(node_u, base_sel, signal_sel, truth_denom_sel, wcol, hmodel_sig, hmodel_bkg));
            }

            booked_weight_systs.push_back(std::move(booked));
        }

        std::vector<BookedDetVarSystematic> booked_detvars;
        booked_detvars.reserve(detvar_systs.size());
        for (const auto &cfg : detvar_systs)
        {
            if (cfg.alt_event_list_path.empty())
                continue;
            if (!looks_like_event_list_root(cfg.alt_event_list_path))
            {
                std::cout << "[scan_first_score_cut_xsec_systs] skip " << cfg.label
                          << ": invalid alt detvar file '" << cfg.alt_event_list_path << "'.\n";
                continue;
            }

            BookedDetVarSystematic booked;
            booked.cfg = cfg;
            booked.alt_io = std::make_unique<EventListIO>(cfg.alt_event_list_path);
            ROOT::RDF::RNode node_alt = prepare_mc_node(*booked.alt_io, mc_weight);
            booked.alt = book_yields(node_alt, base_sel, signal_sel, truth_denom_sel, "__w_cv__", hmodel_sig, hmodel_bkg);

            if (!cfg.cv_reference_event_list_path.empty() && looks_like_event_list_root(cfg.cv_reference_event_list_path))
            {
                booked.cv_ref_io = std::make_unique<EventListIO>(cfg.cv_reference_event_list_path);
                ROOT::RDF::RNode node_ref = prepare_mc_node(*booked.cv_ref_io, mc_weight);
                booked.cv_ref = book_yields(node_ref, base_sel, signal_sel, truth_denom_sel, "__w_cv__", hmodel_sig, hmodel_bkg);
                booked.has_cv_ref_override = true;
            }

            booked_detvars.push_back(std::move(booked));
        }

        std::vector<BookedAltFileSystematic> booked_altfile_systs;
        booked_altfile_systs.reserve(altfile_systs.size());
        for (const auto &cfg : altfile_systs)
        {
            if (cfg.universe_event_list_paths.empty())
                continue;

            BookedAltFileSystematic booked;
            booked.cfg = cfg;
            booked.ios.reserve(cfg.universe_event_list_paths.size());
            booked.universes.reserve(cfg.universe_event_list_paths.size());

            for (const auto &path : cfg.universe_event_list_paths)
            {
                if (!looks_like_event_list_root(path))
                    continue;

                booked.ios.push_back(std::make_unique<EventListIO>(path));
                ROOT::RDF::RNode node_u = prepare_mc_node(*booked.ios.back(), mc_weight);
                booked.universes.push_back(book_yields(node_u, base_sel, signal_sel, truth_denom_sel, "__w_cv__", hmodel_sig, hmodel_bkg));
            }

            if (!booked.universes.empty())
                booked_altfile_systs.push_back(std::move(booked));
        }

        // Force evaluation of the nominal result first. RDataFrame lazy actions will then execute.
        const TH1D &h_sig_cv = *cv_booked.h_sig;
        const TH1D &h_bkg_cv = *cv_booked.h_bkg;
        (void)h_sig_cv;
        (void)h_bkg_cv;

        std::vector<double> cuts;
        cuts.reserve(static_cast<std::size_t>(nbins + 1));
        for (int b = 1; b <= nbins + 1; ++b)
            cuts.push_back(h_sig_cv.GetXaxis()->GetBinLowEdge(b));

        std::vector<double> sigma_cv(cuts.size(), std::numeric_limits<double>::quiet_NaN());
        std::vector<double> eff_cv(cuts.size(), std::numeric_limits<double>::quiet_NaN());
        std::vector<double> purity_cv(cuts.size(), std::numeric_limits<double>::quiet_NaN());
        std::map<std::string, std::vector<double>> abs_var_map;

        auto ensure_curve = [&](const std::string &name) -> std::vector<double> & {
            auto it = abs_var_map.find(name);
            if (it == abs_var_map.end())
                it = abs_var_map.emplace(name, std::vector<double>(cuts.size(), 0.0)).first;
            return it->second;
        };

        for (std::size_t ic = 0; ic < cuts.size(); ++ic)
        {
            const double cut = cuts[ic];
            const EvaluatedYields ycv = evaluate_yields(cv_booked, cut, keep_greater_than);
            const double sig_cv = compute_sigma_hat(ycv, ycv, true, TruthDenomMode::kUseCVTruthDenom, flux_times_targets_cv);
            sigma_cv[ic] = sig_cv;
            eff_cv[ic] = safe_div(ycv.s_sel, ycv.t_sig, std::numeric_limits<double>::quiet_NaN());
            purity_cv[ic] = safe_div(ycv.s_sel, ycv.s_sel + ycv.b_sel, std::numeric_limits<double>::quiet_NaN());

            ensure_curve("MCstats")[ic] = compute_mcstat_variance_nominal(ycv, fixed_cv_asimov, flux_times_targets_cv);

            for (const auto &cfg : fullcorr_systs)
                ensure_curve(cfg.label)[ic] = (cfg.frac * sig_cv) * (cfg.frac * sig_cv);

            for (const auto &booked : booked_weight_systs)
            {
                double var = 0.0;
                int nused = 0;
                for (const auto &u : booked.universes)
                {
                    const EvaluatedYields yu = evaluate_yields(u, cut, keep_greater_than);
                    const double sig_u = compute_sigma_hat(ycv, yu, fixed_cv_asimov,
                                                           booked.cfg.truth_denom_mode,
                                                           flux_times_targets_var_default);
                    if (!std::isfinite(sig_u) || !std::isfinite(sig_cv))
                        continue;
                    const double d = sig_cv - sig_u;
                    var += d * d;
                    ++nused;
                }
                if (booked.cfg.average_over_universes && nused > 0)
                    var /= static_cast<double>(nused);
                ensure_curve(booked.cfg.label)[ic] = var;
            }

            for (const auto &booked : booked_detvars)
            {
                const EvaluatedYields yref = booked.has_cv_ref_override
                                                 ? evaluate_yields(booked.cv_ref, cut, keep_greater_than)
                                                 : ycv;
                const EvaluatedYields yalt = evaluate_yields(booked.alt, cut, keep_greater_than);
                const double sig_ref = booked.has_cv_ref_override
                                           ? compute_sigma_hat(yref, yref, true, TruthDenomMode::kUseCVTruthDenom, flux_times_targets_cv)
                                           : sig_cv;
                const double sig_alt = compute_sigma_hat(yref, yalt, fixed_cv_asimov,
                                                         booked.cfg.truth_denom_mode,
                                                         flux_times_targets_var_default);
                const double d = sig_ref - sig_alt;
                ensure_curve(booked.cfg.label)[ic] = (std::isfinite(d) ? d * d : 0.0);
            }

            for (const auto &booked : booked_altfile_systs)
            {
                double var = 0.0;
                int nused = 0;
                for (const auto &u : booked.universes)
                {
                    const EvaluatedYields yu = evaluate_yields(u, cut, keep_greater_than);
                    const double sig_u = compute_sigma_hat(ycv, yu, fixed_cv_asimov,
                                                           booked.cfg.truth_denom_mode,
                                                           flux_times_targets_var_default);
                    if (!std::isfinite(sig_u) || !std::isfinite(sig_cv))
                        continue;
                    const double d = sig_cv - sig_u;
                    var += d * d;
                    ++nused;
                }
                if (booked.cfg.average_over_universes && nused > 0)
                    var /= static_cast<double>(nused);
                ensure_curve(booked.cfg.label)[ic] = var;
            }
        }

        const std::vector<std::string> xsec_unisim_names = {
            "xsec_AxFFCCQEshape",
            "xsec_DecayAngMEC",
            "xsec_NormCCCOH",
            "xsec_NormNCCOH",
            "xsec_RPA_CCQE",
            "xsec_ThetaDelta2NRad",
            "xsec_Theta_Delta2Npi",
            "xsec_VecFFCCQEshape",
            "xsec_XSecShape_CCMEC",
        };

        const std::vector<std::string> detvar_names = {
            "detVarLYatten",
            "detVarLYdown",
            "detVarLYrayl",
            "detVarRecomb2",
            "detVarSCE",
            "detVarWMAngleXZ",
            "detVarWMAngleYZ",
            "detVarWMdEdx",
            "detVarWMX",
            "detVarWMYZ",
        };

        auto collect_existing = [&](const std::vector<std::string> &names) {
            std::vector<std::vector<double>> curves;
            curves.reserve(names.size());
            for (const auto &n : names)
            {
                auto it = abs_var_map.find(n);
                if (it != abs_var_map.end())
                    curves.push_back(it->second);
            }
            return curves;
        };

        abs_var_map["xsec_unisim"] = sum_curves(collect_existing(xsec_unisim_names));
        abs_var_map["xsec_total"] = sum_curves({
            abs_var_map.count("xsec_multi") ? abs_var_map.at("xsec_multi") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("xsec_unisim") ? abs_var_map.at("xsec_unisim") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("xsec_xsr_scc_Fa3_SCC") ? abs_var_map.at("xsec_xsr_scc_Fa3_SCC") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("xsec_xsr_scc_Fv3_SCC") ? abs_var_map.at("xsec_xsr_scc_Fv3_SCC") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("NuWroGenie") ? abs_var_map.at("NuWroGenie") : std::vector<double>(cuts.size(), 0.0),
        });
        abs_var_map["detVar_total"] = sum_curves(collect_existing(detvar_names));
        abs_var_map["PredTotal"] = sum_curves({
            abs_var_map.count("detVar_total") ? abs_var_map.at("detVar_total") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("flux") ? abs_var_map.at("flux") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("reint") ? abs_var_map.at("reint") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("xsec_total") ? abs_var_map.at("xsec_total") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("POT") ? abs_var_map.at("POT") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("numTargets") ? abs_var_map.at("numTargets") : std::vector<double>(cuts.size(), 0.0),
            abs_var_map.count("MCstats") ? abs_var_map.at("MCstats") : std::vector<double>(cuts.size(), 0.0),
        });
        abs_var_map["total"] = abs_var_map["PredTotal"];

        std::map<std::string, std::vector<double>> rel_map;
        for (const auto &kv : abs_var_map)
            rel_map[kv.first] = rel_from_abs_var(kv.second, sigma_cv);

        std::size_t best_idx = 0;
        double best_rel = std::numeric_limits<double>::infinity();
        auto it_total = rel_map.find("total");
        if (it_total != rel_map.end())
        {
            for (std::size_t i = 0; i < it_total->second.size(); ++i)
            {
                const double v = it_total->second[i];
                if (std::isfinite(v) && v < best_rel)
                {
                    best_rel = v;
                    best_idx = i;
                }
            }
        }

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[scan_first_score_cut_xsec_systs] best cut = " << cuts[best_idx]
                  << (keep_greater_than ? "  for score >= cut\n" : "  for score < cut\n");
        std::cout << "  eff(CV)    = " << eff_cv[best_idx] << "\n";
        std::cout << "  purity(CV) = " << purity_cv[best_idx] << "\n";
        if (it_total != rel_map.end())
            std::cout << "  rel total  = " << it_total->second[best_idx] << "\n";

        for (const auto &name : {std::string("flux"), std::string("reint"), std::string("xsec_total"), std::string("detVar_total"), std::string("MCstats"), std::string("POT"), std::string("numTargets")})
        {
            auto it = rel_map.find(name);
            if (it != rel_map.end())
                std::cout << "  rel " << std::setw(12) << std::left << name << " = " << it->second[best_idx] << "\n";
        }

        const auto csv_path = plot_output_file(output_stem + "_scan").replace_extension(".csv");
        {
            std::ofstream csv(csv_path);
            csv << "cut,sigma_cv,eff_cv,purity_cv";
            for (const auto &kv : rel_map)
                csv << ",rel_" << kv.first;
            csv << "\n";

            for (std::size_t i = 0; i < cuts.size(); ++i)
            {
                csv << std::setprecision(12)
                    << cuts[i] << "," << sigma_cv[i] << "," << eff_cv[i] << "," << purity_cv[i];
                for (const auto &kv : rel_map)
                {
                    double v = std::numeric_limits<double>::quiet_NaN();
                    if (i < kv.second.size())
                        v = kv.second[i];
                    csv << "," << v;
                }
                csv << "\n";
            }
        }
        std::cout << "[scan_first_score_cut_xsec_systs] wrote: " << csv_path << "\n";

        Plotter plotter;
        plotter.set_global_style();
        gStyle->SetOptStat(0);

        {
            TCanvas c("c_cut_syst_total", "Cut scan total/systematic groups", 1200, 850);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            double ymax = 0.0;
            for (const auto &name : {std::string("total"), std::string("flux"), std::string("reint"), std::string("xsec_total"), std::string("detVar_total"), std::string("MCstats")})
            {
                auto it = rel_map.find(name);
                if (it == rel_map.end())
                    continue;
                for (double v : it->second)
                    if (std::isfinite(v))
                        ymax = std::max(ymax, v);
            }
            if (!(ymax > 0.0))
                ymax = 1.0;

            TH1D frame("hframe_cut_syst_total",
                       keep_greater_than ? ";cut on inf_scores[0]  (keep score #geq cut);relative uncertainty on #hat{#sigma}"
                                         : ";cut on inf_scores[0]  (keep score < cut);relative uncertainty on #hat{#sigma}",
                       100, xmin, xmax);
            frame.SetMinimum(0.0);
            frame.SetMaximum(1.15 * ymax);
            frame.Draw("AXIS");

            auto g_total = make_graph(cuts, rel_map["total"], "g_total");
            auto g_flux = rel_map.count("flux") ? make_graph(cuts, rel_map["flux"], "g_flux") : nullptr;
            auto g_reint = rel_map.count("reint") ? make_graph(cuts, rel_map["reint"], "g_reint") : nullptr;
            auto g_xsec = rel_map.count("xsec_total") ? make_graph(cuts, rel_map["xsec_total"], "g_xsec") : nullptr;
            auto g_det = rel_map.count("detVar_total") ? make_graph(cuts, rel_map["detVar_total"], "g_det") : nullptr;
            auto g_mc = rel_map.count("MCstats") ? make_graph(cuts, rel_map["MCstats"], "g_mc") : nullptr;

            g_total->SetLineWidth(3);
            g_total->Draw("L SAME");
            if (g_flux)
            {
                g_flux->SetLineStyle(2);
                g_flux->Draw("L SAME");
            }
            if (g_reint)
            {
                g_reint->SetLineStyle(3);
                g_reint->Draw("L SAME");
            }
            if (g_xsec)
            {
                g_xsec->SetLineStyle(4);
                g_xsec->Draw("L SAME");
            }
            if (g_det)
            {
                g_det->SetLineStyle(5);
                g_det->Draw("L SAME");
            }
            if (g_mc)
            {
                g_mc->SetLineStyle(6);
                g_mc->Draw("L SAME");
            }

            TLine lbest(cuts[best_idx], 0.0, cuts[best_idx], 1.15 * ymax);
            lbest.SetLineStyle(7);
            lbest.SetLineWidth(2);
            lbest.Draw("SAME");

            TLegend leg(0.12, 0.62, 0.54, 0.90);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_total.get(), "total", "l");
            if (g_flux)
                leg.AddEntry(g_flux.get(), "flux", "l");
            if (g_reint)
                leg.AddEntry(g_reint.get(), "reint", "l");
            if (g_xsec)
                leg.AddEntry(g_xsec.get(), "xsec_total", "l");
            if (g_det)
                leg.AddEntry(g_det.get(), "detVar_total", "l");
            if (g_mc)
                leg.AddEntry(g_mc.get(), "MCstats", "l");
            leg.AddEntry(&lbest, "best cut", "l");
            leg.Draw();

            c.RedrawAxis();
            const auto out = plot_output_file(output_stem + "_rel_unc_groups").string();
            c.SaveAs(out.c_str());
            std::cout << "[scan_first_score_cut_xsec_systs] saved: " << out << "\n";
        }

        {
            TCanvas c("c_cut_eff_purity", "Cut scan efficiency/purity", 1200, 850);
            c.SetLeftMargin(0.11);
            c.SetRightMargin(0.04);
            c.SetBottomMargin(0.12);

            TH1D frame("hframe_cut_perf",
                       keep_greater_than ? ";cut on inf_scores[0]  (keep score #geq cut);CV efficiency / purity"
                                         : ";cut on inf_scores[0]  (keep score < cut);CV efficiency / purity",
                       100, xmin, xmax);
            frame.SetMinimum(0.0);
            frame.SetMaximum(1.05);
            frame.Draw("AXIS");

            auto g_eff = make_graph(cuts, eff_cv, "g_eff");
            auto g_pur = make_graph(cuts, purity_cv, "g_pur");
            g_eff->SetLineWidth(3);
            g_pur->SetLineStyle(2);
            g_eff->Draw("L SAME");
            g_pur->Draw("L SAME");

            TLine lbest(cuts[best_idx], 0.0, cuts[best_idx], 1.05);
            lbest.SetLineStyle(7);
            lbest.SetLineWidth(2);
            lbest.Draw("SAME");

            TLegend leg(0.12, 0.75, 0.42, 0.90);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_eff.get(), "efficiency (CV)", "l");
            leg.AddEntry(g_pur.get(), "purity (CV)", "l");
            leg.AddEntry(&lbest, "best cut", "l");
            leg.Draw();

            c.RedrawAxis();
            const auto out = plot_output_file(output_stem + "_eff_purity").string();
            c.SaveAs(out.c_str());
            std::cout << "[scan_first_score_cut_xsec_systs] saved: " << out << "\n";
        }

        std::cout << "[scan_first_score_cut_xsec_systs] done\n";
        return 0;
    });
}

int scan_score_cut_xsec_systs(const std::string &event_list_path = "",
                              const std::string &base_sel = "sel_muon",
                              const std::string &signal_sel = "is_signal",
                              const std::string &truth_denom_sel = "is_signal",
                              const std::string &mc_weight = "w_nominal",
                              int nbins = 60,
                              double xmin = -15.0,
                              double xmax = 15.0,
                              bool keep_greater_than = true,
                              bool fixed_cv_asimov = true,
                              double flux_times_targets_cv = 1.0,
                              double flux_times_targets_var_default = 1.0,
                              const std::string &output_stem = "scan_score_cut_xsec_systs")
{
    return scan_first_score_cut_xsec_systs(event_list_path,
                                           base_sel,
                                           signal_sel,
                                           truth_denom_sel,
                                           mc_weight,
                                           nbins,
                                           xmin,
                                           xmax,
                                           keep_greater_than,
                                           fixed_cv_asimov,
                                           flux_times_targets_cv,
                                           flux_times_targets_var_default,
                                           output_stem);
}
