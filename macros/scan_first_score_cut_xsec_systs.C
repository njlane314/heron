// scan_first_score_cut_xsec_systs.C
//
// Single-cut expected systematic uncertainty study for a cross-section measurement
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
//     'scan_first_score_cut_xsec_systs("./scratch/out/event_list_myana.root", "sel_muon", "is_signal", "is_signal", "w_nominal", 7.0, true, true, 1.0, 1.0, "score0_cut7_xsec")'

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

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
#include "PlotEnv.hh"
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

double abs_unc_from_var(double abs_var)
{
    if (!(abs_var >= 0.0) || !std::isfinite(abs_var))
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(abs_var);
}

double rel_from_abs_var(double abs_var, double sigma_cv)
{
    if (!(abs_var >= 0.0) || !std::isfinite(abs_var) || !(std::abs(sigma_cv) > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(abs_var) / std::abs(sigma_cv);
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
    ROOT::RDF::RResultPtr<double> sum_sig;
    ROOT::RDF::RResultPtr<double> sum_bkg;
    ROOT::RDF::RResultPtr<double> sum_truth_sig;
    ROOT::RDF::RResultPtr<double> sumw2_sig;
    ROOT::RDF::RResultPtr<double> sumw2_bkg;
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
                         double cut_value,
                         bool keep_greater_than)
{
    ROOT::RDF::RNode node_base = base_sel.empty() ? node : node.Filter(base_sel);

    const std::string pass_col = unique_name("__pass_cut");
    ROOT::RDF::RNode node_pass = node_base
                                     .Define(pass_col,
                                             [cut_value, keep_greater_than](double score) {
                                                 return keep_greater_than ? score >= cut_value : score < cut_value;
                                             },
                                             {"inf_score_0"})
                                     .Filter(pass_col);

    ROOT::RDF::RNode node_sig = node_pass.Filter(signal_sel);
    ROOT::RDF::RNode node_bkg = node_pass.Filter("!(" + signal_sel + ")");
    ROOT::RDF::RNode node_truth_sig = node.Filter(truth_denom_sel);

    const std::string w2_sig_col = unique_name("__w2_sig");
    const std::string w2_bkg_col = unique_name("__w2_bkg");
    const std::string w2_truth_col = unique_name("__w2_truth_sig");

    BookedYields out;
    out.sum_sig = node_sig.Sum<double>(weight_col);
    out.sum_bkg = node_bkg.Sum<double>(weight_col);
    out.sum_truth_sig = node_truth_sig.Sum<double>(weight_col);
    out.sumw2_sig = node_sig.Define(w2_sig_col, weight_col + " * " + weight_col).Sum<double>(w2_sig_col);
    out.sumw2_bkg = node_bkg.Define(w2_bkg_col, weight_col + " * " + weight_col).Sum<double>(w2_bkg_col);
    out.sumw2_truth_sig = node_truth_sig.Define(w2_truth_col, weight_col + " * " + weight_col).Sum<double>(w2_truth_col);
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

EvaluatedYields evaluate_yields(const BookedYields &booked)
{
    EvaluatedYields out;
    out.s_sel = booked.sum_sig.GetValue();
    out.b_sel = booked.sum_bkg.GetValue();
    out.var_s_sel = booked.sumw2_sig.GetValue();
    out.var_b_sel = booked.sumw2_bkg.GetValue();
    out.t_sig = booked.sum_truth_sig.GetValue();
    out.var_t_sig = booked.sumw2_truth_sig.GetValue();
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

double sum_existing_variances(const std::map<std::string, double> &abs_var_map,
                              const std::vector<std::string> &names)
{
    double out = 0.0;
    for (const auto &n : names)
    {
        auto it = abs_var_map.find(n);
        if (it != abs_var_map.end())
            out += it->second;
    }
    return out;
}

} // namespace

int scan_first_score_cut_xsec_systs(const std::string &event_list_path = "",
                                    const std::string &base_sel = "sel_muon",
                                    const std::string &signal_sel = "is_signal",
                                    const std::string &truth_denom_sel = "is_signal",
                                    const std::string &mc_weight = "w_nominal",
                                    double cut_value = 7.0,
                                    bool keep_greater_than = true,
                                    bool fixed_cv_asimov = true,
                                    double flux_times_targets_cv = 1.0,
                                    double flux_times_targets_var_default = 1.0,
                                    const std::string &output_stem = "first_score_cut_xsec_systs_eval")
{
    return heron::macro::run_with_guard("scan_first_score_cut_xsec_systs", [&]() -> int {
        if (implicit_mt_enabled())
            ROOT::EnableImplicitMT();

        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[scan_first_score_cut_xsec_systs] input=" << input_path << "\n";
        std::cout << "[scan_first_score_cut_xsec_systs] cut=" << cut_value
                  << (keep_greater_than ? "  (keep score >= cut)\n" : "  (keep score < cut)\n");

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
            {"xsec_AxFFCCQEshape", "weight_AxFFCCQEshape_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_DecayAngMEC", "weight_DecayAngMEC_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_NormCCCOH", "weight_NormCCCOH_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_NormNCCOH", "weight_NormNCCOH_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_RPA_CCQE", "weight_RPA_CCQE_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_ThetaDelta2NRad", "weight_ThetaDelta2NRad_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_Theta_Delta2Npi", "weight_Theta_Delta2Npi_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_VecFFCCQEshape", "weight_VecFFCCQEshape_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
            {"xsec_XSecShape_CCMEC", "weight_XSecShape_CCMEC_UBGenie", kGenieCvComponent, false, TruthDenomMode::kUseVariedTruthDenom},
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

        if (!has_column(node_cv, "inf_scores"))
        {
            std::cerr << "[scan_first_score_cut_xsec_systs] missing inf_scores branch.\n";
            return 1;
        }

        BookedYields cv_booked = book_yields(node_cv,
                                             base_sel,
                                             signal_sel,
                                             truth_denom_sel,
                                             "__w_cv__",
                                             cut_value,
                                             keep_greater_than);

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
                booked.universes.push_back(book_yields(node_u,
                                                       base_sel,
                                                       signal_sel,
                                                       truth_denom_sel,
                                                       wcol,
                                                       cut_value,
                                                       keep_greater_than));
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
            booked.alt = book_yields(node_alt,
                                     base_sel,
                                     signal_sel,
                                     truth_denom_sel,
                                     "__w_cv__",
                                     cut_value,
                                     keep_greater_than);

            if (!cfg.cv_reference_event_list_path.empty() && looks_like_event_list_root(cfg.cv_reference_event_list_path))
            {
                booked.cv_ref_io = std::make_unique<EventListIO>(cfg.cv_reference_event_list_path);
                ROOT::RDF::RNode node_ref = prepare_mc_node(*booked.cv_ref_io, mc_weight);
                booked.cv_ref = book_yields(node_ref,
                                            base_sel,
                                            signal_sel,
                                            truth_denom_sel,
                                            "__w_cv__",
                                            cut_value,
                                            keep_greater_than);
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
                booked.universes.push_back(book_yields(node_u,
                                                       base_sel,
                                                       signal_sel,
                                                       truth_denom_sel,
                                                       "__w_cv__",
                                                       cut_value,
                                                       keep_greater_than));
            }

            if (!booked.universes.empty())
                booked_altfile_systs.push_back(std::move(booked));
        }

        const EvaluatedYields ycv = evaluate_yields(cv_booked);
        const double sigma_cv = compute_sigma_hat(ycv, ycv, true, TruthDenomMode::kUseCVTruthDenom, flux_times_targets_cv);
        const double eff_cv = safe_div(ycv.s_sel, ycv.t_sig, std::numeric_limits<double>::quiet_NaN());
        const double purity_cv = safe_div(ycv.s_sel, ycv.s_sel + ycv.b_sel, std::numeric_limits<double>::quiet_NaN());

        std::map<std::string, double> abs_var_map;
        abs_var_map["MCstats"] = compute_mcstat_variance_nominal(ycv, fixed_cv_asimov, flux_times_targets_cv);

        for (const auto &cfg : fullcorr_systs)
            abs_var_map[cfg.label] = (cfg.frac * sigma_cv) * (cfg.frac * sigma_cv);

        for (const auto &booked : booked_weight_systs)
        {
            double var = 0.0;
            int nused = 0;
            for (const auto &u : booked.universes)
            {
                const EvaluatedYields yu = evaluate_yields(u);
                const double sig_u = compute_sigma_hat(ycv, yu, fixed_cv_asimov,
                                                       booked.cfg.truth_denom_mode,
                                                       flux_times_targets_var_default);
                if (!std::isfinite(sig_u) || !std::isfinite(sigma_cv))
                    continue;
                const double d = sigma_cv - sig_u;
                var += d * d;
                ++nused;
            }
            if (booked.cfg.average_over_universes && nused > 0)
                var /= static_cast<double>(nused);
            abs_var_map[booked.cfg.label] = var;
        }

        for (const auto &booked : booked_detvars)
        {
            const EvaluatedYields yref = booked.has_cv_ref_override ? evaluate_yields(booked.cv_ref) : ycv;
            const EvaluatedYields yalt = evaluate_yields(booked.alt);
            const double sig_ref = booked.has_cv_ref_override
                                       ? compute_sigma_hat(yref, yref, true, TruthDenomMode::kUseCVTruthDenom, flux_times_targets_cv)
                                       : sigma_cv;
            const double sig_alt = compute_sigma_hat(yref, yalt, fixed_cv_asimov,
                                                     booked.cfg.truth_denom_mode,
                                                     flux_times_targets_var_default);
            const double d = sig_ref - sig_alt;
            abs_var_map[booked.cfg.label] = (std::isfinite(d) ? d * d : 0.0);
        }

        for (const auto &booked : booked_altfile_systs)
        {
            double var = 0.0;
            int nused = 0;
            for (const auto &u : booked.universes)
            {
                const EvaluatedYields yu = evaluate_yields(u);
                const double sig_u = compute_sigma_hat(ycv, yu, fixed_cv_asimov,
                                                       booked.cfg.truth_denom_mode,
                                                       flux_times_targets_var_default);
                if (!std::isfinite(sig_u) || !std::isfinite(sigma_cv))
                    continue;
                const double d = sigma_cv - sig_u;
                var += d * d;
                ++nused;
            }
            if (booked.cfg.average_over_universes && nused > 0)
                var /= static_cast<double>(nused);
            abs_var_map[booked.cfg.label] = var;
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

        abs_var_map["xsec_unisim"] = sum_existing_variances(abs_var_map, xsec_unisim_names);
        abs_var_map["xsec_total"] = sum_existing_variances(abs_var_map, {
            "xsec_multi",
            "xsec_unisim",
            "xsec_xsr_scc_Fa3_SCC",
            "xsec_xsr_scc_Fv3_SCC",
            "NuWroGenie",
        });
        abs_var_map["detVar_total"] = sum_existing_variances(abs_var_map, detvar_names);
        abs_var_map["PredTotal"] = sum_existing_variances(abs_var_map, {
            "detVar_total",
            "flux",
            "reint",
            "xsec_total",
            "POT",
            "numTargets",
            "MCstats",
        });
        abs_var_map["total"] = abs_var_map["PredTotal"];

        std::map<std::string, double> abs_unc_map;
        std::map<std::string, double> rel_map;
        for (const auto &kv : abs_var_map)
        {
            abs_unc_map[kv.first] = abs_unc_from_var(kv.second);
            rel_map[kv.first] = rel_from_abs_var(kv.second, sigma_cv);
        }

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  S_sel(CV)     = " << ycv.s_sel << "\n";
        std::cout << "  B_sel(CV)     = " << ycv.b_sel << "\n";
        std::cout << "  T_sig(CV)     = " << ycv.t_sig << "\n";
        std::cout << "  sigma_hat(CV) = " << sigma_cv << "\n";
        std::cout << "  eff(CV)       = " << eff_cv << "\n";
        std::cout << "  purity(CV)    = " << purity_cv << "\n";

        for (const auto &name : {std::string("flux"),
                                 std::string("reint"),
                                 std::string("xsec_total"),
                                 std::string("detVar_total"),
                                 std::string("MCstats"),
                                 std::string("POT"),
                                 std::string("numTargets"),
                                 std::string("total")})
        {
            auto it_abs = abs_unc_map.find(name);
            auto it_rel = rel_map.find(name);
            if (it_abs == abs_unc_map.end() || it_rel == rel_map.end())
                continue;

            std::cout << "  abs " << std::setw(12) << std::left << name << " = " << it_abs->second << "\n";
            std::cout << "  rel " << std::setw(12) << std::left << name << " = " << it_rel->second << "\n";
        }

        const auto csv_path = plot_output_file(output_stem + "_cut_eval").replace_extension(".csv");
        {
            std::ofstream csv(csv_path);
            csv << "cut,keep_greater_than,s_sel_cv,b_sel_cv,t_sig_cv,sigma_cv,eff_cv,purity_cv";
            for (const auto &kv : abs_unc_map)
                csv << ",abs_" << kv.first;
            for (const auto &kv : rel_map)
                csv << ",rel_" << kv.first;
            csv << "\n";

            csv << std::setprecision(12)
                << cut_value << ","
                << (keep_greater_than ? 1 : 0) << ","
                << ycv.s_sel << ","
                << ycv.b_sel << ","
                << ycv.t_sig << ","
                << sigma_cv << ","
                << eff_cv << ","
                << purity_cv;

            for (const auto &kv : abs_unc_map)
                csv << "," << kv.second;
            for (const auto &kv : rel_map)
                csv << "," << kv.second;
            csv << "\n";
        }
        std::cout << "[scan_first_score_cut_xsec_systs] wrote: " << csv_path << "\n";

        std::cout << "[scan_first_score_cut_xsec_systs] done\n";
        return 0;
    });
}

int scan_score_cut_xsec_systs(const std::string &event_list_path = "",
                              const std::string &base_sel = "sel_muon",
                              const std::string &signal_sel = "is_signal",
                              const std::string &truth_denom_sel = "is_signal",
                              const std::string &mc_weight = "w_nominal",
                              double cut_value = 7.0,
                              bool keep_greater_than = true,
                              bool fixed_cv_asimov = true,
                              double flux_times_targets_cv = 1.0,
                              double flux_times_targets_var_default = 1.0,
                              const std::string &output_stem = "score_cut_xsec_systs_eval")
{
    return scan_first_score_cut_xsec_systs(event_list_path,
                                           base_sel,
                                           signal_sel,
                                           truth_denom_sel,
                                           mc_weight,
                                           cut_value,
                                           keep_greater_than,
                                           fixed_cv_asimov,
                                           flux_times_targets_cv,
                                           flux_times_targets_var_default,
                                           output_stem);
}
