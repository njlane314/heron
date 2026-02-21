/* -- C++ -- */
/**
 *  @file  ana/src/DefaultAnalysisModel.cpp
 *
 *  @brief Default analysis model that provides concrete service logic.
 */

#include "DefaultAnalysisModel.hh"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RVec.hxx>

#include "AnalysisChannels.hh"

namespace
{

constexpr float min_x = 5.f;
constexpr float max_x = 251.f;
constexpr float min_y = -110.f;
constexpr float max_y = 110.f;
constexpr float min_z = 20.f;
constexpr float max_z = 986.f;

constexpr float reco_gap_min_z = 675.f;
constexpr float reco_gap_max_z = 775.f;

inline ROOT::RDF::RNode filter_on(ROOT::RDF::RNode node, const char *col)
{
    return node.Filter([](bool pass) { return pass; }, {col});
}

template <typename T>
bool is_within(const T &value, float low, float high)
{
    return value > low && value < high;
}

template <typename X, typename Y, typename Z>
bool is_in_active_volume(const X &x, const Y &y, const Z &z)
{
    return is_within(x, min_x, max_x) &&
           is_within(y, min_y, max_y) &&
           is_within(z, min_z, max_z);
}

} // namespace

const DefaultAnalysisModel &DefaultAnalysisModel::instance()
{
    static const DefaultAnalysisModel model{};
    return model;
}

DefaultAnalysisModel::DefaultAnalysisModel()
{
    m_name = "heron";
    m_tree_name = "nuselection/EventSelectionFilter";
}

const std::string &DefaultAnalysisModel::name() const noexcept
{
    return m_name;
}

const std::string &DefaultAnalysisModel::tree_name() const noexcept
{
    return m_tree_name;
}

ProcessorEntry DefaultAnalysisModel::make_processor(const SampleIO::Sample &sample) const noexcept
{
    ProcessorEntry proc_entry;

    switch (sample.origin)
    {
    case SampleIO::SampleOrigin::kData:
        proc_entry.source = Type::kData;
        break;
    case SampleIO::SampleOrigin::kEXT:
        proc_entry.source = Type::kExt;
        proc_entry.trig_nom = sample.db_tor101_pot_sum;
        proc_entry.trig_eqv = sample.subrun_pot_sum;
        break;
    case SampleIO::SampleOrigin::kOverlay:
    case SampleIO::SampleOrigin::kDirt:
    case SampleIO::SampleOrigin::kStrangeness:
        proc_entry.source = Type::kMC;
        proc_entry.pot_nom = sample.db_tortgt_pot_sum;
        proc_entry.pot_eqv = sample.subrun_pot_sum;
        break;
    default:
        proc_entry.source = Type::kUnknown;
        break;
    }

    return proc_entry;
}

ROOT::RDF::RNode DefaultAnalysisModel::define(ROOT::RDF::RNode node, const ProcessorEntry &rec) const
{
    const bool is_data = (rec.source == Type::kData);
    const bool is_ext = (rec.source == Type::kExt);
    const bool is_mc = (rec.source == Type::kMC);

    const double scale_mc =
        (is_mc && rec.pot_nom > 0.0 && rec.pot_eqv > 0.0) ? (rec.pot_nom / rec.pot_eqv) : 1.0;
    const double scale_ext =
        (is_ext && rec.trig_nom > 0.0 && rec.trig_eqv > 0.0) ? (rec.trig_nom / rec.trig_eqv) : 1.0;

    node = node.Define("w_base", [is_mc, is_ext, scale_mc, scale_ext]() -> double {
        const double scale = is_mc ? scale_mc : (is_ext ? scale_ext : 1.0);
        return scale;
    });

    {
        const auto cnames = node.GetColumnNames();
        auto has = [&](const std::string &name) {
            return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
        };

        if (!has("ppfx_cv"))
            node = node.Define("ppfx_cv", [] { return 1.0f; });
        if (!has("weightSpline"))
            node = node.Define("weightSpline", [] { return 1.0f; });
        if (!has("weightTune"))
            node = node.Define("weightTune", [] { return 1.0f; });
        if (!has("RootinoFix"))
            node = node.Define("RootinoFix", [] { return 1.0; });
    }

    if (is_mc)
    {
        node = node.Define(
            "w_nominal",
            [](double w_base, float w_spline, float w_tune, float w_flux_cv, double w_root) -> double {
                auto sanitise_weight = [](double w) {
                    if (!std::isfinite(w) || w <= 0.0)
                        return 1.0;
                    return w;
                };
                const double out = w_base *
                                   sanitise_weight(w_spline) *
                                   sanitise_weight(w_tune) *
                                   sanitise_weight(w_flux_cv) *
                                   sanitise_weight(w_root);
                if (!std::isfinite(out))
                    return 0.0;
                if (out < 0.0)
                    return 0.0;
                return out;
            },
            {"w_base", "weightSpline", "weightTune", "ppfx_cv", "RootinoFix"});
    }
    else
    {
        node = node.Define("w_nominal", [](double w) -> double { return w; }, {"w_base"});
    }

    if (is_mc)
    {
        node = node.Define(
            "in_fiducial",
            [this](float x, float y, float z) {
                return is_in_truth_volume(x, y, z);
            },
            {"nu_vtx_x", "nu_vtx_y", "nu_vtx_z"});

        node = node.Define(
            "count_strange",
            [](int kplus, int kminus, int kzero, int lambda0, int sigplus, int sigzero, int sigminus) {
                return kplus + kminus + kzero + lambda0 + sigplus + sigzero + sigminus;
            },
            {"n_K_plus", "n_K_minus", "n_K0", "n_lambda", "n_sigma_plus", "n_sigma0", "n_sigma_minus"});

        node = node.Define("is_strange", [](int strange) { return strange > 0; }, {"count_strange"});

        {
            const auto mc_cnames = node.GetColumnNames();
            auto has_mc = [&](const std::string &name) {
                return std::find(mc_cnames.begin(), mc_cnames.end(), name) != mc_cnames.end();
            };

            if (!has_mc("interaction_mode"))
            {
                if (has_mc("int_mode"))
                    node = node.Define("interaction_mode", [](int m) { return m; }, {"int_mode"});
                else
                    node = node.Define("interaction_mode", [] { return -1; });
            }

            if (has_mc("interaction_type"))
            {
            }
            else if (has_mc("int_type"))
            {
                node = node.Define("interaction_type", [](int t) { return t; }, {"int_type"});
            }
            else if (has_mc("interaction_mode"))
            {
                node = node.Define("interaction_type", [](int m) { return m; }, {"interaction_mode"});
            }
            else if (has_mc("int_mode"))
            {
                node = node.Define("interaction_type", [](int m) { return m; }, {"int_mode"});
            }
            else
            {
                node = node.Define("interaction_type", [] { return -1; });
            }
        }

        node = node.Define(
            "analysis_channels",
            [](bool in_fiducial,
               int nu_pdg,
               int ccnc,
               int interaction_type,
               int n_p,
               int n_pi_minus,
               int n_pi_plus,
               int n_pi0,
               int n_gamma,
               bool is_nu_mu_cc,
               int lam_pdg,
               float mu_p,
               float p_p,
               float pi_p,
               float lam_decay_sep) {
                return AnalysisChannels::to_int(
                    AnalysisChannels::classify_analysis_channel(
                        in_fiducial,
                        nu_pdg,
                        ccnc,
                        interaction_type,
                        n_p,
                        n_pi_minus,
                        n_pi_plus,
                        n_pi0,
                        n_gamma,
                        is_nu_mu_cc,
                        lam_pdg,
                        mu_p,
                        p_p,
                        pi_p,
                        lam_decay_sep));
            },
            {"in_fiducial",
             "nu_pdg",
             "int_ccnc",
             "interaction_type",
             "n_p",
             "n_pi_minus",
             "n_pi_plus",
             "n_pi0",
             "n_gamma",
             "is_nu_mu_cc",
             "lam_pdg",
             "mu_p",
             "p_p",
             "pi_p",
             "lam_decay_sep"});

        node = node.Define(
            "is_signal",
            [](bool is_nu_mu_cc, int ccnc, bool in_fiducial, int lam_pdg, float mu_p, float p_p, float pi_p, float lam_decay_sep) {
                return AnalysisChannels::is_signal(
                    is_nu_mu_cc,
                    ccnc,
                    in_fiducial,
                    lam_pdg,
                    mu_p,
                    p_p,
                    pi_p,
                    lam_decay_sep);
            },
            {"is_nu_mu_cc", "int_ccnc", "in_fiducial", "lam_pdg", "mu_p", "p_p", "pi_p", "lam_decay_sep"});
    }
    else
    {
        const int nonmc_channel = is_ext ? static_cast<int>(AnalysisChannels::AnalysisChannel::External)
                                : (is_data ? static_cast<int>(AnalysisChannels::AnalysisChannel::DataInclusive)
                                : static_cast<int>(AnalysisChannels::AnalysisChannel::Unknown));

        auto has_nonmc = [&node](const std::string &name) {
            const auto cnames = node.GetColumnNames();
            return std::find(cnames.begin(), cnames.end(), name) != cnames.end();
        };

        if (!has_nonmc("nu_vtx_x"))
            node = node.Define("nu_vtx_x", [] { return -9999.0f; });
        if (!has_nonmc("nu_vtx_y"))
            node = node.Define("nu_vtx_y", [] { return -9999.0f; });
        if (!has_nonmc("nu_vtx_z"))
            node = node.Define("nu_vtx_z", [] { return -9999.0f; });

        if (!has_nonmc("in_fiducial"))
            node = node.Define("in_fiducial", [] { return false; });
        if (!has_nonmc("is_strange"))
            node = node.Define("is_strange", [] { return false; });
        if (!has_nonmc("analysis_channels"))
            node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        if (!has_nonmc("interaction_mode"))
            node = node.Define("interaction_mode", [] { return -1; });
        if (!has_nonmc("interaction_type"))
            node = node.Define("interaction_type", [] { return -1; });
        if (!has_nonmc("is_signal"))
            node = node.Define("is_signal", [] { return false; });
        if (!has_nonmc("recognised_signal"))
            node = node.Define("recognised_signal", [] { return false; });
    }

    node = node.Define(
        "in_reco_fiducial",
        [this](float x, float y, float z) {
            return is_in_reco_volume(x, y, z);
        },
        {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y", "reco_neutrino_vertex_sce_z"});

    node = decorate(node);

    return node;
}

const char *DefaultAnalysisModel::filter_stage(SampleIO::SampleOrigin origin) const
{
    using SampleOrigin = SampleIO::SampleOrigin;

    if (origin == SampleOrigin::kOverlay)
        return "filter_overlay";
    if (origin == SampleOrigin::kStrangeness)
        return "filter_strangeness";
    return NULL;
}

ROOT::RDF::RNode DefaultAnalysisModel::apply(ROOT::RDF::RNode node, SampleIO::SampleOrigin origin) const
{
    using SampleOrigin = SampleIO::SampleOrigin;

    if (origin == SampleOrigin::kOverlay)
        return node.Filter([](int strange) { return strange == 0; }, {"count_strange"});
    if (origin == SampleOrigin::kStrangeness)
        return node.Filter([](int strange) { return strange > 0; }, {"count_strange"});
    return node;
}

ROOT::RDataFrame DefaultAnalysisModel::load_sample(const SampleIO::Sample &sample,
                                                   const std::string &tree_name) const
{
    std::vector<std::string> files = SampleIO::resolve_root_files(sample);
    return ROOT::RDataFrame(tree_name, files);
}

ROOT::RDF::RNode DefaultAnalysisModel::define_variables(ROOT::RDF::RNode node,
                                                        const std::vector<Column> &definitions) const
{
    ROOT::RDF::RNode updated_node = std::move(node);
    for (const Column &definition : definitions)
    {
        updated_node = updated_node.Define(definition.name, definition.expression);
    }
    return updated_node;
}

int DefaultAnalysisModel::slice_required_count() const noexcept { return 1; }
float DefaultAnalysisModel::slice_min_topology_score() const noexcept { return 0.06f; }
float DefaultAnalysisModel::muon_min_track_score() const noexcept { return 0.5f; }
float DefaultAnalysisModel::muon_min_track_length() const noexcept { return 10.0f; }
float DefaultAnalysisModel::muon_max_track_distance() const noexcept { return 4.0f; }
unsigned DefaultAnalysisModel::muon_required_generation() const noexcept { return 2u; }

ROOT::RDF::RNode DefaultAnalysisModel::apply(ROOT::RDF::RNode node, Preset p, SelectionEntry *selection) const
{
    node = decorate(node);

    ROOT::RDF::RNode filtered = node;
    switch (p)
    {
    case Preset::Empty:
        break;
    case Preset::Trigger:
        filtered = filter_on(node, "sel_trigger");
        break;
    case Preset::Slice:
        filtered = filter_on(node, "sel_slice");
        break;
    case Preset::Fiducial:
        filtered = filter_on(node, "sel_fiducial");
        break;
    case Preset::Muon:
        filtered = filter_on(node, "sel_muon");
        break;
    default:
        break;
    }

    if (selection != NULL)
    {
        selection->nominal.node = filtered;
    }

    return filtered;
}

ROOT::RDF::RNode DefaultAnalysisModel::decorate(ROOT::RDF::RNode node) const
{
    std::vector<std::string> names = node.GetColumnNames();
    auto has = [&](const std::string &name) {
        return std::find(names.begin(), names.end(), name) != names.end();
    };

    auto define_if_missing = [&](const char *name, auto &&f, std::initializer_list<const char *> deps) {
        if (has(name))
            return;
        const std::vector<std::string> columns{deps.begin(), deps.end()};
        node = node.Define(name, std::forward<decltype(f)>(f), columns);
        names.emplace_back(name);
    };

    if (has("beam_mode") && has("run") && has("software_trigger") && has("software_trigger_pre") && has("software_trigger_post"))
    {
        define_if_missing(
            "sel_trigger",
            [](int beam_mode, int run, int sw, int sw_pre, int sw_post) {
                const int numi_beam_mode = static_cast<int>(SampleIO::BeamMode::kNuMI);
                if (beam_mode == numi_beam_mode)
                {
                    constexpr int numi_run_boundary = 16880;
                    if (run < numi_run_boundary)
                        return sw_pre > 0;
                    return sw_post > 0;
                }
                return sw > 0;
            },
            {"beam_mode", "run", "software_trigger", "software_trigger_pre", "software_trigger_post"});
    }
    else
    {
        define_if_missing("sel_trigger", [](int sw) { return sw > 0; }, {"software_trigger"});
    }

    define_if_missing(
        "sel_slice",
        [this](int ns, float topo) {
            return ns == slice_required_count() && topo > slice_min_topology_score();
        },
        {"num_slices", "topological_score"});

    define_if_missing(
        "sel_fiducial",
        [](bool slice, bool fv) { return slice && fv; },
        {"sel_slice", "in_reco_fiducial"});

    define_if_missing(
        "sel_topology",
        [](bool fid) { return fid; },
        {"sel_fiducial"});

    define_if_missing(
        "sel_muon",
        [this](bool topo,
               const ROOT::RVec<float> &scores,
               const ROOT::RVec<float> &lengths,
               const ROOT::RVec<float> &distances,
               const ROOT::RVec<unsigned> &generations) {
            if (!topo)
                return false;

            const auto n = scores.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                if (i >= lengths.size() || i >= distances.size() || i >= generations.size())
                    return false;
                const bool ok = scores[i] > muon_min_track_score() &&
                                lengths[i] > muon_min_track_length() &&
                                distances[i] < muon_max_track_distance() &&
                                generations[i] == muon_required_generation();
                if (ok)
                    return true;
            }
            return false;
        },
        {"sel_topology",
         "track_shower_scores",
         "track_length",
         "track_distance_to_vertex",
         "pfp_generations"});

    return node;
}

std::string DefaultAnalysisModel::selection_label(Preset p) const
{
    switch (p)
    {
    case Preset::Trigger:
        return "Trigger Selection";
    case Preset::Slice:
        return "Slice Selection";
    case Preset::Fiducial:
        return "Fiducial Selection";
    case Preset::Muon:
        return "Muon Selection";
    case Preset::Empty:
    default:
        return "Empty Selection";
    }
}

bool DefaultAnalysisModel::is_in_truth_volume(float x, float y, float z) const noexcept
{
    return is_in_active_volume(x, y, z);
}

bool DefaultAnalysisModel::is_in_reco_volume(float x, float y, float z) const noexcept
{
    return is_in_active_volume(x, y, z) && (z < reco_gap_min_z || z > reco_gap_max_z);
}
