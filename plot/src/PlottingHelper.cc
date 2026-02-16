/* -- C++ -- */
/**
 *  @file  plot/src/PlottingHelper.cc
 *
 *  @brief Plotting helper utilities and plotStackedHist implementation.
 */

#include "PlottingHelper.hh"

#include <string>

#include <TSystem.h>

#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"


namespace nu
{

std::string env_or(const char *key, const std::string &fallback)
{
    const char *v = gSystem->Getenv(key);
    if (v && *v)
    {
        return std::string(v);
    }
    return fallback;
}

std::string default_samples_tsv()
{
    const std::string repo_root = env_or("HERON_REPO_ROOT", ".");
    const std::string set_name = env_or("HERON_SET", "template");
    const std::string out_base = env_or("HERON_OUT_BASE", repo_root + "/scratch/out");
    return out_base + "/" + set_name + "/sample/samples.tsv";
}

std::string default_event_list_root()
{
    const std::string user = env_or("USER", "");
    const std::string out_base = env_or("HERON_OUT_BASE", "/exp/uboone/data/users/" + user);
    const auto &analysis = AnalysisConfigService::instance();
    return out_base + "/event_list_" + analysis.name() + ".root";
}

bool is_data_origin(SampleIO::SampleOrigin o)
{
    return o == SampleIO::SampleOrigin::kData;
}

double pick_pot_nom(const SampleIO::Sample &s)
{
    if (s.normalised_pot_sum > 0.0)
        return s.normalised_pot_sum;
    if (s.subrun_pot_sum > 0.0)
        return s.subrun_pot_sum;
    if (s.db_tortgt_pot_sum > 0.0)
        return s.db_tortgt_pot_sum;
    return 0.0;
}

Entry make_entry(ROOT::RDF::RNode node, const ProcessorEntry &proc_entry)
{
    SelectionEntry selection{proc_entry.source, Frame{std::move(node)}};
    return Entry{std::move(selection), 0.0, 0.0, std::string(), std::string()};
}

TH1DModel make_spec(const std::string &expr,
                    int nbins,
                    double xmin,
                    double xmax,
                    const std::string &weight)
{
    TH1DModel spec;
    spec.id = "h_" + expr;
    spec.name = expr;
    spec.expr = expr;
    spec.weight = weight;
    spec.nbins = nbins;
    spec.xmin = xmin;
    spec.xmax = xmax;
    return spec;
}

} // namespace nu
