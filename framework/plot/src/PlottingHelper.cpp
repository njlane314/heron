/* -- C++ -- */
/**
 *  @file  framework/plot/src/PlottingHelper.cpp
 *
 *  @brief Plotting helper utilities and plotStackedHist implementation.
 */

#include "PlottingHelper.hh"

#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TSystem.h>

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
    const std::string set_name = env_or("HERON_SET", "out");
    const std::string out_base = env_or("HERON_OUT_BASE", repo_root + "/scratch/out");
    return out_base + "/" + set_name + "/sample/samples.tsv";
}

std::string default_event_list_root()
{
    return "/exp/uboone/data/users/nlane/heron/out/event/events.root";
}

bool looks_like_event_list_root(const std::string &path)
{
    const auto n = path.size();
    if (n < 5 || path.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie())
        return false;

    const bool has_refs = (f->Get("sample_refs") != nullptr);
    const bool has_events_tree = (f->Get("events") != nullptr);
    const bool has_event_tree_key = (f->Get("event_tree") != nullptr);

    return has_refs && (has_events_tree || has_event_tree_key);
}

ROOT::RDF::RNode filter_by_sample_mask(ROOT::RDF::RNode node,
                                       std::shared_ptr<const std::vector<char>> mask,
                                       const std::string &sample_id_column)
{
    if (!mask)
        return node;

    return node.Filter(
        [mask](int sid) {
            return sid >= 0
                   && sid < static_cast<int>(mask->size())
                   && (*mask)[static_cast<size_t>(sid)];
        },
        {sample_id_column});
}

ROOT::RDF::RNode filter_not_sample_mask(ROOT::RDF::RNode node,
                                        std::shared_ptr<const std::vector<char>> mask,
                                        const std::string &sample_id_column)
{
    if (!mask)
        return node;

    return node.Filter(
        [mask](int sid) {
            return !(sid >= 0
                     && sid < static_cast<int>(mask->size())
                     && (*mask)[static_cast<size_t>(sid)]);
        },
        {sample_id_column});
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
