// macros/inspect_event_weight_branches.C
//
// Print which event-weight branches are present in the event-list file and,
// when available, inspect the first MC entry to show vector sizes and map keys.

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include "EventListIO.hh"
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

void print_first_vector_size(ROOT::RDF::RNode node, const std::string &branch)
{
    if (!has_column(node, branch))
        return;

    const std::string size_col = "__size_" + branch;
    auto with_size = node.Define(size_col, "static_cast<int>(" + branch + ".size())");

    auto vals = with_size.Range(1).Take<int>(size_col);
    if (!vals || vals->empty())
    {
        std::cout << "  first MC entry: " << branch << " unavailable\n";
        return;
    }

    std::cout << "  first MC entry: " << branch << ".size() = "
              << (*vals)[0] << "\n";
}

} // namespace

int inspect_event_weight_branches(const std::string &event_list_path = "")
{
    return heron::macro::run_with_guard("inspect_event_weight_branches", [&]() -> int {
        const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
        std::cout << "[inspect_event_weight_branches] input=" << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[inspect_event_weight_branches] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        EventListIO el(input_path);
        ROOT::RDF::RNode rdf = SelectionService::decorate(el.rdf());

        const std::vector<std::string> interesting = {
            "weights",
            "weightsGenie",
            "weightsGenieUp",
            "weightsGenieDn",
            "weightsFlux",
            "weightsReint",
            "weightsPPFX",
            "weightSpline",
            "weightTune",
            "weightSplineTimesTune",
            "ppfx_cv",
            "RootinoFix",
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
            "knobxsr_scc_Fa3dn"};

        std::cout << "[inspect_event_weight_branches] column presence\n";
        for (const auto &name : interesting)
            std::cout << "  " << name << " : " << (has_column(rdf, name) ? "yes" : "no") << "\n";

        auto mask_ext = el.mask_for_ext();
        auto mask_mc_like = el.mask_for_mc_like();
        ROOT::RDF::RNode node_all = filter_by_sample_mask(rdf, mask_mc_like, "sample_id");
        ROOT::RDF::RNode node_mc = filter_not_sample_mask(node_all, mask_ext, "sample_id");

        std::cout << "[inspect_event_weight_branches] first MC-entry vector sizes\n";
        print_first_vector_size(node_mc, "weightsGenie");
        print_first_vector_size(node_mc, "weightsGenieUp");
        print_first_vector_size(node_mc, "weightsGenieDn");
        print_first_vector_size(node_mc, "weightsFlux");
        print_first_vector_size(node_mc, "weightsReint");
        print_first_vector_size(node_mc, "weightsPPFX");

        if (has_column(node_mc, "weights"))
        {
            auto maps = node_mc.Range(1).Take<std::map<std::string, std::vector<double>>>("weights");
            if (maps && !maps->empty())
            {
                const auto &wmap = (*maps)[0];
                std::vector<std::string> keys;
                keys.reserve(wmap.size());
                for (const auto &kv : wmap)
                    keys.push_back(kv.first);
                std::sort(keys.begin(), keys.end());

                std::cout << "[inspect_event_weight_branches] first MC-entry weights-map keys\n";
                for (const auto &k : keys)
                {
                    auto it = wmap.find(k);
                    std::cout << "  " << k << " : size=" << it->second.size() << "\n";
                }
            }
            else
            {
                std::cout << "[inspect_event_weight_branches] weights map exists but first MC entry was empty\n";
            }
        }

        std::cout << "[inspect_event_weight_branches] done\n";
        return 0;
    });
}
