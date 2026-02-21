/* -- C++ -- */
/// \file evd/macro/plotSelMuonInferenceThresholdEventDisplays.C
/// \brief Render detector and semantic event displays for events passing sel_muon and inf_scores[0] > 6.9, split by signal/background.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AnalysisConfigService.hh"

#include "../../../framework/evd/include/EventDisplay.hh"
#include "MacroGuard.hh"

using namespace heron::evd;

namespace
{

bool file_exists(const std::string &path)
{
    std::ifstream f(path.c_str());
    return f.good();
}

std::string find_default_event_list_path()
{
    const auto &analysis = AnalysisConfigService::instance();
    std::ostringstream analysis_path;

    const char *out_base = std::getenv("HERON_OUT_BASE");
    if (out_base && *out_base)
    {
        analysis_path << out_base << "/event_list_" << analysis.name() << ".root";
    }
    else
    {
        const char *user = std::getenv("USER");
        if (user && *user)
        {
            analysis_path << "/exp/uboone/data/users/" << user
                          << "/event_list_" << analysis.name() << ".root";
        }
    }

    const std::vector<std::string> candidates{
        "/exp/uboone/data/users/nlane/heron/out/event/events.root",
        analysis_path.str()};

    for (const auto &path : candidates)
    {
        if (file_exists(path))
            return path;
    }

    return "";
}

void render_category_displays(ROOT::RDF::RNode node,
                              const std::string &label,
                              EventDisplay::Mode mode,
                              unsigned long long n_events,
                              const std::string &out_dir)
{
    const auto n_pass = node.Count().GetValue();
    if (n_pass == 0)
    {
        std::cout << "[plotSelMuonInferenceThresholdEventDisplays] mode="
                  << (mode == EventDisplay::Mode::Detector ? "detector" : "semantic")
                  << ", category=" << label << ": no matching events.\n";
        return;
    }

    EventDisplay::BatchOptions opt;
    opt.mode = mode;
    opt.out_dir = out_dir;
    opt.n_events = n_events;
    opt.planes = {"U", "V", "W"};

    std::cout << "[plotSelMuonInferenceThresholdEventDisplays] mode="
              << (mode == EventDisplay::Mode::Detector ? "detector" : "semantic")
              << ", category=" << label
              << ", matched_events=" << n_pass
              << ", rendering up to " << n_events << " event(s).\n";

    EventDisplay::render_from_rdf(node, opt);
}

} // namespace

void plotSelMuonInferenceThresholdEventDisplays(const std::string &input_file,
                                                unsigned long long n_events = 25,
                                                double threshold = 6.9,
                                                const std::string &out_dir = "./plots/evd_sel_muon_inf69",
                                                const std::string &tree_name = "events",
                                                const std::string &signal_expr = "is_signal")
{
  heron::macro::run_with_guard("plotSelMuonInferenceThresholdEventDisplays", [&]() {
    ROOT::RDataFrame df(tree_name, input_file);

    auto base = df.Define("inf_score_0",
                          [](const ROOT::RVec<float> &scores) {
                              if (scores.empty())
                                  return -1.0f;
                              return scores[0];
                          },
                          {"inf_scores"})
                    .Filter("sel_muon && inf_score_0 > " + std::to_string(threshold));

    const std::string signal_selection = "(" + signal_expr + ")";
    const std::string background_selection = "!(" + signal_expr + ")";

    auto signal = base.Filter(signal_selection);
    auto background = base.Filter(background_selection);

    render_category_displays(signal,
                             "signal",
                             EventDisplay::Mode::Detector,
                             n_events,
                             out_dir + "/detector_signal");
    render_category_displays(background,
                             "background",
                             EventDisplay::Mode::Detector,
                             n_events,
                             out_dir + "/detector_background");

    render_category_displays(signal,
                             "signal",
                             EventDisplay::Mode::Semantic,
                             n_events,
                             out_dir + "/semantic_signal");
    render_category_displays(background,
                             "background",
                             EventDisplay::Mode::Semantic,
                             n_events,
                             out_dir + "/semantic_background");

  });
}

void plotSelMuonInferenceThresholdEventDisplays(unsigned long long n_events = 25,
                                                double threshold = 6.9,
                                                const std::string &out_dir = "./plots/evd_sel_muon_inf69",
                                                const std::string &tree_name = "events",
                                                const std::string &signal_expr = "is_signal")
{
    const auto input_file = find_default_event_list_path();
    if (input_file.empty())
    {
        std::cerr
            << "No default event list file found.\n"
            << "Usage: plotSelMuonInferenceThresholdEventDisplays(\"/path/to/event_list.root\""
            << "[, n_events[, threshold[, out_dir[, tree_name[, signal_expr]]]]])\n";
        return;
    }

    std::cout << "[plotSelMuonInferenceThresholdEventDisplays] Using default event list: " << input_file << "\n";
    plotSelMuonInferenceThresholdEventDisplays(input_file,
                                               n_events,
                                               threshold,
                                               out_dir,
                                               tree_name,
                                               signal_expr);
}
