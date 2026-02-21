/* -- C++ -- */
/// \file evd/macro/plotSignalTriggerFailDetectorImages.C
/// \brief Render detector images for signal events that fail trigger selection.

#include <ROOT/RDataFrame.hxx>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AnalysisConfigService.hh"

#include "../../../framework/modules/evd/include/EventDisplay.hh"
#include "include/MacroGuard.hh"
#include "include/MacroEnv.hh"

using namespace heron::evd;

namespace
{


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
        if (heron::macro::file_exists(path))
            return path;
    }

    return "";
}

} // namespace

void plotSignalTriggerFailDetectorImages(const std::string &input_file,
                                         unsigned long long n_events = 25,
                                         const std::string &out_dir = "./plots/evd_detector_signal_trigger_fail",
                                         const std::string &tree_name = "events",
                                         const std::string &signal_expr = "is_signal",
                                         const std::string &trigger_expr = "sel_trigger")
{
  heron::macro::run_with_guard("plotSignalTriggerFailDetectorImages", [&]() {
    ROOT::RDataFrame df(tree_name, input_file);

    const std::string selection = "(" + signal_expr + ") && !(" + trigger_expr + ")";
    auto filtered = df.Filter(selection);

    const auto n_pass = filtered.Count().GetValue();
    if (n_pass == 0)
    {
        std::cerr << "[plotSignalTriggerFailDetectorImages] No matching events found in '"
                  << input_file << "' for selection: " << selection << ".\n";
        return;
    }

    EventDisplay::BatchOptions opt;
    opt.mode = EventDisplay::Mode::Detector;
    opt.out_dir = out_dir;
    opt.n_events = n_events;
    opt.planes = {"U", "V", "W"};

    std::cout << "[plotSignalTriggerFailDetectorImages] matched_events=" << n_pass
              << ", rendering up to " << n_events << " event(s).\n"
              << "[plotSignalTriggerFailDetectorImages] selection=" << selection << "\n";

    EventDisplay::render_from_rdf(filtered, opt);

  });
}

void plotSignalTriggerFailDetectorImages(unsigned long long n_events = 25,
                                         const std::string &out_dir = "./plots/evd_detector_signal_trigger_fail",
                                         const std::string &tree_name = "events",
                                         const std::string &signal_expr = "is_signal",
                                         const std::string &trigger_expr = "sel_trigger")
{
    const auto input_file = find_default_event_list_path();
    if (input_file.empty())
    {
        std::cerr
            << "No default event list file found.\n"
            << "Usage: plotSignalTriggerFailDetectorImages(\"/path/to/event_list.root\""
            << "[, n_events[, out_dir[, tree_name[, signal_expr[, trigger_expr]]]]])\n";
        return;
    }

    std::cout << "[plotSignalTriggerFailDetectorImages] Using default event list: " << input_file << "\n";
    plotSignalTriggerFailDetectorImages(input_file,
                                        n_events,
                                        out_dir,
                                        tree_name,
                                        signal_expr,
                                        trigger_expr);
}
