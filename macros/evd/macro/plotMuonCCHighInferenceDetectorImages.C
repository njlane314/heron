/* -- C++ -- */
/// \file evd/macro/plotMuonCCHighInferenceDetectorImages.C
/// \brief Render detector images for events passing inclusive \f$\nu_\mu\f$ CC with inf_scores[0] > 8.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AnalysisConfigService.hh"

#include "../../../mod/evd/include/EventDisplay.hh"

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

} // namespace

void plotMuonCCHighInferenceDetectorImages(const std::string &input_file,
                                           unsigned long long n_events = 25,
                                           const std::string &out_dir = "./plots/evd_detector_mucc_inf8",
                                           const std::string &tree_name = "events")
{
    ROOT::RDataFrame df(tree_name, input_file);

    auto filtered = df.Define("inf_score_0",
                              [](const ROOT::RVec<float> &scores) {
                                  if (scores.empty())
                                      return -1.0f;
                                  return scores[0];
                              },
                              {"inf_scores"})
                        .Filter("sel_triggered_muon && inf_score_0 > 8.0");

    const auto n_pass = filtered.Count().GetValue();
    if (n_pass == 0)
    {
        std::cerr << "[plotMuonCCHighInferenceDetectorImages] No matching events found in '"
                  << input_file << "'.\n";
        return;
    }

    EventDisplay::BatchOptions opt;
    opt.mode = EventDisplay::Mode::Detector;
    opt.out_dir = out_dir;
    opt.n_events = n_events;
    opt.planes = {"U", "V", "W"};

    std::cout << "[plotMuonCCHighInferenceDetectorImages] matched_events=" << n_pass
              << ", rendering up to " << n_events << " event(s).\n";

    EventDisplay::render_from_rdf(filtered, opt);
}

void plotMuonCCHighInferenceDetectorImages(unsigned long long n_events = 25,
                                           const std::string &out_dir = "./plots/evd_detector_mucc_inf8",
                                           const std::string &tree_name = "events")
{
    const auto input_file = find_default_event_list_path();
    if (input_file.empty())
    {
        std::cerr
            << "No default event list file found.\n"
            << "Usage: plotMuonCCHighInferenceDetectorImages(\"/path/to/event_list.root\""
            << "[, n_events[, out_dir[, tree_name]]])\n";
        return;
    }

    std::cout << "[plotMuonCCHighInferenceDetectorImages] Using default event list: " << input_file << "\n";
    plotMuonCCHighInferenceDetectorImages(input_file, n_events, out_dir, tree_name);
}
