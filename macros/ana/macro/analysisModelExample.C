/* -- C++ -- */

#include <iostream>

#include "AnalysisModel.hh"

class MuonSelectionModel : public heron::AnalysisModel
{
  public:
    void define() override
    {
        clear();

        const auto p_muon_p = var("muon_p", []() { return 0.0; }, {"reco_muon_p"});
        const auto c_fiducial = cut("fiducial", []() { return true; },
                                    {"reco_vertex_x", "reco_vertex_y", "reco_vertex_z"});
        const auto w_cv = weight("cv", []() { return 1.0; }, {"event_weight_cv"});
        const auto s_nominal = selection("nominal", c_fiducial, w_cv);

        hist1d("h_muon_p", p_muon_p.name, 40, 0.0, 2.0,
               "Muon momentum;p [GeV];Events", s_nominal.name, w_cv.name);
        snapshot("selected_events", {"run", "subrun", "event", p_muon_p.name}, s_nominal.name);
    }
};

void analysisModelExample()
{
    MuonSelectionModel model;
    model.define();

    std::cout << "declared vars: " << model.vars().size() << "\n";
    std::cout << "declared cuts: " << model.cuts().size() << "\n";
    std::cout << "declared weights: " << model.weights().size() << "\n";
    std::cout << "declared selections: " << model.selections().size() << "\n";
    std::cout << "declared h1: " << model.h1().size() << "\n";
    std::cout << "declared snapshots: " << model.snapshots().size() << "\n";

    for (std::vector<heron::Hist1DSpec>::const_iterator it = model.h1().begin(); it != model.h1().end(); ++it)
        std::cout << "hist: " << it->name << " variable=" << it->variable << " selection=" << it->selection << "\n";
}
