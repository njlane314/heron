#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "SampleCLI.hh"
#include "SelectionService.hh"
#include "EventListIO.hh"
#include "PlotEnv.hh"
#include "PlottingHelper.hh"

using namespace nu;

namespace {

std::shared_ptr<const std::vector<char>> build_truth_mc_mask(const EventListIO& event_list) {
    const auto& refs = event_list.sample_refs();
    int max_sample_id = -1;
    for (const auto& kv : refs) {
        if (kv.first > max_sample_id) {
            max_sample_id = kv.first;
        }
    }

    auto mask = std::make_shared<std::vector<char>>(static_cast<size_t>(max_sample_id + 1), 0);

    const int overlay = static_cast<int>(SampleIO::SampleOrigin::kOverlay);
    const int dirt = static_cast<int>(SampleIO::SampleOrigin::kDirt);
    const int strangeness = static_cast<int>(SampleIO::SampleOrigin::kStrangeness);

    for (const auto& kv : refs) {
        const int sid = kv.first;
        const int origin = kv.second.sample_origin;

        if (sid < 0 || sid > max_sample_id) {
            continue;
        }

        (*mask)[static_cast<size_t>(sid)] =
            (origin == overlay || origin == dirt || origin == strangeness) ? 1 : 0;
    }

    return mask;
}

double asimov_significance(const double signal, const double background) {
    if (signal <= 0.0 || background <= 0.0) {
        return 0.0;
    }

    const double ratio = signal / background;
    const double value = 2.0 * ((signal + background) * std::log(1.0 + ratio) - signal);

    return value > 0.0 ? std::sqrt(value) : 0.0;
}

}

int plot_asimov_significance_distribution(
    const std::string& variable,
    const int nbins = 20,
    const double xmin = 0.0,
    const double xmax = 5.0,
    const std::string& event_list_path = "",
    const std::string& signal_sel = "is_signal",
    const std::string& stage_sel = "sel_trigger && sel_slice && sel_fiducial && sel_muon",
    const std::string& mc_weight = "w_nominal",
    const std::string& output_stem = "asimov_significance_distribution") {
    ROOT::EnableImplicitMT();

    const std::string variable_expr = variable.empty() ? "inf_score_0" : variable;

    const std::string input_path = event_list_path.empty() ? default_event_list_root() : event_list_path;
    if (!looks_like_event_list_root(input_path)) {
        std::cerr << "[plot_asimov_significance_distribution] input is not an event-list root file: " << input_path << "\n";
        return 1;
    }

    EventListIO event_list(input_path);
    ROOT::RDF::RNode rdf = SelectionService::decorate(event_list.rdf())
                              .Define(
                                  "inf_score_0",
                                  [](const ROOT::RVec<float>& scores) {
                                      if (scores.empty()) {
                                          return -1.0f;
                                      }
                                      return scores[0];
                                  },
                                  {"inf_scores"});

    auto mask_mc = build_truth_mc_mask(event_list);
    auto mask_ext = event_list.mask_for_ext();

    auto filter_by_mask = [](ROOT::RDF::RNode node, std::shared_ptr<const std::vector<char>> mask) {
        return node.Filter(
            [mask](int sid) {
                return sid >= 0
                    && sid < static_cast<int>(mask->size())
                    && (*mask)[static_cast<size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode node_mc = filter_by_mask(rdf, mask_mc).Define("__w__", mc_weight);
    ROOT::RDF::RNode node_ext = filter_by_mask(rdf, mask_ext).Define("__w__", mc_weight);

    const std::string selected = "(" + stage_sel + ")";
    const std::string signal_selected = "(" + signal_sel + ") && " + selected;
    const std::string background_selected = "!(" + signal_sel + ") && " + selected;

    const auto model = ROOT::RDF::TH1DModel("h_tmp", ";;", nbins, xmin, xmax);

    auto h_signal_mc = node_mc.Filter(signal_selected).Histo1D(model, variable_expr, "__w__");
    auto h_background_mc = node_mc.Filter(background_selected).Histo1D(model, variable_expr, "__w__");
    auto h_background_ext = node_ext.Filter(selected).Histo1D(model, variable_expr, "__w__");

    TH1D h_asimov("h_asimov", ("Asimov significance;" + variable_expr + ";Z_{A} per bin").c_str(), nbins, xmin, xmax);

    std::cout << "\n[plot_asimov_significance_distribution] Variable: " << variable_expr << "\n";
    std::cout << "bin\txlow\txhigh\tS\tB\tZ_A\n";

    for (int i = 1; i <= nbins; ++i) {
        const double s = h_signal_mc->GetBinContent(i);
        const double b = h_background_mc->GetBinContent(i) + h_background_ext->GetBinContent(i);
        const double z = asimov_significance(s, b);

        h_asimov.SetBinContent(i, z);

        std::cout << i << "\t"
                    << h_asimov.GetXaxis()->GetBinLowEdge(i) << "\t"
                    << h_asimov.GetXaxis()->GetBinUpEdge(i) << "\t"
                    << s << "\t" << b << "\t" << z << "\n";
    }

    TCanvas c("c_asimov_significance", "Asimov significance distribution", 1000, 700);
    gStyle->SetOptStat(0);
    c.SetBottomMargin(0.13);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.07);

    h_asimov.SetLineColor(kBlue + 1);
    h_asimov.SetMarkerColor(kBlue + 1);
    h_asimov.SetLineWidth(2);
    h_asimov.SetMarkerStyle(20);
    h_asimov.Draw("HIST P");

    c.RedrawAxis();

    const auto out = plot_output_file(output_stem).string();
    c.SaveAs(out.c_str());

    std::cout << "\n[plot_asimov_significance_distribution] saved plot: " << out << "\n";

    return 0;
}

int plot_asimov_significance() {
    return plot_asimov_significance_distribution("inf_score_0");
}
