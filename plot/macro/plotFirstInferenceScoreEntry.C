// plot/macro/plotFirstInferenceScoreEntry.C
//
// Stacked distribution of the first entry in the inf_scores vector.
//
// Run with:
//   ./nuxsec macro plotFirstInferenceScoreEntry.C
//   ./nuxsec macro plotFirstInferenceScoreEntry.C \
//     'plotFirstInferenceScoreEntry("./scratch/out/event_list_myana.root")'
//   ./nuxsec macro plotFirstInferenceScoreEntry.C \
//     'plotFirstInferenceScoreEntry("./scratch/out/event_list_myana.root","sel_triggered_muon",false,true)'

#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &path)
{
    const auto n = path.size();
    if (n < 5 || path.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> input_file(TFile::Open(path.c_str(), "READ"));
    if (!input_file || input_file->IsZombie())
        return false;

    const bool has_refs = (input_file->Get("sample_refs") != nullptr);
    const bool has_events_tree = (input_file->Get("events") != nullptr);
    const bool has_event_tree_key = (input_file->Get("event_tree") != nullptr);
    return has_refs && (has_events_tree || has_event_tree_key);
}

bool implicit_mt_enabled()
{
    const char *env = std::getenv("NUXSEC_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

struct RocCurve
{
    std::vector<double> fpr;
    std::vector<double> tpr;
    double auc = 0.0;
};

RocCurve make_roc_curve(const std::vector<float> &scores, const std::vector<int> &labels)
{
    RocCurve out;
    if (scores.size() != labels.size() || scores.empty())
        return out;

    size_t n_pos = 0;
    size_t n_neg = 0;
    std::vector<std::pair<float, int>> ranked;
    ranked.reserve(scores.size());

    for (size_t i = 0; i < scores.size(); ++i)
    {
        const int y = labels[i] ? 1 : 0;
        ranked.emplace_back(scores[i], y);
        if (y)
            ++n_pos;
        else
            ++n_neg;
    }

    if (n_pos == 0 || n_neg == 0)
        return out;

    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<float, int> &a, const std::pair<float, int> &b) {
                  return a.first > b.first;
              });

    double tp = 0.0;
    double fp = 0.0;
    out.fpr.push_back(0.0);
    out.tpr.push_back(0.0);

    size_t i = 0;
    while (i < ranked.size())
    {
        const float thr = ranked[i].first;
        while (i < ranked.size() && ranked[i].first == thr)
        {
            if (ranked[i].second)
                tp += 1.0;
            else
                fp += 1.0;
            ++i;
        }
        out.tpr.push_back(tp / static_cast<double>(n_pos));
        out.fpr.push_back(fp / static_cast<double>(n_neg));
    }

    for (size_t k = 1; k < out.fpr.size(); ++k)
    {
        const double dx = out.fpr[k] - out.fpr[k - 1];
        const double yavg = 0.5 * (out.tpr[k] + out.tpr[k - 1]);
        out.auc += dx * yavg;
    }

    return out;
}
}

int plotFirstInferenceScoreEntry(const std::string &samples_tsv = "",
                                 const std::string &extra_sel = "sel_triggered_muon",
                                 bool use_logy = true,
                                 bool include_data = false)
{
    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotFirstInferenceScoreEntry] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotFirstInferenceScoreEntry] input is not an event list ROOT file: "
                  << list_path << "\n";
        return 2;
    }

    if (implicit_mt_enabled())
        ROOT::EnableImplicitMT();

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    auto filter_by_mask = [](ROOT::RDF::RNode node, std::shared_ptr<const std::vector<char>> mask) {
        return node.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<std::size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf.Define("inf_score_0",
                                       [](const ROOT::RVec<float> &scores) {
                                           if (scores.empty())
                                               return -1.0f;
                                           return scores[0];
                                       },
                                       {"inf_scores"})
                               .Define("inf_score_0_sigmoid",
                                       [](float x) {
                                           return 1.0f / (1.0f + std::exp(-x));
                                       },
                                       {"inf_score_0"});

    if (rdf.HasColumn("is_signal"))
    {
        base = base.Define("is_signal_label", [](bool is_signal) { return is_signal ? 1 : 0; }, {"is_signal"});
    }
    else
    {
        base = base.Define("is_signal_label",
                           [](int analysis_channel) { return analysis_channel == 15 ? 1 : 0; },
                           {"analysis_channels"});
    }

    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<std::size_t>(sid)]);
                                   },
                                   {"sample_id"});
    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

    std::vector<Entry> entries;
    entries.reserve(include_data ? 3 : 2);

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;

    ProcessorEntry rec_mc;
    rec_mc.source = Type::kMC;
    entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
    Entry &e_mc = entries.back();
    mc.push_back(&e_mc);

    ProcessorEntry rec_ext;
    rec_ext.source = Type::kExt;
    entries.emplace_back(make_entry(std::move(node_ext), rec_ext));
    Entry &e_ext = entries.back();
    mc.push_back(&e_ext);

    Entry *p_data = nullptr;
    if (include_data)
    {
        ProcessorEntry rec_data;
        rec_data.source = Type::kData;
        entries.emplace_back(make_entry(std::move(node_data), rec_data));
        p_data = &entries.back();
        data.push_back(p_data);
    }

    if (!extra_sel.empty())
    {
        e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(extra_sel);
        e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(extra_sel);
        if (p_data != nullptr)
            p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(extra_sel);
    }

    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.overlay_signal = true;
    opt.show_ratio = include_data;
    opt.show_ratio_band = include_data;
    // Use the standard signal-channel definition (#nu_{#mu}CC single+multi-strange).
    opt.signal_channels = Channels::signal_keys();
    opt.y_title = "Events";
    opt.run_numbers = {"1"};
    opt.image_format = "pdf";

    const double pot_data = el.total_pot_data();
    const double pot_mc = el.total_pot_mc();
    opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
    opt.beamline = el.beamline_label();

    TH1DModel spec = make_spec("inf_score_0", 50, -10.0, 10.0, "w_nominal");
    spec.sel = Preset::Empty;

    opt.x_title = "Inference score [0]";

    if (include_data)
        plotter.draw_stack(spec, mc, data);
    else
        plotter.draw_stack(spec, mc);

    TH1DModel spec_sigmoid = make_spec("inf_score_0_sigmoid", 50, 0.0, 1.0, "w_nominal");
    spec_sigmoid.sel = Preset::Empty;
    opt.x_title = "Sigmoid(inference score [0])";

    if (include_data)
        plotter.draw_stack(spec_sigmoid, mc, data);
    else
        plotter.draw_stack(spec_sigmoid, mc);

    ROOT::RDF::RNode auc_node = filter_by_mask(base, mask_mc);
    if (!extra_sel.empty())
        auc_node = auc_node.Filter(extra_sel);

    auto v_logit = auc_node.Take<float>("inf_score_0");
    auto v_sigmoid = auc_node.Take<float>("inf_score_0_sigmoid");
    auto v_label = auc_node.Take<int>("is_signal_label");

    const RocCurve roc_logit = make_roc_curve(*v_logit, *v_label);
    const RocCurve roc_sigmoid = make_roc_curve(*v_sigmoid, *v_label);

    if (!roc_logit.fpr.empty() && !roc_sigmoid.fpr.empty())
    {
        gStyle->SetOptStat(0);
        TCanvas c("c_infscore_auc", "c_infscore_auc", 800, 700);

        TGraph g_logit(static_cast<int>(roc_logit.fpr.size()), roc_logit.fpr.data(), roc_logit.tpr.data());
        g_logit.SetTitle(";False positive rate;True positive rate");
        g_logit.SetLineColor(kBlue + 1);
        g_logit.SetLineWidth(3);
        g_logit.Draw("AL");

        TGraph g_sigmoid(static_cast<int>(roc_sigmoid.fpr.size()), roc_sigmoid.fpr.data(), roc_sigmoid.tpr.data());
        g_sigmoid.SetLineColor(kRed + 1);
        g_sigmoid.SetLineStyle(2);
        g_sigmoid.SetLineWidth(3);
        g_sigmoid.Draw("L SAME");

        TLine diagonal(0.0, 0.0, 1.0, 1.0);
        diagonal.SetLineStyle(3);
        diagonal.SetLineColor(kGray + 2);
        diagonal.Draw("SAME");

        TLegend leg(0.16, 0.70, 0.60, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        char label_logit[64];
        char label_sigmoid[64];
        std::snprintf(label_logit, sizeof(label_logit), "Raw logit (AUC = %.4f)", roc_logit.auc);
        std::snprintf(label_sigmoid, sizeof(label_sigmoid), "Sigmoid(logit) (AUC = %.4f)", roc_sigmoid.auc);
        leg.AddEntry(&g_logit, label_logit, "l");
        leg.AddEntry(&g_sigmoid, label_sigmoid, "l");
        leg.Draw();

        const std::string out_dir = env_or("NUXSEC_PLOT_OUT_DIR", "./scratch/out");
        const std::string out_fmt = env_or("NUXSEC_PLOT_OUT_FMT", "pdf");
        gSystem->mkdir(out_dir.c_str(), true);
        const std::string out_path = out_dir + "/first_inference_score_roc_auc." + out_fmt;
        c.SaveAs(out_path.c_str());
        std::cout << "[plotFirstInferenceScoreEntry] wrote " << out_path << "\n";
    }
    else
    {
        std::cerr << "[plotFirstInferenceScoreEntry] unable to compute ROC/AUC (missing signal/background entries).\n";
    }

    std::cout << "[plotFirstInferenceScoreEntry] done\n";
    return 0;
}
