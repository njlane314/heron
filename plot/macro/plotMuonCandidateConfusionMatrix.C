// plot/macro/plotMuonCandidateConfusionMatrix.C
//
// Build a 2x2 confusion matrix for muon-candidate selection:
//   prediction: sel_muon (muon candidate selected or not)
//   truth:      is_nu_mu_cc (inclusive true νμ CC or not)
//
// Default phase-space is after software trigger, neutrino-slice, and reco-fiducial cuts.
//
// Run with:
//   ./heron macro plotMuonCandidateConfusionMatrix.C
//   ./heron macro plotMuonCandidateConfusionMatrix.C \
//     'plotMuonCandidateConfusionMatrix("./scratch/out/event_list_myana.root")'
//   ./heron macro plotMuonCandidateConfusionMatrix.C \
//     'plotMuonCandidateConfusionMatrix("./scratch/out/event_list_myana.root","true","w_nominal",true,0.5)'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>

#include "EventListIO.hh"
#include "PlottingHelper.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &p)
{
    const auto n = p.size();
    if (n < 5 || p.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
    if (!f || f->IsZombie())
        return false;

    const bool has_refs = (f->Get("sample_refs") != nullptr);
    const bool has_events_tree = (f->Get("events") != nullptr);
    const bool has_event_tree_key = (f->Get("event_tree") != nullptr);

    return has_refs && (has_events_tree || has_event_tree_key);
}

std::string getenv_or(const char *key, const std::string &fallback)
{
    const char *v = std::getenv(key);
    if (!v || std::string(v).empty())
        return fallback;
    return std::string(v);
}

ROOT::RDF::RNode filter_by_mask(ROOT::RDF::RNode n,
                                std::shared_ptr<const std::vector<char>> mask)
{
    if (!mask)
        return n;

    return n.Filter(
        [mask](int sid) {
            return sid >= 0
                   && sid < static_cast<int>(mask->size())
                   && (*mask)[static_cast<size_t>(sid)];
        },
        {"sample_id"});
}

bool passes_mipness_discriminant(const ROOT::RVec<float> &scores, float threshold)
{
    for (std::size_t i = 0; i < scores.size(); ++i)
    {
        if (scores[i] > threshold)
            return true;
    }
    return false;
}

bool has_column(const std::vector<std::string> &cols, const std::string &name)
{
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

ROOT::RVec<float> derive_track_mipness_median_plane_score(const ROOT::RVec<float> &T70u,
                                                           const ROOT::RVec<float> &Q50u,
                                                           const ROOT::RVec<float> &Q90u,
                                                           const ROOT::RVec<float> &T70v,
                                                           const ROOT::RVec<float> &Q50v,
                                                           const ROOT::RVec<float> &Q90v,
                                                           const ROOT::RVec<float> &T70y,
                                                           const ROOT::RVec<float> &Q50y,
                                                           const ROOT::RVec<float> &Q90y)
{
    constexpr double Mref_mev_per_cm = 1.8;
    constexpr double tail_weight_w = 0.5;

    auto D_from_stats = [](float T70, float Q50, float Q90) -> double {
        if (!(std::isfinite(T70) && std::isfinite(Q50) && std::isfinite(Q90)))
            return std::numeric_limits<double>::quiet_NaN();
        if (!(T70 > 0.0f && Q50 > 0.0f && Q90 > 0.0f))
            return std::numeric_limits<double>::quiet_NaN();
        return std::log(static_cast<double>(T70) / Mref_mev_per_cm) +
               tail_weight_w * std::log(static_cast<double>(Q90) / static_cast<double>(Q50));
    };

    const std::size_t n = std::min({T70u.size(), Q50u.size(), Q90u.size(),
                                    T70v.size(), Q50v.size(), Q90v.size(),
                                    T70y.size(), Q50y.size(), Q90y.size()});

    ROOT::RVec<float> out(n, std::numeric_limits<float>::quiet_NaN());
    for (std::size_t i = 0; i < n; ++i)
    {
        std::vector<double> Ds;
        Ds.reserve(3);

        const double Du = D_from_stats(T70u[i], Q50u[i], Q90u[i]);
        const double Dv = D_from_stats(T70v[i], Q50v[i], Q90v[i]);
        const double Dy = D_from_stats(T70y[i], Q50y[i], Q90y[i]);

        if (std::isfinite(Du))
            Ds.push_back(Du);
        if (std::isfinite(Dv))
            Ds.push_back(Dv);
        if (std::isfinite(Dy))
            Ds.push_back(Dy);

        if (Ds.empty())
            continue;

        std::sort(Ds.begin(), Ds.end());
        const double Dmed = Ds[Ds.size() / 2];
        out[i] = static_cast<float>(std::exp(-Dmed));
    }

    return out;
}

void draw_cell_text(const TH2D &h_count,
                    const TH2D &h_row_frac)
{
    TLatex latex;
    latex.SetTextAlign(22);
    latex.SetTextFont(42);
    latex.SetTextSize(0.028f);
    latex.SetNDC(false);

    for (int by = 1; by <= h_count.GetNbinsY(); ++by)
    {
        for (int bx = 1; bx <= h_count.GetNbinsX(); ++bx)
        {
            const double x = h_count.GetXaxis()->GetBinCenter(bx);
            const double y = h_count.GetYaxis()->GetBinCenter(by);
            const double c = h_count.GetBinContent(bx, by);
            const double f = h_row_frac.GetBinContent(bx, by);

            std::ostringstream os;
            os << "#splitline{" << std::fixed << std::setprecision(0) << c << "}";
            os << "{" << std::setprecision(1) << (100.0 * f) << "%}";
            latex.DrawLatex(x, y, os.str().c_str());
        }
    }
}

} // namespace

int plotMuonCandidateConfusionMatrix(const std::string &input = "",
                                     const std::string &extra_sel = "true",
                                     const std::string &weight_expr = "w_nominal",
                                     bool row_normalise = true,
                                     double mipness_threshold = 0.5)
{
    ROOT::EnableImplicitMT();

    const std::string list_path = input.empty() ? default_event_list_root() : input;
    std::cout << "[plotMuonCandidateConfusionMatrix] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotMuonCandidateConfusionMatrix] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();
    const std::vector<std::string> cols = rdf.GetColumnNames();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();

    ROOT::RDF::RNode node_mc = filter_by_mask(rdf, mask_mc);
    if (mask_ext)
    {
        node_mc = node_mc.Filter(
            [mask_ext](int sid) {
                return !(sid >= 0
                         && sid < static_cast<int>(mask_ext->size())
                         && (*mask_ext)[static_cast<size_t>(sid)]);
            },
            {"sample_id"});
    }

    ROOT::RDF::RNode node_base = node_mc;

    const bool has_mipness_scores = has_column(cols, "track_mipness_median_plane_score");
    const bool has_mipness_stats =
        has_column(cols, "track_dedx_T70_u") && has_column(cols, "track_dedx_Q50_u") && has_column(cols, "track_dedx_Q90_u") &&
        has_column(cols, "track_dedx_T70_v") && has_column(cols, "track_dedx_Q50_v") && has_column(cols, "track_dedx_Q90_v") &&
        has_column(cols, "track_dedx_T70_y") && has_column(cols, "track_dedx_Q50_y") && has_column(cols, "track_dedx_Q90_y");

    if (!has_mipness_scores && has_mipness_stats)
    {
        std::cout << "[plotMuonCandidateConfusionMatrix] deriving track_mipness_median_plane_score from dE/dx quantiles\n";
        node_base = node_base.Define(
            "track_mipness_median_plane_score",
            derive_track_mipness_median_plane_score,
            {"track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
             "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
             "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});
    }
    else if (!has_mipness_scores)
    {
        std::cerr << "[plotMuonCandidateConfusionMatrix] missing required MIPness columns. Need either '\n"
                  << "  - track_mipness_median_plane_score\n"
                  << "or dE/dx quantiles track_dedx_{T70,Q50,Q90}_{u,v,y}\n";
        return 1;
    }

    const std::string selection = "(" + (extra_sel.empty() ? std::string("true") : extra_sel) + ")"
                                  " && software_trigger > 0"
                                  " && sel_slice"
                                  " && in_reco_fiducial";

    auto node = node_base.Filter(selection)
                    .Define("cm_truth", [](bool is_numu_cc) { return is_numu_cc ? 1 : 0; }, {"is_nu_mu_cc"})
                    .Define("cm_pred", [](bool sel) { return sel ? 1 : 0; }, {"sel_muon"})
                    .Define("cm_pred_mipness", [mipness_threshold](const ROOT::RVec<float> &mipness_scores) {
                        return passes_mipness_discriminant(mipness_scores, static_cast<float>(mipness_threshold)) ? 1 : 0;
                    }, {"track_mipness_median_plane_score"});

    auto h = node.Histo2D(
        {"h_muon_candidate_confusion", ";Predicted label;True label", 2, -0.5, 1.5, 2, -0.5, 1.5},
        "cm_pred", "cm_truth", weight_expr);

    auto h_mipness = node.Histo2D(
        {"h_muon_candidate_confusion_mipness", ";Predicted label;True label", 2, -0.5, 1.5, 2, -0.5, 1.5},
        "cm_pred_mipness", "cm_truth", weight_expr);

    TH2D h_count = *h;
    h_count.SetDirectory(nullptr);
    h_count.GetXaxis()->SetBinLabel(1, "No muon candidate");
    h_count.GetXaxis()->SetBinLabel(2, "Muon candidate");
    h_count.GetYaxis()->SetBinLabel(1, "True not #nu_{#mu} CC");
    h_count.GetYaxis()->SetBinLabel(2, "True #nu_{#mu} CC");

    TH2D h_row_frac = h_count;
    for (int by = 1; by <= h_row_frac.GetNbinsY(); ++by)
    {
        const double row_sum = h_row_frac.GetBinContent(1, by) + h_row_frac.GetBinContent(2, by);
        if (row_sum <= 0.0)
            continue;
        for (int bx = 1; bx <= h_row_frac.GetNbinsX(); ++bx)
            h_row_frac.SetBinContent(bx, by, h_row_frac.GetBinContent(bx, by) / row_sum);
    }

    const double tn = h_count.GetBinContent(1, 1);
    const double fn = h_count.GetBinContent(1, 2);
    const double fp = h_count.GetBinContent(2, 1);
    const double tp = h_count.GetBinContent(2, 2);
    const double total = tn + fn + fp + tp;

    const double accuracy = total > 0.0 ? (tp + tn) / total : 0.0;
    const double efficiency = (tp + fn) > 0.0 ? tp / (tp + fn) : 0.0;
    const double purity = (tp + fp) > 0.0 ? tp / (tp + fp) : 0.0;

    std::cout << "[plotMuonCandidateConfusionMatrix] selection=" << selection << "\n";
    std::cout << "[plotMuonCandidateConfusionMatrix] TN=" << tn
              << ", FP=" << fp
              << ", FN=" << fn
              << ", TP=" << tp << "\n";
    std::cout << std::fixed << std::setprecision(4)
              << "[plotMuonCandidateConfusionMatrix] accuracy=" << accuracy
              << ", efficiency=" << efficiency
              << ", purity=" << purity << "\n";

    if (gROOT)
        gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".1f");

    const std::string out_dir = getenv_or("HERON_PLOT_DIR", "./scratch/plots");
    const std::string out_fmt = getenv_or("HERON_PLOT_FORMAT", "pdf");
    gSystem->mkdir(out_dir.c_str(), true);

    TCanvas c("c_muon_candidate_confusion", "Muon candidate confusion matrix", 900, 760);
    c.SetLeftMargin(0.18f);
    c.SetRightMargin(0.15f);
    c.SetBottomMargin(0.15f);
    c.SetTopMargin(0.1f);

    if (row_normalise)
    {
        h_row_frac.SetTitle("");
        h_row_frac.GetZaxis()->SetTitle("Fraction per truth row");
        h_row_frac.SetMinimum(0.0);
        h_row_frac.SetMaximum(1.0);
        h_row_frac.Draw("COLZ");
        draw_cell_text(h_count, h_row_frac);
    }
    else
    {
        h_count.SetTitle("");
        h_count.GetZaxis()->SetTitle("Weighted events");
        h_count.Draw("COLZ TEXT");
    }

    std::ostringstream ptxt;
    ptxt << std::fixed << std::setprecision(1)
         << "Purity=" << (100.0 * purity) << "%, Efficiency=" << (100.0 * efficiency) << "%";

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize(0.03f);
    latex.DrawLatex(0.18f, 0.93f, ptxt.str().c_str());

    const std::string out_path = out_dir + "/muon_candidate_confusion_matrix." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plotMuonCandidateConfusionMatrix] wrote " << out_path << "\n";

    TH2D h_count_mipness = *h_mipness;
    h_count_mipness.SetDirectory(nullptr);
    h_count_mipness.GetXaxis()->SetBinLabel(1, "No muon candidate");
    h_count_mipness.GetXaxis()->SetBinLabel(2, "Muon candidate");
    h_count_mipness.GetYaxis()->SetBinLabel(1, "True not #nu_{#mu} CC");
    h_count_mipness.GetYaxis()->SetBinLabel(2, "True #nu_{#mu} CC");

    TH2D h_row_frac_mipness = h_count_mipness;
    for (int by = 1; by <= h_row_frac_mipness.GetNbinsY(); ++by)
    {
        const double row_sum = h_row_frac_mipness.GetBinContent(1, by) + h_row_frac_mipness.GetBinContent(2, by);
        if (row_sum <= 0.0)
            continue;
        for (int bx = 1; bx <= h_row_frac_mipness.GetNbinsX(); ++bx)
            h_row_frac_mipness.SetBinContent(bx, by, h_row_frac_mipness.GetBinContent(bx, by) / row_sum);
    }

    const double tn_mipness = h_count_mipness.GetBinContent(1, 1);
    const double fn_mipness = h_count_mipness.GetBinContent(1, 2);
    const double fp_mipness = h_count_mipness.GetBinContent(2, 1);
    const double tp_mipness = h_count_mipness.GetBinContent(2, 2);
    const double efficiency_mipness = (tp_mipness + fn_mipness) > 0.0 ? tp_mipness / (tp_mipness + fn_mipness) : 0.0;
    const double purity_mipness = (tp_mipness + fp_mipness) > 0.0 ? tp_mipness / (tp_mipness + fp_mipness) : 0.0;

    std::cout << "[plotMuonCandidateConfusionMatrix] mipness_threshold=" << mipness_threshold
              << ", TN=" << tn_mipness
              << ", FP=" << fp_mipness
              << ", FN=" << fn_mipness
              << ", TP=" << tp_mipness << "\n";
    std::cout << std::fixed << std::setprecision(4)
              << "[plotMuonCandidateConfusionMatrix] mipness_purity=" << purity_mipness
              << ", mipness_efficiency=" << efficiency_mipness << "\n";

    TCanvas c_mipness("c_muon_candidate_confusion_mipness", "Muon candidate confusion matrix (mipness)", 900, 760);
    c_mipness.SetLeftMargin(0.18f);
    c_mipness.SetRightMargin(0.15f);
    c_mipness.SetBottomMargin(0.15f);
    c_mipness.SetTopMargin(0.1f);

    if (row_normalise)
    {
        h_row_frac_mipness.SetTitle("");
        h_row_frac_mipness.GetZaxis()->SetTitle("Fraction per truth row");
        h_row_frac_mipness.SetMinimum(0.0);
        h_row_frac_mipness.SetMaximum(1.0);
        h_row_frac_mipness.Draw("COLZ");
        draw_cell_text(h_count_mipness, h_row_frac_mipness);
    }
    else
    {
        h_count_mipness.SetTitle("");
        h_count_mipness.GetZaxis()->SetTitle("Weighted events");
        h_count_mipness.Draw("COLZ TEXT");
    }

    std::ostringstream ptxt_mipness;
    ptxt_mipness << std::fixed << std::setprecision(1)
                 << "Purity=" << (100.0 * purity_mipness) << "%, Efficiency=" << (100.0 * efficiency_mipness) << "%"
                 << ", MIPness > " << mipness_threshold;

    TLatex latex_mipness;
    latex_mipness.SetNDC(true);
    latex_mipness.SetTextFont(42);
    latex_mipness.SetTextSize(0.03f);
    latex_mipness.DrawLatex(0.18f, 0.93f, ptxt_mipness.str().c_str());

    const std::string out_path_mipness = out_dir + "/muon_candidate_confusion_matrix_mipness." + out_fmt;
    c_mipness.SaveAs(out_path_mipness.c_str());
    std::cout << "[plotMuonCandidateConfusionMatrix] wrote " << out_path_mipness << "\n";
    return 0;
}
