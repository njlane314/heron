// plot/macro/plotMuonCandidateConfusionMatrix.C
//
// Build a 2x2 confusion matrix for muon-candidate selection:
//   prediction: sel_muon (muon candidate selected or not)
//   truth:      is_nu_mu_cc (inclusive true νμ CC or not)
//
// Default phase-space is after software trigger, neutrino-slice, and reco-fiducial cuts.
//
// Run with:
//   ./nuxsec macro plotMuonCandidateConfusionMatrix.C
//   ./nuxsec macro plotMuonCandidateConfusionMatrix.C \
//     'plotMuonCandidateConfusionMatrix("./scratch/out/event_list_myana.root")'
//   ./nuxsec macro plotMuonCandidateConfusionMatrix.C \
//     'plotMuonCandidateConfusionMatrix("./scratch/out/event_list_myana.root","true","w_nominal",true)'

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

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
            os << std::fixed << std::setprecision(0) << c << "\\n(";
            os << std::setprecision(1) << (100.0 * f) << "%)";
            latex.DrawLatex(x, y, os.str().c_str());
        }
    }
}

} // namespace

int plotMuonCandidateConfusionMatrix(const std::string &input = "",
                                     const std::string &extra_sel = "true",
                                     const std::string &weight_expr = "w_nominal",
                                     bool row_normalise = true)
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

    const std::string selection = "(" + (extra_sel.empty() ? std::string("true") : extra_sel) + ")"
                                  " && software_trigger > 0"
                                  " && sel_slice"
                                  " && in_reco_fiducial";

    auto node = node_mc.Filter(selection)
                    .Define("cm_truth", [](bool is_numu_cc) { return is_numu_cc ? 1 : 0; }, {"is_nu_mu_cc"})
                    .Define("cm_pred", [](bool sel) { return sel ? 1 : 0; }, {"sel_muon"});

    auto h = node.Histo2D(
        {"h_muon_candidate_confusion", ";Predicted muon-candidate label;True inclusive #nu_{#mu}CC label", 2, -0.5, 1.5, 2, -0.5, 1.5},
        "cm_pred", "cm_truth", weight_expr);

    TH2D h_count = *h;
    h_count.SetDirectory(nullptr);
    h_count.GetXaxis()->SetBinLabel(1, "not selected");
    h_count.GetXaxis()->SetBinLabel(2, "selected");
    h_count.GetYaxis()->SetBinLabel(1, "not #nu_{#mu}CC");
    h_count.GetYaxis()->SetBinLabel(2, "#nu_{#mu}CC");

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

    const std::string out_dir = getenv_or("NUXSEC_PLOT_DIR", "./scratch/plots");
    const std::string out_fmt = getenv_or("NUXSEC_PLOT_FORMAT", "pdf");
    gSystem->mkdir(out_dir.c_str(), true);

    TCanvas c("c_muon_candidate_confusion", "Muon candidate confusion matrix", 900, 760);
    c.SetLeftMargin(0.18f);
    c.SetRightMargin(0.15f);
    c.SetBottomMargin(0.15f);
    c.SetTopMargin(0.1f);

    if (row_normalise)
    {
        h_row_frac.SetTitle("Muon-candidate confusion matrix (row-normalised)");
        h_row_frac.GetZaxis()->SetTitle("Fraction per truth row");
        h_row_frac.SetMinimum(0.0);
        h_row_frac.SetMaximum(1.0);
        h_row_frac.Draw("COLZ");
        draw_cell_text(h_count, h_row_frac);
    }
    else
    {
        h_count.SetTitle("Muon-candidate confusion matrix (weighted counts)");
        h_count.GetZaxis()->SetTitle("Weighted events");
        h_count.Draw("COLZ TEXT");
    }

    std::ostringstream ptxt;
    ptxt << std::fixed << std::setprecision(1)
         << "Purity=" << (100.0 * purity) << "%, Eff.=" << (100.0 * efficiency) << "%";

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(42);
    latex.SetTextSize(0.035f);
    latex.DrawLatex(0.18f, 0.96f, ptxt.str().c_str());

    const std::string out_path = out_dir + "/muon_candidate_confusion_matrix." + out_fmt;
    c.SaveAs(out_path.c_str());

    std::cout << "[plotMuonCandidateConfusionMatrix] wrote " << out_path << "\n";
    return 0;
}
