// plotLambdaVisibility.C
// Usage: root -l -q plotLambdaVisibility.C

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "../include/PlotEnv.hh"
#include "../include/Plotter.hh"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace phys { // PDG masses [GeV]
  constexpr double mL  = 1.115683;
  constexpr double mp  = 0.9382720813;
  constexpr double mpi = 0.13957039;
}

// --- tiny configuration block (edit here) ---
namespace cfg {
  const char* file     = "mc.root";
  const char* tree     = "events";
  const char* br_pL    = "lambda_p_exit"; // Λ momentum at nuclear exit [GeV/c]
  const char* br_isSig = "isExitSignal";  // bool/int: exit-defined Λ signal?
  const char* br_isSel = "isSelected";    // bool/int: passes full selection?
  constexpr double pmin_p  = 0.25;        // proton threshold [GeV/c]
  constexpr double pmin_pi = 0.10;        // pion   threshold [GeV/c]
  constexpr int    nbins = 30;
  constexpr double xlow  = 0.0, xhigh = 3.0; // p_Λ range [GeV/c]
}

void plotLambdaVisibility() {
  nuxsec::plot::Plotter{}.set_global_style();
  gStyle->SetOptStat(0);

  // --------- Analytic two-body constants (Λ→pπ) ----------
  const double Ep_star  = (phys::mL*phys::mL + phys::mp*phys::mp  - phys::mpi*phys::mpi)/(2.0*phys::mL);
  const double Epi_star = (phys::mL*phys::mL + phys::mpi*phys::mpi - phys::mp*phys::mp )/(2.0*phys::mL);
  const double p_star   = std::sqrt(std::max(0.0, Ep_star*Ep_star - phys::mp*phys::mp));
  const double Ep_min   = std::sqrt(phys::mp*phys::mp  + cfg::pmin_p*cfg::pmin_p);
  const double Epi_min  = std::sqrt(phys::mpi*phys::mpi + cfg::pmin_pi*cfg::pmin_pi);
  const double pvis     = 0.42; // Λ visibility threshold [GeV/c]

  // --------- Analytic visibility fraction F_kin(p_Λ) ----------
  const int N = 600;
  auto gF = new TGraph(N);
  gF->SetLineWidth(3);
  for (int i=0; i<N; ++i) {
    const double p = cfg::xlow + (cfg::xhigh - cfg::xlow)*(i + 0.5)/N;
    const double E = std::sqrt(p*p + phys::mL*phys::mL);
    const double gamma = E/phys::mL;
    const double beta  = (E>0.0) ? p/E : 0.0;

    double frac = 0.0;
    if (beta > 0.0) {
      const double mu_min = (Ep_min/gamma - Ep_star)/(beta*p_star); // require E_p ≥ Ep_min
      const double mu_max = (Epi_star     - Epi_min/gamma)/(beta*p_star); // and E_π ≥ Eπ_min
      const double a = std::max(-1.0, mu_min);
      const double b = std::min(+1.0, mu_max);
      if (b > a) frac = 0.5*(b - a); // isotropic in cosθ*
      if (frac < 0.0) frac = 0.0;
      if (frac > 1.0) frac = 1.0;
    }
    gF->SetPoint(i, p, frac);
  }

  // --------- MC efficiency vs p_Λ (exit-defined denominator) ----------
  TEfficiency* eff = nullptr;
  if (gSystem->AccessPathName(cfg::file)) {
    std::cout << "Warning: MC file '" << cfg::file
              << "' not found; skipping MC efficiency overlay.\n";
  } else {
    TFile f(cfg::file, "READ");
    if (f.IsZombie()) {
      std::cout << "Warning: MC file '" << cfg::file
                << "' could not be opened; skipping MC efficiency overlay.\n";
    } else {
      TTree* T = (TTree*) f.Get(cfg::tree);
      if (!T) {
        std::cout << "Warning: tree '" << cfg::tree
                  << "' not found in '" << cfg::file
                  << "'; skipping MC efficiency overlay.\n";
      } else if (!T->GetBranch(cfg::br_pL)
                 || !T->GetBranch(cfg::br_isSig)
                 || !T->GetBranch(cfg::br_isSel)) {
        std::cout << "Warning: required branches not found in '" << cfg::file
                  << "'; skipping MC efficiency overlay.\n";
      } else {
        double pL = 0.0; int isSig=0, isSel=0;
        T->SetBranchAddress(cfg::br_pL,    &pL);
        T->SetBranchAddress(cfg::br_isSig, &isSig);
        T->SetBranchAddress(cfg::br_isSel, &isSel);

        eff = new TEfficiency("eff",";p_{#Lambda} at Nuclear Exit [GeV/c];#varepsilon_{fid} (MC)",
                              cfg::nbins, cfg::xlow, cfg::xhigh);
        eff->SetStatisticOption(TEfficiency::kFCP);

        const Long64_t n = T->GetEntries();
        for (Long64_t i=0;i<n;++i) {
          T->GetEntry(i);
          if (!isSig) continue;                // denominator: exit-defined Λ (any decay)
          eff->Fill(isSel, pL);                // numerator: selected events
        }
      }
    }
  }

  // --------- Draw ----------
  TCanvas c("c","Λ visibility & MC efficiency", 900, 650);
  auto* frame = c.DrawFrame(cfg::xlow, 0.0, cfg::xhigh, 1.05,
                            ";p_{#Lambda} at Nuclear Exit [GeV/c];Kinematic Visibility Efficiency");
  frame->GetYaxis()->SetTitleOffset(1.15);

  gF->Draw("L");
  if (eff) eff->Draw("pe same");

  TLine L(pvis, 0.0, pvis, 1.05); L.SetLineStyle(2); L.SetLineWidth(2); L.Draw("same");

  TLegend leg(0.50, 0.18, 0.88, 0.38); leg.SetBorderSize(0);
  if (eff) leg.AddEntry(eff, "#varepsilon_{fid}(p_{#Lambda}) (MC)", "pe");
  leg.AddEntry(gF,  "F_{kin}(p_{#Lambda}) (two-body)", "l");
  leg.AddEntry(&L,  Form("p^{vis}_{#Lambda} = %.2f GeV/c", pvis), "l");
  leg.Draw();

  const std::string out = nuxsec::plot::plot_output_file("lambda_visibility_efficiency").string();
  c.SaveAs(out.c_str());
}
