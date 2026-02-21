// plotSigma0LambdaVisibility.C
// Usage: root -l -q plotSigma0LambdaVisibility.C
//
// Computes kinematic visibility for the cascade Σ0→Λγ.
// Now draws TWO analytic curves vs p_Σ:
//   1) F_Λ(p_Σ): fraction with a visible Λ→pπ (γ not required)
//   2) F_{Λ+γ}(p_Σ): fraction with a visible Λ→pπ AND a visible γ (E_γ ≥ Emin_gamma)

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TSystem.h"
#include "../include/PlotEnv.hh"
#include "../include/Plotter.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "MacroGuard.hh"

using namespace nu;

  constexpr double mS0 = 1.192642;     // Σ0
  constexpr double mL  = 1.115683;     // Λ
  constexpr double mp  = 0.9382720813; // p
  constexpr double mpi = 0.13957039;   // π
}

// --- configuration (edit here) ---
  const char* file     = "mc.root";         // optional overlay (if available)
  const char* tree     = "events";
  const char* br_pS    = "sigma0_p_exit";   // Σ0 momentum at nuclear exit [GeV/c]
  const char* br_isSig = "isExitSignal";    // bool/int: exit-defined Σ0 signal?
  const char* br_isSel = "isSelected";      // bool/int: passes full selection (Λ→pπ)
  constexpr double pmin_p     = 0.25;       // proton threshold [GeV/c]
  constexpr double pmin_pi    = 0.10;       // pion   threshold [GeV/c]
  constexpr double Emin_gamma = 0.050;      // photon lab energy threshold [GeV]  <-- set to your detector
  constexpr double vline_pS   = 0.50;       // dashed vertical reference line (set <0 to use analytic pS_thr)
  constexpr int    nbins = 30;
  constexpr double xlow  = 0.0, xhigh = 3.0; // p_Σ range [GeV/c]
  constexpr int    Np  = 600;               // sampling in p_Σ
  constexpr int    Nmu = 1200;              // sampling in μ1 = cosθ*_Σ
}

static inline double clip(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}

void plotSigma0LambdaVisibility() {
  heron::macro::run_with_guard("plotSigma0LambdaVisibility", [&]() {
  Plotter{}.set_global_style();
  gStyle->SetOptStat(0);

  // --------- Two-body constants ----------
  // Λ→pπ at Λ rest
  const double Ep_star  = (mL*mL + mp*mp  - mpi*mpi)/(2.0*mL);
  const double Epi_star = (mL*mL + mpi*mpi - mp*mp )/(2.0*mL);
  const double p_ppi    = std::sqrt(std::max(0.0, Ep_star*Ep_star - mp*mp));

  // Σ0→Λγ at Σ0 rest
  const double EL_star  = (mS0*mS0 + mL*mL)/(2.0*mS0);
  const double pL_star  = (mS0*mS0 - mL*mL)/(2.0*mS0); // = E_γ^*

  // Lab thresholds as energies for tracks
  const double Ep_min  = std::sqrt(mp*mp  + pmin_p*pmin_p);
  const double Epi_min = std::sqrt(mpi*mpi + pmin_pi*pmin_pi);

  // Single-decay Λ→pπ "gamma threshold"
  const double gamma_thr = (Ep_min + Epi_min)/mL;

  // Analytic Σ0 turn-on for Λ-only: solve quadratic for β_Σ with μ1=+1
  double pS_thr = 0.0;
  {
    const double A = pL_star*pL_star + (gamma_thr*mL)*(gamma_thr*mL);
    const double B = 2.0*EL_star*pL_star;
    const double C = EL_star*EL_star - (gamma_thr*mL)*(gamma_thr*mL);
    const double disc = B*B - 4.0*A*C;
    if (disc>0.0) {
      const double betaS = (-B + std::sqrt(disc))/(2.0*A); // physical (+) root
      if (betaS>0.0 && betaS<1.0) {
        const double gammaS = 1.0/std::sqrt(1.0 - betaS*betaS);
        pS_thr = mS0 * betaS * gammaS;
      }
    }
  }

  // --------- Build analytic-numeric curves vs p_Σ ----------
  auto gF     = new TGraph(Np); // Λ-only
  auto gBoth  = new TGraph(Np); // Λ + γ
  gF   ->SetLineWidth(3);
  gBoth->SetLineWidth(3); gBoth->SetLineStyle(7); // dashed
  gF   ->SetTitle(";p_{#Sigma^{0}} at Nuclear Exit [GeV/c];Kinematic Visibility Efficiency");

  for (int i=0; i<Np; ++i) {
    const double pS = xlow + (xhigh - xlow)*(i + 0.5)/Np;
    const double ES = std::sqrt(pS*pS + mS0*mS0);
    const double gammaS = ES/mS0;
    const double betaS  = (ES>0.0) ? pS/ES : 0.0;

    double sumLambda = 0.0;
    double sumBoth   = 0.0;

    // Average over μ1 = cosθ*_Σ (isotropic Σ0→Λγ)
    for (int j=0; j<Nmu; ++j) {
      const double mu1 = -1.0 + 2.0*(j + 0.5)/Nmu;

      // Λ four-momentum after boost from Σ-rest (boost along z)
      const double EL_lab = gammaS * (EL_star + betaS * pL_star * mu1);
      const double pLz    = gammaS * (pL_star * mu1 + betaS * EL_star);
      const double pLperp = pL_star * std::sqrt(std::max(0.0, 1.0 - mu1*mu1));
      const double pL_lab = std::sqrt(pLz*pLz + pLperp*pLperp);

      const double gammaL = EL_lab / mL;
      const double betaL  = (EL_lab>0.0) ? pL_lab / EL_lab : 0.0;

      // Conditional Λ→pπ visibility at this μ1 (analytic in μ2)
      double frac_mu2 = 0.0;
      if (betaL>0.0) {
        const double mu2_min = (Ep_min/gammaL - Ep_star) / (betaL * p_ppi);
        const double mu2_max = (Epi_star     - Epi_min/gammaL) / (betaL * p_ppi);
        const double a = std::max(-1.0, mu2_min);
        const double b = std::min(+1.0, mu2_max);
        if (b>a) frac_mu2 = 0.5*(b - a);
        if (frac_mu2<0.0) frac_mu2 = 0.0;
        if (frac_mu2>1.0) frac_mu2 = 1.0;
      } else {
        const bool pass = (gammaL*Ep_star >= Ep_min) && (gammaL*Epi_star >= Epi_min);
        frac_mu2 = pass ? 1.0 : 0.0;
      }

      sumLambda += frac_mu2;

      // Photon lab energy at this μ1 (γ points opposite Λ in Σ rest)
      const double Egamma_lab = gammaS * pL_star * (1.0 - betaS * mu1);
      if (Egamma_lab >= Emin_gamma) sumBoth += frac_mu2;
    }

    const double fLambda = sumLambda / Nmu;
    const double fBoth   = sumBoth   / Nmu;

    gF   ->SetPoint(i, pS, fLambda);
    gBoth->SetPoint(i, pS, fBoth);
  }

  // --------- Optional TEfficiency overlay from a file (if available) ----------
  TEfficiency* eff = nullptr;
  if (gSystem->AccessPathName(file)) {
    std::cout << "Warning: MC file '" << file
              << "' not found; skipping MC efficiency overlay.\n";
  } else {
    TFile f(file, "READ");
    if (f.IsZombie()) {
      std::cout << "Warning: MC file '" << file
                << "' could not be opened; skipping MC efficiency overlay.\n";
    } else {
      TTree* T = (TTree*) f.Get(tree);
      if (!T) {
        std::cout << "Warning: tree '" << tree
                  << "' not found in '" << file
                  << "'; skipping MC efficiency overlay.\n";
      } else if (!T->GetBranch(br_pS)
                 || !T->GetBranch(br_isSig)
                 || !T->GetBranch(br_isSel)) {
        std::cout << "Warning: required branches not found in '" << file
                  << "'; skipping MC efficiency overlay.\n";
      } else {
        double pS=0.0; int isSig=0, isSel=0;
        T->SetBranchAddress(br_pS,    &pS);
        T->SetBranchAddress(br_isSig, &isSig);
        T->SetBranchAddress(br_isSel, &isSel);

        eff = new TEfficiency("effS",";p_{#Sigma^{0}} at Nuclear Exit [GeV/c];#varepsilon_{fid} (MC)",
                              nbins, xlow, xhigh);
        eff->SetStatisticOption(TEfficiency::kFCP);

        const Long64_t n = T->GetEntries();
        for (Long64_t k=0;k<n;++k) {
          T->GetEntry(k);
          if (!isSig) continue;   // denominator: exit-defined Σ0 decays
          eff->Fill(isSel, pS);   // numerator: selected Λ→pπ events
        }
      }
    }
  }

  // --------- Draw ----------
  TCanvas c("c","Σ^{0}→Λγ: Λ-only vs Λ+γ visibility", 900, 650);
  auto* frame = c.DrawFrame(xlow, 0.0, xhigh, 1.05,
    ";p_{#Sigma^{0}} at Nuclear Exit [GeV/c];Kinematic Visibility Efficiency");
  frame->GetYaxis()->SetTitleOffset(1.15);

  gF   ->Draw("L");
  gBoth->Draw("L same");
  if (eff) eff->Draw("pe same");

  // Vertical line at analytic Σ0 turn-on (Λ-only)
  const double vline_pS = (vline_pS >= 0.0) ? vline_pS : pS_thr;
  const bool usingAnalyticLine = (vline_pS < 0.0);
  TLine L(vline_pS, 0.0, vline_pS, 1.05);
  L.SetLineStyle(2); L.SetLineWidth(2); L.Draw("same");

  // Legend
  TLegend leg(0.48, 0.16, 0.88, 0.38); leg.SetBorderSize(0);
  if (eff) leg.AddEntry(eff, "#varepsilon_{fid}(p_{#Sigma^{0}}) (MC)", "pe");
  leg.AddEntry(gF,    "F_{kin}^{#Lambda only}(p_{#Sigma^{0}})", "l");
  leg.AddEntry(gBoth, Form("F_{kin}^{#Lambda + #gamma}(p_{#Sigma^{0}}), E_{#gamma}^{min}=%.0f MeV",
                           1000.0*Emin_gamma), "l");
  leg.AddEntry(&L,    usingAnalyticLine
                    ? Form("p^{thr}_{#Sigma^{0}} (#Lambda only) = %.3f GeV/c", vline_pS)
                    : Form("Reference p_{#Sigma^{0}} = %.3f GeV/c", vline_pS), "l");
  leg.Draw();

  const std::string out =
      plot_output_file("sigma0_lambda_gamma_visibility_efficiency").string();
  c.SaveAs(out.c_str());

  });
}
