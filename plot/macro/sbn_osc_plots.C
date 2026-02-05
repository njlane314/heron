// sbn_osc_plots.C
// Schematic short-baseline oscillation plots (ROOT macro).
//
// Usage examples (batch mode):
//   root -l -q 'sbn_osc_plots.C("prob_le")'
//   root -l -q 'sbn_osc_plots.C("prob_E")'
//   root -l -q 'sbn_osc_plots.C("oscillogram")'
//   root -l -q 'sbn_osc_plots.C("smear")'
//   root -l -q 'sbn_osc_plots.C("bias")'
//   root -l -q 'sbn_osc_plots.C("nearfar")'
//   root -l -q 'sbn_osc_plots.C("dm2_sin22_template")'
//
// Each mode writes a PDF in the current directory.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TStyle.h"

static double sinsq(double x) {
  const double s = std::sin(x);
  return s*s;
}

// Effective 2-flavour (vacuum) oscillation phase convention:
//   sin^2( 1.267 * (Δm^2[eV^2]) * (L[km]) / (E[GeV]) )
static double phase(double L_km, double E_GeV, double dm2_eV2) {
  return 1.267 * dm2_eV2 * L_km / E_GeV;
}

// Appearance probability (e.g. νμ→νe in a 3+1 SBL limit):
//   P_app ≈ sin^2(2θ) * sin^2(phase)
static double P_app(double L_km, double E_GeV, double dm2_eV2, double sin22th) {
  return sin22th * sinsq( phase(L_km, E_GeV, dm2_eV2) );
}

// Survival probability (toy 2-flavour disappearance):
//   P_surv ≈ 1 - sin^2(2θ) * sin^2(phase)
static double P_surv(double L_km, double E_GeV, double dm2_eV2, double sin22th) {
  return 1.0 - sin22th * sinsq( phase(L_km, E_GeV, dm2_eV2) );
}

// Average P_app over a Gaussian energy resolution around E0.
// This is a schematic “resolution damping” model.
static double P_app_Esmeared_avg(double L_km, double E0_GeV, double dm2_eV2,
                                double sin22th, double frac_sigmaE_over_E)
{
  if (frac_sigmaE_over_E <= 0.0) return P_app(L_km, E0_GeV, dm2_eV2, sin22th);

  const double sigma = frac_sigmaE_over_E * E0_GeV;
  if (sigma <= 0.0) return P_app(L_km, E0_GeV, dm2_eV2, sin22th);

  const int    N     = 600;
  const double Emin  = std::max(0.001, E0_GeV - 5.0*sigma);
  const double Emax  = E0_GeV + 5.0*sigma;

  double num = 0.0, den = 0.0;
  for (int i = 0; i < N; ++i) {
    const double E = Emin + (Emax - Emin) * (i + 0.5) / N;
    const double z = (E - E0_GeV) / sigma;
    const double w = std::exp(-0.5 * z*z);
    num += w * P_app(L_km, E, dm2_eV2, sin22th);
    den += w;
  }
  return (den > 0.0) ? (num / den) : P_app(L_km, E0_GeV, dm2_eV2, sin22th);
}

// Minimal “beam-like” toy spectrum shape (arbitrary units).
static double toy_flux(double E_GeV) {
  const double mu  = 0.75;
  const double sig = 0.25;
  const double z   = (E_GeV - mu) / sig;
  return std::exp(-0.5*z*z);
}

// Naive rising cross-section factor (schematic).
static double toy_xsec(double E_GeV) {
  return std::max(0.0, E_GeV);
}

static void SetNiceStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleOffset(1.10, "X");
  gStyle->SetTitleOffset(1.25, "Y");
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void sbn_prob_vs_LE() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double sin22th = 0.01;            // schematic appearance amplitude
  const std::vector<double> dm2 = {0.3, 1.0, 3.0};  // eV^2
  const std::vector<int> colors = {kBlack, kBlue+1, kRed+1};

  const int N = 1400;
  const double xMin = 0.0, xMax = 5.0;    // L/E in km/GeV

  TCanvas* c = new TCanvas("c_prob_le", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();

  std::vector<TGraph*> gs;
  gs.reserve(dm2.size());

  for (size_t i = 0; i < dm2.size(); ++i) {
    TGraph* g = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      const double LE = xMin + (xMax - xMin) * j / (N - 1.0);
      const double P  = sin22th * sinsq(1.267 * dm2[i] * LE);
      g->SetPoint(j, LE, P);
    }
    g->SetLineColor(colors[i]);
    g->SetLineWidth(3);
    g->SetLineStyle(1 + (int)i);
    mg->Add(g, "L");
    gs.push_back(g);
  }

  mg->SetTitle(";L/E  [km/GeV];P(#nu_{#mu}#rightarrow#nu_{e})  (effective 2-flavour)");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 1.05*sin22th);

  TLegend* leg = new TLegend(0.56, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (size_t i = 0; i < dm2.size(); ++i) {
    leg->AddEntry(gs[i], Form("#Delta m^{2}=%.1f eV^{2}", dm2[i]), "l");
  }
  leg->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.15, 0.88, Form("sin^{2}2#theta = %.3f", sin22th));

  c->SaveAs("sbn_prob_vs_LE.pdf");
}

void sbn_prob_vs_E() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double dm2     = 1.0;     // eV^2
  const double sin22th = 0.01;
  const std::vector<double> L = {0.10, 0.30, 0.60}; // km (e.g. near/intermediate/far)
  const std::vector<int> colors = {kBlack, kBlue+1, kRed+1};

  const int N = 1400;
  const double Emin = 0.2, Emax = 3.0; // GeV

  TCanvas* c = new TCanvas("c_prob_E", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();

  std::vector<TGraph*> gs;
  gs.reserve(L.size());

  for (size_t i = 0; i < L.size(); ++i) {
    TGraph* g = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      const double E = Emin + (Emax - Emin) * j / (N - 1.0);
      const double P = P_app(L[i], E, dm2, sin22th);
      g->SetPoint(j, E, P);
    }
    g->SetLineColor(colors[i]);
    g->SetLineWidth(3);
    g->SetLineStyle(1 + (int)i);
    mg->Add(g, "L");
    gs.push_back(g);
  }

  mg->SetTitle(";E_{#nu}  [GeV];P(#nu_{#mu}#rightarrow#nu_{e})  (effective 2-flavour)");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 1.05*sin22th);

  TLegend* leg = new TLegend(0.56, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (size_t i = 0; i < L.size(); ++i) {
    leg->AddEntry(gs[i], Form("L = %.2f km", L[i]), "l");
  }
  leg->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.15, 0.88, Form("#Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f", dm2, sin22th));

  c->SaveAs("sbn_prob_vs_E.pdf");
}

void sbn_oscillogram() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double dm2     = 1.0;   // eV^2
  const double sin22th = 0.01;

  const int nE = 240;
  const int nL = 220;
  const double Emin = 0.2, Emax = 3.0;    // GeV
  const double Lmin = 0.05, Lmax = 2.0;   // km

  TCanvas* c = new TCanvas("c_oscillogram", "", 950, 700);
  gStyle->SetPadRightMargin(0.16);

  TH2D* h = new TH2D("h_osc",
                     ";E_{#nu}  [GeV];Baseline L  [km];P(#nu_{#mu}#rightarrow#nu_{e})",
                     nE, Emin, Emax, nL, Lmin, Lmax);

  for (int ix = 1; ix <= nE; ++ix) {
    const double E = h->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= nL; ++iy) {
      const double L = h->GetYaxis()->GetBinCenter(iy);
      h->SetBinContent(ix, iy, P_app(L, E, dm2, sin22th));
    }
  }

  h->GetZaxis()->SetTitleOffset(1.25);
  h->SetMinimum(0.0);
  h->SetMaximum(sin22th);
  h->Draw("COLZ");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.14, 0.92, Form("#Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f", dm2, sin22th));

  c->SaveAs("sbn_oscillogram.pdf");
}

void sbn_smear() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double dm2     = 1.0;   // eV^2
  const double sin22th = 0.01;
  const double L       = 0.60;  // km

  const std::vector<double> frac = {0.0, 0.05, 0.15}; // σE/E
  const std::vector<int> colors  = {kBlack, kBlue+1, kRed+1};

  const int N = 1400;
  const double Emin = 0.2, Emax = 3.0;

  TCanvas* c = new TCanvas("c_smear", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();

  std::vector<TGraph*> gs;
  for (size_t i = 0; i < frac.size(); ++i) {
    TGraph* g = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      const double E = Emin + (Emax - Emin) * j / (N - 1.0);
      const double P = P_app_Esmeared_avg(L, E, dm2, sin22th, frac[i]);
      g->SetPoint(j, E, P);
    }
    g->SetLineColor(colors[i]);
    g->SetLineWidth(3);
    g->SetLineStyle(1 + (int)i);
    mg->Add(g, "L");
    gs.push_back(g);
  }

  mg->SetTitle(";E_{#nu}  [GeV];#LT P(#nu_{#mu}#rightarrow#nu_{e}) #GT  (Gaussian energy averaging)");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 1.05*sin22th);

  TLegend* leg = new TLegend(0.50, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (size_t i = 0; i < frac.size(); ++i) {
    leg->AddEntry(gs[i], Form("#sigma_{E}/E = %.0f%%", 100.0*frac[i]), "l");
  }
  leg->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.15, 0.88, Form("L=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f", L, dm2, sin22th));

  c->SaveAs("sbn_smearing_damping.pdf");
}

void sbn_bias() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double dm2     = 1.0;   // eV^2
  const double sin22th = 0.10;  // larger so phase-shift is visually obvious
  const double L       = 0.60;  // km

  // Energy scale bias: E_rec = (1 + b) * E_true
  const std::vector<double> bias = {0.0, +0.02, -0.02}; // ±2%
  const std::vector<int> colors  = {kBlack, kRed+1, kBlue+1};

  const int N = 1400;
  const double Emin = 0.2, Emax = 3.0;

  TCanvas* c = new TCanvas("c_bias", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();

  std::vector<TGraph*> gs;
  for (size_t i = 0; i < bias.size(); ++i) {
    TGraph* g = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      const double Erec  = Emin + (Emax - Emin) * j / (N - 1.0);
      const double Etrue = Erec / (1.0 + bias[i]); // what the neutrino energy really was
      const double P     = P_surv(L, Etrue, dm2, sin22th);
      g->SetPoint(j, Erec, P);
    }
    g->SetLineColor(colors[i]);
    g->SetLineWidth(3);
    g->SetLineStyle(1 + (int)i);
    mg->Add(g, "L");
    gs.push_back(g);
  }

  mg->SetTitle(";E_{rec}  [GeV];P(#nu_{#mu}#rightarrow#nu_{#mu})  (effective 2-flavour)");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(1.0 - 1.2*sin22th, 1.01);

  TLegend* leg = new TLegend(0.52, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(gs[0], "Nominal energy scale", "l");
  leg->AddEntry(gs[1], "+2% energy scale", "l");
  leg->AddEntry(gs[2], "-2% energy scale", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.15, 0.88, Form("L=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.2f", L, dm2, sin22th));

  c->SaveAs("sbn_energy_scale_bias.pdf");
}

void sbn_nearfar() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Disappearance-style toy: show how baseline-dependent oscillations distort spectra.
  const double dm2     = 1.0;   // eV^2
  const double sin22th = 0.10;  // schematic
  const double Lnear   = 0.10;  // km
  const double Lfar    = 0.60;  // km

  const int nb = 120;
  const double Emin = 0.2, Emax = 3.0;

  TH1D* hN = new TH1D("hN", ";E_{#nu}  [GeV];Shape-normalised rate (arb.)", nb, Emin, Emax);
  TH1D* hF = new TH1D("hF", ";E_{#nu}  [GeV];Shape-normalised rate (arb.)", nb, Emin, Emax);

  for (int i = 1; i <= nb; ++i) {
    const double E = hN->GetBinCenter(i);
    const double N0 = toy_flux(E) * toy_xsec(E);

    hN->SetBinContent(i, N0 * P_surv(Lnear, E, dm2, sin22th));
    hF->SetBinContent(i, N0 * P_surv(Lfar,  E, dm2, sin22th));
  }

  // Shape normalise (remove trivial normalisation differences).
  const double IN = hN->Integral("width");
  const double IF = hF->Integral("width");
  if (IN > 0) hN->Scale(1.0/IN);
  if (IF > 0) hF->Scale(1.0/IF);

  TH1D* hR = (TH1D*)hF->Clone("hR");
  hR->SetTitle(";E_{#nu}  [GeV];(Far spectrum)/(Near spectrum)  (shape-normalised)");
  hR->Divide(hN);

  TCanvas* c = new TCanvas("c_nearfar", "", 900, 650);
  hR->SetLineWidth(3);
  hR->SetLineColor(kBlack);
  hR->GetYaxis()->SetRangeUser(0.8, 1.2);
  hR->Draw("HIST");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.15, 0.88, Form("L_{near}=%.2f km, L_{far}=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.2f",
                                 Lnear, Lfar, dm2, sin22th));

  c->SaveAs("sbn_nearfar_ratio.pdf");
}

// Template for a Δm^2 vs sin^2(2θ) parameter-space panel.
// You can overlay contours from two-column text files via TGraph("file.txt").
void sbn_dm2_sin22_template() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  TCanvas* c = new TCanvas("c_dm2_sin22", "", 900, 650);
  c->SetLogx();
  c->SetLogy();

  // Frame ranges are typical for SBL sterile discussions (adjust as needed).
  TH2D* frame = new TH2D("frame",
                         ";#Delta m^{2}  [eV^{2}];sin^{2}2#theta",
                         10, 1e-3, 1e2,
                         10, 1e-4, 1.0);
  frame->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.15, 0.92, "Template: overlay allowed/excluded contours (two-column text files)");

  // Example (uncomment after you have files):
  // TGraph* g_allowed = new TGraph("data/allowed_region.txt"); // columns: dm2  sin22th
  // g_allowed->SetLineColor(kRed+1);
  // g_allowed->SetLineWidth(3);
  // g_allowed->Draw("L SAME");
  //
  // TGraph* g_excl = new TGraph("data/exclusion_curve.txt");
  // g_excl->SetLineColor(kBlue+1);
  // g_excl->SetLineWidth(3);
  // g_excl->SetLineStyle(2);
  // g_excl->Draw("L SAME");
  //
  // TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.88);
  // leg->AddEntry(g_allowed, "Allowed (example)", "l");
  // leg->AddEntry(g_excl,   "Excluded (example)", "l");
  // leg->Draw();

  c->SaveAs("sbn_dm2_sin22_template.pdf");
}

void sbn_osc_plots(const char* which = "prob_le") {
  std::string w(which ? which : "prob_le");

  if      (w == "prob_le")            sbn_prob_vs_LE();
  else if (w == "prob_E")             sbn_prob_vs_E();
  else if (w == "oscillogram")        sbn_oscillogram();
  else if (w == "smear")              sbn_smear();
  else if (w == "bias")               sbn_bias();
  else if (w == "nearfar")            sbn_nearfar();
  else if (w == "dm2_sin22_template") sbn_dm2_sin22_template();
  else {
    std::cout << "Unknown mode: " << w << "\n"
              << "Options: prob_le, prob_E, oscillogram, smear, bias, nearfar, dm2_sin22_template\n";
  }
}
