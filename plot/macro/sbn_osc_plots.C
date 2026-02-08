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
//   root -l -q 'sbn_osc_plots.C("fig1p9")'   // 3+1 multi-panel P(νμ→να) vs L/E, like the example figure
//   root -l -q 'sbn_osc_plots.C("all")'
//
// Each mode writes a PDF in the current directory.

#include <array>
#include <algorithm>
#include <cmath>
#include <complex>
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

// ---------------------------
// 3+1 vacuum oscillations (4-flavour) helper code
//   P_{α→β} = δ_{αβ}
//            -4 Σ_{i>j} Re(U_{αi}U*_{βi}U*_{αj}U_{βj}) sin^2(Δij)
//            +2 Σ_{i>j} Im(...) sin(2Δij),
// where Δij = 1.267 Δm^2_ij [eV^2] * (L/E) [km/GeV]
// ---------------------------
using cplx = std::complex<double>;
using Mat4 = std::array<std::array<cplx,4>,4>;

static double deg2rad(double deg) { return deg * (std::acos(-1.0) / 180.0); }

static double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static Mat4 Identity4() {
  Mat4 U{};
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) U[r][c] = (r == c) ? cplx(1.0,0.0) : cplx(0.0,0.0);
  }
  return U;
}

// Left-multiply U by a complex rotation R_ij(theta, delta) in the (i,j) subspace.
// This reproduces the usual PDG-like build order when applied as:
//   U = R34 * R24 * R14 * R23 * R13 * R12
static void LeftMultiplyRotation(Mat4& U, int i, int j, double theta, double delta) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const cplx e_m = std::exp(cplx(0.0, -delta));
  const cplx e_p = std::exp(cplx(0.0, +delta));

  // Build the 2x2 block action on rows i and j of U: U <- R * U
  for (int col = 0; col < 4; ++col) {
    const cplx Ui = U[i][col];
    const cplx Uj = U[j][col];
    U[i][col] =  c*Ui + s*e_m*Uj;
    U[j][col] = -s*e_p*Ui + c*Uj;
  }
}

struct OscParams3p1 {
  // Standard 3-flavour (NO) parameters
  double th12 = deg2rad(33.44);
  double th13 = deg2rad(8.57);
  double th23 = deg2rad(49.20);
  double d13  = deg2rad(195.0);
  double dm21 = 7.42e-5;   // eV^2
  double dm31 = 2.517e-3;  // eV^2 (NO: m3^2 - m1^2)

  // Sterile (3+1) parameters (example close to LSND-style best-fit quoted in the screenshot caption)
  double dm41 = 1.20;              // eV^2
  double th14 = deg2rad(18.4);     // degrees -> radians
  double th24 = deg2rad(5.2);
  double th34 = deg2rad(0.0);
  double d24  = deg2rad(0.0);
  double d34  = deg2rad(0.0);
};

static Mat4 BuildPMNS_3p1(const OscParams3p1& p) {
  Mat4 U = Identity4();

  // Build U = R34 * R24 * R14 * R23 * R13 * R12 (PDG-like extension).
  LeftMultiplyRotation(U, 0, 1, p.th12, 0.0);
  LeftMultiplyRotation(U, 0, 2, p.th13, p.d13);
  LeftMultiplyRotation(U, 1, 2, p.th23, 0.0);
  LeftMultiplyRotation(U, 0, 3, p.th14, 0.0);
  LeftMultiplyRotation(U, 1, 3, p.th24, p.d24);
  LeftMultiplyRotation(U, 2, 3, p.th34, p.d34);
  return U;
}

// Flavour indices: 0=e, 1=mu, 2=tau, 3=sterile
static double P_vac_3p1(int alpha, int beta, double LE_km_per_GeV,
                        const std::array<double,4>& m2, const Mat4& U)
{
  double P = (alpha == beta) ? 1.0 : 0.0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < i; ++j) {
      const double dm2 = m2[i] - m2[j];
      const double Dij = 1.267 * dm2 * LE_km_per_GeV;
      const cplx X = U[alpha][i] * std::conj(U[beta][i]) * std::conj(U[alpha][j]) * U[beta][j];
      P -= 4.0 * std::real(X) * sinsq(Dij);
      P += 2.0 * std::imag(X) * std::sin(2.0 * Dij);
    }
  }
  return clamp01(P);
}

static TGraph* MakeProbGraphLE(int alpha, int beta,
                               double xMin, double xMax, int N,
                               bool logX,
                               const std::array<double,4>& m2, const Mat4& U)
{
  TGraph* g = new TGraph(N);
  for (int k = 0; k < N; ++k) {
    double x = 0.0;
    if (logX) {
      const double lx0 = std::log10(xMin);
      const double lx1 = std::log10(xMax);
      const double t   = (N == 1) ? 0.0 : (double)k / (double)(N - 1);
      x = std::pow(10.0, lx0 + (lx1 - lx0) * t);
    } else {
      x = xMin + (xMax - xMin) * (double)k / (double)(N - 1);
    }
    const double P = P_vac_3p1(alpha, beta, x, m2, U);
    g->SetPoint(k, x, P);
  }
  return g;
}

static void SetNiceStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleOffset(1.10, "X");
  gStyle->SetTitleOffset(1.35, "Y");
  // Slightly roomier margins: avoids clipped y-titles and gives a clean header band.
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

// Two-line header drawn in the *top margin* (so it doesn't collide with top ticks).
static void DrawHeader(const char* line1, const char* line2 = nullptr) {
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.15, 0.965, line1);
  if (line2 && std::string(line2).size() > 0) {
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.15, 0.925, line2);
  }
}

// Multi-panel (a)-(f) P(νμ→να) vs L/E figure (3+1 vacuum), similar to the screenshot layout:
//   (a) 0..40000, (b) 0..4000, (c) 0..200, (d) 0..20 (y-zoom),
//   (e) logx + logy (0.1..1e5), (f) logx (0.1..1e5).
void sbn_fig1p9() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Physics inputs (edit as desired).
  const OscParams3p1 p; // defaults set above
  const Mat4 U_3p1 = BuildPMNS_3p1(p);
  const std::array<double,4> m2 = {0.0, p.dm21, p.dm31, p.dm41};

  // “No sterile” comparison: same 3-flavour parameters, but θ14=θ24=θ34=0.
  OscParams3p1 p0 = p;
  p0.th14 = 0.0; p0.th24 = 0.0; p0.th34 = 0.0;
  p0.d24  = 0.0; p0.d34  = 0.0;
  const Mat4 U_3fl = BuildPMNS_3p1(p0);

  struct PanelSpec {
    double xMin, xMax;
    double yMin, yMax;
    int    N;
    bool   logX;
    bool   logY;
    bool   drawNoSterile;
    const char* tag;
  };

  const PanelSpec panels[6] = {
    {0.0,   40000.0, 0.0,   1.0,   40000, false, false, false, "(a)"},
    {0.0,    4000.0, 0.0,   1.0,   16000, false, false, false, "(b)"},
    {0.0,     200.0, 0.0,   1.0,    4000, false, false, false, "(c)"},
    {0.0,      20.0, 0.0,   0.020,  4000, false, false, false, "(d)"},
    {1e-1,   1e5,    1e-9,  1.0,    6000, true,  true,  true,  "(e)"},
    {1e-1,   1e5,    0.0,   1.0,    6000, true,  false, true,  "(f)"},
  };

  // Colors roughly matching the screenshot: blue (μ→e), orange (μ→μ), green (μ→τ)
  const int col_mue  = kBlue+1;
  const int col_mumu = kOrange+7;
  const int col_mut  = kGreen+2;

  TCanvas* c = new TCanvas("c_fig1p9", "", 1200, 1500);
  c->Divide(2, 3, 0.01, 0.02);

  for (int ip = 0; ip < 6; ++ip) {
    c->cd(ip + 1);
    TPad* pad = (TPad*)gPad;
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.05);
    pad->SetTopMargin(0.06);
    pad->SetBottomMargin(0.18);
    pad->SetLogx(panels[ip].logX);
    pad->SetLogy(panels[ip].logY);

    TMultiGraph* mg = new TMultiGraph();

    // νμ -> να, with α = {e, μ, τ}
    TGraph* g_mue  = MakeProbGraphLE(1, 0, panels[ip].xMin, panels[ip].xMax, panels[ip].N, panels[ip].logX, m2, U_3p1);
    TGraph* g_mumu = MakeProbGraphLE(1, 1, panels[ip].xMin, panels[ip].xMax, panels[ip].N, panels[ip].logX, m2, U_3p1);
    TGraph* g_mut  = MakeProbGraphLE(1, 2, panels[ip].xMin, panels[ip].xMax, panels[ip].N, panels[ip].logX, m2, U_3p1);

    g_mue->SetLineColor(col_mue);
    g_mue->SetLineWidth(2);

    g_mumu->SetLineColor(col_mumu);
    g_mumu->SetLineWidth(2);

    g_mut->SetLineColor(col_mut);
    g_mut->SetLineWidth(2);

    mg->Add(g_mue,  "L");
    mg->Add(g_mumu, "L");
    mg->Add(g_mut,  "L");

    // Optional “no sterile” dashed reference (matches the screenshot: shown in (e) and (f))
    TGraph* g_mue_3fl = nullptr;
    if (panels[ip].drawNoSterile) {
      g_mue_3fl = MakeProbGraphLE(1, 0, panels[ip].xMin, panels[ip].xMax, panels[ip].N, panels[ip].logX, m2, U_3fl);
      g_mue_3fl->SetLineColor(kBlack);
      g_mue_3fl->SetLineStyle(2);
      g_mue_3fl->SetLineWidth(2);
      mg->Add(g_mue_3fl, "L");
    }

    mg->SetTitle(";L / E (km / GeV);P(#nu_{#mu} #rightarrow #nu)");
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(panels[ip].xMin, panels[ip].xMax);
    mg->GetYaxis()->SetRangeUser(panels[ip].yMin, panels[ip].yMax);
    mg->GetXaxis()->SetNdivisions(510);

    if (panels[ip].logX) {
      mg->GetXaxis()->SetMoreLogLabels(true);
      mg->GetXaxis()->SetNoExponent(true);
    }
    if (panels[ip].logY) {
      mg->GetYaxis()->SetMoreLogLabels(true);
      mg->GetYaxis()->SetNoExponent(true);
    }

    // Legend placement (tuned for readability; adjust to taste)
    TLegend* leg = nullptr;
    if (ip < 4)      leg = new TLegend(0.70, 0.74, 0.93, 0.93);
    else             leg = new TLegend(0.60, 0.18, 0.93, 0.45);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->SetFillStyle(0);
    leg->AddEntry(g_mue,  "#nu_{#mu} #rightarrow #nu_{e}",   "l");
    leg->AddEntry(g_mumu, "#nu_{#mu} #rightarrow #nu_{#mu}", "l");
    leg->AddEntry(g_mut,  "#nu_{#mu} #rightarrow #nu_{#tau}","l");
    if (g_mue_3fl) leg->AddEntry(g_mue_3fl, "Oscillations Without A Sterile Neutrino", "l");
    leg->Draw();

    // Panel tag (a)-(f) placed in the bottom margin band.
    TLatex tag;
    tag.SetNDC(true);
    tag.SetTextFont(42);
    tag.SetTextSize(0.060);
    tag.DrawLatex(0.47, 0.03, panels[ip].tag);
  }

  c->SaveAs("sbn_fig1p9_mu_LE_3p1.pdf");
}

void sbn_prob_vs_LE() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double sin22th = 0.01;            // schematic appearance amplitude
  const std::vector<double> dm2 = {0.3, 1.0, 3.0};  // eV^2
  const std::vector<int> colors = {kBlack, kBlue+1, kRed+1};

  const int N = 1400;
  // Extend L/E so even the smallest Δm² curve shows a full max/min within the frame.
  const double xMin = 0.0, xMax = 10.0;    // L/E in km/GeV

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

  mg->SetTitle(";L/E  [km/GeV];P_{#mu e}");
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

  DrawHeader("Vacuum SBL: oscillation phase #propto #Delta m^{2}#times(L/E)",
             Form("P_{#mu e}=sin^{2}2#theta#times sin^{2}(1.267#Delta m^{2}L/E),  sin^{2}2#theta=%.3f", sin22th));

  c->SaveAs("sbn_prob_vs_LE.pdf");
}

void sbn_prob_vs_E() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double dm2     = 1.0;     // eV^2
  const double sin22th = 0.01;
  // Choose baselines so the first oscillation maximum sits inside ~0.2–3 GeV for Δm²=1 eV².
  const std::vector<double> L = {0.30, 0.60, 1.20}; // km
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

  mg->SetTitle(";E_{#nu}  [GeV];P_{#mu e}");
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

  DrawHeader("Vacuum SBL: changing L shifts oscillation features in E_{#nu}",
             Form("#Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f  (curves: L)", dm2, sin22th));

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
  // Give the color palette breathing room without mutating global style.
  c->SetRightMargin(0.18);

  TH2D* h = new TH2D("h_osc",
                     ";E_{#nu}  [GeV];Baseline L  [km];P_{#mu e}",
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

  DrawHeader("Oscillogram: constant L/E bands (vacuum, 2-flavour)",
             Form("#Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f", dm2, sin22th));

  c->SaveAs("sbn_oscillogram.pdf");
}

void sbn_smear() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Use a faster-oscillation benchmark so resolution damping is visually obvious.
  const double dm2     = 3.0;   // eV^2
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

  mg->SetTitle(";E_{#nu}  [GeV];#LT P_{#mu e} #GT");
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

  DrawHeader("Energy resolution averages rapid oscillations #rightarrow damping",
             Form("Gaussian smearing in E:  L=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f", L, dm2, sin22th));

  c->SaveAs("sbn_smearing_damping.pdf");
}

void sbn_bias() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Faster-oscillation benchmark so a small scale bias makes a visible phase shift.
  const double dm2     = 3.0;   // eV^2
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

  mg->SetTitle(";E_{rec}  [GeV];P_{#mu#mu}");
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

  DrawHeader("Energy-scale bias shifts oscillation phase in E_{rec}",
             Form("E_{rec}=(1+b)E_{true}:  L=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.2f", L, dm2, sin22th));

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

  // Use a smooth graph instead of a coarse bin-by-bin ratio to avoid jagged steps.
  const int N = 2000;
  const double Emin = 0.2, Emax = 3.0;

  std::vector<double> E(N), rN(N), rF(N);
  for (int i = 0; i < N; ++i) {
    const double t = (N == 1) ? 0.0 : (double)i / (double)(N - 1);
    E[i] = Emin + (Emax - Emin) * t;
    const double N0 = toy_flux(E[i]) * toy_xsec(E[i]);
    rN[i] = N0 * P_surv(Lnear, E[i], dm2, sin22th);
    rF[i] = N0 * P_surv(Lfar,  E[i], dm2, sin22th);
  }

  auto trapz = [&](const std::vector<double>& y) -> double {
    double sum = 0.0;
    for (int i = 0; i < N - 1; ++i) {
      const double dx = E[i+1] - E[i];
      sum += 0.5 * (y[i] + y[i+1]) * dx;
    }
    return sum;
  };

  const double IN = trapz(rN);
  const double IF = trapz(rF);

  TGraph* gR = new TGraph(N);
  for (int i = 0; i < N; ++i) {
    const double n = (IN > 0.0) ? (rN[i] / IN) : rN[i];
    const double f = (IF > 0.0) ? (rF[i] / IF) : rF[i];
    const double R = (n > 0.0)  ? (f / n)      : 0.0;
    gR->SetPoint(i, E[i], R);
  }

  TCanvas* c = new TCanvas("c_nearfar", "", 900, 650);
  gR->SetLineWidth(3);
  gR->SetLineColor(kBlack);
  gR->SetTitle(";E_{#nu}  [GeV];Far/Near (shape norm.)");
  gR->Draw("AL");
  gR->GetXaxis()->SetNdivisions(510);
  gR->GetYaxis()->SetRangeUser(0.90, 1.10);

  // Reference line at unity.
  TLine* one = new TLine(Emin, 1.0, Emax, 1.0);
  one->SetLineStyle(2);
  one->SetLineWidth(2);
  one->SetLineColor(kGray+2);
  one->Draw("SAME");

  DrawHeader("Near/Far shape ratio isolates oscillation-driven distortions",
             Form("Toy flux#timesxsec (shape-normalised):  L_{near}=%.2f km, L_{far}=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.2f",
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

  frame->GetXaxis()->SetMoreLogLabels(true);
  frame->GetXaxis()->SetNoExponent(true);
  frame->GetYaxis()->SetMoreLogLabels(true);
  frame->GetYaxis()->SetNoExponent(true);

  DrawHeader("Sterile-parameter plane (vacuum SBL limit)",
             "Overlay allowed/excluded contours (TGraph from two-column text files)");

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

void sbn_plot_all() {
  sbn_prob_vs_LE();
  sbn_prob_vs_E();
  sbn_oscillogram();
  sbn_smear();
  sbn_bias();
  sbn_nearfar();
  sbn_dm2_sin22_template();
  // sbn_fig1p9(); // uncomment if you want the multi-panel 3+1 figure in the "all" bundle
}

void sbn_osc_plots(const char* which = "all") {
  std::string w(which ? which : "all");

  if      (w == "prob_le")            sbn_prob_vs_LE();
  else if (w == "prob_E")             sbn_prob_vs_E();
  else if (w == "oscillogram")        sbn_oscillogram();
  else if (w == "smear")              sbn_smear();
  else if (w == "bias")               sbn_bias();
  else if (w == "nearfar")            sbn_nearfar();
  else if (w == "dm2_sin22_template") sbn_dm2_sin22_template();
  else if (w == "fig1p9")             sbn_fig1p9();
  else if (w == "all")                sbn_plot_all();
  else {
    std::cout << "Unknown mode: " << w << "\n"
              << "Options: prob_le, prob_E, oscillogram, smear, bias, nearfar, dm2_sin22_template, fig1p9, all\n";
  }
}
