// plotOscillation.C
// Schematic short-baseline oscillation plots (ROOT macro).
//
// Usage examples (batch mode):
//   root -l -q 'plotOscillation.C("prob_le")'
//   root -l -q 'plotOscillation.C("prob_E")'
//   root -l -q 'plotOscillation.C("oscillogram")'
//   root -l -q 'plotOscillation.C("smear")'
//   root -l -q 'plotOscillation.C("bias")'
//   root -l -q 'plotOscillation.C("nearfar")'
//   root -l -q 'plotOscillation.C("dm2_sin22_template")'
//
//   root -l -q 'plotOscillation.C("fig1p9_a")'
//   root -l -q 'plotOscillation.C("fig1p9_b")'
//   root -l -q 'plotOscillation.C("fig1p9_c")'
//   root -l -q 'plotOscillation.C("fig1p9_d")'
//   root -l -q 'plotOscillation.C("fig1p9_e")'
//   root -l -q 'plotOscillation.C("fig1p9_f")'
//   root -l -q 'plotOscillation.C("fig1p9_all")'
//
//   root -l -q 'plotOscillation.C("3fl_LE_overview")'
//   root -l -q 'plotOscillation.C("3fl_biprob_SBND")'
//   root -l -q 'plotOscillation.C("3fl_biprob_MicroBooNE")'
//   root -l -q 'plotOscillation.C("3fl_biprob_ICARUS")'
//   root -l -q 'plotOscillation.C("3fl_biprob_SBN")'
//   root -l -q 'plotOscillation.C("3fl_all")'
//
//   root -l -q 'plotOscillation.C("all")'
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
#include "TGraph.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"
#include "include/MacroGuard.hh"

static double sinsq(double x) {
  const double s = std::sin(x);
  return s * s;
}

// Effective 2-flavour (vacuum) oscillation phase convention:
//   sin^2( 1.267 * (Δm^2[eV^2]) * (L[km]) / (E[GeV]) )
static double phase(double L_km, double E_GeV, double dm2_eV2) {
  return 1.267 * dm2_eV2 * L_km / E_GeV;
}

// Appearance probability (toy 2-flavour):
//   P_app ≈ sin^2(2θ) * sin^2(phase)
static double P_app(double L_km, double E_GeV, double dm2_eV2, double sin22th) {
  return sin22th * sinsq(phase(L_km, E_GeV, dm2_eV2));
}

// Survival probability (toy 2-flavour disappearance):
//   P_surv ≈ 1 - sin^2(2θ) * sin^2(phase)
static double P_surv(double L_km, double E_GeV, double dm2_eV2, double sin22th) {
  return 1.0 - sin22th * sinsq(phase(L_km, E_GeV, dm2_eV2));
}

// Average P_app over a Gaussian energy resolution around E0.
// Schematic “resolution damping” model.
static double P_app_Esmeared_avg(double L_km, double E0_GeV, double dm2_eV2,
                                 double sin22th, double frac_sigmaE_over_E)
{
  if (frac_sigmaE_over_E <= 0.0) return P_app(L_km, E0_GeV, dm2_eV2, sin22th);

  const double sigma = frac_sigmaE_over_E * E0_GeV;
  if (sigma <= 0.0) return P_app(L_km, E0_GeV, dm2_eV2, sin22th);

  const int    N    = 600;
  const double Emin = std::max(0.001, E0_GeV - 5.0 * sigma);
  const double Emax = E0_GeV + 5.0 * sigma;

  double num = 0.0, den = 0.0;
  for (int i = 0; i < N; ++i) {
    const double E = Emin + (Emax - Emin) * (i + 0.5) / N;
    const double z = (E - E0_GeV) / sigma;
    const double w = std::exp(-0.5 * z * z);
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
  return std::exp(-0.5 * z * z);
}

// Naive rising cross-section factor (schematic).
static double toy_xsec(double E_GeV) { return std::max(0.0, E_GeV); }

// ---------------------------
// 3-flavour vacuum oscillations helper code
//   P_{α→β} = δ_{αβ}
//            -4 Σ_{i>j} Re(U_{αi}U*_{βi}U*_{αj}U_{βj}) sin^2(Δij)
//            +2 Σ_{i>j} Im(...) sin(2Δij),
// where Δij = 1.267 Δm^2_ij [eV^2] * (L/E) [km/GeV]
// ---------------------------
using cplx = std::complex<double>;
using Mat3 = std::array<std::array<cplx, 3>, 3>;

static double deg2rad(double deg) { return deg * (std::acos(-1.0) / 180.0); }

static double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static Mat3 Identity3() {
  Mat3 U{};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      U[r][c] = (r == c) ? cplx(1.0, 0.0) : cplx(0.0, 0.0);
  return U;
}

// Left-multiply U by a complex rotation R_ij(theta, delta) in the (i,j) subspace.
// PDG-like 3-flavour build order:
//   U = R23 * R13(δ) * R12
static void LeftMultiplyRotation(Mat3& U, int i, int j, double theta, double delta) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const cplx e_m = std::exp(cplx(0.0, -delta));
  const cplx e_p = std::exp(cplx(0.0, +delta));

  for (int col = 0; col < 3; ++col) {
    const cplx Ui = U[i][col];
    const cplx Uj = U[j][col];
    U[i][col] =  c * Ui + s * e_m * Uj;
    U[j][col] = -s * e_p * Ui + c * Uj;
  }
}

struct OscParams3fl {
  double th12 = deg2rad(33.44);
  double th13 = deg2rad(8.57);
  double th23 = deg2rad(49.20);
  double d13  = deg2rad(195.0);
  double dm21 = 7.42e-5;   // eV^2
  double dm31 = 2.517e-3;  // eV^2 (NO: m3^2 - m1^2)
};

static OscParams3fl ParamsNO(double delta_deg = -90.0) {
  OscParams3fl p;
  p.d13  = deg2rad(delta_deg);
  p.dm31 = +std::abs(p.dm31);
  return p;
}
static OscParams3fl ParamsIO(double delta_deg = -90.0) {
  OscParams3fl p = ParamsNO(delta_deg);
  p.dm31 = -std::abs(p.dm31);
  return p;
}

static Mat3 BuildPMNS_3fl(const OscParams3fl& p) {
  Mat3 U = Identity3();
  LeftMultiplyRotation(U, 0, 1, p.th12, 0.0);
  LeftMultiplyRotation(U, 0, 2, p.th13, p.d13);
  LeftMultiplyRotation(U, 1, 2, p.th23, 0.0);
  return U;
}

static void BuildVacuumObjects(const OscParams3fl& p, std::array<double, 3>& m2, Mat3& U) {
  U  = BuildPMNS_3fl(p);
  m2 = {0.0, p.dm21, p.dm31};
}

// Flavour indices: 0=e, 1=mu, 2=tau
static double P_vac_3fl(int alpha, int beta, double LE_km_per_GeV,
                        const std::array<double, 3>& m2, const Mat3& U)
{
  double P = (alpha == beta) ? 1.0 : 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < i; ++j) {
      const double dm2ij = m2[i] - m2[j];
      const double Dij   = 1.267 * dm2ij * LE_km_per_GeV;
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
                               const std::array<double, 3>& m2, const Mat3& U,
                               double yFloor = 0.0,
                               bool oneMinus = false)
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
    double P = P_vac_3fl(alpha, beta, x, m2, U);
    double y = oneMinus ? clamp01(1.0 - P) : P;
    if (yFloor > 0.0 && y < yFloor) y = yFloor;
    g->SetPoint(k, x, y);
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
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNumberContours(60);
  TGaxis::SetMaxDigits(4);
}

// Two-line header drawn in the top margin.
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

// --------------------------
// 3-flavour SBL fig1p9 panels
// --------------------------

// Plots (μ-flavour initial state):
//   blue : P(νμ→νe)
//   green: P(νμ→ντ)
//   red  : 1 - P(νμ→νμ)
static void sbn_fig1p9_panel(const char* cname,
                             const char* outPdf,
                             double xMin, double xMax,
                             double yMin, double yMax,
                             int N,
                             bool logX,
                             bool logY)
{
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const OscParams3fl p = ParamsNO(-90.0);
  Mat3 U;
  std::array<double, 3> m2;
  BuildVacuumObjects(p, m2, U);

  const double yFloor = (logY ? std::max(1e-14, 0.1 * yMin) : 0.0);

  TCanvas* c = new TCanvas(cname, "", 900, 650);
  c->SetLogx(logX);
  c->SetLogy(logY);

  TMultiGraph* mg = new TMultiGraph();

  TGraph* g_mue = MakeProbGraphLE(1, 0, xMin, xMax, N, logX, m2, U, yFloor, false);
  TGraph* g_mut = MakeProbGraphLE(1, 2, xMin, xMax, N, logX, m2, U, yFloor, false);
  TGraph* g_dis = MakeProbGraphLE(1, 1, xMin, xMax, N, logX, m2, U, yFloor, true);

  g_mue->SetLineColor(kBlue + 1);
  g_mue->SetLineWidth(3);

  g_mut->SetLineColor(kGreen + 2);
  g_mut->SetLineWidth(3);

  g_dis->SetLineColor(kRed + 1);
  g_dis->SetLineWidth(3);

  mg->Add(g_mue, "L");
  mg->Add(g_mut, "L");
  mg->Add(g_dis, "L");

  mg->SetTitle(";L/E  [km/GeV];Probability");
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(xMin, xMax);
  mg->GetYaxis()->SetRangeUser(yMin, yMax);
  mg->GetXaxis()->SetNdivisions(510);

  if (logX) {
    mg->GetXaxis()->SetMoreLogLabels(true);
    mg->GetXaxis()->SetNoExponent(true);
  }
  if (logY) {
    mg->GetYaxis()->SetMoreLogLabels(false);
    mg->GetYaxis()->SetNoExponent(false);
    mg->GetYaxis()->SetTitleOffset(1.55);
  }

  TLegend* leg = (!logX && !logY) ? new TLegend(0.56, 0.72, 0.88, 0.88)
                                  : new TLegend(0.52, 0.18, 0.88, 0.42);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(g_mue, "P_{#mu e}", "l");
  leg->AddEntry(g_mut, "P_{#mu #tau}", "l");
  leg->AddEntry(g_dis, "1 - P_{#mu#mu} (=P_{#mu e}+P_{#mu #tau})", "l");
  leg->Draw();

  const double rad2deg = 180.0 / std::acos(-1.0);
  DrawHeader("3-flavour vacuum: short-baseline limit (standard oscillations are small)",
             Form("#Delta m^{2}_{31}=%.3e eV^{2},  #theta_{13}=%.2f^{#circ},  #theta_{23}=%.2f^{#circ}   (note: #mu-disappearance dominated by #nu_{#tau})",
                  p.dm31, p.th13 * rad2deg, p.th23 * rad2deg));

  c->SaveAs(outPdf);
}

void sbn_fig1p9_a() { sbn_fig1p9_panel("c_fig1p9_a", "sbn_3fl_sbl_LE_0to10.pdf",     0.0, 10.0, 0.0,   1.2e-3, 5000, false, false); }
void sbn_fig1p9_b() { sbn_fig1p9_panel("c_fig1p9_b", "sbn_3fl_sbl_LE_0to5.pdf",      0.0,  5.0, 0.0,   3.2e-4, 5000, false, false); }
void sbn_fig1p9_c() { sbn_fig1p9_panel("c_fig1p9_c", "sbn_3fl_sbl_LE_0to2.pdf",      0.0,  2.0, 0.0,   6.0e-5, 5000, false, false); }
void sbn_fig1p9_d() { sbn_fig1p9_panel("c_fig1p9_d", "sbn_3fl_sbl_LE_0to1.pdf",      0.0,  1.0, 0.0,   1.6e-5, 5000, false, false); }
void sbn_fig1p9_e() { sbn_fig1p9_panel("c_fig1p9_e", "sbn_3fl_sbl_LE_logy.pdf",      0.05, 10.0, 1e-10, 2.0e-3, 6000, false, true ); }
void sbn_fig1p9_f() { sbn_fig1p9_panel("c_fig1p9_f", "sbn_3fl_sbl_LE_logx_logy.pdf", 0.05, 10.0, 1e-10, 2.0e-3, 6000, true,  true ); }

void sbn_fig1p9_all() {
  sbn_fig1p9_a();
  sbn_fig1p9_b();
  sbn_fig1p9_c();
  sbn_fig1p9_d();
  sbn_fig1p9_e();
  sbn_fig1p9_f();
}

// -------------------------------------------
// 3-flavour wide-range L/E overview (vacuum)
// -------------------------------------------

static double YUserAtFrac(double frac) {
  const double y1 = gPad->GetUymin();
  const double y2 = gPad->GetUymax();
  if (gPad->GetLogy()) {
    const double ly1 = std::log10(std::max(1e-300, y1));
    const double ly2 = std::log10(std::max(1e-300, y2));
    return std::pow(10.0, ly1 + frac * (ly2 - ly1));
  }
  return y1 + frac * (y2 - y1);
}

static void DrawLEMarker(double x, const char* label, double yFrac = 0.86) {
  gPad->Update();
  const double y1 = gPad->GetUymin();
  const double y2 = gPad->GetUymax();

  TLine* ln = new TLine(x, y1, x, y2);
  ln->SetLineStyle(2);
  ln->SetLineWidth(2);
  ln->SetLineColor(kGray + 2);
  ln->Draw("SAME");

  TLatex lat;
  lat.SetNDC(false);
  lat.SetTextFont(42);
  lat.SetTextSize(0.030);
  lat.SetTextAngle(90);
  lat.SetTextAlign(12);
  lat.DrawLatex(x, YUserAtFrac(yFrac), label);
}

void sbn_3fl_LE_overview() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const OscParams3fl p = ParamsNO(-90.0);
  Mat3 U;
  std::array<double, 3> m2;
  BuildVacuumObjects(p, m2, U);

  const int    N    = 8000;
  const double xMin = 1e-2;
  const double xMax = 3e4;

  TCanvas* c = new TCanvas("c_3fl_LE_overview", "", 950, 650);
  c->SetLogx();

  TMultiGraph* mg = new TMultiGraph();
  TGraph* g_mue  = MakeProbGraphLE(1, 0, xMin, xMax, N, true, m2, U, 0.0, false);
  TGraph* g_mumu = MakeProbGraphLE(1, 1, xMin, xMax, N, true, m2, U, 0.0, false);
  TGraph* g_mut  = MakeProbGraphLE(1, 2, xMin, xMax, N, true, m2, U, 0.0, false);

  g_mue->SetLineColor(kBlue + 1);   g_mue->SetLineWidth(3);
  g_mumu->SetLineColor(kRed + 1);   g_mumu->SetLineWidth(3);
  g_mut->SetLineColor(kGreen + 2);  g_mut->SetLineWidth(3);

  mg->Add(g_mue,  "L");
  mg->Add(g_mumu, "L");
  mg->Add(g_mut,  "L");

  mg->SetTitle(";L/E  [km/GeV];P(#nu_{#mu} #rightarrow #nu_{#alpha})");
  mg->Draw("A");
  mg->GetXaxis()->SetMoreLogLabels(false);
  mg->GetXaxis()->SetNoExponent(false);
  mg->GetXaxis()->SetNdivisions(505);
  mg->GetYaxis()->SetRangeUser(0.0, 1.0);

  TLegend* leg = new TLegend(0.56, 0.70, 0.86, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(g_mue,  "#nu_{#mu}#rightarrow#nu_{e}",   "l");
  leg->AddEntry(g_mumu, "#nu_{#mu}#rightarrow#nu_{#mu}", "l");
  leg->AddEntry(g_mut,  "#nu_{#mu}#rightarrow#nu_{#tau}","l");
  leg->Draw();

  DrawLEMarker(0.7,        "SBN (~0.5 km / 0.7 GeV)", 0.84);
  DrawLEMarker(295.0/0.60, "T2K (295/0.6)",           0.86);
  DrawLEMarker(810.0/2.0,  "NOvA (810/2.0)",          0.88);
  DrawLEMarker(1300.0/2.5, "DUNE (1300/2.5)",         0.90);

  DrawHeader("3-flavour vacuum: wide-range L/E overview",
             "Colours = final flavour; markers = typical experimental L/E (SBL and LBL)");

  c->SaveAs("sbn_3fl_LE_overview.pdf");
}

// -------------------------------------------
// Bi-probability (vacuum) for SBN detectors
// -------------------------------------------

static void BiProbPointVac(double L_km, double E_GeV, bool inverted,
                           double delta_deg, double& Pnu, double& Pnub)
{
  const double LE = (E_GeV > 0.0) ? (L_km / E_GeV) : 0.0;

  const OscParams3fl pnu = inverted ? ParamsIO(delta_deg) : ParamsNO(delta_deg);
  Mat3 U_nu;
  std::array<double, 3> m2_nu;
  BuildVacuumObjects(pnu, m2_nu, U_nu);
  Pnu = P_vac_3fl(1, 0, LE, m2_nu, U_nu);

  const OscParams3fl pnub = inverted ? ParamsIO(-delta_deg) : ParamsNO(-delta_deg);
  Mat3 U_nb;
  std::array<double, 3> m2_nb;
  BuildVacuumObjects(pnub, m2_nb, U_nb);
  Pnub = P_vac_3fl(1, 0, LE, m2_nb, U_nb);
}

static TGraph* MakeBiProbCurveVac(double L_km, double E_GeV, bool inverted,
                                  int N = 361,
                                  double dmin_deg = -180.0,
                                  double dmax_deg = +180.0)
{
  TGraph* g = new TGraph(N);
  for (int i = 0; i < N; ++i) {
    const double t     = (N == 1) ? 0.0 : (double)i / (double)(N - 1);
    const double delta = dmin_deg + (dmax_deg - dmin_deg) * t;
    double Pnu = 0.0, Pnub = 0.0;
    BiProbPointVac(L_km, E_GeV, inverted, delta, Pnu, Pnub);
    g->SetPoint(i, Pnu, Pnub);
  }
  return g;
}

static void DrawBiProbVac(const char* cname,
                          const char* outPdf,
                          const char* subtitle,
                          double L_km, double E_GeV,
                          double axisMax)
{
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  TCanvas* c = new TCanvas(cname, "", 900, 700);

  TH2D* frame = new TH2D(Form("frame_%s", cname),
                         ";P(#nu_{#mu}#rightarrow#nu_{e});P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e})",
                         10, 0.0, axisMax,
                         10, 0.0, axisMax);
  frame->Draw();
  frame->GetXaxis()->SetNdivisions(510);
  frame->GetYaxis()->SetNdivisions(510);

  TGraph* gNO = MakeBiProbCurveVac(L_km, E_GeV, false);
  TGraph* gIO = MakeBiProbCurveVac(L_km, E_GeV, true);

  gNO->SetLineColor(kBlue + 1);
  gNO->SetLineWidth(3);

  gIO->SetLineColor(kRed + 1);
  gIO->SetLineWidth(3);
  gIO->SetLineStyle(2);

  gNO->Draw("L SAME");
  gIO->Draw("L SAME");

  TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.034);
  leg->SetFillStyle(0);
  leg->AddEntry(gNO, "NO (solid): #delta_{CP} sweep", "l");
  leg->AddEntry(gIO, "IO (dashed): #delta_{CP} sweep", "l");
  leg->Draw();

  DrawHeader("Bi-probability (vacuum): P_{#mu e}(#nu) vs P_{#mu e}(#bar{#nu})",
             Form("%s   L=%.2f km,  E=%.2f GeV", subtitle, L_km, E_GeV));

  c->SaveAs(outPdf);
}

void sbn_3fl_biprob_SBND() {
  DrawBiProbVac("c_3fl_biprob_SBND",
                "sbn_3fl_biprob_SBND.pdf",
                "SBN: SBND",
                0.11, 0.70,
                1.2e-6);
}
void sbn_3fl_biprob_MicroBooNE() {
  DrawBiProbVac("c_3fl_biprob_MicroBooNE",
                "sbn_3fl_biprob_MicroBooNE.pdf",
                "SBN: MicroBooNE",
                0.47, 0.70,
                1.2e-6);
}
void sbn_3fl_biprob_ICARUS() {
  DrawBiProbVac("c_3fl_biprob_ICARUS",
                "sbn_3fl_biprob_ICARUS.pdf",
                "SBN: ICARUS",
                0.60, 0.70,
                1.2e-6);
}
void sbn_3fl_biprob_SBN() {
  sbn_3fl_biprob_SBND();
  sbn_3fl_biprob_MicroBooNE();
  sbn_3fl_biprob_ICARUS();
}

// -------------------------------------------
// 2-flavour schematic plots
// -------------------------------------------

void sbn_prob_vs_LE() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double sin22th = 0.01;
  const std::vector<double> dm2 = {0.3, 1.0, 3.0};
  const std::vector<int> colors = {kBlack, kBlue + 1, kRed + 1};

  const int    N    = 1400;
  const double xMin = 0.0, xMax = 10.0;

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
  mg->GetYaxis()->SetRangeUser(0.0, 1.05 * sin22th);

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

  const double dm2     = 1.0;
  const double sin22th = 0.01;
  const std::vector<double> L = {0.30, 0.60, 1.20};
  const std::vector<int> colors = {kBlack, kBlue + 1, kRed + 1};

  const int    N    = 1400;
  const double Emin = 0.2, Emax = 3.0;

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
  mg->GetYaxis()->SetRangeUser(0.0, 1.05 * sin22th);

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

  const double dm2     = 1.0;
  const double sin22th = 0.01;

  const int    nE   = 240;
  const int    nL   = 220;
  const double Emin = 0.2, Emax = 3.0;
  const double Lmin = 0.05, Lmax = 2.0;

  TCanvas* c = new TCanvas("c_oscillogram", "", 950, 700);
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

  const double dm2     = 3.0;
  const double sin22th = 0.01;
  const double L       = 0.60;

  const std::vector<double> frac   = {0.0, 0.05, 0.15};
  const std::vector<int>    colors = {kBlack, kBlue + 1, kRed + 1};

  const int    N    = 1400;
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
  mg->GetYaxis()->SetRangeUser(0.0, 1.05 * sin22th);

  TLegend* leg = new TLegend(0.50, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (size_t i = 0; i < frac.size(); ++i) {
    leg->AddEntry(gs[i], Form("#sigma_{E}/E = %.0f%%", 100.0 * frac[i]), "l");
  }
  leg->Draw();

  DrawHeader("Energy resolution averages rapid oscillations #rightarrow damping",
             Form("Gaussian smearing in E:  L=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.3f", L, dm2, sin22th));

  c->SaveAs("sbn_smearing_damping.pdf");
}

void sbn_bias() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const double dm2     = 3.0;
  const double sin22th = 0.10;
  const double L       = 0.60;

  const std::vector<double> bias   = {0.0, +0.02, -0.02};
  const std::vector<int>    colors = {kBlack, kRed + 1, kBlue + 1};

  const int    N    = 1400;
  const double Emin = 0.2, Emax = 3.0;

  TCanvas* c = new TCanvas("c_bias", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();

  std::vector<TGraph*> gs;
  for (size_t i = 0; i < bias.size(); ++i) {
    TGraph* g = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      const double Erec  = Emin + (Emax - Emin) * j / (N - 1.0);
      const double Etrue = Erec / (1.0 + bias[i]);
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
  mg->GetYaxis()->SetRangeUser(1.0 - 1.2 * sin22th, 1.01);

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

  const double dm2     = 1.0;
  const double sin22th = 0.10;
  const double Lnear   = 0.10;
  const double Lfar    = 0.60;

  const int    N    = 2000;
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
      const double dx = E[i + 1] - E[i];
      sum += 0.5 * (y[i] + y[i + 1]) * dx;
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

  TLine* one = new TLine(Emin, 1.0, Emax, 1.0);
  one->SetLineStyle(2);
  one->SetLineWidth(2);
  one->SetLineColor(kGray + 2);
  one->Draw("SAME");

  DrawHeader("Near/Far shape ratio isolates oscillation-driven distortions",
             Form("Toy flux#timesxsec (shape-normalised):  L_{near}=%.2f km, L_{far}=%.2f km,  #Delta m^{2}=%.2f eV^{2},  sin^{2}2#theta=%.2f",
                  Lnear, Lfar, dm2, sin22th));

  c->SaveAs("sbn_nearfar_ratio.pdf");
}

void sbn_dm2_sin22_template() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  TCanvas* c = new TCanvas("c_dm2_sin22", "", 900, 650);
  c->SetLogx();
  c->SetLogy();

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

  c->SaveAs("sbn_dm2_sin22_template.pdf");
}

// -------------------------------------------
// Aggregates
// -------------------------------------------

void sbn_3fl_all() {
  sbn_fig1p9_all();
  sbn_3fl_LE_overview();
  sbn_3fl_biprob_SBN();
}

void sbn_plot_all() {
  sbn_prob_vs_LE();
  sbn_prob_vs_E();
  sbn_oscillogram();
  sbn_smear();
  sbn_bias();
  sbn_nearfar();
  sbn_dm2_sin22_template();
  sbn_3fl_all();
}

// -------------------------------------------
// Entry point (matches filename)
// -------------------------------------------

void plotOscillation(const char* which = "all") {
  heron::macro::run_with_guard("plotOscillation", [&]() {
  std::string w(which ? which : "all");

  if      (w == "prob_le")            sbn_prob_vs_LE();
  else if (w == "prob_E")             sbn_prob_vs_E();
  else if (w == "oscillogram")        sbn_oscillogram();
  else if (w == "smear")              sbn_smear();
  else if (w == "bias")               sbn_bias();
  else if (w == "nearfar")            sbn_nearfar();
  else if (w == "dm2_sin22_template") sbn_dm2_sin22_template();

  else if (w == "fig1p9_a")           sbn_fig1p9_a();
  else if (w == "fig1p9_b")           sbn_fig1p9_b();
  else if (w == "fig1p9_c")           sbn_fig1p9_c();
  else if (w == "fig1p9_d")           sbn_fig1p9_d();
  else if (w == "fig1p9_e")           sbn_fig1p9_e();
  else if (w == "fig1p9_f")           sbn_fig1p9_f();
  else if (w == "fig1p9_all")         sbn_fig1p9_all();

  else if (w == "3fl_LE_overview")        sbn_3fl_LE_overview();
  else if (w == "3fl_biprob_SBND")        sbn_3fl_biprob_SBND();
  else if (w == "3fl_biprob_MicroBooNE")  sbn_3fl_biprob_MicroBooNE();
  else if (w == "3fl_biprob_ICARUS")      sbn_3fl_biprob_ICARUS();
  else if (w == "3fl_biprob_SBN")         sbn_3fl_biprob_SBN();
  else if (w == "3fl_all")                sbn_3fl_all();

  else if (w == "all")                    sbn_plot_all();
  else {
    std::cout << "Unknown mode: " << w << "\n"
              << "Options: prob_le, prob_E, oscillogram, smear, bias, nearfar, dm2_sin22_template, "
              << "fig1p9_a, fig1p9_b, fig1p9_c, fig1p9_d, fig1p9_e, fig1p9_f, fig1p9_all, "
              << "3fl_LE_overview, 3fl_biprob_SBND, 3fl_biprob_MicroBooNE, 3fl_biprob_ICARUS, 3fl_biprob_SBN, 3fl_all, "
              << "all\n";
  }

  });
}
