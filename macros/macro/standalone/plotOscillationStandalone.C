// plotOscillation.C
// Schematic short-baseline oscillation plots (ROOT macro).
// Supersedes older standalone oscillation plotting macros.
//
// Usage examples (batch mode):
//   root -l -q 'plotOscillation.C("prob_le")'
//   root -l -q 'plotOscillation.C("prob_E")'
//   root -l -q 'plotOscillation.C("oscillogram")'
//   root -l -q 'plotOscillation.C("smear")'
//   root -l -q 'plotOscillation.C("bias")'
//   root -l -q 'plotOscillation.C("nearfar")'
//   root -l -q 'plotOscillation.C("dm2_sin22_template")'
//   root -l -q 'plotOscillation.C("fig1p9_a")' // 3-flavour SBL: P(νμ→νe) and (1-P(νμ→νμ)) vs L/E (separate PDFs)
//   root -l -q 'plotOscillation.C("fig1p9_b")'
//   root -l -q 'plotOscillation.C("fig1p9_c")'
//   root -l -q 'plotOscillation.C("fig1p9_d")'
//   root -l -q 'plotOscillation.C("fig1p9_e")' // log y (SBL)
//   root -l -q 'plotOscillation.C("fig1p9_f")' // log x + log y (SBL)
//   root -l -q 'plotOscillation.C("fig1p9_all")'
//   root -l -q 'plotOscillation.C("3fl_LE_overview")'            // wide-range L/E with experiment markers
//   root -l -q 'plotOscillation.C("3fl_LE_long_mu_channels")'     // long-baseline regime: P(νμ→νe,νμ,ντ) vs L/E
//   root -l -q 'plotOscillation.C("3fl_LE_ordering_vac")'         // NO vs IO comparison (vacuum)
//   root -l -q 'plotOscillation.C("3fl_T2K_app_CP")'              // Pμe(E) at T2K baseline, δCP scan (ν and ν̄)
//   root -l -q 'plotOscillation.C("3fl_NOvA_app_CP")'             // Pμe(E) at NOvA baseline, δCP scan (ν and ν̄)
//   root -l -q 'plotOscillation.C("3fl_DUNE_app_CP")'             // Pμe(E) at DUNE baseline, δCP scan (ν and ν̄)
//   root -l -q 'plotOscillation.C("3fl_DUNE_matter_ordering")'     // constant-density matter: ordering + ν/ν̄
//   root -l -q 'plotOscillation.C("3fl_lbl_oscillogram")'          // 3-flavour (vacuum) oscillogram (Pμe vs L,E)
//   root -l -q 'plotOscillation.C("3fl_biprob_DUNE")'              // bi-probability (CP ellipse): DUNE-like
//   root -l -q 'plotOscillation.C("3fl_biprob_SBND")'              // bi-probability: SBND (SBN)
//   root -l -q 'plotOscillation.C("3fl_biprob_MicroBooNE")'         // bi-probability: MicroBooNE (SBN)
//   root -l -q 'plotOscillation.C("3fl_biprob_ICARUS")'            // bi-probability: ICARUS (SBN)
//   root -l -q 'plotOscillation.C("3fl_biprob_SBN")'               // all SBN bi-probability plots
//   root -l -q 'plotOscillation.C("3fl_all")'
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
#include "TColor.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "include/MacroGuard.hh"

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

// Phase factor for exp(-i m^2 L / (2E)) with m^2[eV^2], L[km], E[GeV]
static double phase_mass2(double m2_eV2, double L_km, double E_GeV) { return 2.534 * m2_eV2 * L_km / E_GeV; }

// ---------------------------
// 3-flavour vacuum oscillations helper code
//   P_{α→β} = δ_{αβ}
//            -4 Σ_{i>j} Re(U_{αi}U*_{βi}U*_{αj}U_{βj}) sin^2(Δij)
//            +2 Σ_{i>j} Im(...) sin(2Δij),
// where Δij = 1.267 Δm^2_ij [eV^2] * (L/E) [km/GeV]
// ---------------------------
using cplx = std::complex<double>;
using Mat3 = std::array<std::array<cplx,3>,3>;

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
      U[r][c] = (r == c) ? cplx(1.0,0.0) : cplx(0.0,0.0);
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

  // Build the 2x2 block action on rows i and j of U: U <- R * U
  for (int col = 0; col < 3; ++col) {
    const cplx Ui = U[i][col];
    const cplx Uj = U[j][col];
    U[i][col] =  c*Ui + s*e_m*Uj;
    U[j][col] = -s*e_p*Ui + c*Uj;
  }
}

struct OscParams3fl {
  // Standard 3-flavour vacuum parameters (NO defaults; tweak as desired)
  double th12 = deg2rad(33.44);
  double th13 = deg2rad(8.57);
  double th23 = deg2rad(49.20);
  double d13  = deg2rad(195.0);
  double dm21 = 7.42e-5;   // eV^2
  double dm31 = 2.517e-3;  // eV^2 (NO: m3^2 - m1^2)
};

enum class MassOrdering { kNO, kIO };

static OscParams3fl ParamsNO(double delta_deg = -90.0) {
  OscParams3fl p;
  p.d13 = deg2rad(delta_deg);
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

static void BuildVacuumObjects(const OscParams3fl& p, std::array<double,3>& m2, Mat3& U) {
  U  = BuildPMNS_3fl(p);
  m2 = {0.0, p.dm21, p.dm31};
}

// Flavour indices: 0=e, 1=mu, 2=tau
static double P_vac_3fl(int alpha, int beta, double LE_km_per_GeV,
                        const std::array<double,3>& m2, const Mat3& U)
{
  double P = (alpha == beta) ? 1.0 : 0.0;
  for (int i = 0; i < 3; ++i) {
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

// Constant-density matter potential a = 2*sqrt(2)*G_F*N_e*E in eV^2:
//   a[eV^2] ≈ 1.52e-4 * rho[g/cm^3] * Ye * E[GeV]
static double MatterA_eV2(double E_GeV, double rho_gcm3, double Ye = 0.5) {
  return 1.52e-4 * rho_gcm3 * Ye * E_GeV;
}

// 3-flavour constant-density matter probability (δCP forced to 0 so U is real and we can use a symmetric diagonaliser).
// This is intended for *illustrative* long-baseline phenomenology (ordering + ν/ν̄ matter asymmetry).
static double P_matter_3fl_const(int alpha, int beta,
                                 double L_km, double E_GeV,
                                 const OscParams3fl& pin,
                                 double rho_gcm3, double Ye,
                                 bool antineutrino)
{
  if (E_GeV <= 0.0) return (alpha == beta) ? 1.0 : 0.0;

  // Force δ=0 for a real mixing matrix (keeps the matter treatment simple + robust in a ROOT macro).
  OscParams3fl p = pin;
  p.d13 = 0.0;

  Mat3 Uc;
  std::array<double,3> m2;
  BuildVacuumObjects(p, m2, Uc);

  // Build flavor-basis (mass^2) matrix: M = U diag(m^2) U^T  (U real when δ=0)
  TMatrixDSym M(3);
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b <= a; ++b) {
      double sum = 0.0;
      for (int i = 0; i < 3; ++i) sum += std::real(Uc[a][i]) * std::real(Uc[b][i]) * m2[i];
      M(a,b) = sum;
      M(b,a) = sum;
    }
  }

  const double A = MatterA_eV2(E_GeV, rho_gcm3, Ye) * (antineutrino ? -1.0 : +1.0);
  M(0,0) += A;

  // Diagonalise the symmetric matrix M = V diag(λ) V^T
  TMatrixDSymEigen eig(M);
  const TVectorD eval = eig.GetEigenValues();
  const TMatrixD V    = eig.GetEigenVectors(); // columns are eigenvectors in flavour basis

  // Evolution operator: S = V diag(exp(-i φ_k)) V^T, with φ_k = (m^2_eff)_k * L/(2E)
  const double L_over_E = L_km / E_GeV;
  cplx amp(0.0, 0.0);
  for (int k = 0; k < 3; ++k) {
    const double phi = 2.534 * eval[k] * L_over_E;
    const cplx ex    = std::exp(cplx(0.0, -phi));
    amp += V(beta, k) * V(alpha, k) * ex;
  }
  return clamp01(std::norm(amp));
}

static TGraph* MakeProbGraphLE(int alpha, int beta,
                               double xMin, double xMax, int N,
                               bool logX,
                               const std::array<double,3>& m2, const Mat3& U,
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
    if (yFloor > 0.0 && y < yFloor) y = yFloor; // avoid y<=0 points on log-y pads
    g->SetPoint(k, x, y);
  }
  return g;
}

static TGraph* MakeProbGraphE_vac(int alpha, int beta,
                                  double L_km,
                                  double Emin, double Emax, int N,
                                  const std::array<double,3>& m2, const Mat3& U,
                                  double yFloor = 0.0,
                                  bool oneMinus = false)
{
  TGraph* g = new TGraph(N);
  for (int i = 0; i < N; ++i) {
    const double t = (N == 1) ? 0.0 : (double)i / (double)(N - 1);
    const double E = Emin + (Emax - Emin) * t;
    const double LE = (E > 0.0) ? (L_km / E) : 0.0;
    const double P = P_vac_3fl(alpha, beta, LE, m2, U);
    double y = oneMinus ? clamp01(1.0 - P) : P;
    if (yFloor > 0.0 && y < yFloor) y = yFloor;
    g->SetPoint(i, E, y);
  }
  return g;
}

static TGraph* MakeProbGraphE_matter(int alpha, int beta,
                                     double L_km,
                                     double Emin, double Emax, int N,
                                     const OscParams3fl& p,
                                     double rho_gcm3, double Ye,
                                     bool antineutrino,
                                     double yFloor = 0.0,
                                     bool oneMinus = false)
{
  TGraph* g = new TGraph(N);
  for (int i = 0; i < N; ++i) {
    const double t = (N == 1) ? 0.0 : (double)i / (double)(N - 1);
    const double E = Emin + (Emax - Emin) * t;
    const double P = P_matter_3fl_const(alpha, beta, L_km, E, p, rho_gcm3, Ye, antineutrino);
    double y = oneMinus ? clamp01(1.0 - P) : P;
    if (yFloor > 0.0 && y < yFloor) y = yFloor;
    g->SetPoint(i, E, y);
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
  // Smoother 2D colour plots (oscillograms).
  gStyle->SetNumberContours(60);
  // Help avoid unreadable "0" labels when axes span many decades.
  TGaxis::SetMaxDigits(4);
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

// 3-flavour SBL-focused plot as a standalone PDF (same styling as the other single-canvas outputs).
// Plots (μ-flavour initial state):
//   blue : P(νμ→νe)
//   green: P(νμ→ντ)
//   red  : 1 - P(νμ→νμ)   (disappearance probability, equals Pμe+Pμτ by unitarity)
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

  // SBL: ordering/δ choices barely matter numerically; keep a simple NO benchmark.
  const OscParams3fl p = ParamsNO(-90.0);
  Mat3 U;
  std::array<double,3> m2;
  BuildVacuumObjects(p, m2, U);

  // For log-y pads: ensure we never feed y<=0 points to TGraph.
  const double yFloor = (logY ? std::max(1e-14, 0.1*yMin) : 0.0);

  TCanvas* c = new TCanvas(cname, "", 900, 650);
  c->SetLogx(logX);
  c->SetLogy(logY);

  TMultiGraph* mg = new TMultiGraph();

  // Curves: Pμe, Pμτ, and 1-Pμμ
  TGraph* g_mue  = MakeProbGraphLE(1, 0, xMin, xMax, N, logX, m2, U, yFloor, false);
  TGraph* g_mut  = MakeProbGraphLE(1, 2, xMin, xMax, N, logX, m2, U, yFloor, false);
  TGraph* g_dis  = MakeProbGraphLE(1, 1, xMin, xMax, N, logX, m2, U, yFloor, true);

  g_mue->SetLineColor(kBlue+1);
  g_mue->SetLineWidth(3);

  g_mut->SetLineColor(kGreen+2);
  g_mut->SetLineWidth(3);

  g_dis->SetLineColor(kRed+1);
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
    // For tiny probabilities, forcing fixed-point labels produces "0" due to rounding.
    // Use scientific/exponent-style labels on log-y.
    mg->GetYaxis()->SetMoreLogLabels(false);
    mg->GetYaxis()->SetNoExponent(false);
    mg->GetYaxis()->SetTitleOffset(1.55);
  }

  // Legend placement similar to other single-canvas plots.
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
                  p.dm31, p.th13*rad2deg, p.th23*rad2deg));

  c->SaveAs(outPdf);
}

// SBL-focused ranges (separate PDFs).
// Typical accelerator SBL (SBN-like): L ~ 0.1–1 km, E ~ 0.2–3 GeV  =>  L/E ~ O(0.05–5) km/GeV.
// We extend to 10 km/GeV to show the (still small) rise of standard 3-flavour effects.
void sbn_fig1p9_a() { sbn_fig1p9_panel("c_fig1p9_a", "sbn_3fl_sbl_LE_0to10.pdf",        0.0, 10.0, 0.0,   1.2e-3, 5000, false, false); }
void sbn_fig1p9_b() { sbn_fig1p9_panel("c_fig1p9_b", "sbn_3fl_sbl_LE_0to5.pdf",         0.0,  5.0, 0.0,   3.2e-4, 5000, false, false); }
void sbn_fig1p9_c() { sbn_fig1p9_panel("c_fig1p9_c", "sbn_3fl_sbl_LE_0to2.pdf",         0.0,  2.0, 0.0,   6.0e-5, 5000, false, false); }
void sbn_fig1p9_d() { sbn_fig1p9_panel("c_fig1p9_d", "sbn_3fl_sbl_LE_0to1.pdf",         0.0,  1.0, 0.0,   1.6e-5, 5000, false, false); }
void sbn_fig1p9_e() { sbn_fig1p9_panel("c_fig1p9_e", "sbn_3fl_sbl_LE_logy.pdf",         0.05, 10.0, 1e-10, 2.0e-3, 6000, false, true ); }
void sbn_fig1p9_f() { sbn_fig1p9_panel("c_fig1p9_f", "sbn_3fl_sbl_LE_logx_logy.pdf",    0.05, 10.0, 1e-10, 2.0e-3, 6000, true,  true ); }

void sbn_fig1p9_all() {
  sbn_fig1p9_a();
  sbn_fig1p9_b();
  sbn_fig1p9_c();
  sbn_fig1p9_d();
  sbn_fig1p9_e();
  sbn_fig1p9_f();
}

// ------------------------------------------------------------
// More informative 3-flavour phenomenology plots
// ------------------------------------------------------------

static double YUserAtFrac(double frac) {
  // frac in [0,1] of the visible y-range; handles log-y and linear-y pads.
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
  // Draw vertical marker and label in *user coordinates* (robust for log-x).
  gPad->Update();
  const double y1 = gPad->GetUymin();
  const double y2 = gPad->GetUymax();

  TLine* ln = new TLine(x, y1, x, y2);
  ln->SetLineStyle(2);
  ln->SetLineWidth(2);
  ln->SetLineColor(kGray+2);
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

  // Wide-range L/E overview: show both the SBL regime and the atmospheric/solar-scale structure.
  const OscParams3fl p = ParamsNO(-90.0);
  Mat3 U;
  std::array<double,3> m2;
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

  g_mue->SetLineColor(kBlue+1);   g_mue->SetLineWidth(3);
  g_mumu->SetLineColor(kRed+1);   g_mumu->SetLineWidth(3);
  g_mut->SetLineColor(kGreen+2);  g_mut->SetLineWidth(3);

  mg->Add(g_mue,  "L");
  mg->Add(g_mumu, "L");
  mg->Add(g_mut,  "L");

  mg->SetTitle(";L/E  [km/GeV];P(#nu_{#mu} #rightarrow #nu_{#alpha})");
  mg->Draw("A");
  // Very wide log-x ranges become unreadable if minor labels are printed as decimals.
  // Prefer decade-style labels (10^{n}) here.
  mg->GetXaxis()->SetMoreLogLabels(false);
  mg->GetXaxis()->SetNoExponent(false);
  mg->GetXaxis()->SetNdivisions(505);
  mg->GetYaxis()->SetRangeUser(0.0, 1.0);

  // Pull legend left a bit to avoid any clipping on narrow canvases/PDF viewers.
  TLegend* leg = new TLegend(0.56, 0.70, 0.86, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(g_mue,  "#nu_{#mu}#rightarrow#nu_{e}",   "l");
  leg->AddEntry(g_mumu, "#nu_{#mu}#rightarrow#nu_{#mu}", "l");
  leg->AddEntry(g_mut,  "#nu_{#mu}#rightarrow#nu_{#tau}","l");
  leg->Draw();

  // Context markers (typical L/E for well-known setups)
  // SBN-like: L~0.5 km, E~0.7 GeV -> ~0.7
  DrawLEMarker(0.7,        "SBN (~0.5 km / 0.7 GeV)", 0.84);
  // T2K: 295 km / 0.6 GeV -> ~492
  DrawLEMarker(295.0/0.60, "T2K (295/0.6)", 0.86);
  // NOvA: 810 km / 2.0 GeV -> 405
  DrawLEMarker(810.0/2.0,  "NOvA (810/2.0)", 0.88);
  // DUNE: 1300 km / 2.5 GeV -> 520
  DrawLEMarker(1300.0/2.5, "DUNE (1300/2.5)", 0.90);

  DrawHeader("3-flavour vacuum: wide-range L/E overview",
             "Colours = final flavour; markers = typical experimental L/E (SBL and LBL)");

  c->SaveAs("sbn_3fl_LE_overview.pdf");
}

void sbn_3fl_LE_long_mu_channels() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Long-baseline regime in L/E: show full three-channel unitarity dance.
  const OscParams3fl p = ParamsNO(-90.0);
  Mat3 U;
  std::array<double,3> m2;
  BuildVacuumObjects(p, m2, U);

  const int    N    = 20000;
  const double xMin = 0.0;
  const double xMax = 2000.0; // covers the atmospheric-scale maxima around ~500 km/GeV

  TCanvas* c = new TCanvas("c_3fl_LE_long_mu_channels", "", 950, 650);
  TMultiGraph* mg = new TMultiGraph();

  TGraph* g_mue  = MakeProbGraphLE(1, 0, xMin, xMax, N, false, m2, U, 0.0, false);
  TGraph* g_mumu = MakeProbGraphLE(1, 1, xMin, xMax, N, false, m2, U, 0.0, false);
  TGraph* g_mut  = MakeProbGraphLE(1, 2, xMin, xMax, N, false, m2, U, 0.0, false);

  g_mue->SetLineColor(kBlue+1);   g_mue->SetLineWidth(3);
  g_mumu->SetLineColor(kRed+1);   g_mumu->SetLineWidth(3);
  g_mut->SetLineColor(kGreen+2);  g_mut->SetLineWidth(3);

  mg->Add(g_mue,  "L");
  mg->Add(g_mumu, "L");
  mg->Add(g_mut,  "L");

  mg->SetTitle(";L/E  [km/GeV];P(#nu_{#mu} #rightarrow #nu_{#alpha})");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 1.0);

  TLegend* leg = new TLegend(0.58, 0.70, 0.90, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(g_mue,  "#nu_{#mu}#rightarrow#nu_{e}",   "l");
  leg->AddEntry(g_mumu, "#nu_{#mu}#rightarrow#nu_{#mu}", "l");
  leg->AddEntry(g_mut,  "#nu_{#mu}#rightarrow#nu_{#tau}","l");
  leg->Draw();

  // Mark the atmospheric first maximum in L/E for |Δm^2_31|
  const double LE_atm1 = (0.5 * std::acos(-1.0)) / (1.267 * std::abs(p.dm31)); // π/2 / (1.267*dm31)
  TLine* lmax = new TLine(LE_atm1, 0.0, LE_atm1, 1.0);
  lmax->SetLineStyle(2);
  lmax->SetLineWidth(2);
  lmax->SetLineColor(kGray+2);
  lmax->Draw("SAME");
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.032);
  lat.DrawLatex(0.16, 0.84, Form("Atmospheric 1st max: L/E #approx %.0f km/GeV", LE_atm1));

  DrawHeader("3-flavour vacuum: long-baseline regime in L/E (NO, #delta_{CP}=-90^{#circ})",
             "Shows the full μ-flavour unitarity (Pμe+Pμμ+Pμτ=1)");

  c->SaveAs("sbn_3fl_LE_long_mu_channels.pdf");
}

void sbn_3fl_LE_ordering_vac() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Vacuum ordering comparison (illustrative; in matter the separation becomes much stronger).
  const OscParams3fl pNO = ParamsNO(-90.0);
  const OscParams3fl pIO = ParamsIO(-90.0);

  Mat3 U_NO, U_IO;
  std::array<double,3> m2_NO, m2_IO;
  BuildVacuumObjects(pNO, m2_NO, U_NO);
  BuildVacuumObjects(pIO, m2_IO, U_IO);

  const int    N    = 20000;
  const double xMin = 0.0;
  const double xMax = 2000.0;

  TCanvas* c = new TCanvas("c_3fl_LE_ordering_vac", "", 950, 650);
  TMultiGraph* mg = new TMultiGraph();

  // Focus on the appearance channel where hierarchy/δ interplay is most visible.
  TGraph* gNO = MakeProbGraphLE(1, 0, xMin, xMax, N, false, m2_NO, U_NO, 0.0, false);
  TGraph* gIO = MakeProbGraphLE(1, 0, xMin, xMax, N, false, m2_IO, U_IO, 0.0, false);

  gNO->SetLineColor(kBlue+1); gNO->SetLineWidth(3); gNO->SetLineStyle(1);
  gIO->SetLineColor(kRed+1);  gIO->SetLineWidth(3); gIO->SetLineStyle(2);

  mg->Add(gNO, "L");
  mg->Add(gIO, "L");

  mg->SetTitle(";L/E  [km/GeV];P_{#mu e}");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 0.12);

  // Shift left so legend text cannot clip on the right border.
  TLegend* leg = new TLegend(0.50, 0.74, 0.86, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(gNO, "Normal ordering (NO), #delta_{CP}=-90^{#circ}", "l");
  leg->AddEntry(gIO, "Inverted ordering (IO), #delta_{CP}=-90^{#circ}", "l");
  leg->Draw();

  DrawHeader("3-flavour vacuum: ordering comparison in P_{#mu e}",
             "Line style encodes ordering (solid=NO, dashed=IO); matter effects amplify separation");

  c->SaveAs("sbn_3fl_LE_ordering_vac.pdf");
}

static void StyleDeltaCurves(std::vector<TGraph*>& gs, const std::vector<int>& cols, bool dashed = false) {
  for (size_t i = 0; i < gs.size(); ++i) {
    gs[i]->SetLineColor(cols[i]);
    gs[i]->SetLineWidth(3);
    gs[i]->SetLineStyle(dashed ? 2 : 1);
  }
}

static void PlotAppCP_vsE(const char* cname, const char* outPdf,
                          double L_km, double Emin, double Emax,
                          const char* titleLine,
                          bool showAntinu)
{
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Use vacuum here to isolate “geometry”: oscillation peaks and CP interference in a clean way.
  // (Matter effects are handled in a separate dedicated plot below.)
  const std::vector<double> deltas = {-90.0, 0.0, +90.0};
  const std::vector<int>    cols   = {kBlue+1, kBlack, kRed+1};

  const int N = 2000;

  TCanvas* c = new TCanvas(cname, "", 950, 650);
  TMultiGraph* mg = new TMultiGraph();

  std::vector<TGraph*> g_nu;
  std::vector<TGraph*> g_nub;
  g_nu.reserve(deltas.size());
  g_nub.reserve(deltas.size());

  for (size_t i = 0; i < deltas.size(); ++i) {
    const OscParams3fl pnu  = ParamsNO(deltas[i]);
    const OscParams3fl pnub = ParamsNO(-deltas[i]); // vacuum: ν̄ equivalent to δ -> -δ

    Mat3 U_nu, U_nub;
    std::array<double,3> m2_nu, m2_nub;
    BuildVacuumObjects(pnu,  m2_nu,  U_nu);
    BuildVacuumObjects(pnub, m2_nub, U_nub);

    TGraph* g1 = MakeProbGraphE_vac(1, 0, L_km, Emin, Emax, N, m2_nu,  U_nu);
    g_nu.push_back(g1);
    mg->Add(g1, "L");

    if (showAntinu) {
      TGraph* g2 = MakeProbGraphE_vac(1, 0, L_km, Emin, Emax, N, m2_nub, U_nub);
      g_nub.push_back(g2);
      mg->Add(g2, "L");
    }
  }

  StyleDeltaCurves(g_nu, cols, false);
  if (showAntinu) StyleDeltaCurves(g_nub, cols, true);

  mg->SetTitle(";E_{#nu}  [GeV];P_{#mu e}");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 0.20);

  TLegend* leg = new TLegend(0.52, 0.62, 0.90, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.033);
  leg->SetFillStyle(0);
  leg->AddEntry(g_nu[0],  "#nu,  #delta_{CP}=-90^{#circ}", "l");
  leg->AddEntry(g_nu[1],  "#nu,  #delta_{CP}=0",           "l");
  leg->AddEntry(g_nu[2],  "#nu,  #delta_{CP}=+90^{#circ}", "l");
  if (showAntinu) {
    leg->AddEntry(g_nub[0], "#bar{#nu},  #delta_{CP}=-90^{#circ}  (dashed)", "l");
  }
  leg->Draw();

  DrawHeader(titleLine,
             showAntinu ? "Vacuum; dashed = #bar{#nu} (implemented as #delta_{CP}#rightarrow-#delta_{CP})"
                        : "Vacuum; shows CP-phase modulation of the appearance peak structure");

  c->SaveAs(outPdf);
}

void sbn_3fl_T2K_app_CP()  { PlotAppCP_vsE("c_3fl_T2K_app_CP",  "sbn_3fl_T2K_app_CP.pdf",  295.0, 0.2, 1.5, "3-flavour vacuum: T2K-like baseline (L=295 km) — P_{#mu e}(E) vs #delta_{CP}", true); }
void sbn_3fl_NOvA_app_CP() { PlotAppCP_vsE("c_3fl_NOvA_app_CP", "sbn_3fl_NOvA_app_CP.pdf", 810.0, 0.5, 4.0, "3-flavour vacuum: NOvA-like baseline (L=810 km) — P_{#mu e}(E) vs #delta_{CP}", true); }
void sbn_3fl_DUNE_app_CP() { PlotAppCP_vsE("c_3fl_DUNE_app_CP", "sbn_3fl_DUNE_app_CP.pdf", 1300.0,0.5, 6.0, "3-flavour vacuum: DUNE-like baseline (L=1300 km) — P_{#mu e}(E) vs #delta_{CP}", true); }

void sbn_3fl_DUNE_matter_ordering() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Constant-density matter comparison at DUNE baseline.
  // Here we *intentionally* set δ=0 in the matter calculation to keep the treatment exact and stable
  // with a real-symmetric diagonalisation. This isolates the ordering + ν/ν̄ matter asymmetry.
  const double L_km = 1300.0;
  const double Emin = 0.5;
  const double Emax = 6.0;
  const int    N    = 2000;

  const double rho = 2.8; // g/cm^3 (schematic mantle/crust scale)
  const double Ye  = 0.5;

  OscParams3fl pNO = ParamsNO(0.0);
  OscParams3fl pIO = ParamsIO(0.0);

  TCanvas* c = new TCanvas("c_3fl_DUNE_matter_ordering", "", 950, 650);
  TMultiGraph* mg = new TMultiGraph();

  // Colour encodes ordering; line style encodes ν vs ν̄.
  TGraph* gNO_nu  = MakeProbGraphE_matter(1, 0, L_km, Emin, Emax, N, pNO, rho, Ye, false);
  TGraph* gNO_nub = MakeProbGraphE_matter(1, 0, L_km, Emin, Emax, N, pNO, rho, Ye, true );
  TGraph* gIO_nu  = MakeProbGraphE_matter(1, 0, L_km, Emin, Emax, N, pIO, rho, Ye, false);
  TGraph* gIO_nub = MakeProbGraphE_matter(1, 0, L_km, Emin, Emax, N, pIO, rho, Ye, true );

  gNO_nu ->SetLineColor(kBlue+1); gNO_nu ->SetLineWidth(3); gNO_nu ->SetLineStyle(1);
  gNO_nub->SetLineColor(kBlue+1); gNO_nub->SetLineWidth(3); gNO_nub->SetLineStyle(2);
  gIO_nu ->SetLineColor(kRed+1);  gIO_nu ->SetLineWidth(3); gIO_nu ->SetLineStyle(1);
  gIO_nub->SetLineColor(kRed+1);  gIO_nub->SetLineWidth(3); gIO_nub->SetLineStyle(2);

  mg->Add(gNO_nu,  "L");
  mg->Add(gNO_nub, "L");
  mg->Add(gIO_nu,  "L");
  mg->Add(gIO_nub, "L");

  mg->SetTitle(";E_{#nu}  [GeV];P_{#mu e}");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 0.25);

  TLegend* leg = new TLegend(0.50, 0.62, 0.90, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.033);
  leg->SetFillStyle(0);
  leg->AddEntry(gNO_nu,  "NO: #nu (solid),  #delta_{CP}=0", "l");
  leg->AddEntry(gNO_nub, "NO: #bar{#nu} (dashed)",          "l");
  leg->AddEntry(gIO_nu,  "IO: #nu (solid),  #delta_{CP}=0", "l");
  leg->AddEntry(gIO_nub, "IO: #bar{#nu} (dashed)",          "l");
  leg->Draw();

  DrawHeader("3-flavour constant-density matter: DUNE-like baseline (L=1300 km)",
             Form("rho=%.1f g/cm^{3}, Ye=%.2f; illustrates ordering + #nu/#bar{#nu} matter asymmetry (δ forced to 0 here)", rho, Ye));

  c->SaveAs("sbn_3fl_DUNE_matter_ordering.pdf");
}

void sbn_3fl_lbl_oscillogram() {
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  // Vacuum oscillogram for Pμe over a long-baseline-like L and E range.
  const OscParams3fl p = ParamsNO(-90.0);
  Mat3 U;
  std::array<double,3> m2;
  BuildVacuumObjects(p, m2, U);

  const int nE = 260;
  const int nL = 260;
  const double Emin = 0.2, Emax = 6.0;   // GeV
  const double Lmin = 50.0, Lmax = 2000.0; // km

  TCanvas* c = new TCanvas("c_3fl_lbl_oscillogram", "", 950, 700);
  c->SetRightMargin(0.18);

  TH2D* h = new TH2D("h_3fl_lbl_osc",
                     ";E_{#nu}  [GeV];Baseline L  [km];P_{#mu e}",
                     nE, Emin, Emax,
                     nL, Lmin, Lmax);

  for (int ix = 1; ix <= nE; ++ix) {
    const double E = h->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= nL; ++iy) {
      const double L = h->GetYaxis()->GetBinCenter(iy);
      const double LE = L / E;
      h->SetBinContent(ix, iy, P_vac_3fl(1, 0, LE, m2, U));
    }
  }

  h->SetContour(60);
  h->SetMinimum(0.0);
  h->SetMaximum(0.25);
  h->GetZaxis()->SetTitleOffset(1.25);
  h->Draw("COLZ");

  DrawHeader("3-flavour vacuum oscillogram: P_{#mu e}(L,E)",
             "Vacuum, NO, #delta_{CP}=-90^{#circ}; constant L/E structures appear as diagonal bands");

  c->SaveAs("sbn_3fl_lbl_oscillogram.pdf");
}

// -------------------------------------------------------------------
// Bi-probability (CP ellipse) plots: P(νμ→νe) vs P(ν̄μ→ν̄e), δCP swept
//   - Vacuum version (exact with complex U): ν̄ implemented as δ → -δ.
//   - Shows ordering (NO vs IO) and the δCP “ellipse geometry”.
// -------------------------------------------------------------------

static void BiProbPointVac(double L_km, double E_GeV, bool inverted,
                           double delta_deg, double& Pnu, double& Pnub)
{
  const double LE = (E_GeV > 0.0) ? (L_km / E_GeV) : 0.0;

  // ν
  const OscParams3fl pnu = inverted ? ParamsIO(delta_deg) : ParamsNO(delta_deg);
  Mat3 U_nu;
  std::array<double,3> m2_nu;
  BuildVacuumObjects(pnu, m2_nu, U_nu);
  Pnu = P_vac_3fl(1, 0, LE, m2_nu, U_nu);

  // ν̄ in vacuum: U -> U*  ⇔  δ -> -δ (same ordering)
  const OscParams3fl pnub = inverted ? ParamsIO(-delta_deg) : ParamsNO(-delta_deg);
  Mat3 U_nb;
  std::array<double,3> m2_nb;
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
    const double t = (N == 1) ? 0.0 : (double)i / (double)(N - 1);
    const double delta = dmin_deg + (dmax_deg - dmin_deg) * t;
    double Pnu = 0.0, Pnub = 0.0;
    BiProbPointVac(L_km, E_GeV, inverted, delta, Pnu, Pnub);
    g->SetPoint(i, Pnu, Pnub);
  }
  return g;
}

static TGraph* MakeBiProbMarkerVac(double L_km, double E_GeV, bool inverted,
                                   double delta_deg, int color, int mstyle)
{
  double Pnu = 0.0, Pnub = 0.0;
  BiProbPointVac(L_km, E_GeV, inverted, delta_deg, Pnu, Pnub);
  TGraph* g = new TGraph(1);
  g->SetPoint(0, Pnu, Pnub);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(mstyle);
  g->SetMarkerSize(1.2);
  return g;
}

static void DrawBiProbVac(const char* cname,
                          const char* outPdf,
                          const char* subtitle,
                          double L_km, double E_GeV,
                          double axisMax,
                          bool annotateDeltaMarkers)
{
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  TCanvas* c = new TCanvas(cname, "", 900, 700);

  // Frame
  TH2D* frame = new TH2D(Form("frame_%s", cname),
                         ";P(#nu_{#mu}#rightarrow#nu_{e});P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e})",
                         10, 0.0, axisMax,
                         10, 0.0, axisMax);
  frame->Draw();
  frame->GetXaxis()->SetNdivisions(510);
  frame->GetYaxis()->SetNdivisions(510);
  frame->GetXaxis()->SetNoExponent(false);
  frame->GetYaxis()->SetNoExponent(false);

  // Curves: ordering
  TGraph* gNO = MakeBiProbCurveVac(L_km, E_GeV, false);
  TGraph* gIO = MakeBiProbCurveVac(L_km, E_GeV, true);

  gNO->SetLineColor(kBlue+1);
  gNO->SetLineWidth(3);

  gIO->SetLineColor(kRed+1);
  gIO->SetLineWidth(3);
  gIO->SetLineStyle(2);

  gNO->Draw("L SAME");
  gIO->Draw("L SAME");

  // Optional δCP markers (use NO markers only to keep the plot readable).
  TGraph* m_d90 = nullptr;
  TGraph* m_d0  = nullptr;
  TGraph* m_p90 = nullptr;
  TGraph* m_p180 = nullptr;
  if (annotateDeltaMarkers) {
    m_d90  = MakeBiProbMarkerVac(L_km, E_GeV, false, -90.0, kBlue+1, 20); // circle
    m_d0   = MakeBiProbMarkerVac(L_km, E_GeV, false,   0.0, kBlue+1, 21); // square
    m_p90  = MakeBiProbMarkerVac(L_km, E_GeV, false, +90.0, kBlue+1, 22); // triangle
    m_p180 = MakeBiProbMarkerVac(L_km, E_GeV, false, 180.0, kBlue+1, 29); // star
    m_d90 ->Draw("P SAME");
    m_d0  ->Draw("P SAME");
    m_p90 ->Draw("P SAME");
    m_p180->Draw("P SAME");
  }

  TLegend* leg = new TLegend(0.52, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.034);
  leg->SetFillStyle(0);
  leg->AddEntry(gNO, "NO (solid): #delta_{CP} sweep", "l");
  leg->AddEntry(gIO, "IO (dashed): #delta_{CP} sweep", "l");
  if (annotateDeltaMarkers) {
    leg->AddEntry(m_d90,  "Markers (NO): #delta=-90^{#circ}", "p");
    leg->AddEntry(m_d0,   "Markers (NO): #delta=0",          "p");
    leg->AddEntry(m_p90,  "Markers (NO): #delta=+90^{#circ}","p");
    leg->AddEntry(m_p180, "Markers (NO): #delta=180^{#circ}","p");
  }
  leg->Draw();

  DrawHeader("Bi-probability (CP ellipse): P_{#mu e}(#nu) vs P_{#mu e}(#bar{#nu})",
             Form("%s   (vacuum)   L=%.2f km,  E=%.2f GeV", subtitle, L_km, E_GeV));

  c->SaveAs(outPdf);
}

// DUNE-like: choose E near the first maximum (~2.5 GeV for 1300 km).
void sbn_3fl_biprob_DUNE() {
  DrawBiProbVac("c_3fl_biprob_DUNE",
                "sbn_3fl_biprob_DUNE.pdf",
                "DUNE-like",
                1300.0, 2.50,
                0.20,
                true);
}

// SBN: BNB peak energy ~0.7 GeV; baselines from target to detectors (approx).
// Standard 3-flavour effects are tiny here; the ellipse collapses very near the origin.
void sbn_3fl_biprob_SBND() {
  DrawBiProbVac("c_3fl_biprob_SBND",
                "sbn_3fl_biprob_SBND.pdf",
                "SBN: SBND",
                0.11, 0.70,
                1.2e-6,
                false);
}

void sbn_3fl_biprob_MicroBooNE() {
  DrawBiProbVac("c_3fl_biprob_MicroBooNE",
                "sbn_3fl_biprob_MicroBooNE.pdf",
                "SBN: MicroBooNE",
                0.47, 0.70,
                1.2e-6,
                false);
}

void sbn_3fl_biprob_ICARUS() {
  DrawBiProbVac("c_3fl_biprob_ICARUS",
                "sbn_3fl_biprob_ICARUS.pdf",
                "SBN: ICARUS",
                0.60, 0.70,
                1.2e-6,
                false);
}

void sbn_3fl_biprob_SBN() {
  sbn_3fl_biprob_SBND();
  sbn_3fl_biprob_MicroBooNE();
  sbn_3fl_biprob_ICARUS();
}

void sbn_3fl_all() {
  // Keep the original “SBL zoom” set (now with Pμτ added).
  sbn_fig1p9_all();

  // Add richer, more global phenomenology views.
  sbn_3fl_LE_overview();
  sbn_3fl_LE_long_mu_channels();
  sbn_3fl_LE_ordering_vac();
  sbn_3fl_T2K_app_CP();
  sbn_3fl_NOvA_app_CP();
  sbn_3fl_DUNE_app_CP();
  sbn_3fl_DUNE_matter_ordering();
  sbn_3fl_lbl_oscillogram();

  // Advanced: CP-ellipse (bi-probability) views for DUNE and SBN.
  sbn_3fl_biprob_DUNE();
  sbn_3fl_biprob_SBN();
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
  sbn_3fl_all();
}

void sbn_osc_plots(const char* which = "all") {
  heron::macro::run_with_guard("sbn_osc_plots", [&]() {
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
  else if (w == "3fl_LE_overview")            sbn_3fl_LE_overview();
  else if (w == "3fl_LE_long_mu_channels")     sbn_3fl_LE_long_mu_channels();
  else if (w == "3fl_LE_ordering_vac")         sbn_3fl_LE_ordering_vac();
  else if (w == "3fl_T2K_app_CP")              sbn_3fl_T2K_app_CP();
  else if (w == "3fl_NOvA_app_CP")             sbn_3fl_NOvA_app_CP();
  else if (w == "3fl_DUNE_app_CP")             sbn_3fl_DUNE_app_CP();
  else if (w == "3fl_DUNE_matter_ordering")    sbn_3fl_DUNE_matter_ordering();
  else if (w == "3fl_lbl_oscillogram")         sbn_3fl_lbl_oscillogram();
  else if (w == "3fl_biprob_DUNE")             sbn_3fl_biprob_DUNE();
  else if (w == "3fl_biprob_SBND")             sbn_3fl_biprob_SBND();
  else if (w == "3fl_biprob_MicroBooNE")       sbn_3fl_biprob_MicroBooNE();
  else if (w == "3fl_biprob_ICARUS")           sbn_3fl_biprob_ICARUS();
  else if (w == "3fl_biprob_SBN")              sbn_3fl_biprob_SBN();
  else if (w == "3fl_all")                     sbn_3fl_all();
  else if (w == "all")                sbn_plot_all();
  else {
    std::cout << "Unknown mode: " << w << "\n"
              << "Options: prob_le, prob_E, oscillogram, smear, bias, nearfar, dm2_sin22_template, "
              << "fig1p9_a, fig1p9_b, fig1p9_c, fig1p9_d, fig1p9_e, fig1p9_f, fig1p9_all, "
              << "3fl_LE_overview, 3fl_LE_long_mu_channels, 3fl_LE_ordering_vac, "
              << "3fl_T2K_app_CP, 3fl_NOvA_app_CP, 3fl_DUNE_app_CP, 3fl_DUNE_matter_ordering, 3fl_lbl_oscillogram, 3fl_all, "
              << "3fl_biprob_DUNE, 3fl_biprob_SBND, 3fl_biprob_MicroBooNE, 3fl_biprob_ICARUS, 3fl_biprob_SBN, "
              << "all\n";
  }

  });
}
