// LArTPC_PhysicsPlots.C
//
// ROOT macro: schematic “detector physics” plots for LArTPCs.
//
// Combines a few Geant4 Physics Reference Manual–style formulas already in your macro
// (restricted Bethe–Bloch, Bohr straggling, Highland multiple scattering)
// with LArTPC-specific charge transport / response relations that you typically need
// to go from Edep -> collected charge vs drift time:
//
//  - electron drift velocity v_d(E,T) (parameterization copied from LArSoft DetectorPropertiesStandard)
//  - recombination (Birks / Modified Box; form copied from LArSoft ISCalcSeparate)
//  - attachment / lifetime attenuation: Q(t)=Q0 exp(-t/tau_e)
//  - diffusion: sigma^2 = 2 D t (longitudinal/transverse)
//  - Shockley–Ramo (toy 1D parallel-plate) + simple CR-RC electronics shaping
//
// Cherenkov radiation intentionally omitted.
//
// Usage (batch):
//   root -l -q 'LArTPC_PhysicsPlots.C("all","LArTPC")'
//   root -l -q 'LArTPC_PhysicsPlots.C("drift","LArTPC")'
//
// Modes: all, dedx, straggle, mcs, drift, recomb, lifetime, diffusion, ramo, scint
// Output:
//   - if which="all":  <outPrefix>_larTPC_physics.pdf (multi-page)
//   - else:            <outPrefix>_<which>.pdf
//
// ------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include "TString.h"

// --------------------------
// Style (match sbn_osc_plots)
// --------------------------
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
}

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

static void StyleGraph(TGraph* g, int color, int style = 1, int width = 3) {
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}

// --------------------------
// Material + constants
// --------------------------
struct Material {
  double Z;        // atomic number
  double A_gmol;   // g/mol
  double rho_gcm3; // g/cm^3
  double I_MeV;    // mean excitation energy [MeV]
  double X0_cm;    // radiation length [cm]
};

static const double kRe_cm  = 2.8179403227e-13; // classical electron radius [cm]
static const double kMe_MeV = 0.51099895;       // electron mass-energy [MeV]
static const double kNa     = 6.02214076e23;    // Avogadro [1/mol]
static const double kPi     = TMath::Pi();

static const double kQe_C   = 1.602176634e-19;  // |electron charge| [C]
static const double kFC_per_e = kQe_C / 1.0e-15; // fC per electron

// Default-ish LAr numbers (edit to match your detector/temperature assumptions)
static Material MakeLAr() {
  Material m;
  m.Z        = 18.0;
  m.A_gmol   = 39.948;
  m.rho_gcm3 = 1.396;          // ~87–89 K; update if you prefer a specific density
  m.I_MeV    = 188e-6;         // 188 eV
  m.X0_cm    = 14.0;           // ~14 cm (rough)
  return m;
}

// --------------------------
// Kinematics helpers
// --------------------------
static double BetaFromP(double p_MeVc, double M_MeV) {
  const double E = std::sqrt(p_MeVc*p_MeVc + M_MeV*M_MeV);
  return (E > 0.0) ? (p_MeVc / E) : 0.0;
}

static double GammaFromP(double p_MeVc, double M_MeV) {
  const double E = std::sqrt(p_MeVc*p_MeVc + M_MeV*M_MeV);
  return (M_MeV > 0.0) ? (E / M_MeV) : 1.0;
}

static double ElectronDensity_cm3(double Z, double A_gmol, double rho_gcm3) {
  // n_el = Z * N_A * rho / A  [electrons/cm^3]
  return Z * kNa * rho_gcm3 / A_gmol;
}

// --------------------------
// Geant4 manual-style pieces
// --------------------------

// Eq. (13.1): Tmax (max kinetic energy transferable to an electron)
static double Tmax_MeV(double p_MeVc, double M_MeV) {
  const double gamma = GammaFromP(p_MeVc, M_MeV);
  const double r = kMe_MeV / M_MeV;
  const double num = 2.0 * kMe_MeV * (gamma*gamma - 1.0);
  const double den = 1.0 + 2.0*gamma*r + r*r;
  return (den > 0.0) ? (num / den) : 0.0;
}

// Eq. (13.2): restricted Bethe–Bloch dE/dx (with optional correction terms)
static double BetheBlochRestricted_MeVpercm(double p_MeVc,
                                            double M_MeV,
                                            int z,
                                            const Material &mat,
                                            double Tcut_MeV,
                                            double delta = 0.0,
                                            double Ce    = 0.0,
                                            double S     = 0.0,
                                            double F     = 0.0)
{
  const double beta  = BetaFromP(p_MeVc, M_MeV);
  const double gamma = GammaFromP(p_MeVc, M_MeV);
  if (beta <= 0.0) return 0.0;

  const double Tmax = Tmax_MeV(p_MeVc, M_MeV);
  const double Tup  = (Tcut_MeV < Tmax) ? Tcut_MeV : Tmax;

  const double n_el = ElectronDensity_cm3(mat.Z, mat.A_gmol, mat.rho_gcm3);

  // 2π r_e^2 m_e c^2 n_el  [MeV/cm]
  const double pref = 2.0 * kPi * (kRe_cm*kRe_cm) * kMe_MeV * n_el;

  const double beta2 = beta*beta;

  double arg = 2.0 * kMe_MeV * beta2 * gamma*gamma * Tup / (mat.I_MeV * mat.I_MeV);
  if (arg <= 1.0) arg = 1.0000001;

  const double bracket =
      std::log(arg)
    - beta2 * (1.0 + Tup/Tmax)
    - delta
    - 2.0 * Ce / mat.Z
    + S
    + F;

  return pref * ( (double)(z*z) / beta2 ) * bracket; // MeV/cm
}

// Eq. (8.8) (as coded in your original macro): Bohr straggling RMS for a step of length s
static double BohrSigma_MeV(double p_MeVc,
                            double M_MeV,
                            double step_cm,
                            int z,
                            const Material &mat,
                            double Tcut_MeV)
{
  const double beta = BetaFromP(p_MeVc, M_MeV);
  if (beta <= 0.0) return 0.0;
  const double beta2 = beta*beta;

  const double Tmax = Tmax_MeV(p_MeVc, M_MeV);
  const double Tc   = Tcut_MeV;

  const double n_el = ElectronDensity_cm3(mat.Z, mat.A_gmol, mat.rho_gcm3);
  const double pref = 2.0 * kPi * (kRe_cm*kRe_cm) * kMe_MeV * n_el; // MeV/cm

  double factor = 1.0 - (beta2/2.0) * (Tc / Tmax);
  if (factor < 0.0) factor = 0.0;

  double Omega2 = pref * ( (double)(z*z) / beta2 ) * Tmax * step_cm * factor;
  if (Omega2 < 0.0) Omega2 = 0.0;

  return std::sqrt(Omega2); // MeV
}

// Highland–Lynch–Dahl multiple scattering width theta0 (hc = 0.038)
static double Theta0_Highland_rad(double p_MeVc,
                                  double M_MeV,
                                  double t_cm,
                                  double X0_cm,
                                  int z,
                                  double hc = 0.038)
{
  const double beta = BetaFromP(p_MeVc, M_MeV);
  if (beta <= 0.0) return 0.0;

  const double tx0 = t_cm / X0_cm;
  if (tx0 <= 0.0) return 0.0;

  const double corr = 1.0 + hc * std::log(tx0);
  return (13.6 / (beta * p_MeVc)) * z * std::sqrt(tx0) * corr; // rad
}

// --------------------------
// LArTPC-specific pieces (LArSoft-like)
// --------------------------

// Drift velocity parameterization copied from LArSoft DetectorPropertiesStandard::DriftVelocity()
// Returns v_d in cm/us.  Efield in kV/cm, temperature in K.
// Two models:
//  - "Walkowiak + ICARUS + transition" (useIcarusMicroBooNEModel=false)
//  - "MicroBooNE+ICARUS polynomial" (useIcarusMicroBooNEModel=true)
static double DriftVelocity_LArSoft_cm_per_us(double efield_kVcm,
                                              double temperature_K,
                                              bool useIcarusMicroBooNEModel,
                                              double driftVelFudgeFactor = 1.0)
{
  if (efield_kVcm <= 0.0) return 0.0;

  double vd = 0.0;

  if (!useIcarusMicroBooNEModel) {

    const double tshift = -87.203 + temperature_K;
    const double xFit   = 0.0938163 - 0.0052563 * tshift - 0.0001470 * tshift * tshift;
    const double uFit   = 5.18406   + 0.01448  * tshift - 0.003497  * tshift * tshift
                        - 0.000516  * tshift * tshift * tshift;

    // ICARUS parameter set
    const double P1 = -0.04640;
    const double P2 =  0.01712;
    const double P3 =  1.88125;
    const double P4 =  0.99408;
    const double P5 =  0.01172;
    const double P6 =  4.20214;
    const double T0 = 105.749;

    // Walkowiak parameter set
    const double P1W = -0.01481;
    const double P2W = -0.0075;
    const double P3W =  0.141;
    const double P4W = 12.4;
    const double P5W =  1.627;
    const double P6W =  0.317;
    const double T0W = 90.371;

    auto Icarus = [&](double E) -> double {
      return ( (P1 * (temperature_K - T0) + 1.0)
               * (P3 * E * std::log(1.0 + P4 / E) + P5 * std::pow(E, P6))
               + P2 * (temperature_K - T0) );
    };

    auto Walkowiak = [&](double E) -> double {
      return ( (P1W * (temperature_K - T0W) + 1.0)
               * (P3W * E * std::log(1.0 + P4W / E) + P5W * std::pow(E, P6W))
               + P2W * (temperature_K - T0W) );
    };

    // Smooth transition (Craig Thorne note in LArSoft)
    if (efield_kVcm < xFit) {
      vd = efield_kVcm * uFit;
    }
    else if (efield_kVcm < 0.619) {
      vd = Icarus(efield_kVcm);
    }
    else if (efield_kVcm < 0.699) {
      // Linear blend between ICARUS and Walkowiak in this band:
      // weights: 12.5*(E-0.619) and 12.5*(0.699-E) -> sum to 1
      vd = 12.5 * (efield_kVcm - 0.619) * Walkowiak(efield_kVcm)
         + 12.5 * (0.699 - efield_kVcm) * Icarus(efield_kVcm);
    }
    else {
      vd = Walkowiak(efield_kVcm);
    }
  }
  else {
    // MicroBooNE+ICARUS model polynomial (arXiv:2008.09765) with temperature dependence
    const double P0 = 0.0;
    const double P1 =  5.53416;
    const double P2 = -6.53093;
    const double P3 =  3.20752;
    const double P4 =  0.389696;
    const double P5 = -0.556184;

    const double E = efield_kVcm;
    const double poly = P0
                      + P1*E
                      + P2*E*E
                      + P3*E*E*E
                      + P4*E*E*E*E
                      + P5*E*E*E*E*E;

    vd = (1.0 - 0.0184 * (temperature_K - 89.0)) * poly;
  }

  // LArSoft returns cm/us after this scaling.
  vd *= driftVelFudgeFactor / 10.0;

  return vd; // cm/us
}

// Recombination models in the same *forward* direction used in simulation:
// given dE/dx and E field, compute recombination factor R (fraction of electrons surviving recombination)
// Forms match LArSoft ISCalcSeparate.
struct RecombParams {
  double BirksA   = 0.800; // "RecombA"  (example default; experiment-specific)
  double Birksk   = 0.0486; // "Recombk"  (example default; experiment-specific)
  double ModBoxA  = 0.930; // "ModBoxA"  (example default; experiment-specific)
  double ModBoxB  = 0.212; // "ModBoxB"  (example default; experiment-specific)
};

static double Recomb_Birks(double dEdx_MeVcm,
                           double E_kVcm,
                           double rho_gcm3,
                           const RecombParams& p)
{
  if (E_kVcm <= 0.0) return 0.0;
  // LArSoft: fRecombk = Recombk / Density()
  const double k_over_rho = p.Birksk / rho_gcm3;
  return p.BirksA / (1.0 + dEdx_MeVcm * k_over_rho / E_kVcm);
}

static double Recomb_ModBox(double dEdx_MeVcm,
                            double E_kVcm,
                            double rho_gcm3,
                            const RecombParams& p)
{
  if (E_kVcm <= 0.0) return 0.0;
  // LArSoft: fModBoxB = ModBoxB / Density()
  const double B_over_rho = p.ModBoxB / rho_gcm3;
  const double Xi = B_over_rho * dEdx_MeVcm / E_kVcm;
  if (Xi <= 0.0) return 0.0;
  return std::log(p.ModBoxA + Xi) / Xi;
}

// W_ion conversion: mean ionization work function (use 23.6 eV unless you override)
static double dQdx_electrons_per_cm(double dEdx_MeVcm, double recombFactor, double Wion_eV = 23.6) {
  const double e_per_MeV = 1.0e6 / Wion_eV; // electrons per MeV (mean)
  return dEdx_MeVcm * e_per_MeV * recombFactor;
}

// Lifetime attenuation: Q(t)=Q0*exp(-t/tau)
static double SurvivalFraction(double drift_cm, double vd_cm_per_us, double tau_ms) {
  if (vd_cm_per_us <= 0.0 || tau_ms <= 0.0) return 0.0;
  const double t_us   = drift_cm / vd_cm_per_us;
  const double tau_us = tau_ms * 1000.0;
  return std::exp(-t_us / tau_us);
}

// Diffusion: sigma^2 = 2 D t
static double DiffusionSigma_cm(double drift_cm,
                                double vd_cm_per_us,
                                double D_cm2_per_s)
{
  if (vd_cm_per_us <= 0.0 || D_cm2_per_s <= 0.0) return 0.0;
  const double t_s = (drift_cm / vd_cm_per_us) * 1.0e-6; // us -> s
  return std::sqrt(2.0 * D_cm2_per_s * t_s);
}

// Shockley–Ramo (toy parallel plate): i = q * v / d  (|E_w| = 1/d)
static double RamoCurrent_nA(double Ne, double vd_cm_per_us, double gap_cm) {
  if (gap_cm <= 0.0 || vd_cm_per_us <= 0.0) return 0.0;
  const double q_C = std::abs(Ne) * kQe_C;
  const double v_cm_per_s = vd_cm_per_us * 1.0e6;
  const double i_A = q_C * (v_cm_per_s / gap_cm);
  return i_A * 1.0e9; // nA
}

// Simple CR-RC impulse response (unit area not enforced): h(t) = (t/tau^2) exp(-t/tau) for t>=0
static std::vector<double> CRRCImpulse(const std::vector<double>& t_us, double tau_us) {
  std::vector<double> h(t_us.size(), 0.0);
  if (tau_us <= 0.0) return h;
  for (size_t i = 0; i < t_us.size(); ++i) {
    const double t = t_us[i];
    if (t < 0.0) { h[i] = 0.0; continue; }
    h[i] = (t / (tau_us * tau_us)) * std::exp(-t / tau_us);
  }
  return h;
}

// Naive O(N^2) discrete convolution y(t)=∫ x(t') h(t-t') dt' (Riemann sum)
static std::vector<double> ConvolveSameSize(const std::vector<double>& x,
                                            const std::vector<double>& h,
                                            double dt_us)
{
  const size_t N = x.size();
  std::vector<double> y(N, 0.0);
  for (size_t i = 0; i < N; ++i) {
    double sum = 0.0;
    for (size_t j = 0; j <= i; ++j) {
      sum += x[j] * h[i - j];
    }
    y[i] = sum * dt_us;
  }
  return y;
}

// MicroBooNE example electronics impulse response (time-domain) as coded in LArSoft SimWire::SetElectResponse()
// Here: t_us and To_us in microseconds; Ao is an overall scale.
static double ElectResponse_MicroBooNE(double t_us, double Ao, double To_us) {
  if (To_us <= 0.0) return 0.0;
  const double x = t_us / To_us;
  const double e1 = std::exp(-2.94809 * x);
  const double e2 = std::exp(-2.82833 * x);
  const double e3 = std::exp(-2.40318 * x);

  const double c1 = std::cos(1.19361 * x);
  const double s1 = std::sin(1.19361 * x);
  const double c2 = std::cos(2.38722 * x);
  const double s2 = std::sin(2.38722 * x);
  const double c3 = std::cos(2.5928  * x);
  const double s3 = std::sin(2.5928  * x);
  const double c5 = std::cos(5.18561 * x);
  const double s5 = std::sin(5.18561 * x);

  double r = 0.0;
  r +=  4.31054  * e1 * Ao;
  r += -2.6202   * e2 * c1 * Ao;
  r += -2.6202   * e2 * c1 * c2 * Ao;
  r +=  0.464924 * e3 * c3 * Ao;
  r +=  0.464924 * e3 * c3 * c5 * Ao;
  r +=  0.762456 * e2 * s1 * Ao;
  r += -0.762456 * e2 * c2 * s1 * Ao;
  r +=  0.762456 * e2 * c1 * s2 * Ao;
  r += -2.6202   * e2 * s1 * s2 * Ao;
  r += -0.327684 * e3 * s3 * Ao;
  r +=  0.327684 * e3 * c5 * s3 * Ao;
  r += -0.327684 * e3 * c3 * s5 * Ao;
  r +=  0.464924 * e3 * s3 * s5 * Ao;

  return r;
}

// --------------------------
// Plot builders
// --------------------------

static TCanvas* Plot_dEdx_and_dQdx(const Material& lar,
                                  double Tcut_MeV,
                                  double Efield_kVcm,
                                  double temperature_K,
                                  bool useIcarusMicroBooNEModel,
                                  const RecombParams& rec,
                                  double Wion_eV)
{
  // Particle masses [MeV]
  const double m_mu = 105.6583745;
  const double m_p  = 938.2720813;
  const double m_pi = 139.57039;

  const int Np = 600;
  const double pmin_GeV = 0.05, pmax_GeV = 10.0;

  TGraph* g_dedx_mu = new TGraph(Np);
  TGraph* g_dedx_p  = new TGraph(Np);
  TGraph* g_dedx_pi = new TGraph(Np);

  TGraph* g_dqdx_mu = new TGraph(Np);
  TGraph* g_dqdx_p  = new TGraph(Np);

  for (int i = 0; i < Np; ++i) {
    const double logp = std::log10(pmin_GeV)
                      + (std::log10(pmax_GeV) - std::log10(pmin_GeV)) * i / (Np - 1.0);
    const double p_GeV = std::pow(10.0, logp);
    const double p_MeV = 1000.0 * p_GeV;

    const double dedx_mu = BetheBlochRestricted_MeVpercm(p_MeV, m_mu, 1, lar, Tcut_MeV);
    const double dedx_p  = BetheBlochRestricted_MeVpercm(p_MeV, m_p , 1, lar, Tcut_MeV);
    const double dedx_pi = BetheBlochRestricted_MeVpercm(p_MeV, m_pi, 1, lar, Tcut_MeV);

    const double R_mu = Recomb_ModBox(dedx_mu, Efield_kVcm, lar.rho_gcm3, rec);
    const double R_p  = Recomb_ModBox(dedx_p , Efield_kVcm, lar.rho_gcm3, rec);

    const double dqdx_mu = dQdx_electrons_per_cm(dedx_mu, R_mu, Wion_eV);
    const double dqdx_p  = dQdx_electrons_per_cm(dedx_p , R_p , Wion_eV);

    g_dedx_mu->SetPoint(i, p_GeV, dedx_mu);
    g_dedx_p ->SetPoint(i, p_GeV, dedx_p);
    g_dedx_pi->SetPoint(i, p_GeV, dedx_pi);

    g_dqdx_mu->SetPoint(i, p_GeV, dqdx_mu);
    g_dqdx_p ->SetPoint(i, p_GeV, dqdx_p);
  }

  StyleGraph(g_dedx_mu, kBlack, 1, 3);
  StyleGraph(g_dedx_p , kBlue+1, 1, 3);
  StyleGraph(g_dedx_pi, kRed+1, 2, 3);

  StyleGraph(g_dqdx_mu, kBlack, 1, 3);
  StyleGraph(g_dqdx_p , kBlue+1, 1, 3);

  TCanvas* c = new TCanvas("c_dedx", "", 900, 750);
  c->SetLogx(true);
  c->Divide(1,2);

  // Top: dE/dx
  c->cd(1);
  gPad->SetLogx(true);
  TMultiGraph* mg1 = new TMultiGraph();
  mg1->Add(g_dedx_mu, "L");
  mg1->Add(g_dedx_p , "L");
  mg1->Add(g_dedx_pi, "L");
  mg1->SetTitle(";p  [GeV/c];Restricted dE/dx  [MeV/cm]");
  mg1->Draw("A");
  mg1->GetXaxis()->SetNdivisions(510);
  mg1->GetYaxis()->SetRangeUser(0.0, 40.0);

  TLegend* leg1 = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(g_dedx_mu, "#mu^{#pm}", "l");
  leg1->AddEntry(g_dedx_p , "p", "l");
  leg1->AddEntry(g_dedx_pi, "#pi^{#pm}", "l");
  leg1->Draw();

  DrawHeader("Energy loss and ionization yield in LAr (schematic)",
             Form("Restricted Bethe#minusBloch (T_{cut}=%.2f MeV) and Modified Box recombination @ E=%.0f V/cm, T=%.1f K",
                  Tcut_MeV, 1000.0*Efield_kVcm, temperature_K));

  // Bottom: dQ/dx (after recombination, mean Wion)
  c->cd(2);
  gPad->SetLogx(true);
  TMultiGraph* mg2 = new TMultiGraph();
  mg2->Add(g_dqdx_mu, "L");
  mg2->Add(g_dqdx_p , "L");
  mg2->SetTitle(";p  [GeV/c];#LT dQ/dx #GT  [e^{-}/cm]");
  mg2->Draw("A");
  mg2->GetXaxis()->SetNdivisions(510);
  mg2->GetYaxis()->SetRangeUser(0.0, 8.0e5);

  TLegend* leg2 = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(g_dqdx_mu, "#mu^{#pm} (ModBox)", "l");
  leg2->AddEntry(g_dqdx_p , "p (ModBox)", "l");
  leg2->Draw();

  // Tiny annotation in-plot (kept minimal to avoid clutter)
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.033);
  lat.DrawLatex(0.16, 0.83, Form("v_{d}(E,T)=%.3f cm/#mus (LArSoft-like)", DriftVelocity_LArSoft_cm_per_us(Efield_kVcm, temperature_K, useIcarusMicroBooNEModel)));

  return c;
}

static TCanvas* Plot_Straggling(const Material& lar, double Tcut_MeV)
{
  const double m_mu = 105.6583745;
  const double p0_GeV = 1.0;
  const double p0_MeV = 1000.0 * p0_GeV;

  const int Ns = 350;
  const double smin_cm = 0.01, smax_cm = 50.0;

  TGraph* g_sig = new TGraph(Ns);
  TGraph* g_rel = new TGraph(Ns);

  // Use the same restricted mean dE/dx to define a mean step loss for the relative plot.
  const double dedx0 = BetheBlochRestricted_MeVpercm(p0_MeV, m_mu, 1, lar, Tcut_MeV);

  for (int i = 0; i < Ns; ++i) {
    const double logs = std::log10(smin_cm)
                      + (std::log10(smax_cm) - std::log10(smin_cm)) * i / (Ns - 1.0);
    const double s_cm = std::pow(10.0, logs);

    const double sigma = BohrSigma_MeV(p0_MeV, m_mu, s_cm, 1, lar, Tcut_MeV);
    const double meanE = dedx0 * s_cm;
    const double rel   = (meanE > 0.0) ? (sigma / meanE) : 0.0;

    g_sig->SetPoint(i, s_cm, sigma);
    g_rel->SetPoint(i, s_cm, rel);
  }

  StyleGraph(g_sig, kBlack, 1, 3);
  StyleGraph(g_rel, kBlue+1, 1, 3);

  TCanvas* c = new TCanvas("c_straggle", "", 900, 750);
  c->SetLogx(true);
  c->Divide(1,2);

  c->cd(1);
  gPad->SetLogx(true);
  g_sig->SetTitle(";Step length s  [cm];Bohr straggling RMS  #sigma_{#DeltaE}  [MeV]");
  g_sig->Draw("AL");
  g_sig->GetXaxis()->SetNdivisions(510);

  DrawHeader("Energy-loss fluctuations (Bohr straggling; schematic)",
             Form("1 GeV/c #mu, restricted dE/dx step model; T_{cut}=%.2f MeV", Tcut_MeV));

  c->cd(2);
  gPad->SetLogx(true);
  g_rel->SetTitle(";Step length s  [cm];Relative RMS  #sigma_{#DeltaE}/#LT#DeltaE#GT");
  g_rel->Draw("AL");
  g_rel->GetXaxis()->SetNdivisions(510);
  g_rel->GetYaxis()->SetRangeUser(0.0, 1.0);

  return c;
}

static TCanvas* Plot_MultipleScattering(const Material& lar)
{
  const double m_mu = 105.6583745;
  const double m_p  = 938.2720813;
  const double m_e  = 0.51099895;

  const int Np = 600;
  const double pmin_GeV = 0.05, pmax_GeV = 10.0;

  const std::vector<double> thickness_cm = {1.0, 10.0, 100.0};
  const std::vector<int> colors = {kBlack, kBlue+1, kRed+1};
  const std::vector<int> styles = {1, 2, 3};

  std::vector<TGraph*> g_mu;
  g_mu.reserve(thickness_cm.size());

  for (size_t it = 0; it < thickness_cm.size(); ++it) {
    TGraph* g = new TGraph(Np);
    for (int i = 0; i < Np; ++i) {
      const double logp = std::log10(pmin_GeV)
                        + (std::log10(pmax_GeV) - std::log10(pmin_GeV)) * i / (Np - 1.0);
      const double p_GeV = std::pow(10.0, logp);
      const double p_MeV = 1000.0 * p_GeV;

      const double th = Theta0_Highland_rad(p_MeV, m_mu, thickness_cm[it], lar.X0_cm, 1) * 1e3; // mrad
      g->SetPoint(i, p_GeV, th);
    }
    StyleGraph(g, colors[it], styles[it], 3);
    g_mu.push_back(g);
  }

  // Also show electron vs p at a representative thickness.
  const double t_ref = 10.0;
  TGraph* g_e = new TGraph(Np);
  for (int i = 0; i < Np; ++i) {
    const double logp = std::log10(pmin_GeV)
                      + (std::log10(pmax_GeV) - std::log10(pmin_GeV)) * i / (Np - 1.0);
    const double p_GeV = std::pow(10.0, logp);
    const double p_MeV = 1000.0 * p_GeV;

    const double th = Theta0_Highland_rad(p_MeV, m_e, t_ref, lar.X0_cm, 1) * 1e3; // mrad
    g_e->SetPoint(i, p_GeV, th);
  }
  StyleGraph(g_e, kGreen+2, 1, 3);

  TCanvas* c = new TCanvas("c_mcs", "", 900, 650);
  c->SetLogx(true);

  TMultiGraph* mg = new TMultiGraph();
  for (auto* g : g_mu) mg->Add(g, "L");
  mg->Add(g_e, "L");

  mg->SetTitle(";p  [GeV/c];#theta_{0}  [mrad]");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 600.0);

  TLegend* leg = new TLegend(0.50, 0.66, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (size_t it = 0; it < thickness_cm.size(); ++it) {
    leg->AddEntry(g_mu[it], Form("#mu: t=%.0f cm", thickness_cm[it]), "l");
  }
  leg->AddEntry(g_e, Form("e: t=%.0f cm", t_ref), "l");
  leg->Draw();

  DrawHeader("Multiple scattering in LAr (Highland–Lynch–Dahl; schematic)",
             Form("#theta_{0}#propto (13.6 MeV)/(#beta p) #sqrt{t/X_{0}} (1+h_{c}ln(t/X_{0})),  h_{c}=0.038,  X_{0}#approx%.1f cm", lar.X0_cm));

  return c;
}

static TCanvas* Plot_DriftVelocityAndMapping(const Material& lar,
                                             double temperature_K)
{
  (void)lar; // not used directly here, but kept to match signature pattern

  const int N = 600;
  const double Emin = 0.05, Emax = 2.0; // kV/cm

  // 2 models x 3 temperatures (color encodes temperature; line style encodes model)
  const std::vector<double> temps = {87.0, 89.0, 91.0};
  const std::vector<int> tcolors  = {kBlue+1, kBlack, kRed+1};

  std::vector<TGraph*> gs;
  std::vector<std::string> labels;

  for (size_t it = 0; it < temps.size(); ++it) {
    // model A: Walkowiak+ICARUS transition
    TGraph* gA = new TGraph(N);
    // model B: MicroBooNE+ICARUS polynomial
    TGraph* gB = new TGraph(N);

    for (int i = 0; i < N; ++i) {
      const double E = Emin + (Emax - Emin) * i / (N - 1.0);
      const double vA = DriftVelocity_LArSoft_cm_per_us(E, temps[it], false);
      const double vB = DriftVelocity_LArSoft_cm_per_us(E, temps[it], true);
      gA->SetPoint(i, E, vA);
      gB->SetPoint(i, E, vB);
    }
    StyleGraph(gA, tcolors[it], 2, 3);
    StyleGraph(gB, tcolors[it], 1, 3);

    gs.push_back(gA);
    labels.push_back(Form("T=%.0f K, Walkowiak/ICARUS", temps[it]));

    gs.push_back(gB);
    labels.push_back(Form("T=%.0f K, MicroBooNE+ICARUS", temps[it]));
  }

  // Mapping panel: t(x) for a few fields at T=temperature_K using MicroBooNE+ICARUS
  const std::vector<double> fields = {0.3, 0.5, 0.8}; // kV/cm
  const std::vector<int> fcolors   = {kBlue+1, kBlack, kRed+1};

  const int Nx = 300;
  const double xmax_cm = 400.0;

  std::vector<TGraph*> g_tofx;
  for (size_t iE = 0; iE < fields.size(); ++iE) {
    TGraph* g = new TGraph(Nx);
    const double vd = DriftVelocity_LArSoft_cm_per_us(fields[iE], temperature_K, true);
    for (int ix = 0; ix < Nx; ++ix) {
      const double x_cm = xmax_cm * ix / (Nx - 1.0);
      const double t_ms = (vd > 0.0) ? ( (x_cm / vd) * 1e-3 ) : 0.0; // us -> ms
      g->SetPoint(ix, x_cm, t_ms);
    }
    StyleGraph(g, fcolors[iE], 1, 3);
    g_tofx.push_back(g);
  }

  TCanvas* c = new TCanvas("c_drift", "", 900, 750);
  c->Divide(1,2);

  c->cd(1);
  TMultiGraph* mg1 = new TMultiGraph();
  for (auto* g : gs) mg1->Add(g, "L");
  mg1->SetTitle(";E  [kV/cm];v_{d}  [cm/#mus]");
  mg1->Draw("A");
  mg1->GetXaxis()->SetNdivisions(510);
  mg1->GetYaxis()->SetRangeUser(0.0, 0.35);

  TLegend* leg1 = new TLegend(0.40, 0.55, 0.88, 0.88);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.030);
  for (size_t i = 0; i < gs.size(); ++i) leg1->AddEntry(gs[i], labels[i].c_str(), "l");
  leg1->Draw();

  DrawHeader("Electron drift velocity in LAr (LArSoft-like parameterizations)",
             "Solid: MicroBooNE+ICARUS poly; Dashed: Walkowiak+ICARUS transition (both return cm/#mus)");

  c->cd(2);
  TMultiGraph* mg2 = new TMultiGraph();
  for (auto* g : g_tofx) mg2->Add(g, "L");
  mg2->SetTitle(";Drift distance x  [cm];Drift time t  [ms]");
  mg2->Draw("A");
  mg2->GetXaxis()->SetNdivisions(510);
  mg2->GetYaxis()->SetRangeUser(0.0, 4.0);

  TLegend* leg2 = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.035);
  for (size_t iE = 0; iE < fields.size(); ++iE) {
    leg2->AddEntry(g_tofx[iE], Form("E=%.1f kV/cm", fields[iE]), "l");
  }
  leg2->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.034);
  lat.DrawLatex(0.16, 0.83, Form("t(x)=x/v_{d}(E,T),  T=%.1f K", temperature_K));

  return c;
}

static TCanvas* Plot_Recombination(const Material& lar,
                                   const RecombParams& rec,
                                   double Wion_eV)
{
  const int N = 700;
  const double dEdxMin = 0.3;  // MeV/cm
  const double dEdxMax = 60.0; // MeV/cm

  const std::vector<double> fields = {0.2, 0.5, 1.0}; // kV/cm
  const std::vector<int> colors = {kBlue+1, kBlack, kRed+1};

  std::vector<TGraph*> gR_birks, gR_box;
  std::vector<TGraph*> gQ_birks, gQ_box;

  for (size_t iE = 0; iE < fields.size(); ++iE) {
    TGraph* gRb = new TGraph(N);
    TGraph* gRx = new TGraph(N);
    TGraph* gQb = new TGraph(N);
    TGraph* gQx = new TGraph(N);

    for (int i = 0; i < N; ++i) {
      const double u = (double)i / (double)(N - 1);
      // log spacing in dEdx for readability
      const double dEdx = std::pow(10.0, std::log10(dEdxMin) + u*(std::log10(dEdxMax) - std::log10(dEdxMin)));

      const double Rb = Recomb_Birks (dEdx, fields[iE], lar.rho_gcm3, rec);
      const double Rx = Recomb_ModBox(dEdx, fields[iE], lar.rho_gcm3, rec);

      const double Qb = dQdx_electrons_per_cm(dEdx, Rb, Wion_eV);
      const double Qx = dQdx_electrons_per_cm(dEdx, Rx, Wion_eV);

      gRb->SetPoint(i, dEdx, Rb);
      gRx->SetPoint(i, dEdx, Rx);
      gQb->SetPoint(i, dEdx, Qb);
      gQx->SetPoint(i, dEdx, Qx);
    }

    StyleGraph(gRb, colors[iE], 2, 3);
    StyleGraph(gRx, colors[iE], 1, 3);
    StyleGraph(gQb, colors[iE], 2, 3);
    StyleGraph(gQx, colors[iE], 1, 3);

    gR_birks.push_back(gRb);
    gR_box  .push_back(gRx);
    gQ_birks.push_back(gQb);
    gQ_box  .push_back(gQx);
  }

  TCanvas* c = new TCanvas("c_recomb", "", 900, 750);
  c->Divide(1,2);

  c->cd(1);
  gPad->SetLogx(true);
  TMultiGraph* mg1 = new TMultiGraph();
  for (size_t iE = 0; iE < fields.size(); ++iE) {
    mg1->Add(gR_box[iE],   "L");
    mg1->Add(gR_birks[iE], "L");
  }
  mg1->SetTitle(";dE/dx  [MeV/cm];Recombination factor R");
  mg1->Draw("A");
  mg1->GetXaxis()->SetNdivisions(510);
  mg1->GetYaxis()->SetRangeUser(0.0, 1.05);

  TLegend* leg1 = new TLegend(0.45, 0.62, 0.88, 0.88);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.030);
  leg1->AddEntry((TObject*)nullptr, "Solid = Modified Box, dashed = Birks", "");
  for (size_t iE = 0; iE < fields.size(); ++iE) {
    leg1->AddEntry(gR_box[iE], Form("E=%.1f kV/cm", fields[iE]), "l");
  }
  leg1->Draw();

  DrawHeader("Recombination in LAr: Birks vs Modified Box (LArSoft-like forms)",
             "R(dE/dx,E) multiplies the mean ionization yield: dQ/dx = (dE/dx)/W_{ion} #times R");

  c->cd(2);
  gPad->SetLogx(true);
  TMultiGraph* mg2 = new TMultiGraph();
  for (size_t iE = 0; iE < fields.size(); ++iE) {
    mg2->Add(gQ_box[iE],   "L");
    mg2->Add(gQ_birks[iE], "L");
  }
  mg2->SetTitle(";dE/dx  [MeV/cm];#LT dQ/dx #GT  [e^{-}/cm]");
  mg2->Draw("A");
  mg2->GetXaxis()->SetNdivisions(510);
  mg2->GetYaxis()->SetRangeUser(0.0, 1.1e6);

  return c;
}

static TCanvas* Plot_Lifetime(const Material& lar,
                              double Efield_kVcm,
                              double temperature_K,
                              bool useIcarusMicroBooNEModel)
{
  const double vd = DriftVelocity_LArSoft_cm_per_us(Efield_kVcm, temperature_K, useIcarusMicroBooNEModel);

  const std::vector<double> taus_ms = {1.0, 3.0, 10.0};
  const std::vector<int> colors = {kRed+1, kBlack, kBlue+1};

  const int N = 400;
  const double xmax_cm = 400.0;

  std::vector<TGraph*> gs;
  for (size_t i = 0; i < taus_ms.size(); ++i) {
    TGraph* g = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      const double x = xmax_cm * j / (N - 1.0);
      const double f = SurvivalFraction(x, vd, taus_ms[i]);
      g->SetPoint(j, x, f);
    }
    StyleGraph(g, colors[i], 1, 3);
    gs.push_back(g);
  }

  TCanvas* c = new TCanvas("c_lifetime", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();
  for (auto* g : gs) mg->Add(g, "L");
  mg->SetTitle(";Drift distance x  [cm];Charge survival Q/Q_{0}");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 1.05);

  TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (size_t i = 0; i < taus_ms.size(); ++i) {
    leg->AddEntry(gs[i], Form("#tau_{e}=%.0f ms", taus_ms[i]), "l");
  }
  leg->Draw();

  DrawHeader("Attachment / electron lifetime attenuation",
             Form("Q(x)=Q_{0} exp(-t/ #tau_{e}),  t=x/v_{d};  E=%.0f V/cm, T=%.1f K, v_{d}=%.3f cm/#mus",
                  1000.0*Efield_kVcm, temperature_K, vd));

  return c;
}

static TCanvas* Plot_Diffusion(const Material& lar,
                               double Efield_kVcm,
                               double temperature_K,
                               bool useIcarusMicroBooNEModel)
{
  // Toy diffusion constants in cm^2/s (typical order-of-magnitude at ~500 V/cm).
  // Replace with your experiment’s numbers (or pull from detProp in LArSoft).
  const double DL = 4.8;  // longitudinal diffusion coefficient [cm^2/s]
  const double DT = 13.0; // transverse diffusion coefficient [cm^2/s]

  const double vd = DriftVelocity_LArSoft_cm_per_us(Efield_kVcm, temperature_K, useIcarusMicroBooNEModel);

  const int N = 400;
  const double xmax_cm = 400.0;

  TGraph* gL = new TGraph(N);
  TGraph* gT = new TGraph(N);

  for (int i = 0; i < N; ++i) {
    const double x = xmax_cm * i / (N - 1.0);
    const double sL_cm = DiffusionSigma_cm(x, vd, DL);
    const double sT_cm = DiffusionSigma_cm(x, vd, DT);
    gL->SetPoint(i, x, 10.0*sL_cm); // mm
    gT->SetPoint(i, x, 10.0*sT_cm); // mm
  }

  StyleGraph(gL, kBlue+1, 1, 3);
  StyleGraph(gT, kRed+1,  1, 3);

  TCanvas* c = new TCanvas("c_diff", "", 900, 650);
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gT, "L");
  mg->Add(gL, "L");
  mg->SetTitle(";Drift distance x  [cm];Diffusion RMS  #sigma  [mm]");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(0.0, 6.0);

  TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(gT, Form("#sigma_{T} (D_{T}=%.1f cm^{2}/s)", DT), "l");
  leg->AddEntry(gL, Form("#sigma_{L} (D_{L}=%.1f cm^{2}/s)", DL), "l");
  leg->Draw();

  DrawHeader("Electron cloud diffusion during drift (schematic)",
             Form("#sigma^{2}=2Dt,  t=x/v_{d};  E=%.0f V/cm, T=%.1f K, v_{d}=%.3f cm/#mus",
                  1000.0*Efield_kVcm, temperature_K, vd));

  return c;
}

static TCanvas* Plot_RamoAndShapingToy(double Efield_kVcm,
                                      double temperature_K,
                                      bool useIcarusMicroBooNEModel,
                                      double Wion_eV,
                                      const RecombParams& rec,
                                      const Material& lar)
{
  // A toy “charge packet” corresponding to a small energy deposit, with recombination applied.
  // Pick something like a 1 MeV local deposit.
  const double Edep_MeV = 1.0;

  // Use a representative stopping power to pick a recombination factor (purely illustrative).
  const double dEdx_rep = 2.1; // MeV/cm ~ MIP-ish
  const double R = Recomb_ModBox(dEdx_rep, Efield_kVcm, lar.rho_gcm3, rec);

  const double Ne = (Edep_MeV * 1.0e6 / Wion_eV) * R; // electrons

  // Drift velocity and a toy gap for weighting field.
  const double vd = DriftVelocity_LArSoft_cm_per_us(Efield_kVcm, temperature_K, useIcarusMicroBooNEModel);

  const double gap_cm = 0.30; // 3 mm (illustrative)
  const double tdrift_us = (vd > 0.0) ? (gap_cm / vd) : 0.0;

  // Time axis for waveform (us)
  const double tmax_us = std::max(10.0, 6.0 * tdrift_us);
  const double dt_us   = 0.01; // 10 ns
  const int N = (int)(tmax_us / dt_us) + 1;

  std::vector<double> t_us(N, 0.0), i_nA(N, 0.0);
  const double I0 = RamoCurrent_nA(Ne, vd, gap_cm);

  for (int i = 0; i < N; ++i) {
    t_us[i] = i * dt_us;
    // Rectangular induced current while drifting across gap (parallel-plate toy)
    i_nA[i] = (t_us[i] >= 0.0 && t_us[i] <= tdrift_us) ? I0 : 0.0;
  }

  // Two shaping times for comparison
  const std::vector<double> taus_us = {0.5, 2.0};
  const std::vector<int> colors = {kBlack, kRed+1};

  std::vector<TGraph*> gIshape;
  for (size_t k = 0; k < taus_us.size(); ++k) {
    auto h = CRRCImpulse(t_us, taus_us[k]);
    auto y = ConvolveSameSize(i_nA, h, dt_us); // output in nA * (dimensionless) * us

    // Normalize for easy visual comparison
    double ymax = 0.0;
    for (auto v : y) ymax = std::max(ymax, std::abs(v));
    if (ymax > 0.0) for (auto& v : y) v /= ymax;

    TGraph* g = new TGraph(N);
    for (int i = 0; i < N; ++i) g->SetPoint(i, t_us[i], y[i]);
    StyleGraph(g, colors[k], 1, 3);
    gIshape.push_back(g);
  }

  // Also: show MicroBooNE SimWire electronics response function shape (normalized) for To=0.5 us.
  const double To_us = 0.5;
  std::vector<double> resp(N, 0.0);
  double rmax = 0.0;
  for (int i = 0; i < N; ++i) {
    const double r = ElectResponse_MicroBooNE(t_us[i], 1.0, To_us);
    resp[i] = r;
    rmax = std::max(rmax, std::abs(r));
  }
  if (rmax > 0.0) for (auto& v : resp) v /= rmax;

  TGraph* gResp = new TGraph(N);
  for (int i = 0; i < N; ++i) gResp->SetPoint(i, t_us[i], resp[i]);
  StyleGraph(gResp, kBlue+1, 2, 3);

  // Current graph
  TGraph* gI = new TGraph(N);
  for (int i = 0; i < N; ++i) gI->SetPoint(i, t_us[i], i_nA[i]);
  StyleGraph(gI, kBlack, 1, 3);

  TCanvas* c = new TCanvas("c_ramo", "", 900, 750);
  c->Divide(1,2);

  c->cd(1);
  gI->SetTitle(";t  [#mus];i_{ind}  [nA]");
  gI->Draw("AL");
  gI->GetXaxis()->SetNdivisions(510);

  DrawHeader("Toy Shockley–Ramo induced current and electronics shaping",
             Form("Parallel-plate toy: i=q v/d,  gap=%.1f mm;  E=%.0f V/cm, v_{d}=%.3f cm/#mus;  E_{dep}=%.1f MeV -> N_{e}#approx%.0f",
                  10.0*gap_cm, 1000.0*Efield_kVcm, vd, Edep_MeV, Ne));

  // Draw drift-time marker
  TLine* l = new TLine(tdrift_us, 0.0, tdrift_us, 1.05*I0);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  l->SetLineColor(kGray+2);
  l->Draw("SAME");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.033);
  lat.DrawLatex(0.17, 0.80, Form("t_{drift}=%.2f #mus", tdrift_us));

  c->cd(2);
  TMultiGraph* mg = new TMultiGraph();
  for (auto* g : gIshape) mg->Add(g, "L");
  mg->Add(gResp, "L");
  mg->SetTitle(";t  [#mus];Normalized output");
  mg->Draw("A");
  mg->GetXaxis()->SetNdivisions(510);
  mg->GetYaxis()->SetRangeUser(-0.2, 1.2);

  TLegend* leg = new TLegend(0.48, 0.68, 0.88, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);
  leg->AddEntry(gIshape[0], "CR-RC shape, #tau=0.5 #mus", "l");
  leg->AddEntry(gIshape[1], "CR-RC shape, #tau=2.0 #mus", "l");
  leg->AddEntry(gResp, Form("MicroBooNE SimWire electronics resp (To=%.1f #mus)", To_us), "l");
  leg->Draw();

  return c;
}

static TCanvas* Plot_ScintillationToy()
{
  // Keep it simple: mean yield + Poisson toy
  const double Y_ph_per_MeV = 40000.0; // placeholder; field/particle dependence not included here
  TGraph* g_sc = new TGraph();
  const int Ne = 250;
  for (int i = 0; i < Ne; ++i) {
    const double E = 0.0 + 10.0 * i/(Ne-1.0); // MeV
    g_sc->SetPoint(i, E, Y_ph_per_MeV * E);
  }
  StyleGraph(g_sc, kBlack, 1, 3);

  TRandom3 rng(12345);
  const double Etoy_MeV = 0.01; // 10 keV
  const double meanN = Y_ph_per_MeV * Etoy_MeV;

  TH1D* hN = new TH1D("hN",
      Form(";N_{#gamma};Entries  (E_{dep}=%.0f keV,  #LTN#GT=%.1f)", 1000.0*Etoy_MeV, meanN),
      120, std::max(0.0, meanN - 8.0*std::sqrt(std::max(1.0, meanN))),
           meanN + 8.0*std::sqrt(std::max(1.0, meanN)));

  for (int i = 0; i < 30000; ++i) hN->Fill(rng.PoissonD(meanN));
  hN->SetLineWidth(3);

  TCanvas* c = new TCanvas("c_scint", "", 900, 750);
  c->Divide(1,2);

  c->cd(1);
  g_sc->SetTitle(";E_{dep}  [MeV];#LT N_{#gamma} #GT");
  g_sc->Draw("AL");

  DrawHeader("Scintillation yield (toy model)",
             Form("#LTN_{#gamma}#GT = Y #times E_{dep},  Y=%.0f ph/MeV (placeholder)", Y_ph_per_MeV));

  c->cd(2);
  hN->Draw("HIST");

  return c;
}

// --------------------------
// Driver / dispatcher
// --------------------------
static void PrintCanvasesToPDF(const std::vector<TCanvas*>& pages, const TString& pdf) {
  if (pages.empty()) return;
  pages.front()->Print((pdf + "(").Data());
  for (size_t i = 1; i + 1 < pages.size(); ++i) pages[i]->Print(pdf.Data());
  pages.back()->Print((pdf + ")").Data());
}

void LArTPC_PhysicsPlots(const char* which = "all", const char* outPrefix = "LArTPC")
{
  SetNiceStyle();
  gROOT->SetBatch(kTRUE);

  const std::string w(which ? which : "all");

  // Defaults (edit these to match your detector/FHiCL)
  const Material lar = MakeLAr();

  const double Tcut_MeV = 0.1;      // restricted dE/dx delta-ray threshold
  const double Efield_kVcm = 0.5;   // 500 V/cm
  const double temperature_K = 89.0;

  const bool useIcarusMicroBooNEModel = true; // match your LArSoft config choice

  RecombParams rec; // defaults above (example)
  const double Wion_eV = 23.6;      // mean work function (ion pairs)

  std::vector<TCanvas*> pages;

  auto want = [&](const std::string& key) { return (w == "all" || w == key); };

  if (want("dedx"))      pages.push_back(Plot_dEdx_and_dQdx(lar, Tcut_MeV, Efield_kVcm, temperature_K, useIcarusMicroBooNEModel, rec, Wion_eV));
  if (want("straggle"))  pages.push_back(Plot_Straggling(lar, Tcut_MeV));
  if (want("mcs"))       pages.push_back(Plot_MultipleScattering(lar));
  if (want("drift"))     pages.push_back(Plot_DriftVelocityAndMapping(lar, temperature_K));
  if (want("recomb"))    pages.push_back(Plot_Recombination(lar, rec, Wion_eV));
  if (want("lifetime"))  pages.push_back(Plot_Lifetime(lar, Efield_kVcm, temperature_K, useIcarusMicroBooNEModel));
  if (want("diffusion")) pages.push_back(Plot_Diffusion(lar, Efield_kVcm, temperature_K, useIcarusMicroBooNEModel));
  if (want("ramo"))      pages.push_back(Plot_RamoAndShapingToy(Efield_kVcm, temperature_K, useIcarusMicroBooNEModel, Wion_eV, rec, lar));
  if (want("scint"))     pages.push_back(Plot_ScintillationToy());

  if (pages.empty()) {
    std::cout << "Unknown mode: " << w << "\n"
              << "Options: all, dedx, straggle, mcs, drift, recomb, lifetime, diffusion, ramo, scint\n";
    return;
  }

  const TString base = TString::Format("%s", outPrefix);

  if (w == "all") {
    const TString pdf = base + "_larTPC_physics.pdf";
    PrintCanvasesToPDF(pages, pdf);
    std::cout << "Wrote " << pdf << "\n";
  }
  else {
    const TString pdf = base + "_" + w + ".pdf";
    pages.front()->SaveAs(pdf.Data());
    std::cout << "Wrote " << pdf << "\n";
  }
}
