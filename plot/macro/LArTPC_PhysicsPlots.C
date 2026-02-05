// LArTPC_PhysicsPlots.C
// Plots a few Geant4 Physics Reference Manual equations commonly relevant to LArTPC detector physics:
//  - Tmax (Eq. 13.1) and restricted Bethe-Bloch dE/dx (Eq. 13.2)
//  - Bohr straggling variance (Eq. 8.8)
//  - Highland-Lynch-Dahl multiple scattering width theta0 (hc = 0.038)
//  - Cerenkov angle + photon yield
//  - Scintillation: <N_gamma>=Y*Edep and a Poisson toy

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TString.h"

struct Material {
  double Z;        // atomic number
  double A_gmol;   // g/mol
  double rho_gcm3; // g/cm^3
  double I_MeV;    // mean excitation energy [MeV]
  double X0_cm;    // radiation length [cm]
  double n;        // refractive index (for Cerenkov plots)
};

// Constants
static const double kRe_cm   = 2.8179403227e-13; // classical electron radius [cm]
static const double kMe_MeV  = 0.51099895;       // electron mass-energy [MeV]
static const double kNa      = 6.02214076e23;    // Avogadro [1/mol]
static const double kPi      = TMath::Pi();

// --- Kinematics helpers
double BetaFromP(double p_MeVc, double M_MeV) {
  const double E = TMath::Sqrt(p_MeVc*p_MeVc + M_MeV*M_MeV);
  return (E > 0.0) ? (p_MeVc / E) : 0.0;
}

double GammaFromP(double p_MeVc, double M_MeV) {
  const double E = TMath::Sqrt(p_MeVc*p_MeVc + M_MeV*M_MeV);
  return (M_MeV > 0.0) ? (E / M_MeV) : 1.0;
}

double ElectronDensity_cm3(double Z, double A_gmol, double rho_gcm3) {
  // n_el = Z * N_A * rho / A [electrons/cm^3]
  return Z * kNa * rho_gcm3 / A_gmol;
}

// --- Eq. (13.1): Tmax
double Tmax_MeV(double p_MeVc, double M_MeV) {
  const double gamma = GammaFromP(p_MeVc, M_MeV);
  const double r = kMe_MeV / M_MeV;
  const double num = 2.0 * kMe_MeV * (gamma*gamma - 1.0);
  const double den = 1.0 + 2.0*gamma*r + r*r;
  return (den > 0.0) ? (num / den) : 0.0;
}

// --- Eq. (13.2): restricted Bethe-Bloch dE/dx (correction terms default to 0 here)
double BetheBlochRestricted_MeVpercm(double p_MeVc,
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
  const double Tup  = (Tcut_MeV < Tmax) ? Tcut_MeV : Tmax; // min(Tcut, Tmax)

  const double n_el = ElectronDensity_cm3(mat.Z, mat.A_gmol, mat.rho_gcm3);

  // 2π r_e^2 m_e c^2 n_el has units MeV/cm (with r_e in cm, m_e c^2 in MeV, n_el in 1/cm^3)
  const double pref = 2.0 * kPi * (kRe_cm*kRe_cm) * kMe_MeV * n_el;

  const double beta2 = beta*beta;

  double arg = 2.0 * kMe_MeV * beta2 * gamma*gamma * Tup / (mat.I_MeV * mat.I_MeV);
  if (arg <= 1.0) arg = 1.0000001; // avoid log(<=1) pathologies for very low beta

  const double bracket =
      TMath::Log(arg)
    - beta2 * (1.0 + Tup/Tmax)
    - delta
    - 2.0 * Ce / mat.Z
    + S
    + F;

  return pref * ( (double)(z*z) / beta2 ) * bracket; // MeV/cm
}

// --- Eq. (8.8): Bohr straggling sigma = sqrt(Omega^2) for a step of length s
double BohrSigma_MeV(double p_MeVc,
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

  return TMath::Sqrt(Omega2); // MeV
}

// --- Highland-Lynch-Dahl theta0 (original form, with hc=0.038)
double Theta0_Highland_rad(double p_MeVc,
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

  const double corr = 1.0 + hc * TMath::Log(tx0);
  return (13.6 / (beta * p_MeVc)) * z * TMath::Sqrt(tx0) * corr; // rad
}

// --- Cerenkov: angle and yield (constant n, constant band)
double CerenkovTheta_deg(double beta, double n) {
  if (beta <= 0.0) return 0.0;
  if (beta*n <= 1.0) return 0.0;
  const double costh = 1.0 / (beta*n);
  return TMath::ACos(costh) * 180.0 / kPi;
}

double Cerenkov_dNdx_phot_per_cm(double beta, double n, double epsMin_eV, double epsMax_eV, int z) {
  if (beta*n <= 1.0) return 0.0;
  const double term = 1.0 - 1.0/(n*n*beta*beta);
  double dE = epsMax_eV - epsMin_eV;
  if (dE < 0.0) dE = 0.0;
  return 370.0 * z * z * dE * term; // photons/cm
}

// --- Default-ish LAr numbers (edit to match your assumptions)
Material MakeLAr() {
  Material m;
  m.Z        = 18.0;
  m.A_gmol   = 39.948;
  m.rho_gcm3 = 1.396;
  m.I_MeV    = 188e-6; // 188 eV -> 188e-6 MeV
  m.X0_cm    = 14.0;   // ~ (19.55 g/cm^2)/(1.396 g/cm^3) ≈ 14 cm
  m.n        = 1.23;   // placeholder constant n for Cerenkov plots
  return m;
}

void LArTPC_PhysicsPlots(const char* outPrefix="LArTPC")
{
  gStyle->SetOptStat(0);

  const Material lar = MakeLAr();

  // Particle masses [MeV]
  const double m_mu = 105.6583745;
  const double m_p  = 938.2720813;
  const double m_e  = 0.51099895;

  // Plot output
  const TString pdf = TString::Format("%s_physics_plots.pdf", outPrefix);

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.035);

  // -------------------------
  // 1) Bethe-Bloch dE/dx vs p
  // -------------------------
  const int Np = 400;
  const double pmin_GeV = 0.05, pmax_GeV = 10.0;
  const double Tcut_MeV = 0.1; // delta-ray threshold (edit)

  TGraph *g_mu = new TGraph();
  TGraph *g_pr = new TGraph();

  for (int i=0; i<Np; ++i) {
    const double logp = TMath::Log10(pmin_GeV)
                      + (TMath::Log10(pmax_GeV)-TMath::Log10(pmin_GeV)) * i/(Np-1.0);
    const double p_GeV = TMath::Power(10.0, logp);
    const double p_MeV = 1000.0 * p_GeV;

    const double dedx_mu = BetheBlochRestricted_MeVpercm(p_MeV, m_mu, 1, lar, Tcut_MeV);
    const double dedx_pr = BetheBlochRestricted_MeVpercm(p_MeV, m_p , 1, lar, Tcut_MeV);

    g_mu->SetPoint(i, p_GeV, dedx_mu);
    g_pr->SetPoint(i, p_GeV, dedx_pr);
  }

  TCanvas *c1 = new TCanvas("c1","Bethe-Bloch",900,700);
  c1->SetLogx(true);

  g_mu->SetLineWidth(2);
  g_pr->SetLineWidth(2);

  g_mu->SetTitle("Restricted Bethe-Bloch dE/dx in LAr; p [GeV/c]; dE/dx [MeV/cm]");
  g_mu->Draw("AL");
  g_pr->Draw("L SAME");

  TLegend *leg1 = new TLegend(0.60,0.70,0.88,0.88);
  leg1->AddEntry(g_mu,"#mu^{#pm} (z=1)","l");
  leg1->AddEntry(g_pr,"p (z=1)","l");
  leg1->Draw();

  lat.DrawLatex(0.14,0.86,"Eq. (13.1, 13.2)");
  c1->Print((pdf+"(").Data());

  // -----------------------------------------
  // 2) Bohr straggling sigma vs step length s
  // -----------------------------------------
  const double p0_GeV = 1.0;
  const double p0_MeV = 1000.0 * p0_GeV;

  TGraph *g_sig = new TGraph();
  const int Ns = 250;
  const double smin_cm = 0.01, smax_cm = 20.0;

  for (int i=0; i<Ns; ++i) {
    const double logs = TMath::Log10(smin_cm)
                      + (TMath::Log10(smax_cm)-TMath::Log10(smin_cm)) * i/(Ns-1.0);
    const double s_cm = TMath::Power(10.0, logs);

    const double sigma = BohrSigma_MeV(p0_MeV, m_mu, s_cm, 1, lar, Tcut_MeV);
    g_sig->SetPoint(i, s_cm, sigma);
  }

  TCanvas *c2 = new TCanvas("c2","Straggling",900,700);
  c2->SetLogx(true);

  g_sig->SetLineWidth(2);
  g_sig->SetTitle("Bohr straggling RMS (#sigma) for a 1 GeV/c #mu in LAr; step length s [cm]; #sigma_{#DeltaE} [MeV]");
  g_sig->Draw("AL");

  lat.DrawLatex(0.14,0.86,"Eq. (8.8)");
  c2->Print(pdf.Data());

  // ----------------------------------------
  // 3) Multiple scattering theta0 vs momentum
  // ----------------------------------------
  const double t_cm = 10.0; // thickness
  TGraph *g_th_mu = new TGraph();
  TGraph *g_th_e  = new TGraph();

  for (int i=0; i<Np; ++i) {
    const double logp = TMath::Log10(pmin_GeV)
                      + (TMath::Log10(pmax_GeV)-TMath::Log10(pmin_GeV)) * i/(Np-1.0);
    const double p_GeV = TMath::Power(10.0, logp);
    const double p_MeV = 1000.0 * p_GeV;

    const double th_mu = Theta0_Highland_rad(p_MeV, m_mu, t_cm, lar.X0_cm, 1) * 1e3; // mrad
    const double th_e  = Theta0_Highland_rad(p_MeV, m_e , t_cm, lar.X0_cm, 1) * 1e3;

    g_th_mu->SetPoint(i, p_GeV, th_mu);
    g_th_e ->SetPoint(i, p_GeV, th_e);
  }

  TCanvas *c3 = new TCanvas("c3","Multiple scattering",900,700);
  c3->SetLogx(true);

  g_th_mu->SetLineWidth(2);
  g_th_e ->SetLineWidth(2);

  g_th_mu->SetTitle(Form("Highland-Lynch-Dahl #theta_{0} in LAr (t=%.1f cm); p [GeV/c]; #theta_{0} [mrad]", t_cm));
  g_th_mu->Draw("AL");
  g_th_e->Draw("L SAME");

  TLegend *leg3 = new TLegend(0.60,0.70,0.88,0.88);
  leg3->AddEntry(g_th_mu,"#mu^{#pm}","l");
  leg3->AddEntry(g_th_e ,"e^{#pm}","l");
  leg3->Draw();

  lat.DrawLatex(0.14,0.86,"#theta_{0} with h_{c}=0.038");
  c3->Print(pdf.Data());

  // -----------------------------------------
  // 4) Cerenkov angle and photon yield vs beta
  // -----------------------------------------
  const int Nb = 400;
  const double epsMin_eV = 2.0, epsMax_eV = 6.0; // example band
  TGraph *g_ck_th = new TGraph();
  TGraph *g_ck_N  = new TGraph();

  for (int i=0; i<Nb; ++i) {
    const double beta = 0.0 + (0.9999-0.0) * i/(Nb-1.0);
    const double th   = CerenkovTheta_deg(beta, lar.n);
    const double dn   = Cerenkov_dNdx_phot_per_cm(beta, lar.n, epsMin_eV, epsMax_eV, 1);
    g_ck_th->SetPoint(i, beta, th);
    g_ck_N ->SetPoint(i, beta, dn);
  }

  TCanvas *c4 = new TCanvas("c4","Cerenkov",900,700);
  c4->Divide(1,2);

  c4->cd(1);
  g_ck_th->SetLineWidth(2);
  g_ck_th->SetTitle(Form("Cerenkov angle (constant n=%.3f); #beta; #theta_{C} [deg]", lar.n));
  g_ck_th->Draw("AL");
  lat.DrawLatex(0.14,0.83,"cos#theta = 1/(#beta n)");

  c4->cd(2);
  g_ck_N->SetLineWidth(2);
  g_ck_N->SetTitle(Form("Cerenkov photons per cm (%.1f-%.1f eV band); #beta; dN/dx [photons/cm]", epsMin_eV, epsMax_eV));
  g_ck_N->Draw("AL");
  lat.DrawLatex(0.14,0.83,"dN/dx #approx 370 z^{2} #int d#epsilon (1 - 1/(n^{2}#beta^{2}))");

  c4->Print(pdf.Data());

  // ----------------------------------------------
  // 5) Scintillation: mean yield + Poisson toy demo
  // ----------------------------------------------
  const double Y_ph_per_MeV = 40000.0; // user-set
  TGraph *g_sc = new TGraph();

  const int Ne = 200;
  for (int i=0; i<Ne; ++i) {
    const double E = 0.0 + 10.0 * i/(Ne-1.0); // MeV
    g_sc->SetPoint(i, E, Y_ph_per_MeV * E);
  }

  // Poisson toy for a small deposit (so the distribution is visible)
  TRandom3 rng(12345);
  const double Etoy_MeV = 0.01; // 10 keV
  const double meanN = Y_ph_per_MeV * Etoy_MeV;

  TH1D *hN = new TH1D("hN",
      Form("Poisson toy: E_{dep}=%.3f MeV, Y=%.0f ph/MeV; N_{#gamma}; Entries", Etoy_MeV, Y_ph_per_MeV),
      200,
      meanN - 8.0*TMath::Sqrt(meanN),
      meanN + 8.0*TMath::Sqrt(meanN));

  for (int i=0; i<20000; ++i) hN->Fill(rng.PoissonD(meanN));

  TCanvas *c5 = new TCanvas("c5","Scintillation",900,700);
  c5->Divide(1,2);

  c5->cd(1);
  g_sc->SetLineWidth(2);
  g_sc->SetTitle("Scintillation mean yield; E_{dep} [MeV]; <N_{#gamma}>");
  g_sc->Draw("AL");
  lat.DrawLatex(0.14,0.83,"<N_{#gamma}> = Y #times E_{dep}");

  c5->cd(2);
  hN->SetLineWidth(2);
  hN->Draw("HIST");

  c5->Print((pdf+")").Data());
}
