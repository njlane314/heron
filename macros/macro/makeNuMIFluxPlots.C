// makeNuMIFluxPlots.C
//
// Produces a multi-page PDF with “interesting” flux plots from the NuMIFlux_dk2nu_*.root
// files you listed (CV spectra, mode comparisons, wrong-sign, nue contamination,
// parentage stacks/fractions, region breakdown, a couple of 2D maps, and PPFX multisim
// fractional uncertainties + covariance/correlation matrices).
//
// Run (ROOT):
//   root -l -b -q 'makeNuMIFluxPlots.C()'
// or (heron wrapper, like you used before):
//   heron macro makeNuMIFluxPlots.C
//
// Optional:
//   root -l -b -q 'makeNuMIFluxPlots.C("FHC.root","RHC.root","myplots.pdf")'

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TMath.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "MacroGuard.hh"

namespace {

void SetNiceStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  TH1::SetDefaultSumw2();
}

void DrawHeader(const TString& text) {
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.04);
  lat.DrawLatex(0.12, 0.94, text.Data());
}

TString NormalizeOutDir(const TString& out) {
  TString d(out);
  if (d.EndsWith(".pdf")) d.Resize(d.Length() - 4);
  if (d.IsNull()) d = ".";
  return d;
}

void EnsureOutDir(const TString& outdir) {
  if (outdir.IsNull()) return;
  gSystem->mkdir(outdir, kTRUE);
}

void SaveCanvas(TCanvas& c, const TString& outdir) {
  static int idx = 0;
  EnsureOutDir(outdir);
  TString p(outdir);
  if (!p.EndsWith("/")) p += "/";
  p += TString::Format("plot_%03d.pdf", idx++);
  c.Print(p);
}

double SumPOT(TFile* f) {
  if (!f) return 0.0;
  TTree* t = dynamic_cast<TTree*>(f->Get("POT"));
  if (!t) return 0.0;

  double pot = 0.0;
  double x = 0.0;
  t->SetBranchAddress("POT", &x);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    pot += x;
  }
  return pot;
}

TH1D* GetH1Clone(TFile* f, const TString& path, const TString& newname) {
  if (!f) return nullptr;
  TH1* h = dynamic_cast<TH1*>(f->Get(path));
  if (!h) {
    std::cerr << "[WARN] Missing TH1 at: " << path << "\n";
    return nullptr;
  }
  TH1D* c = dynamic_cast<TH1D*>(h->Clone(newname));
  if (!c) {
    // If the stored type isn't TH1D, clone to TH1 then cast via copy into TH1D
    // (rare for your file, but keep it robust).
    TH1* hc = dynamic_cast<TH1*>(h->Clone(newname));
    if (!hc) return nullptr;
    // Create a TH1D with same binning and copy contents.
    const int nb = hc->GetNbinsX();
    TH1D* tmp = new TH1D(newname, hc->GetTitle(),
                         nb, hc->GetXaxis()->GetXmin(), hc->GetXaxis()->GetXmax());
    for (int i = 1; i <= nb; ++i) {
      tmp->SetBinContent(i, hc->GetBinContent(i));
      tmp->SetBinError(i, hc->GetBinError(i));
    }
    delete hc;
    tmp->SetDirectory(0);
    return tmp;
  }
  c->SetDirectory(0);
  return c;
}

TH2D* GetH2Clone(TFile* f, const TString& path, const TString& newname) {
  if (!f) return nullptr;
  TH2* h = dynamic_cast<TH2*>(f->Get(path));
  if (!h) {
    std::cerr << "[WARN] Missing TH2 at: " << path << "\n";
    return nullptr;
  }
  TH2D* c = dynamic_cast<TH2D*>(h->Clone(newname));
  if (!c) {
    std::cerr << "[WARN] Found TH2 but not TH2D at: " << path << "\n";
    return nullptr;
  }
  c->SetDirectory(0);
  return c;
}

TH1D* MakeRatio(const TH1D* num, const TH1D* den, const TString& name, const TString& title) {
  if (!num || !den) return nullptr;
  TH1D* r = dynamic_cast<TH1D*>(num->Clone(name));
  if (!r) return nullptr;
  r->SetDirectory(0);
  r->SetTitle(title);
  r->Reset("ICES");

  const int nb = std::min(num->GetNbinsX(), den->GetNbinsX());
  for (int i = 1; i <= nb; ++i) {
    const double a = num->GetBinContent(i);
    const double b = den->GetBinContent(i);
    double v = 0.0;
    if (b != 0.0) v = a / b;
    r->SetBinContent(i, v);
    r->SetBinError(i, 0.0);
  }
  return r;
}

struct CovPack {
  TH2D* cov = nullptr;
  TH2D* corr = nullptr;
  TH1D* frac = nullptr;       // sqrt(diag(cov))/CV
};

TH2D* MakeEmptyCovLike(const TH1D* cv, const TString& name, const TString& title) {
  if (!cv) return nullptr;
  const int nb = cv->GetNbinsX();
  const double xmin = cv->GetXaxis()->GetXmin();
  const double xmax = cv->GetXaxis()->GetXmax();
  TH2D* h2 = new TH2D(name, title, nb, xmin, xmax, nb, xmin, xmax);
  h2->SetDirectory(0);
  h2->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  h2->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
  return h2;
}

TH2D* CovToCorr(const TH2D* cov, const TString& name, const TString& title) {
  if (!cov) return nullptr;
  TH2D* corr = dynamic_cast<TH2D*>(cov->Clone(name));
  if (!corr) return nullptr;
  corr->SetDirectory(0);
  corr->SetTitle(title);

  const int nb = corr->GetNbinsX();
  for (int i = 1; i <= nb; ++i) {
    const double vii = cov->GetBinContent(i, i);
    const double sii = (vii > 0) ? std::sqrt(vii) : 0.0;
    for (int j = 1; j <= nb; ++j) {
      const double vjj = cov->GetBinContent(j, j);
      const double sjj = (vjj > 0) ? std::sqrt(vjj) : 0.0;
      const double vij = cov->GetBinContent(i, j);
      double rho = 0.0;
      if (sii > 0 && sjj > 0) rho = vij / (sii * sjj);
      corr->SetBinContent(i, j, rho);
    }
  }
  corr->GetZaxis()->SetRangeUser(-1.0, 1.0);
  return corr;
}

double MaxAbsBinContent(const TH2D* h) {
  if (!h) return 0.0;
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  double m = 0.0;
  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      const double v = std::abs(h->GetBinContent(ix, iy));
      if (v > m) m = v;
    }
  }
  return m;
}

void SetSymmetricZRange(TH2D* h, const double maxAbsOverride = -1.0) {
  if (!h) return;
  double m = maxAbsOverride;
  if (m < 0.0) m = MaxAbsBinContent(h);
  if (m <= 0.0) m = 1.0;
  h->GetZaxis()->SetRangeUser(-m, +m);
}

TH1D* DiagFracUnc(const TH2D* cov, const TH1D* cv, const TString& name, const TString& title) {
  if (!cov || !cv) return nullptr;
  TH1D* f = dynamic_cast<TH1D*>(cv->Clone(name));
  if (!f) return nullptr;
  f->SetDirectory(0);
  f->Reset("ICES");
  f->SetTitle(title);
  f->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  f->GetYaxis()->SetTitle("Fractional uncertainty");

  const int nb = std::min(cov->GetNbinsX(), cv->GetNbinsX());
  for (int i = 1; i <= nb; ++i) {
    const double v = cov->GetBinContent(i, i);
    const double s = (v > 0) ? std::sqrt(v) : 0.0;
    const double c = cv->GetBinContent(i);
    const double frac = (c != 0.0) ? (s / c) : 0.0;
    f->SetBinContent(i, frac);
    f->SetBinError(i, 0.0);
  }
  return f;
}

// Compute covariance from a single multisim group for a given flavor.
// If shapeOnly=true, each universe is scaled to match CV integral before covariance.
CovPack ComputeMultisimCov(
  TFile* f,
  const TString& flavor,
  const TString& group,   // e.g. "ppfx_mipppi_PPFXMIPPPion"
  const TH1D* cv,
  const int nuniv = 200,
  const bool shapeOnly = false,
  const bool useCVasMean = false
) {
  CovPack out;
  if (!f || !cv) return out;

  // Collect universes.
  const int nb = cv->GetNbinsX();
  const double cvInt = cv->Integral(1, nb);

  std::vector<std::vector<double>> vals;
  vals.reserve(nuniv);

  int found = 0;
  for (int u = 0; u < nuniv; ++u) {
    TString path;
    path.Form("%s/Multisims/%s_%s_Uni_%d_AV_TPC", flavor.Data(), flavor.Data(), group.Data(), u);

    TH1* hu = dynamic_cast<TH1*>(f->Get(path));
    if (!hu) {
      // Stop early if universes are contiguous and we hit a missing one.
      if (u == 0) {
        std::cerr << "[WARN] No universes found for group " << group << " at " << path << "\n";
      }
      break;
    }

    const double uniInt = hu->Integral(1, nb);
    double scale = 1.0;
    if (shapeOnly && uniInt != 0.0) scale = cvInt / uniInt;

    std::vector<double> v(nb, 0.0);
    for (int i = 1; i <= nb; ++i) {
      v[i - 1] = hu->GetBinContent(i) * scale;
    }
    vals.push_back(v);
    found++;
  }

  if (found < 2) {
    std::cerr << "[WARN] Too few universes for covariance: " << flavor << " " << group << "\n";
    return out;
  }

  // Mean per bin.
  std::vector<double> mean(nb, 0.0);
  if (useCVasMean) {
    for (int i = 1; i <= nb; ++i) mean[i - 1] = cv->GetBinContent(i);
  } else {
    for (int i = 0; i < nb; ++i) {
      double s = 0.0;
      for (int u = 0; u < found; ++u) s += vals[u][i];
      mean[i] = s / found;
    }
  }

  // Covariance.
  TH2D* cov = MakeEmptyCovLike(cv, TString::Format("cov_%s_%s%s", flavor.Data(), group.Data(), shapeOnly ? "_shape" : ""),
                              TString::Format("%s %s covariance%s", flavor.Data(), group.Data(), shapeOnly ? " (shape-only)" : ""));
  if (!cov) return out;

  for (int u = 0; u < found; ++u) {
    for (int i = 0; i < nb; ++i) {
      const double di = vals[u][i] - mean[i];
      for (int j = 0; j < nb; ++j) {
        const double dj = vals[u][j] - mean[j];
        cov->AddBinContent(cov->GetBin(i + 1, j + 1), di * dj);
      }
    }
  }
  cov->Scale(1.0 / (found - 1.0));

  out.cov = cov;
  out.corr = CovToCorr(cov,
                       TString::Format("corr_%s_%s%s", flavor.Data(), group.Data(), shapeOnly ? "_shape" : ""),
                       TString::Format("%s %s correlation%s", flavor.Data(), group.Data(), shapeOnly ? " (shape-only)" : ""));
  out.frac = DiagFracUnc(cov, cv,
                         TString::Format("frac_%s_%s%s", flavor.Data(), group.Data(), shapeOnly ? "_shape" : ""),
                         TString::Format("%s %s fractional uncertainty%s", flavor.Data(), group.Data(), shapeOnly ? " (shape-only)" : ""));
  return out;
}

void PlotTotalFluxByFlavor(TCanvas& c, TFile* f, const TString& mode, const TString& outdir,
                           const bool logy = true,
                           const int rebinFine = 50,
                           const double xmax = 10.0) {
  const std::vector<TString> flavors = {"numu", "numubar", "nue", "nuebar"};
  const std::vector<int> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};

  std::vector<TH1D*> hs;
  for (size_t i = 0; i < flavors.size(); ++i) {
    const TString& fl = flavors[i];
    TString path;
    path.Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", fl.Data(), fl.Data());
    TH1D* h = GetH1Clone(f, path, TString::Format("h_%s_%s", mode.Data(), fl.Data()));
    if (!h) continue;

    if (rebinFine > 1) h->Rebin(rebinFine);
    h->SetLineColor(colors[i]);
    h->SetLineWidth(2);
    h->SetTitle(TString::Format("%s: CV flux spectra (Detsmear/*_CV_AV_TPC_5MeV_bin)", mode.Data()));
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    h->GetYaxis()->SetTitle("Flux (a.u.)");
    h->GetXaxis()->SetRangeUser(0.0, xmax);
    hs.push_back(h);
  }

  if (hs.empty()) return;

  // Set plotting range.
  double ymax = 0.0;
  double ymin_pos = 1e99;
  for (auto* h : hs) {
    ymax = std::max(ymax, h->GetMaximum());
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      const double v = h->GetBinContent(i);
      if (v > 0) ymin_pos = std::min(ymin_pos, v);
    }
  }
  if (!(ymin_pos < 1e90)) ymin_pos = 1e-12;

  c.Clear();
  c.SetLogy(logy);
  hs[0]->SetMaximum(logy ? ymax * 20.0 : ymax * 1.25);
  hs[0]->SetMinimum(logy ? ymin_pos * 0.2 : 0.0);
  hs[0]->Draw("hist");
  for (size_t i = 1; i < hs.size(); ++i) hs[i]->Draw("hist same");

  TLegend leg(0.60, 0.60, 0.92, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  for (auto* h : hs) leg.AddEntry(h, h->GetName(), "l");
  // Replace internal names with clean flavor labels:
  leg.Clear();
  for (size_t i = 0; i < flavors.size(); ++i) {
    // Find the hist we actually kept (some could be missing).
    for (auto* h : hs) {
      if (TString(h->GetName()).Contains(flavors[i])) {
        leg.AddEntry(h, flavors[i], "l");
        break;
      }
    }
  }
  leg.Draw();

  DrawHeader(TString::Format("%s: flux spectra by flavor", mode.Data()));
  SaveCanvas(c, outdir);
}

void PlotFHCvsRHC_WithRatio(TCanvas& c, TFile* fFHC, TFile* fRHC,
                           const TString& flavor, const TString& outdir,
                           const int rebinFine = 50,
                           const double xmax = 10.0,
                           const bool logy = true) {
  TString pF, pR;
  pF.Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", flavor.Data(), flavor.Data());
  pR = pF;

  TH1D* hF = GetH1Clone(fFHC, pF, TString::Format("hF_%s", flavor.Data()));
  TH1D* hR = GetH1Clone(fRHC, pR, TString::Format("hR_%s", flavor.Data()));
  if (!hF || !hR) return;

  if (rebinFine > 1) { hF->Rebin(rebinFine); hR->Rebin(rebinFine); }

  hF->SetLineColor(kBlue + 1);
  hR->SetLineColor(kRed + 1);
  hF->SetLineWidth(2);
  hR->SetLineWidth(2);

  hF->SetTitle(TString::Format("%s: FHC vs RHC (CV, 5MeV bin, rebinned)", flavor.Data()));
  hF->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  hF->GetYaxis()->SetTitle("Flux (a.u.)");
  hF->GetXaxis()->SetRangeUser(0.0, xmax);

  // Ratio = FHC/RHC
  TH1D* ratio = MakeRatio(hF, hR, TString::Format("ratio_%s", flavor.Data()),
                          TString::Format("%s: FHC/RHC", flavor.Data()));
  if (!ratio) return;
  ratio->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  ratio->GetYaxis()->SetTitle("FHC / RHC");
  ratio->GetXaxis()->SetRangeUser(0.0, xmax);
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);

  c.Clear();

  TPad* p1 = new TPad("p1", "p1", 0, 0.30, 1, 1);
  TPad* p2 = new TPad("p2", "p2", 0, 0.00, 1, 0.30);
  p1->SetBottomMargin(0.02);
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.30);
  p1->Draw();
  p2->Draw();

  p1->cd();
  p1->SetLogy(logy);

  double ymax = std::max(hF->GetMaximum(), hR->GetMaximum());
  double ymin_pos = 1e99;
  for (int i = 1; i <= hF->GetNbinsX(); ++i) {
    double v = hF->GetBinContent(i);
    if (v > 0) ymin_pos = std::min(ymin_pos, v);
    v = hR->GetBinContent(i);
    if (v > 0) ymin_pos = std::min(ymin_pos, v);
  }
  if (!(ymin_pos < 1e90)) ymin_pos = 1e-12;

  hF->SetMaximum(logy ? ymax * 20.0 : ymax * 1.25);
  hF->SetMinimum(logy ? ymin_pos * 0.2 : 0.0);
  hF->Draw("hist");
  hR->Draw("hist same");

  TLegend leg(0.60, 0.70, 0.92, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hF, "FHC", "l");
  leg.AddEntry(hR, "RHC", "l");
  leg.Draw();
  DrawHeader(TString::Format("%s: FHC vs RHC", flavor.Data()));

  p2->cd();
  ratio->SetMinimum(0.0);
  ratio->SetMaximum(2.0);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetLabelSize(0.09);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetLabelSize(0.10);
  ratio->Draw("hist");

  TLine line;
  line.SetLineStyle(2);
  line.DrawLine(0.0, 1.0, xmax, 1.0);

  SaveCanvas(c, outdir);

  delete p1;
  delete p2;
  delete ratio;
}

void PlotWrongSignAndNueContamination(TCanvas& c, TFile* fFHC, TFile* fRHC,
                                     const TString& outdir,
                                     const int rebinFine = 50,
                                     const double xmax = 10.0) {
  // Wrong-sign (mode quality):
  //   FHC: numubar/numu
  //   RHC: numu/numubar
  // Use fine 5MeV hist so both go out to 20 GeV; plot up to xmax.
  auto getFine = [&](TFile* f, const TString& fl) -> TH1D* {
    TString p; p.Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", fl.Data(), fl.Data());
    TH1D* h = GetH1Clone(f, p, TString::Format("fine_%s", fl.Data()));
    if (!h) return nullptr;
    if (rebinFine > 1) h->Rebin(rebinFine);
    h->GetXaxis()->SetRangeUser(0.0, xmax);
    return h;
  };

  TH1D* numuF = getFine(fFHC, "numu");
  TH1D* numubF = getFine(fFHC, "numubar");
  TH1D* numuR = getFine(fRHC, "numu");
  TH1D* numubR = getFine(fRHC, "numubar");
  TH1D* nueF = getFine(fFHC, "nue");
  TH1D* nuebR = getFine(fRHC, "nuebar");

  if (numuF && numubF) {
    TH1D* rF = MakeRatio(numubF, numuF, "ws_FHC", "FHC wrong-sign: #bar{#nu}_{#mu} / #nu_{#mu}");
    rF->SetLineColor(kBlue + 1);
    rF->SetLineWidth(2);
    rF->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    rF->GetYaxis()->SetTitle("Wrong-sign fraction");
    rF->GetXaxis()->SetRangeUser(0.0, xmax);

    TH1D* rR = nullptr;
    if (numuR && numubR) {
      rR = MakeRatio(numuR, numubR, "ws_RHC", "RHC wrong-sign: #nu_{#mu} / #bar{#nu}_{#mu}");
      rR->SetLineColor(kRed + 1);
      rR->SetLineWidth(2);
      rR->GetXaxis()->SetRangeUser(0.0, xmax);
    }

    c.Clear();
    c.SetLogy(false);
    rF->SetMinimum(0.0);
    rF->SetMaximum(1.0);
    rF->Draw("hist");
    if (rR) rR->Draw("hist same");

    TLegend leg(0.55, 0.72, 0.92, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(rF, "FHC: numubar/numu", "l");
    if (rR) leg.AddEntry(rR, "RHC: numu/numubar", "l");
    leg.Draw();

    DrawHeader("Wrong-sign fraction vs energy");
    SaveCanvas(c, outdir);

    delete rF;
    if (rR) delete rR;
  }

  // nue contamination:
  //   FHC: nue/numu
  //   RHC: nuebar/numubar
  if (nueF && numuF) {
    TH1D* cF = MakeRatio(nueF, numuF, "cont_FHC", "FHC intrinsic #nu_{e}: #nu_{e}/#nu_{#mu}");
    cF->SetLineColor(kBlue + 1);
    cF->SetLineWidth(2);
    cF->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    cF->GetYaxis()->SetTitle("Contamination ratio");
    cF->GetXaxis()->SetRangeUser(0.0, xmax);

    TH1D* cR = nullptr;
    if (nuebR && numubR) {
      cR = MakeRatio(nuebR, numubR, "cont_RHC", "RHC intrinsic #bar{#nu}_{e}: #bar{#nu}_{e}/#bar{#nu}_{#mu}");
      cR->SetLineColor(kRed + 1);
      cR->SetLineWidth(2);
      cR->GetXaxis()->SetRangeUser(0.0, xmax);
    }

    c.Clear();
    c.SetLogy(true);
    cF->SetMinimum(1e-6);
    cF->SetMaximum(1e-1);
    cF->Draw("hist");
    if (cR) cR->Draw("hist same");

    TLegend leg(0.55, 0.72, 0.92, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(cF, "FHC: nue/numu", "l");
    if (cR) leg.AddEntry(cR, "RHC: nuebar/numubar", "l");
    leg.Draw();

    DrawHeader("Intrinsic #nu_{e} contamination vs energy");
    SaveCanvas(c, outdir);

    delete cF;
    if (cR) delete cR;
  }
}

void PlotParentageStackAndFraction(TCanvas& c, TFile* f, const TString& mode,
                                  const TString& flavor, const TString& outdir,
                                  const int rebinFine = 50,
                                  const double xmax = 10.0) {
  // Total (CV, fine binning)
  TString pTot;
  pTot.Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", flavor.Data(), flavor.Data());
  TH1D* total = GetH1Clone(f, pTot, TString::Format("total_%s_%s", mode.Data(), flavor.Data()));
  if (!total) return;
  if (rebinFine > 1) total->Rebin(rebinFine);
  total->SetLineColor(kBlack);
  total->SetLineWidth(3);
  total->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  total->GetYaxis()->SetTitle("Flux (a.u.)");
  total->GetXaxis()->SetRangeUser(0.0, xmax);

  const std::vector<TString> parents = {"PI_Plus", "PI_Minus", "Kaon_Plus", "Kaon_Minus", "K0L", "Mu_Plus", "Mu_Minus"};
  const std::vector<int> fillColors = {kAzure - 9, kAzure - 4, kOrange - 2, kOrange + 1, kGreen + 1, kMagenta - 7, kMagenta - 2};

  std::vector<TH1D*> comps;
  comps.reserve(parents.size());

  for (size_t i = 0; i < parents.size(); ++i) {
    TString p;
    p.Form("%s/%s/Enu_%s_%s_AV_TPC", flavor.Data(), parents[i].Data(), flavor.Data(), parents[i].Data());
    TH1D* h = GetH1Clone(f, p, TString::Format("parent_%s_%s_%s", mode.Data(), flavor.Data(), parents[i].Data()));
    if (!h) continue;
    if (rebinFine > 1) h->Rebin(rebinFine);
    h->SetFillColor(fillColors[i % fillColors.size()]);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetRangeUser(0.0, xmax);
    comps.push_back(h);
  }

  // Stack plot
  {
    THStack st(TString::Format("st_%s_%s", mode.Data(), flavor.Data()),
               TString::Format("%s %s parentage stack;Neutrino Energy (GeV);Flux (a.u.)", mode.Data(), flavor.Data()));
    for (auto* h : comps) st.Add(h);

    c.Clear();
    c.SetLogy(true);
    st.Draw("hist");
    st.GetXaxis()->SetRangeUser(0.0, xmax);

    // Overlay total
    total->Draw("hist same");

    // Legend
    TLegend leg(0.55, 0.50, 0.92, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(total, "Total (CV detsmear)", "l");
    for (size_t i = 0; i < comps.size(); ++i) {
      // Recover parent name from hist name.
      TString nm = comps[i]->GetName();
      TString lab = nm;
      lab.ReplaceAll("parent_", "");
      lab.ReplaceAll(TString::Format("%s_", mode.Data()), "");
      lab.ReplaceAll(TString::Format("%s_", flavor.Data()), "");
      leg.AddEntry(comps[i], lab, "f");
    }
    leg.Draw();

    DrawHeader(TString::Format("%s %s: flux parentage (stack)", mode.Data(), flavor.Data()));
    SaveCanvas(c, outdir);
  }

  // Fraction (stacked fractions sum to ~1)
  {
    THStack stf(TString::Format("stfrac_%s_%s", mode.Data(), flavor.Data()),
                TString::Format("%s %s parent fraction;Neutrino Energy (GeV);Fraction", mode.Data(), flavor.Data()));

    std::vector<TH1D*> fracs;
    fracs.reserve(comps.size());
    for (auto* h : comps) {
      TH1D* r = dynamic_cast<TH1D*>(h->Clone(TString::Format("%s_frac", h->GetName())));
      if (!r) continue;
      r->SetDirectory(0);
      r->Divide(total);
      r->SetFillColor(h->GetFillColor());
      r->SetLineColor(kBlack);
      r->SetLineWidth(1);
      r->GetXaxis()->SetRangeUser(0.0, xmax);
      fracs.push_back(r);
      stf.Add(r);
    }

    c.Clear();
    c.SetLogy(false);
    stf.Draw("hist");
    stf.GetXaxis()->SetRangeUser(0.0, xmax);
    stf.SetMinimum(0.0);
    stf.SetMaximum(1.05);

    TLegend leg(0.55, 0.50, 0.92, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for (auto* r : fracs) {
      TString nm = r->GetName();
      nm.ReplaceAll("parent_", "");
      nm.ReplaceAll(TString::Format("%s_", mode.Data()), "");
      nm.ReplaceAll(TString::Format("%s_", flavor.Data()), "");
      nm.ReplaceAll("_frac", "");
      leg.AddEntry(r, nm, "f");
    }
    leg.Draw();

    DrawHeader(TString::Format("%s %s: parent fraction vs energy", mode.Data(), flavor.Data()));
    SaveCanvas(c, outdir);

    for (auto* r : fracs) delete r;
  }

  delete total;
  for (auto* h : comps) delete h;
}

void PlotRegionBreakdown(TCanvas& c, TFile* f, const TString& mode,
                         const TString& flavor, const TString& outdir,
                         const int rebinFine = 50,
                         const double xmax = 10.0) {
  // Uses:
  //   flavor/OtherPlots/flavor_flux_targ
  //   flavor/OtherPlots/flavor_flux_pipe
  //   flavor/OtherPlots/flavor_flux_dump
  TString pT, pP, pD;
  pT.Form("%s/OtherPlots/%s_flux_targ", flavor.Data(), flavor.Data());
  pP.Form("%s/OtherPlots/%s_flux_pipe", flavor.Data(), flavor.Data());
  pD.Form("%s/OtherPlots/%s_flux_dump", flavor.Data(), flavor.Data());

  TH1D* hT = GetH1Clone(f, pT, TString::Format("targ_%s_%s", mode.Data(), flavor.Data()));
  TH1D* hP = GetH1Clone(f, pP, TString::Format("pipe_%s_%s", mode.Data(), flavor.Data()));
  TH1D* hD = GetH1Clone(f, pD, TString::Format("dump_%s_%s", mode.Data(), flavor.Data()));
  if (!hT || !hP || !hD) { delete hT; delete hP; delete hD; return; }

  if (rebinFine > 1) { hT->Rebin(rebinFine); hP->Rebin(rebinFine); hD->Rebin(rebinFine); }

  hT->SetLineColor(kBlue + 1);
  hP->SetLineColor(kGreen + 2);
  hD->SetLineColor(kRed + 1);
  hT->SetLineWidth(2);
  hP->SetLineWidth(2);
  hD->SetLineWidth(2);

  hT->SetTitle(TString::Format("%s %s: region breakdown (OtherPlots/*_flux_*)", mode.Data(), flavor.Data()));
  hT->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  hT->GetYaxis()->SetTitle("Flux (a.u.)");
  hT->GetXaxis()->SetRangeUser(0.0, xmax);

  const double ymax = std::max({hT->GetMaximum(), hP->GetMaximum(), hD->GetMaximum()});
  c.Clear();
  c.SetLogy(true);
  hT->SetMaximum(ymax * 20.0);
  hT->SetMinimum(std::max(1e-12, 0.000000000001));
  hT->Draw("hist");
  hP->Draw("hist same");
  hD->Draw("hist same");

  TLegend leg(0.55, 0.70, 0.92, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hT, "targ", "l");
  leg.AddEntry(hP, "pipe", "l");
  leg.AddEntry(hD, "dump", "l");
  leg.Draw();

  DrawHeader(TString::Format("%s %s: flux by production region", mode.Data(), flavor.Data()));
  SaveCanvas(c, outdir);

  delete hT;
  delete hP;
  delete hD;
}

void Plot2DExamples(TCanvas& c, TFile* f, const TString& mode,
                    const TString& flavor, const TString& outdir) {
  // Two 2D examples you listed:
  //   flavor/OtherPlots/flavor_parent_zpos_angle
  //   flavor/OtherPlots/flavor_parent_zpos_angle_energy
  const std::vector<TString> h2names = {
    TString::Format("%s/OtherPlots/%s_parent_zpos_angle", flavor.Data(), flavor.Data()),
    TString::Format("%s/OtherPlots/%s_parent_zpos_angle_energy", flavor.Data(), flavor.Data())
  };

  for (size_t i = 0; i < h2names.size(); ++i) {
    TH2D* h2 = GetH2Clone(f, h2names[i], TString::Format("h2_%s_%s_%zu", mode.Data(), flavor.Data(), i));
    if (!h2) continue;

    c.Clear();
    c.SetLogy(false);
    c.SetLogz(false);
    c.SetRightMargin(0.14);
    h2->SetTitle(TString::Format("%s %s: %s", mode.Data(), flavor.Data(), h2names[i].Data()));
    // Axis titles are blank in file; set something readable based on name.
    h2->GetXaxis()->SetTitle("z position (a.u.)");
    h2->GetYaxis()->SetTitle("angle (deg)");
    h2->Draw("COLZ");

    DrawHeader(TString::Format("%s %s: 2D map (%zu/2)", mode.Data(), flavor.Data(), i + 1));
    SaveCanvas(c, outdir);

    delete h2;
  }
}

void PlotPPFXUncertainties(TCanvas& c, TFile* f, const TString& mode,
                           const TString& flavor, const TString& outdir,
                           const int nuniv = 200) {
  // CV coarse histogram (matches multisim binning):
  //   numu/Detsmear/numu_CV_AV_TPC  (20 bins, 0-10)
  //   nue/Detsmear/nue_CV_AV_TPC    (18 bins, 0-5)
  TString p;
  p.Form("%s/Detsmear/%s_CV_AV_TPC", flavor.Data(), flavor.Data());
  TH1D* cv = GetH1Clone(f, p, TString::Format("cv_%s_%s", mode.Data(), flavor.Data()));
  if (!cv) return;

  cv->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  cv->GetYaxis()->SetTitle("Flux (a.u.)");

  // Groups present in your file listing.
  const std::vector<TString> groups = {
    "ppfx_mipppi_PPFXMIPPPion",
    "ppfx_mippk_PPFXMIPPKaon",
    "ppfx_thinpi_PPFXThinPion",
    "ppfx_thinmes_PPFXThinMeson",
    "ppfx_thinnpi_PPFXThinNeutronPion",
    "ppfx_totabs_PPFXTotAbsorp",
    "ppfx_targatt_PPFXTargAtten",
    "ppfx_other_PPFXOther",
    "ppfx_ms_UBPPFX",
    "ppfx_think_PPFXThinKaon",
    "thinn_PPFXThinNuc",
    "thinna_PPFXThinNucA"
  };

  // We'll explicitly plot a few components + total:
  //   - pion MIPP
  //   - kaon MIPP
  //   - TOTAL (sum of all covariances; assumes groups are independent)
  // Also compute TOTAL shape-only.
  TH2D* covTotal = MakeEmptyCovLike(cv, TString::Format("covTotal_%s_%s", mode.Data(), flavor.Data()),
                                   TString::Format("%s %s total PPFX covariance (sum of groups)", mode.Data(), flavor.Data()));
  TH2D* covTotalShape = MakeEmptyCovLike(cv, TString::Format("covTotalShape_%s_%s", mode.Data(), flavor.Data()),
                                        TString::Format("%s %s total PPFX covariance (shape-only, sum of groups)", mode.Data(), flavor.Data()));

  TH2D* covPion = nullptr;
  TH2D* covKaon = nullptr;

  for (const auto& g : groups) {
    CovPack pack = ComputeMultisimCov(f, flavor, g, cv, nuniv, /*shapeOnly=*/false, /*useCVasMean=*/false);
    CovPack packS = ComputeMultisimCov(f, flavor, g, cv, nuniv, /*shapeOnly=*/true,  /*useCVasMean=*/false);

    if (pack.cov) covTotal->Add(pack.cov);
    if (packS.cov) covTotalShape->Add(packS.cov);

    if (g == "ppfx_mipppi_PPFXMIPPPion" && pack.cov) covPion = dynamic_cast<TH2D*>(pack.cov->Clone("covPion"));
    if (g == "ppfx_mippk_PPFXMIPPKaon" && pack.cov) covKaon = dynamic_cast<TH2D*>(pack.cov->Clone("covKaon"));

    // Clean up pack objects we don't keep.
    if (pack.corr) delete pack.corr;
    if (pack.frac) delete pack.frac;
    if (pack.cov) delete pack.cov;

    if (packS.corr) delete packS.corr;
    if (packS.frac) delete packS.frac;
    if (packS.cov) delete packS.cov;
  }

  if (covPion) covPion->SetDirectory(0);
  if (covKaon) covKaon->SetDirectory(0);

  TH1D* fracTotal = DiagFracUnc(covTotal, cv, TString::Format("fracTotal_%s_%s", mode.Data(), flavor.Data()),
                               TString::Format("%s %s total PPFX fractional uncertainty", mode.Data(), flavor.Data()));
  TH1D* fracTotalShape = DiagFracUnc(covTotalShape, cv, TString::Format("fracTotalShape_%s_%s", mode.Data(), flavor.Data()),
                                    TString::Format("%s %s total PPFX fractional uncertainty (shape-only)", mode.Data(), flavor.Data()));
  TH1D* fracPion = covPion ? DiagFracUnc(covPion, cv, TString::Format("fracPion_%s_%s", mode.Data(), flavor.Data()),
                                        TString::Format("%s %s pion MIPP fractional uncertainty", mode.Data(), flavor.Data())) : nullptr;
  TH1D* fracKaon = covKaon ? DiagFracUnc(covKaon, cv, TString::Format("fracKaon_%s_%s", mode.Data(), flavor.Data()),
                                        TString::Format("%s %s kaon MIPP fractional uncertainty", mode.Data(), flavor.Data())) : nullptr;

  // 1D fractional uncertainty plot
  if (fracTotal) {
    c.Clear();
    c.SetLogy(false);

    fracTotal->SetLineColor(kBlack);
    fracTotal->SetLineWidth(3);
    fracTotal->SetMinimum(0.0);

    double ymax = fracTotal->GetMaximum();
    if (fracTotalShape) ymax = std::max(ymax, fracTotalShape->GetMaximum());
    if (fracPion) ymax = std::max(ymax, fracPion->GetMaximum());
    if (fracKaon) ymax = std::max(ymax, fracKaon->GetMaximum());
    fracTotal->SetMaximum(std::max(0.2, 1.25 * ymax));

    fracTotal->Draw("hist");
    if (fracTotalShape) {
      fracTotalShape->SetLineColor(kGray + 2);
      fracTotalShape->SetLineStyle(2);
      fracTotalShape->SetLineWidth(2);
      fracTotalShape->Draw("hist same");
    }
    if (fracPion) {
      fracPion->SetLineColor(kBlue + 1);
      fracPion->SetLineWidth(2);
      fracPion->Draw("hist same");
    }
    if (fracKaon) {
      fracKaon->SetLineColor(kRed + 1);
      fracKaon->SetLineWidth(2);
      fracKaon->Draw("hist same");
    }

    TLegend leg(0.50, 0.65, 0.92, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(fracTotal, "Total (sum of groups)", "l");
    if (fracTotalShape) leg.AddEntry(fracTotalShape, "Total shape-only", "l");
    if (fracPion) leg.AddEntry(fracPion, "PPFX MIPP pion", "l");
    if (fracKaon) leg.AddEntry(fracKaon, "PPFX MIPP kaon", "l");
    leg.Draw();

    DrawHeader(TString::Format("%s %s: PPFX multisim uncertainty (coarse binning)", mode.Data(), flavor.Data()));
    SaveCanvas(c, outdir);
  }

  // --------------------------------------------------------------------------
  // Covariance matrices (total + a couple of key components + shape-only total)
  // NOTE: covariance can have negative off-diagonals, so do NOT use logz here.
  // We draw with a symmetric z-range around 0 to make sign structure visible.
  // --------------------------------------------------------------------------

  // Total covariance
  if (covTotal) {
    TH2D* covPlot = dynamic_cast<TH2D*>(covTotal->Clone(TString::Format("covPlotTotal_%s_%s", mode.Data(), flavor.Data())));
    if (covPlot) {
      covPlot->SetDirectory(0);
      covPlot->SetTitle(TString::Format("%s %s: total PPFX covariance (sum of groups)", mode.Data(), flavor.Data()));
      covPlot->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      covPlot->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      SetSymmetricZRange(covPlot);

      c.Clear();
      c.SetRightMargin(0.14);
      c.SetLogz(false);
      covPlot->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX covariance matrix (total)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete covPlot;
    }
  }

  // Total shape-only covariance
  if (covTotalShape) {
    TH2D* covPlot = dynamic_cast<TH2D*>(covTotalShape->Clone(TString::Format("covPlotTotalShape_%s_%s", mode.Data(), flavor.Data())));
    if (covPlot) {
      covPlot->SetDirectory(0);
      covPlot->SetTitle(TString::Format("%s %s: total PPFX covariance (shape-only; sum of groups)", mode.Data(), flavor.Data()));
      covPlot->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      covPlot->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      SetSymmetricZRange(covPlot);

      c.Clear();
      c.SetRightMargin(0.14);
      c.SetLogz(false);
      covPlot->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX covariance matrix (shape-only)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete covPlot;
    }
  }

  // Component covariance: MIPP pion
  if (covPion) {
    TH2D* covPlot = dynamic_cast<TH2D*>(covPion->Clone(TString::Format("covPlotPion_%s_%s", mode.Data(), flavor.Data())));
    if (covPlot) {
      covPlot->SetDirectory(0);
      covPlot->SetTitle(TString::Format("%s %s: PPFX MIPP pion covariance", mode.Data(), flavor.Data()));
      covPlot->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      covPlot->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      SetSymmetricZRange(covPlot);

      c.Clear();
      c.SetRightMargin(0.14);
      c.SetLogz(false);
      covPlot->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX covariance (MIPP pion)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete covPlot;
    }
  }

  // Component covariance: MIPP kaon
  if (covKaon) {
    TH2D* covPlot = dynamic_cast<TH2D*>(covKaon->Clone(TString::Format("covPlotKaon_%s_%s", mode.Data(), flavor.Data())));
    if (covPlot) {
      covPlot->SetDirectory(0);
      covPlot->SetTitle(TString::Format("%s %s: PPFX MIPP kaon covariance", mode.Data(), flavor.Data()));
      covPlot->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      covPlot->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      SetSymmetricZRange(covPlot);

      c.Clear();
      c.SetRightMargin(0.14);
      c.SetLogz(false);
      covPlot->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX covariance (MIPP kaon)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete covPlot;
    }
  }

  // Correlation matrix for total
  {
    TH2D* corrTotal = CovToCorr(covTotal, TString::Format("corrTotal_%s_%s", mode.Data(), flavor.Data()),
                                TString::Format("%s %s total PPFX correlation", mode.Data(), flavor.Data()));
    if (corrTotal) {
      c.Clear();
      c.SetRightMargin(0.14);
      corrTotal->SetTitle(TString::Format("%s %s: total PPFX correlation (sum of groups)", mode.Data(), flavor.Data()));
      corrTotal->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      corrTotal->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      corrTotal->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX correlation matrix", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete corrTotal;
    }
  }

  // Correlation matrix for total shape-only
  if (covTotalShape) {
    TH2D* corrTotalShape = CovToCorr(covTotalShape,
                                     TString::Format("corrTotalShape_%s_%s", mode.Data(), flavor.Data()),
                                     TString::Format("%s %s total PPFX correlation (shape-only)", mode.Data(), flavor.Data()));
    if (corrTotalShape) {
      c.Clear();
      c.SetRightMargin(0.14);
      corrTotalShape->SetTitle(TString::Format("%s %s: total PPFX correlation (shape-only; sum of groups)", mode.Data(), flavor.Data()));
      corrTotalShape->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      corrTotalShape->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      corrTotalShape->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX correlation (shape-only)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete corrTotalShape;
    }
  }

  // Correlation matrices for a couple of key components
  if (covPion) {
    TH2D* corrPion = CovToCorr(covPion,
                              TString::Format("corrPion_%s_%s", mode.Data(), flavor.Data()),
                              TString::Format("%s %s PPFX MIPP pion correlation", mode.Data(), flavor.Data()));
    if (corrPion) {
      c.Clear();
      c.SetRightMargin(0.14);
      corrPion->SetTitle(TString::Format("%s %s: PPFX MIPP pion correlation", mode.Data(), flavor.Data()));
      corrPion->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      corrPion->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      corrPion->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX correlation (MIPP pion)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete corrPion;
    }
  }

  if (covKaon) {
    TH2D* corrKaon = CovToCorr(covKaon,
                              TString::Format("corrKaon_%s_%s", mode.Data(), flavor.Data()),
                              TString::Format("%s %s PPFX MIPP kaon correlation", mode.Data(), flavor.Data()));
    if (corrKaon) {
      c.Clear();
      c.SetRightMargin(0.14);
      corrKaon->SetTitle(TString::Format("%s %s: PPFX MIPP kaon correlation", mode.Data(), flavor.Data()));
      corrKaon->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
      corrKaon->GetYaxis()->SetTitle("Neutrino Energy (GeV)");
      corrKaon->Draw("COLZ");
      DrawHeader(TString::Format("%s %s: PPFX correlation (MIPP kaon)", mode.Data(), flavor.Data()));
      SaveCanvas(c, outdir);
      delete corrKaon;
    }
  }

  // Reset canvas margins for callers (avoid leaking a big right margin).
  c.SetRightMargin(0.05);

  delete cv;
  delete covTotal;
  delete covTotalShape;
  if (covPion) delete covPion;
  if (covKaon) delete covKaon;
  if (fracTotal) delete fracTotal;
  if (fracTotalShape) delete fracTotalShape;
  if (fracPion) delete fracPion;
  if (fracKaon) delete fracKaon;
}

} // namespace

void makeNuMIFluxPlots(
  const char* fhc_file = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root",
  const char* rhc_file = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root",
  const char* outdir_in   = "NuMIFlux_dk2nu_plots"
) {
  heron::macro::run_with_guard("makeNuMIFluxPlots", [&]() {
  gROOT->SetBatch(kTRUE);
  SetNiceStyle();

  TFile* fFHC = TFile::Open(fhc_file, "READ");
  TFile* fRHC = TFile::Open(rhc_file, "READ");
  if (!fFHC || fFHC->IsZombie()) {
    std::cerr << "[ERROR] Failed to open FHC file: " << fhc_file << "\n";
    return;
  }
  if (!fRHC || fRHC->IsZombie()) {
    std::cerr << "[ERROR] Failed to open RHC file: " << rhc_file << "\n";
    return;
  }

  const double potF = SumPOT(fFHC);
  const double potR = SumPOT(fRHC);
  std::cout << "[INFO] FHC total POT (sum of POT tree): " << potF << "\n";
  std::cout << "[INFO] RHC total POT (sum of POT tree): " << potR << "\n";

  TString outdir = NormalizeOutDir(outdir_in);
  EnsureOutDir(outdir);

  TCanvas c("c", "c", 1000, 750);


  // 1) Flux by flavor (per mode)
  PlotTotalFluxByFlavor(c, fFHC, "FHC", outdir, /*logy=*/true,  /*rebinFine=*/50, /*xmax=*/10.0);
  PlotTotalFluxByFlavor(c, fRHC, "RHC", outdir, /*logy=*/true,  /*rebinFine=*/50, /*xmax=*/10.0);

  // 2) FHC vs RHC comparison for dominant flavors
  PlotFHCvsRHC_WithRatio(c, fFHC, fRHC, "numu",   outdir, /*rebinFine=*/50, /*xmax=*/10.0, /*logy=*/true);
  PlotFHCvsRHC_WithRatio(c, fFHC, fRHC, "numubar",outdir, /*rebinFine=*/50, /*xmax=*/10.0, /*logy=*/true);

  // 3) Wrong-sign + nue contamination
  PlotWrongSignAndNueContamination(c, fFHC, fRHC, outdir, /*rebinFine=*/50, /*xmax=*/10.0);

  // 4) Parentage stacks + fractions
  PlotParentageStackAndFraction(c, fFHC, "FHC", "numu",    outdir, /*rebinFine=*/50, /*xmax=*/10.0);
  PlotParentageStackAndFraction(c, fRHC, "RHC", "numubar", outdir, /*rebinFine=*/50, /*xmax=*/10.0);

  // 5) Region breakdown (targ/pipe/dump)
  PlotRegionBreakdown(c, fFHC, "FHC", "numu", outdir, /*rebinFine=*/50, /*xmax=*/10.0);
  PlotRegionBreakdown(c, fRHC, "RHC", "numubar", outdir, /*rebinFine=*/50, /*xmax=*/10.0);

  // 6) A couple of 2D “map” plots
  Plot2DExamples(c, fFHC, "FHC", "numu", outdir);
  Plot2DExamples(c, fRHC, "RHC", "numubar", outdir);

  // 7) PPFX multisim uncertainty (coarse binning; matches Multisims hist binning)
  PlotPPFXUncertainties(c, fFHC, "FHC", "numu", outdir, /*nuniv=*/200);
  PlotPPFXUncertainties(c, fRHC, "RHC", "numubar", outdir, /*nuniv=*/200);


  fFHC->Close();
  fRHC->Close();

  std::cout << "[INFO] Wrote PDFs under: " << outdir << "\n";

  });
}
