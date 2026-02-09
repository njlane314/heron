// plotOscPars.C
// ROOT macro to plot NuFIT 6.0 oscillation-parameter constraints (best-fit, 1σ, 3σ)
// Values are taken from the table in your screenshot (IC19 w/o SK-atm, IC24 with SK-atm).
//
// Usage:
//   root -l -q 'plotOscPars.C("IC24","oscpars_IC24.pdf")'
//   root -l -q 'plotOscPars.C("IC19","oscpars_IC19.pdf")'

#include <algorithm>
#include <vector>

#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"
#include "TStyle.h"

struct ParEntry {
  TString name;    // for internal naming
  TString xtitle;  // x-axis title (parameter + units)

  // Normal ordering (NO)
  double no_bf;
  double no_sig_up;   // +1σ
  double no_sig_dn;   // -1σ (stored as positive number)
  double no_3min;     // 3σ min
  double no_3max;     // 3σ max

  // Inverted ordering (IO)
  double io_bf;
  double io_sig_up;
  double io_sig_dn;
  double io_3min;
  double io_3max;
};

static void DrawOnePar(const ParEntry& p, int idx, bool drawLegend=false)
{
  const double yIO = 1.0;
  const double yNO = 2.0;

  // Compute axis range from union of 3σ ranges
  double xmin = std::min(p.no_3min, p.io_3min);
  double xmax = std::max(p.no_3max, p.io_3max);
  double dx = xmax - xmin;
  if (dx <= 0) {
    dx = 1.0;
    xmin -= 0.5;
    xmax += 0.5;
  }
  xmin -= 0.08*dx;
  xmax += 0.08*dx;

  gPad->SetLeftMargin(0.18);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.18);
  gPad->SetTicks(1,0);

  TH2F* frame = new TH2F(Form("frame_%d_%s", idx, p.name.Data()),
                         "", 100, xmin, xmax, 10, 0.5, 2.5);
  frame->SetStats(0);
  frame->GetYaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTickLength(0);
  frame->GetXaxis()->SetTitle(p.xtitle);
  frame->GetXaxis()->SetTitleSize(0.065);
  frame->GetXaxis()->SetLabelSize(0.055);
  frame->GetXaxis()->SetTitleOffset(1.10);
  frame->Draw();

  // Ordering labels (NDC so they don't depend on x-range)
  TLatex lab;
  lab.SetNDC(true);
  lab.SetTextSize(0.07);
  lab.SetTextAlign(12);
  lab.DrawLatex(0.02, 0.70, "NO");
  lab.DrawLatex(0.02, 0.30, "IO");

  // 3σ ranges (thin, gray)
  TLine l3;
  l3.SetLineColor(kGray + 2);
  l3.SetLineWidth(2);
  l3.DrawLine(p.no_3min, yNO, p.no_3max, yNO);
  l3.DrawLine(p.io_3min, yIO, p.io_3max, yIO);

  // 1σ ranges (thick, colored)
  TLine l1;
  l1.SetLineWidth(5);

  l1.SetLineColor(kBlue + 1); // NO
  l1.DrawLine(p.no_bf - p.no_sig_dn, yNO, p.no_bf + p.no_sig_up, yNO);

  l1.SetLineColor(kRed + 1);  // IO
  l1.DrawLine(p.io_bf - p.io_sig_dn, yIO, p.io_bf + p.io_sig_up, yIO);

  // Best-fit markers
  TMarker m;
  m.SetMarkerStyle(20);
  m.SetMarkerSize(1.1);

  m.SetMarkerColor(kBlue + 1);
  m.DrawMarker(p.no_bf, yNO);

  m.SetMarkerColor(kRed + 1);
  m.DrawMarker(p.io_bf, yIO);

  if (drawLegend) {
    // Dummy objects for legend
    TLine* legNO = new TLine(0,0,1,0);
    legNO->SetLineColor(kBlue+1);
    legNO->SetLineWidth(5);

    TLine* legIO = new TLine(0,0,1,0);
    legIO->SetLineColor(kRed+1);
    legIO->SetLineWidth(5);

    TLine* leg3  = new TLine(0,0,1,0);
    leg3->SetLineColor(kGray+2);
    leg3->SetLineWidth(2);

    TMarker* legM = new TMarker(0,0,20);
    legM->SetMarkerSize(1.1);

    TLegend* L = new TLegend(0.45, 0.15, 0.95, 0.42);
    L->SetBorderSize(0);
    L->SetFillStyle(0);
    L->SetTextSize(0.05);
    L->AddEntry(legNO, "NO  1#sigma", "l");
    L->AddEntry(legIO, "IO  1#sigma", "l");
    L->AddEntry(leg3,  "3#sigma range", "l");
    L->AddEntry(legM,  "best fit", "p");
    L->Draw();
  }
}

void plotOscPars(TString variant="IC24", TString out="oscpars.pdf")
{
  gStyle->SetOptStat(0);

  // Choose which block from the table to use
  const bool useIC24 = (variant.Contains("IC24") || variant.Contains("ic24"));

  std::vector<ParEntry> pars;
  pars.reserve(6);

  if (useIC24) {
    // IC24 with SK-atm
    pars.push_back({"s2th12", "#sin^{2}#theta_{12}",
      0.308, 0.012, 0.011, 0.275, 0.345,
      0.308, 0.012, 0.011, 0.275, 0.345});

    pars.push_back({"s2th23", "#sin^{2}#theta_{23}",
      0.470, 0.017, 0.013, 0.435, 0.585,
      0.550, 0.012, 0.015, 0.440, 0.584});

    pars.push_back({"s2th13", "#sin^{2}#theta_{13}",
      0.02215, 0.00056, 0.00058, 0.02030, 0.02388,
      0.02231, 0.00056, 0.00056, 0.02060, 0.02409});

    pars.push_back({"dcp", "#delta_{CP} [deg]",
      212.0, 26.0, 41.0, 124.0, 364.0,
      274.0, 22.0, 25.0, 201.0, 335.0});

    pars.push_back({"dm21", "#Deltam^{2}_{21} [10^{-5} eV^{2}]",
      7.49, 0.19, 0.19, 6.92, 8.05,
      7.49, 0.19, 0.19, 6.92, 8.05});

    pars.push_back({"dm3l", "#Deltam^{2}_{3#ell} [10^{-3} eV^{2}]",
      +2.513, 0.021, 0.019, +2.451, +2.578,
      -2.484, 0.020, 0.020, -2.547, -2.421});
  } else {
    // IC19 w/o SK-atm
    pars.push_back({"s2th12", "#sin^{2}#theta_{12}",
      0.307, 0.012, 0.011, 0.275, 0.345,
      0.308, 0.012, 0.011, 0.275, 0.345});

    pars.push_back({"s2th23", "#sin^{2}#theta_{23}",
      0.561, 0.012, 0.015, 0.430, 0.596,
      0.562, 0.012, 0.015, 0.437, 0.597});

    pars.push_back({"s2th13", "#sin^{2}#theta_{13}",
      0.02195, 0.00054, 0.00058, 0.02023, 0.02376,
      0.02224, 0.00056, 0.00057, 0.02053, 0.02397});

    pars.push_back({"dcp", "#delta_{CP} [deg]",
      177.0, 19.0, 20.0, 96.0, 422.0,
      285.0, 25.0, 28.0, 201.0, 348.0});

    pars.push_back({"dm21", "#Deltam^{2}_{21} [10^{-5} eV^{2}]",
      7.49, 0.19, 0.19, 6.92, 8.05,
      7.49, 0.19, 0.19, 6.92, 8.05});

    pars.push_back({"dm3l", "#Deltam^{2}_{3#ell} [10^{-3} eV^{2}]",
      +2.534, 0.025, 0.023, +2.463, +2.606,
      -2.510, 0.024, 0.025, -2.584, -2.438});
  }

  TCanvas* c = new TCanvas("c", "Oscillation parameters", 1500, 900);
  c->Divide(3, 2, 0.005, 0.005);

  for (int i = 0; i < static_cast<int>(pars.size()); ++i) {
    c->cd(i+1);
    DrawOnePar(pars[i], i, (i==0)); // put legend only in first pad

    if (i==0) {
      TLatex hdr;
      hdr.SetNDC(true);
      hdr.SetTextAlign(13);
      hdr.SetTextSize(0.055);
      hdr.DrawLatex(0.18, 0.98, Form("NuFIT 6.0 (%s)", useIC24 ? "IC24 with SK-atm" : "IC19 w/o SK-atm"));
    }
  }

  c->SaveAs(out);

  // Also write a PNG with the same basename
  TString outPng = out;
  outPng.ReplaceAll(".pdf", ".png");
  c->SaveAs(outPng);
}
