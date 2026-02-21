// osc_evolution_demo.C
// Historical evolution of neutrino oscillation parameters using hard-coded snapshots.
// (Updated aesthetics: cleaner style, colorblind-safe palette, improved label/title placement.)

#include <algorithm>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"

namespace {

struct Bounds {
  double best_fit;
  double lo1;
  double hi1;
  double lo2;
  double hi2;
  double lo3;
  double hi3;
};

TGraphAsymmErrors* BandFromBounds(const std::vector<double>& x,
                                  const std::vector<double>& ylo,
                                  const std::vector<double>& yhi)
{
  const int n = static_cast<int>(x.size());
  TGraphAsymmErrors* g = new TGraphAsymmErrors(n);
  for (int i = 0; i < n; ++i) {
    const double y  = 0.5 * (ylo[i] + yhi[i]);
    const double em = y - ylo[i];
    const double ep = yhi[i] - y;
    g->SetPoint(i, x[i], y);
    g->SetPointError(i, 0.0, 0.0, em, ep);
  }
  g->SetLineWidth(0);
  g->SetMarkerSize(0);
  g->SetFillStyle(1001);
  return g;
}

void DrawBandStack(const std::vector<double>& x,
                   const std::vector<double>& lo3, const std::vector<double>& hi3,
                   const std::vector<double>& lo2, const std::vector<double>& hi2,
                   const std::vector<double>& lo1, const std::vector<double>& hi1,
                   int baseColor)
{
  // Lighter alpha stack (reads better, especially when a band spans the full y-range).
  TGraphAsymmErrors* g3 = BandFromBounds(x, lo3, hi3);
  g3->SetFillColorAlpha(baseColor, 0.14);
  g3->Draw("3 same");

  TGraphAsymmErrors* g2 = BandFromBounds(x, lo2, hi2);
  g2->SetFillColorAlpha(baseColor, 0.25);
  g2->Draw("3 same");

  TGraphAsymmErrors* g1 = BandFromBounds(x, lo1, hi1);
  g1->SetFillColorAlpha(baseColor, 0.40);
  g1->Draw("3 same");
}

void DrawPanel(TPad* p,
               double xmin, double xmax,
               double ymin, double ymax,
               const char* latexLabel,
               int color,
               const std::vector<double>& x,
               const std::vector<double>& lo1, const std::vector<double>& hi1,
               const std::vector<double>& lo2, const std::vector<double>& hi2,
               const std::vector<double>& lo3, const std::vector<double>& hi3,
               const std::vector<double>& bf,
               bool showXLabels)
{
  p->cd();

  TH1* fr = p->DrawFrame(xmin, ymin, xmax, ymax);
  fr->SetTitle("");
  fr->GetXaxis()->SetNdivisions(505);
  fr->GetYaxis()->SetNdivisions(505);

  // Axes cosmetics: consistent fonts + slightly smaller tick labels (more whitespace).
  fr->GetXaxis()->SetLabelFont(42);
  fr->GetYaxis()->SetLabelFont(42);
  fr->GetXaxis()->SetLabelOffset(0.012);
  fr->GetYaxis()->SetLabelOffset(0.010);
  fr->GetXaxis()->SetTickLength(0.030);
  fr->GetYaxis()->SetTickLength(0.020);

  fr->GetYaxis()->SetLabelSize(0.13);
  fr->GetYaxis()->SetTitleSize(0.0);

  if (showXLabels) {
    fr->GetXaxis()->SetLabelSize(0.13);
    fr->GetXaxis()->SetTitleSize(0.0); // draw x-title manually for better placement control
  } else {
    fr->GetXaxis()->SetLabelSize(0.0);
    fr->GetXaxis()->SetTitleSize(0.0);
  }

  DrawBandStack(x, lo3, hi3, lo2, hi2, lo1, hi1, color);

  // Best-fit line + markers (makes the discrete snapshot years obvious).
  TGraph* gBF = new TGraph(static_cast<int>(x.size()));
  for (int i = 0; i < static_cast<int>(x.size()); ++i) {
    gBF->SetPoint(i, x[i], bf[i]);
  }
  gBF->SetLineColor(color);
  gBF->SetLineWidth(3);
  gBF->SetMarkerStyle(20);
  gBF->SetMarkerSize(0.75);
  gBF->SetMarkerColor(color);
  gBF->Draw("LP same");

  // ----- Panel y-label placement -----
  // Place the parameter label centered in the *left margin* and centered vertically
  // in the *inner plotting area* (so it doesn't jump around when pad margins differ).
  const double xLab = 0.50 * p->GetLeftMargin();
  const double yInnerLo = p->GetBottomMargin();
  const double yInnerHi = 1.0 - p->GetTopMargin();
  const double yLab = 0.50 * (yInnerLo + yInnerHi);

  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);
  t.SetTextSize(0.20);
  t.SetTextAlign(22); // centered
  t.DrawLatex(xLab, yLab, latexLabel);

  // ----- Global x-label ("Year") placement -----
  // Draw in bottom margin, centered under the axis.
  if (showXLabels) {
    const double xC = p->GetLeftMargin() + 0.50 * (1.0 - p->GetLeftMargin() - p->GetRightMargin());
    const double yT = 0.30 * p->GetBottomMargin();
    TLatex tx;
    tx.SetNDC(true);
    tx.SetTextFont(42);
    tx.SetTextSize(0.18);
    tx.SetTextAlign(22);
    tx.DrawLatex(xC, yT, "Year");
  }
}

void FillVectors(const std::vector<Bounds>& src,
                 std::vector<double>& bf,
                 std::vector<double>& lo1, std::vector<double>& hi1,
                 std::vector<double>& lo2, std::vector<double>& hi2,
                 std::vector<double>& lo3, std::vector<double>& hi3)
{
  const int n = static_cast<int>(src.size());
  bf.resize(n);
  lo1.resize(n); hi1.resize(n);
  lo2.resize(n); hi2.resize(n);
  lo3.resize(n); hi3.resize(n);

  for (int i = 0; i < n; ++i) {
    bf[i] = src[i].best_fit;
    lo1[i] = src[i].lo1; hi1[i] = src[i].hi1;
    lo2[i] = src[i].lo2; hi2[i] = src[i].hi2;
    lo3[i] = src[i].lo3; hi3[i] = src[i].hi3;
  }
}

} // namespace

void osc_evolution_demo()
{
  // Global style cleanup (keeps plots looking modern and consistent).
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");

  const double xmin = 1998.0;
  const double xmax = 2019.0;

  // Colorblind-safe (Okabe–Ito) palette.
  const int cS23   = TColor::GetColor("#CC79A7"); // purple
  const int cDm31  = TColor::GetColor("#E69F00"); // orange
  const int cS12   = TColor::GetColor("#D55E00"); // vermillion
  const int cDm21  = TColor::GetColor("#009E73"); // bluish green
  const int cS13   = TColor::GetColor("#56B4E9"); // sky blue
  const int cDelta = TColor::GetColor("#0072B2"); // blue

  // Snapshot years hard-coded from literature points in the prompt.
  std::vector<double> yr = {1998, 2001, 2003, 2005, 2006, 2008, 2010, 2011, 2012, 2014, 2018};

  // s23 = sin^2(theta23)
  std::vector<Bounds> s23 = {
    {0.5000, 0.3707, 0.6293, 0.3231, 0.6769, 0.2879, 0.7121},
    {0.5152, 0.4055, 0.6249, 0.3431, 0.6873, 0.2806, 0.7561},
    {0.5000, 0.4333, 0.5667, 0.3667, 0.6333, 0.3000, 0.7000},
    {0.5000, 0.4138, 0.5862, 0.3862, 0.6138, 0.3586, 0.6414},
    {0.5000, 0.4500, 0.5500, 0.4100, 0.5900, 0.3700, 0.6500},
    {0.5000, 0.4400, 0.5700, 0.3900, 0.6300, 0.3600, 0.6700},
    {0.4616, 0.4115, 0.5401, 0.3743, 0.5931, 0.3372, 0.6462},
    {0.4616, 0.4013, 0.5820, 0.3693, 0.6140, 0.3372, 0.6462},
    {0.6130, 0.4000, 0.6350, 0.3800, 0.6600, 0.3600, 0.6800},
    {0.4250, 0.4030, 0.4410, 0.3810, 0.6150, 0.3530, 0.6410},
    {0.5470, 0.5170, 0.5670, 0.4360, 0.6020, 0.4180, 0.6270}
  };

  // |Delta m^2_31| in 10^-3 eV^2
  std::vector<Bounds> dm31 = {
    {2.200, 1.286, 3.114, 0.893, 4.507, 0.500, 6.000},
    {3.500, 2.550, 4.450, 1.600, 5.400, 1.100, 7.300},
    {2.500, 2.067, 2.933, 1.633, 3.367, 1.200, 4.800},
    {2.450, 1.958, 2.943, 1.729, 3.171, 1.500, 3.400},
    {2.500, 2.300, 2.700, 2.100, 2.900, 1.900, 3.100},
    {2.400, 2.290, 2.520, 2.180, 2.640, 2.070, 2.750},
    {2.400, 2.280, 2.520, 2.175, 2.645, 2.070, 2.770},
    {2.400, 2.280, 2.520, 2.175, 2.645, 2.070, 2.770},
    {2.550, 2.460, 2.610, 2.380, 2.680, 2.310, 2.740},
    {2.480, 2.410, 2.530, 2.340, 2.590, 2.270, 2.650},
    {2.500, 2.470, 2.530, 2.440, 2.560, 2.420, 2.580}
  };

  // s12 = sin^2(theta12)
  std::vector<Bounds> s12 = {
    {0.2701, 0.2501, 0.2901, 0.2401, 0.3001, 0.2201, 0.3201},
    {0.2701, 0.2401, 0.3001, 0.2301, 0.3101, 0.2201, 0.3201},
    {0.3151, 0.2850, 0.3451, 0.2550, 0.3751, 0.2248, 0.4595},
    {0.3151, 0.3064, 0.3237, 0.2977, 0.3324, 0.2890, 0.3411},
    {0.3151, 0.3064, 0.3237, 0.2977, 0.3324, 0.2890, 0.3411},
    {0.3040, 0.2880, 0.3260, 0.2700, 0.3500, 0.2500, 0.3700},
    {0.3192, 0.3030, 0.3356, 0.2880, 0.3539, 0.2730, 0.3723},
    {0.3192, 0.3030, 0.3356, 0.2880, 0.3539, 0.2730, 0.3723},
    {0.3200, 0.3030, 0.3360, 0.2900, 0.3500, 0.2700, 0.3700},
    {0.3080, 0.2910, 0.3270, 0.2780, 0.3280, 0.2590, 0.3590},
    {0.3200, 0.3040, 0.3400, 0.2890, 0.3590, 0.2710, 0.3710}
  };

  // Delta m^2_21 in 10^-5 eV^2
  std::vector<Bounds> dm21 = {
    {4.700, 3.700, 5.700, 2.700, 6.700, 1.700, 7.700},
    {4.700, 3.700, 5.700, 2.700, 6.700, 1.700, 7.700},
    {7.200, 6.600, 7.800, 6.000, 8.400, 5.400, 20.000},
    {7.200, 7.042, 7.357, 6.885, 7.515, 6.727, 7.672},
    {7.200, 7.042, 7.357, 6.885, 7.515, 6.727, 7.672},
    {7.650, 7.450, 7.880, 7.250, 8.110, 7.050, 8.340},
    {7.590, 7.390, 7.790, 7.145, 7.995, 6.900, 8.200},
    {7.590, 7.390, 7.790, 7.145, 7.995, 6.900, 8.200},
    {7.620, 7.430, 7.810, 7.270, 8.010, 7.120, 8.200},
    {7.600, 7.420, 7.790, 7.250, 7.990, 7.110, 8.180},
    {7.550, 7.390, 7.750, 7.200, 7.900, 6.990, 8.180}
  };

  // s13 = sin^2(theta13)
  std::vector<Bounds> s13 = {
    {0.0000, 0.0000, 0.0600, 0.0000, 0.1000, 0.0000, 0.1400},
    {0.0000, 0.0000, 0.0412, 0.0000, 0.0576, 0.0000, 0.0741},
    {0.0000, 0.0000, 0.0333, 0.0000, 0.0444, 0.0000, 0.0556},
    {0.0000, 0.0000, 0.0233, 0.0000, 0.0311, 0.0000, 0.0389},
    {0.0000, 0.0000, 0.0467, 0.0000, 0.0933, 0.0000, 0.1400},
    {0.0100, 0.0000, 0.0260, 0.0000, 0.0400, 0.0000, 0.0560},
    {0.0160, 0.0060, 0.0260, 0.0010, 0.0360, 0.0000, 0.0460},
    {0.0283, 0.0076, 0.0757, 0.0038, 0.0927, 0.0000, 0.1098},
    {0.0246, 0.0218, 0.0275, 0.0190, 0.0300, 0.0170, 0.0330},
    {0.0234, 0.0215, 0.0254, 0.0198, 0.0282, 0.0176, 0.0309},
    {0.0216, 0.0209, 0.0224, 0.0200, 0.0232, 0.0192, 0.0241}
  };

  // delta/pi
  std::vector<Bounds> delta = {
    {0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {1.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {1.000, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {0.800, 0.000, 2.000, 0.000, 2.000, 0.000, 2.000},
    {1.340, 0.960, 1.980, 0.000, 2.000, 0.000, 2.000},
    {1.320, 1.170, 1.530, 1.020, 1.970, 0.870, 2.000}
  };

  std::vector<double> bf1, lo1_1, hi1_1, lo2_1, hi2_1, lo3_1, hi3_1;
  std::vector<double> bf2, lo1_2, hi1_2, lo2_2, hi2_2, lo3_2, hi3_2;
  std::vector<double> bf3, lo1_3, hi1_3, lo2_3, hi2_3, lo3_3, hi3_3;
  std::vector<double> bf4, lo1_4, hi1_4, lo2_4, hi2_4, lo3_4, hi3_4;
  std::vector<double> bf5, lo1_5, hi1_5, lo2_5, hi2_5, lo3_5, hi3_5;
  std::vector<double> bf6, lo1_6, hi1_6, lo2_6, hi2_6, lo3_6, hi3_6;

  FillVectors(s23,   bf1, lo1_1, hi1_1, lo2_1, hi2_1, lo3_1, hi3_1);
  FillVectors(dm31,  bf2, lo1_2, hi1_2, lo2_2, hi2_2, lo3_2, hi3_2);
  FillVectors(s12,   bf3, lo1_3, hi1_3, lo2_3, hi2_3, lo3_3, hi3_3);
  FillVectors(dm21,  bf4, lo1_4, hi1_4, lo2_4, hi2_4, lo3_4, hi3_4);
  FillVectors(s13,   bf5, lo1_5, hi1_5, lo2_5, hi2_5, lo3_5, hi3_5);
  FillVectors(delta, bf6, lo1_6, hi1_6, lo2_6, hi2_6, lo3_6, hi3_6);

  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->cd();

  const int nPads = 6;
  const double y0 = 0.08;
  const double y1 = 0.98;
  const double h  = (y1 - y0) / nPads;

  std::vector<TPad*> pads(nPads, nullptr);
  for (int i = 0; i < nPads; ++i) {
    const double ylow  = y0 + (nPads - 1 - i) * h;
    const double yhigh = ylow + h;
    pads[i] = new TPad(Form("p%d", i), "", 0.0, ylow, 1.0, yhigh);

    // More breathing room on the left (parameter labels + y tick labels).
    pads[i]->SetLeftMargin(0.18);
    pads[i]->SetRightMargin(0.03);

    pads[i]->SetTopMargin(i == 0 ? 0.06 : 0.02);
    pads[i]->SetBottomMargin(i == nPads - 1 ? 0.30 : 0.02);

    pads[i]->SetBorderMode(0);
    pads[i]->SetFrameBorderMode(0);
    pads[i]->Draw();
  }

  DrawPanel(pads[0], xmin, xmax, 0.20, 0.85, "s^{2}_{23}",                 cS23,
            yr, lo1_1, hi1_1, lo2_1, hi2_1, lo3_1, hi3_1, bf1, false);

  DrawPanel(pads[1], xmin, xmax, 1.0, 6.2,  "|#Delta m^{2}_{31}|/10^{-3}",  cDm31,
            yr, lo1_2, hi1_2, lo2_2, hi2_2, lo3_2, hi3_2, bf2, false);

  DrawPanel(pads[2], xmin, xmax, 0.15, 0.50, "s^{2}_{12}",                 cS12,
            yr, lo1_3, hi1_3, lo2_3, hi2_3, lo3_3, hi3_3, bf3, false);

  DrawPanel(pads[3], xmin, xmax, 2.0, 10.0,  "#Delta m^{2}_{21}/10^{-5}",   cDm21,
            yr, lo1_4, hi1_4, lo2_4, hi2_4, lo3_4, hi3_4, bf4, false);

  DrawPanel(pads[4], xmin, xmax, 0.0, 0.15,  "s^{2}_{13}",                 cS13,
            yr, lo1_5, hi1_5, lo2_5, hi2_5, lo3_5, hi3_5, bf5, false);

  DrawPanel(pads[5], xmin, xmax, 0.0, 2.0,   "#delta_{CP}/#pi",            cDelta,
            yr, lo1_6, hi1_6, lo2_6, hi2_6, lo3_6, hi3_6, bf6, true);

  c->SaveAs("osc_evolution_demo.pdf");
}
