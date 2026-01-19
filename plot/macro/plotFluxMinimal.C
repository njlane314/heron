#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TString.h"
#include "TColor.h"
#include "TSystem.h"
#include "../include/PlotEnv.hh"
#include "../include/Plotter.hh"
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// Flux unit normalization:
// Most dk2nu/PPFX flux histograms are stored in  # / (1e9 POT) / m^2 / GeV.
// We want to report/plot everything in          # / (6e20 POT) / cm^2 / GeV.
// => multiply hist contents by:
//    (6e20 / 1e9) [POT]  ×  (1 / 1e4) [m^-2 -> cm^-2]  =  6e7.
// If your file is per 1 POT (not 1e9 POT), set kPOT_IN_FILE = 1.0.
// ---------------------------------------------------------------------------
static constexpr double kPOT_TARGET   = 6.0e20;  // what we print/show
static constexpr double kPOT_IN_FILE  = 1.0e9;   // file convention (typical dk2nu/ppfx)
static constexpr double kM2_TO_CM2    = 1.0e-4;  // m^-2 -> cm^-2
static constexpr double kUNIT_SCALE   = (kPOT_TARGET / kPOT_IN_FILE) * kM2_TO_CM2; // = 6e7 by default

static inline void scale_flux_to_release_units(TH1* h){
  if(h) h->Scale(kUNIT_SCALE);
}
static inline void scale_flux_to_release_units(TH1* h1, TH1* h2, TH1* h3, TH1* h4){
  scale_flux_to_release_units(h1);
  scale_flux_to_release_units(h2);
  scale_flux_to_release_units(h3);
  scale_flux_to_release_units(h4);
}
// ---------------------------------------------------------------------------

static void set_global_style(){
  nuxsec::plot::Plotter{}.set_global_style();
}

static TLegend* build_flux_legend_like_stacked(TPad* p_leg, TH1* h_numu, TH1* h_anumu, TH1* h_nue, TH1* h_anue, double split, double s_numu, double s_anumu, double s_nue, double s_anue, double s_tot){
  p_leg->cd();
  TLegend* L=new TLegend(0.12,0.00,0.95,0.75);
  L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42);
  const int n_entries=4; int n_cols=(n_entries>4)?3:2;
  L->SetNColumns(n_cols); L->SetColumnSeparation(0.08); L->SetEntrySeparation(0.00); L->SetMargin(0.25);
  const double s_main=0.045; const double s_leg=s_main*(split/(1.0-split)); L->SetTextSize(s_leg);
  auto pct=[&](double x){ return (s_tot>0?100.0*x/s_tot:0.0); };
  L->AddEntry(h_numu ,Form("#nu_{#mu} (%.1f%%)",      pct(s_numu)) ,"l");
  L->AddEntry(h_anumu,Form("#bar{#nu}_{#mu} (%.1f%%)",pct(s_anumu)),"l");
  L->AddEntry(h_nue  ,Form("#nu_{e} (%.1f%%)",        pct(s_nue))  ,"l");
  L->AddEntry(h_anue ,Form("#bar{#nu}_{e} (%.1f%%)",  pct(s_anue)) ,"l");
  return L;
}

static void style_line(TH1* h,int col,int ls){ h->SetLineColor(col); h->SetLineStyle(ls); h->SetLineWidth(2); h->SetMarkerSize(0); }

static double integral_in(double xmin, double xmax, const TH1* h, bool width=false){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  return width ? h->Integral(bmin, bmax, "width") : h->Integral(bmin, bmax);
}

static double integral_and_error(double xmin, double xmax, const TH1* h,
                                 double& err, bool width=true){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  err = 0.0;
  return h->IntegralAndError(bmin, bmax, err, width ? "width" : "");
}

static void print_flux_window_integrals(const char* tag,
                                        const TH1* h_numu,  const TH1* h_anumu,
                                        const TH1* h_nue,   const TH1* h_anue,
                                        double xmin, double xmax){
  double e_numu, e_anumu, e_nue, e_anue;
  const double i_numu  = integral_and_error(xmin, xmax, h_numu , e_numu , /*width=*/true);
  const double i_anumu = integral_and_error(xmin, xmax, h_anumu, e_anumu, /*width=*/true);
  const double i_nue   = integral_and_error(xmin, xmax, h_nue  , e_nue  , /*width=*/true);
  const double i_anue  = integral_and_error(xmin, xmax, h_anue , e_anue , /*width=*/true);
  const double i_tot   = i_numu + i_anumu + i_nue + i_anue;

  auto pct = [&](double x){ return (i_tot > 0 ? 100.0*x/i_tot : 0.0); };
  const int binwMeV = (int)std::lround(h_numu->GetXaxis()->GetBinWidth(1) * 1000.0);

  printf("[plotFluxMinimal/%s] Energy-window integrals (%.2f–%.2f GeV)\n", tag, xmin, xmax);
  printf("  (TH1::Integral(...,\"width\") used: area = sum(content × binWidth).\n");
  printf("   Units: # / (6×10^20 POT) / cm^2; histogram y-axis is per %d MeV.)\n", binwMeV);

  printf("    nu_mu      : % .6e  ± %.2e   (%.1f%%)\n", i_numu , e_numu , pct(i_numu ));
  printf("    anti-nu_mu : % .6e  ± %.2e   (%.1f%%)\n", i_anumu, e_anumu, pct(i_anumu));
  printf("    nu_e       : % .6e  ± %.2e   (%.1f%%)\n", i_nue  , e_nue  , pct(i_nue  ));
  printf("    anti-nu_e  : % .6e  ± %.2e   (%.1f%%)\n", i_anue , e_anue , pct(i_anue ));
  printf("    TOTAL      : % .6e\n\n", i_tot);
}

static void auto_logy_limits_range(TH1* frame, std::initializer_list<TH1*> hs, double xmin, double xmax){
  double mn=std::numeric_limits<double>::infinity(), mx=0.0;
  for(TH1* h:hs){
    int bmin=std::max(1,h->GetXaxis()->FindFixBin(xmin+1e-9));
    int bmax=std::min(h->GetNbinsX(),h->GetXaxis()->FindFixBin(xmax-1e-9));
    for(int b=bmin;b<=bmax;++b){ double y=h->GetBinContent(b); if(y>0&&y<mn) mn=y; if(y>mx) mx=y; }
  }
  if(!std::isfinite(mn)) mn=1e-18; if(mx<=0.0) mx=1.0;
  frame->SetMinimum(std::max(1e-30,mn*0.8)); frame->SetMaximum(mx*6.0);
}

// -------------------- NEW: CSV table writers for x-sec release --------------------
struct FluxHists {
  TH1D* numu    = nullptr;
  TH1D* numubar = nullptr;
  TH1D* nue     = nullptr;
  TH1D* nuebar  = nullptr;
};

static bool load_flux_hists(const char* file, FluxHists& H, const char* tag){
  TFile f(file,"READ");
  auto fetch = [&](const char* a, const char* b)->TH1D*{
    TH1D* h=(TH1D*)f.Get(a);
    if(!h) h=(TH1D*)f.Get(b);
    return h ? (TH1D*)h->Clone() : nullptr;
  };
  H.numu    = fetch("numu/Detsmear/numu_CV_AV_TPC_5MeV_bin",     "numu/Detsmear/numu_CV_AV_TPC");
  H.numubar = fetch("numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin","numubar/Detsmear/numubar_CV_AV_TPC");
  H.nue     = fetch("nue/Detsmear/nue_CV_AV_TPC_5MeV_bin",       "nue/Detsmear/nue_CV_AV_TPC");
  H.nuebar  = fetch("nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin", "nuebar/Detsmear/nuebar_CV_AV_TPC");
  if(H.numu   ) H.numu   ->SetDirectory(0);
  if(H.numubar) H.numubar->SetDirectory(0);
  if(H.nue    ) H.nue    ->SetDirectory(0);
  if(H.nuebar ) H.nuebar ->SetDirectory(0);
  f.Close();
  // Normalize all histograms to # / (6e20 POT) / cm^2 / GeV
  scale_flux_to_release_units(H.numu, H.numubar, H.nue, H.nuebar);
  if(!H.numu || !H.numubar || !H.nue || !H.nuebar){
    printf("[tables/%s] missing *_CV_AV_TPC in Detsmear\n", tag);
    return false;
  }
  return true;
}

static void write_bin_edges_csv(const TH1* ref, const char* outdir){
  gSystem->mkdir(Form("%s/bins", outdir), true);
  std::ofstream fo(Form("%s/bins/bin_edges.csv", outdir));
  fo << "bin_low_GeV,bin_high_GeV\n";
  fo << std::scientific << std::setprecision(8);
  for(int b=1;b<=ref->GetNbinsX();++b){
    fo << ref->GetXaxis()->GetBinLowEdge(b) << ","
       << ref->GetXaxis()->GetBinUpEdge(b)  << "\n";
  }
}

static void write_flux_spectrum_csv(const char* outdir, const char* mode,
                                    const TH1* a, const TH1* b,
                                    const TH1* c, const TH1* d){
  gSystem->mkdir(Form("%s/flux", outdir), true);
  TString m(mode);
  m.ToLower();
  std::ofstream fo(Form("%s/flux/numi_%s_flux.csv", outdir, m.Data()));
  fo << "# Columns are per energy bin. Values are INTEGRATED per bin "
        "(content × binWidth): # / (6×10^20 POT) / cm^2.\n";
  fo << "E_low_GeV,E_high_GeV,bin_width_GeV,"
        "numu,numu_err,numubar,numubar_err,nue,nue_err,nuebar,nuebar_err\n";
  fo << std::scientific << std::setprecision(8);

  const TH1* hs[4] = {a,b,c,d};
  for(int bidx=1; bidx<=a->GetNbinsX(); ++bidx){
    double lo = a->GetXaxis()->GetBinLowEdge(bidx);
    double hi = a->GetXaxis()->GetBinUpEdge(bidx);
    double w  = hi - lo;

    double v[4], e[4];
    for(int k=0;k<4;++k){
      v[k] = hs[k]->GetBinContent(bidx) * w;  // integrated per bin
      e[k] = hs[k]->GetBinError(bidx)   * w;  // stat error integrated
    }
    fo << lo << "," << hi << "," << w << ","
       << v[0] << "," << e[0] << ","
       << v[1] << "," << e[1] << ","
       << v[2] << "," << e[2] << ","
       << v[3] << "," << e[3] << "\n";
  }
}

static void write_window_integrals_csv(const char* outdir, const char* mode,
                                       const char* tag,
                                       const TH1* a, const TH1* b,
                                       const TH1* c, const TH1* d,
                                       double Emin, double Emax){
  gSystem->mkdir(Form("%s/flux", outdir), true);
  TString m(mode);
  m.ToLower();
  std::ofstream fo(Form("%s/flux/numi_%s_integrals_%s.csv", outdir, m.Data(), tag));
  fo << "# Integrated flux in window (" << std::fixed << std::setprecision(3)
     << Emin << "," << Emax << ") GeV.\n";
  fo << "# Units: # / (6×10^20 POT) / cm^2.\n";
  fo << "mode,window,species,E_min_GeV,E_max_GeV,integral,stat_err,fraction\n";

  double e[4], v[4];
  v[0]=integral_and_error(Emin,Emax,a,e[0],true);
  v[1]=integral_and_error(Emin,Emax,b,e[1],true);
  v[2]=integral_and_error(Emin,Emax,c,e[2],true);
  v[3]=integral_and_error(Emin,Emax,d,e[3],true);
  double tot=v[0]+v[1]+v[2]+v[3];
  const char* nm[4]={"numu","numubar","nue","nuebar"};

  for(int k=0;k<4;++k){
    fo << mode << "," << tag << "," << nm[k] << ","
       << std::fixed << std::setprecision(3) << Emin << ","
       << std::fixed << std::setprecision(3) << Emax << ","
       << std::scientific << std::setprecision(8) << v[k] << ","
       << std::scientific << std::setprecision(8) << e[k] << ","
       << std::fixed << std::setprecision(3) << (tot>0?100.0*v[k]/tot:0.0) << "\n";
  }
  double etot = std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3]); // assumes uncorrelated species
  fo << mode << "," << tag << ",total,"
     << std::fixed << std::setprecision(3) << Emin << ","
     << std::fixed << std::setprecision(3) << Emax << ","
     << std::scientific << std::setprecision(8) << tot << ","
     << std::scientific << std::setprecision(8) << etot << ",100.000\n";
}

static void make_tables_one(const char* file, const char* mode, const char* outdir){
  FluxHists H{};
  if(!load_flux_hists(file, H, mode)) return;

  // write binning and spectra
  write_bin_edges_csv(H.numu, outdir);
  write_flux_spectrum_csv(outdir, mode, H.numu, H.numubar, H.nue, H.nuebar);

  // integrated windows: full and analysis (adjust as needed)
  const double Efull_min=0.0,  Efull_max=10.0;
  const double Eana_min =0.25, Eana_max =10.0;
  write_window_integrals_csv(outdir, mode, "full",     H.numu, H.numubar, H.nue, H.nuebar, Efull_min, Efull_max);
  write_window_integrals_csv(outdir, mode, "analysis", H.numu, H.numubar, H.nue, H.nuebar, Eana_min,  Eana_max);

  delete H.numu; delete H.numubar; delete H.nue; delete H.nuebar;
}
// ------------------ END NEW TABLE WRITERS ------------------

static void draw_one(const char* file,const char* tag,const char* out){
  const double Emin=0.0, Emax=10.0, split=0.85;
  const double EanaMin = 0.25;  // analysis window used everywhere
  TFile f(file,"READ");
  TH1D* a=(TH1D*)f.Get("numu/Detsmear/numu_CV_AV_TPC_5MeV_bin"); if(!a) a=(TH1D*)f.Get("numu/Detsmear/numu_CV_AV_TPC");
  TH1D* b=(TH1D*)f.Get("numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin"); if(!b) b=(TH1D*)f.Get("numubar/Detsmear/numubar_CV_AV_TPC");
  TH1D* c=(TH1D*)f.Get("nue/Detsmear/nue_CV_AV_TPC_5MeV_bin"); if(!c) c=(TH1D*)f.Get("nue/Detsmear/nue_CV_AV_TPC");
  TH1D* d=(TH1D*)f.Get("nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin"); if(!d) d=(TH1D*)f.Get("nuebar/Detsmear/nuebar_CV_AV_TPC");
  if(!a||!b||!c||!d){ printf("[plotFluxMinimal/%s] missing *_CV_AV_TPC in Detsmear\n",tag); return; }
  a=(TH1D*)a->Clone("h_numu"); b=(TH1D*)b->Clone("h_anumu"); c=(TH1D*)c->Clone("h_nue"); d=(TH1D*)d->Clone("h_anue");
  a->SetDirectory(0); b->SetDirectory(0); c->SetDirectory(0); d->SetDirectory(0); f.Close();

  scale_flux_to_release_units(a, b, c, d); // convert to #/(6e20 POT)/cm^2/GeV before any integrals/plots
  print_flux_window_integrals(tag, a, b, c, d, EanaMin, Emax);

  int CR=TColor::GetColor("#e41a1c"), CB=TColor::GetColor("#1f78b4");
  style_line(a,CR,1); style_line(c,CR,2); style_line(b,CB,1); style_line(d,CB,3);

  TCanvas canv(Form("c_%s",tag),Form("%s Mode",tag),
               nuxsec::plot::kCanvasWidth,
               nuxsec::plot::kCanvasHeight);
  TPad* p_main=new TPad("pad_main","pad_main",0.,0.00,1.,split);
  TPad* p_leg =new TPad("pad_legend","pad_legend",0.,split,1.,1.00);
  p_main->SetTopMargin(0.01); p_main->SetBottomMargin(0.12); p_main->SetLeftMargin(0.12); p_main->SetRightMargin(0.05); p_main->SetLogy();
  p_leg->SetTopMargin(0.05); p_leg->SetBottomMargin(0.01); p_leg->SetLeftMargin(0.02); p_leg->SetRightMargin(0.02);
  p_main->Draw(); p_leg->Draw();
  p_main->cd();

  TH1D* frame=new TH1D(Form("frame_%s",tag),"",100,Emin,Emax);
  int binwMeV=(int)std::lround(a->GetXaxis()->GetBinWidth(1)*1000.0);
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  frame->GetYaxis()->SetTitle(Form("#nu / 6 #times 10^{20} POT / %d MeV / cm^{2}",binwMeV));
  a->GetXaxis()->SetRangeUser(Emin,Emax); b->GetXaxis()->SetRangeUser(Emin,Emax); c->GetXaxis()->SetRangeUser(Emin,Emax); d->GetXaxis()->SetRangeUser(Emin,Emax);
  auto_logy_limits_range(frame,{a,b,c,d},Emin,Emax);
  frame->Draw("AXIS");
  a->Draw("HIST SAME"); c->Draw("HIST SAME"); b->Draw("HIST SAME"); d->Draw("HIST SAME");

  p_leg->cd();
  // Use the same area definition ("width") and the same window [0.25, 10]
  double s_numu  = integral_in(EanaMin, Emax, a, /*width=*/true);
  double s_anumu = integral_in(EanaMin, Emax, b, /*width=*/true);
  double s_nue   = integral_in(EanaMin, Emax, c, /*width=*/true);
  double s_anue  = integral_in(EanaMin, Emax, d, /*width=*/true);
  double s_tot=std::max(1e-300,s_numu+s_anumu+s_nue+s_anue);
  TLegend* L=build_flux_legend_like_stacked(p_leg,a,b,c,d,split,s_numu,s_anumu,s_nue,s_anue,s_tot);
  L->Draw();

  canv.cd(); canv.Update(); canv.Print(out);
  printf("[plotFluxMinimal] %s | bin width ≈ %.3f GeV (%d MeV)\n", tag, a->GetXaxis()->GetBinWidth(1), binwMeV);
  delete frame; delete L; delete p_main; delete p_leg; delete a; delete b; delete c; delete d;
}

void plotFluxMinimal(
  const char* fhc_file = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root",
  const char* rhc_file = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root"
){
  set_global_style();
  const std::string out_fhc = nuxsec::plot::plot_output_file("uboone_flux_FHC").string();
  const std::string out_rhc = nuxsec::plot::plot_output_file("uboone_flux_RHC").string();
  draw_one(fhc_file,"FHC",out_fhc.c_str());
  draw_one(rhc_file,"RHC",out_rhc.c_str());

  // Informative print for sanity:
  printf("[units] Histograms scaled by kUNIT_SCALE = %.3e to #/(6e20 POT)/cm^2/GeV "
         "(kPOT_TARGET=%.3e, kPOT_IN_FILE=%.3e, m^2->cm^2=%g)\n",
         kUNIT_SCALE, kPOT_TARGET, kPOT_IN_FILE, kM2_TO_CM2);
  // -------- NEW: write machine-readable flux tables for x-sec release --------
  const std::string outdir = nuxsec::plot::release_dir_path().string();
  make_tables_one(fhc_file,"FHC",outdir.c_str());
  make_tables_one(rhc_file,"RHC",outdir.c_str());
  printf("[tables] Wrote CSVs under %s/{bins,flux}/\n", outdir.c_str());
}
