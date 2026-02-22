// plot/macro/plotMinimal.C
//
// Minimal plotting macro for heron that uses framework helpers.
//
// Usage:
//   heron --set template macro plotMinimal.C
//   heron --set template macro plotMinimal.C 'plotMinimal("./scratch/out/template/event/events.root")'

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"

#include "framework/core/include/EventColumnProvider.hh"
#include "framework/modules/plot/include/PlotEnv.hh"
#include "macros/include/MacroIO.hh"

bool plotMinimal(const char* input_file = "./scratch/out/template/event/events.root",
                 const char* tree_name = "events",
                 const char* branch_name = "reco_nslice",
                 const char* output_name = "plotMinimal",
                 const char* columns_tsv_path = "") {
  if (!input_file || !tree_name || !branch_name || !output_name || !columns_tsv_path) {
    std::cerr << "plotMinimal: invalid NULL argument" << std::endl;
    return false;
  }

  if (!heron::macro::validate_root_input_path(input_file)) {
    std::cerr << "plotMinimal: invalid input ROOT file path" << std::endl;
    return false;
  }

  EventColumnProvider provider(columns_tsv_path);
  bool known_column = (provider.columns().empty() && std::string(branch_name) == "reco_nslice");
  for (const std::string& column : provider.columns()) {
    if (column == branch_name) {
      known_column = true;
      break;
    }
  }

  if (!known_column) {
    std::cerr << "plotMinimal: branch not declared in configured event columns: " << branch_name << std::endl;
    return false;
  }

  std::unique_ptr<TFile> p_file = heron::macro::open_root_file_read(input_file);
  if (!p_file) {
    std::cerr << "plotMinimal: could not open file: " << input_file << std::endl;
    return false;
  }

  TTree* tree = dynamic_cast<TTree*>(p_file->Get(tree_name));
  if (!tree) {
    std::cerr << "plotMinimal: tree not found: " << tree_name << std::endl;
    return false;
  }

  TH1D* hist = new TH1D("h_minimal", branch_name, 40, 0.0, 40.0);
  hist->SetLineWidth(2);

  const std::string draw_expr = std::string(branch_name) + ">>h_minimal";
  tree->Draw(draw_expr.c_str(), "", "goff");

  TCanvas* c = new TCanvas("c_minimal", "c_minimal", nu::kCanvasWidth, nu::kCanvasHeight);
  c->SetGrid();
  hist->GetXaxis()->SetTitle(branch_name);
  hist->GetYaxis()->SetTitle("Entries");
  hist->Draw("HIST");

  const std::filesystem::path output_file = nu::resolve_output_file(output_name, "plotMinimal", "pdf");
  c->SaveAs(output_file.string().c_str());

  delete c;
  delete hist;

  std::cout << "plotMinimal: wrote " << output_file << std::endl;
  return true;
}
