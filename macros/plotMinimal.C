// plot/macro/plotMinimal.C
//
// Minimal plotting macro for muon-neutrino selection stages, using
// Dataset + StackedHist library classes.
//
// Usage:
//   heron --set template macro plotMinimal.C
//   heron --set template macro plotMinimal.C \
//     'plotMinimal("./scratch/out/template/sample/samples.tsv", "./scratch/out/template/event/events.root")'

#include <filesystem>
#include <iostream>
#include <algorithm>
#include <memory>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "framework/core/include/Dataset.hh"
#include "framework/modules/plot/include/PlottingHelper.hh"
#include "framework/modules/plot/include/StackedHist.hh"
#include "macros/include/MacroIO.hh"

bool plotMinimal(const char* sample_list = "./scratch/out/template/sample/samples.tsv",
                 const char* event_list_root = "./scratch/out/template/event/events.root",
                 const char* tree_name = "events",
                 const char* output_name = "plotMinimal") {
  if (!sample_list || !event_list_root || !tree_name || !output_name) {
    std::cerr << "plotMinimal: invalid NULL argument" << std::endl;
    return false;
  }

  if (!heron::macro::validate_root_input_path(event_list_root)) {
    std::cerr << "plotMinimal: invalid input ROOT file path" << std::endl;
    return false;
  }

  Dataset dataset;
  try {
    dataset = Dataset::load(sample_list);
  } catch (const std::exception& e) {
    std::cerr << "plotMinimal: failed to load dataset: " << e.what() << std::endl;
    return false;
  }

  ROOT::RDataFrame df(tree_name, event_list_root);
  const auto base = SelectionService::decorate(df);

  std::vector<char> mc_mask(dataset.inputs().size(), 0);
  std::vector<char> data_mask(dataset.inputs().size(), 0);
  for (size_t i = 0; i < dataset.inputs().size(); ++i) {
    const auto origin = dataset.inputs()[i].sample.origin;
    if (origin == SampleIO::SampleOrigin::kData || origin == SampleIO::SampleOrigin::kEXT) {
      data_mask[i] = 1;
      continue;
    }
    mc_mask[i] = 1;
  }

  const auto p_mc_mask = std::make_shared<const std::vector<char>>(mc_mask);
  const auto p_data_mask = std::make_shared<const std::vector<char>>(data_mask);

  auto mc_node = nu::filter_by_sample_mask(base, p_mc_mask);
  auto data_node = nu::filter_by_sample_mask(base, p_data_mask);

  const auto stage_expr =
      "sel_muon ? 5 : (sel_topology ? 4 : (sel_fiducial ? 3 : (sel_slice ? 2 : (sel_trigger ? 1 : 0))))";
  mc_node = mc_node.Define("sel_stage", stage_expr);
  data_node = data_node.Define("sel_stage", stage_expr);

  std::vector<nu::Entry> mc_entries;
  std::vector<nu::Entry> data_entries;
  std::vector<const nu::Entry*> mc_ptrs;
  std::vector<const nu::Entry*> data_ptrs;

  if (!mc_mask.empty()) {
    const bool have_mc = std::any_of(mc_mask.begin(), mc_mask.end(), [](char v) { return v != 0; });
    if (have_mc) {
      SelectionEntry sel{Type::kMC, Frame{std::move(mc_node)}};
      mc_entries.push_back(nu::Entry{std::move(sel), 0.0, 0.0, std::string(), std::string()});
    }
  }

  if (!data_mask.empty()) {
    const bool have_data = std::any_of(data_mask.begin(), data_mask.end(), [](char v) { return v != 0; });
    if (have_data) {
      SelectionEntry sel{Type::kData, Frame{std::move(data_node)}};
      data_entries.push_back(nu::Entry{std::move(sel), 0.0, 0.0, std::string(), std::string()});
    }
  }

  for (auto& e : mc_entries) mc_ptrs.push_back(&e);
  for (auto& e : data_entries) data_ptrs.push_back(&e);

  nu::TH1DModel spec;
  spec.id = output_name;
  spec.name = "Muon neutrino selection stage";
  spec.title = "Muon neutrino selection stage";
  spec.expr = "sel_stage";
  spec.nbins = 6;
  spec.xmin = -0.5;
  spec.xmax = 5.5;
  spec.sel = Preset::Empty;

  nu::Options opt;
  opt.out_dir = std::filesystem::path(output_name).parent_path().string();
  if (opt.out_dir.empty()) opt.out_dir = ".";
  opt.image_format = "pdf";
  opt.show_ratio = false;
  opt.show_ratio_band = false;
  opt.normalise_by_bin_width = false;
  opt.x_title = "Selection stage (0=All, 1=Trigger, 2=Slice, 3=Fiducial, 4=Topology, 5=Muon)";
  opt.y_title = "Entries";

  nu::StackedHist plot(spec, opt, mc_ptrs, data_ptrs);
  plot.draw_and_save(opt.image_format);

  std::cout << "plotMinimal: wrote " << output_name << ".pdf" << std::endl;
  return true;
}
