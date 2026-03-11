#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
#endif

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <cctype>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "SampleCLI.hh"
#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "PlotEnv.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"

using namespace nu;

namespace {

bool is_simple_identifier(const std::string& expr) {
  if (expr.empty()) return false;
  const unsigned char first = static_cast<unsigned char>(expr.front());
  if (!(std::isalpha(first) || expr.front() == '_')) return false;

  for (size_t i = 1; i < expr.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(expr[i]);
    if (!(std::isalnum(c) || expr[i] == '_')) return false;
  }

  return true;
}

ROOT::RDF::RNode filter_by_mask(ROOT::RDF::RNode node,
                                std::shared_ptr<const std::vector<char>> mask) {
  if (!mask) return node;

  return node.Filter(
      [mask](int sid) {
        return sid >= 0 && sid < static_cast<int>(mask->size()) &&
               (*mask)[static_cast<std::size_t>(sid)];
      },
      {"sample_id"});
}

ROOT::RDF::RNode filter_not_by_mask(ROOT::RDF::RNode node,
                                    std::shared_ptr<const std::vector<char>> mask) {
  if (!mask) return node;

  return node.Filter(
      [mask](int sid) {
        return !(sid >= 0 && sid < static_cast<int>(mask->size()) &&
                 (*mask)[static_cast<std::size_t>(sid)]);
      },
      {"sample_id"});
}

void draw_raw_stack_plot(Plotter& plotter,
                         std::vector<const Entry*>& mc,
                         std::vector<const Entry*>& data,
                         bool include_data) {
  auto& opt = plotter.options();

  TH1DModel spec = make_spec("inf_score_0", 50, -15.0, 15.0, "unit_w");
  spec.id = "h_inf_score_0_raw_counts";
  spec.name = "inf_score_0_raw_counts";
  spec.sel = Preset::Empty;

  opt.x_title = "Inference score [0]";
  opt.y_title = "Entries";
  opt.normalise_by_bin_width = false;

  if (include_data)
    plotter.draw_stack(spec, mc, data);
  else
    plotter.draw_stack(spec, mc);
}

}  // namespace

int plot_model_logit_raw_stack(const std::string& samples_tsv = "",
                               const char* extra_sel = "sel_muon",
                               bool use_logy = true,
                               bool include_data = false) {
  const std::string extra_sel_expr = (extra_sel != nullptr) ? extra_sel : "";
  const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
  std::cout << "[plot_model_logit_raw_stack] input=" << list_path << "\n";

  if (!looks_like_event_list_root(list_path)) {
    std::cerr << "[plot_model_logit_raw_stack] input is not an event list ROOT file: "
              << list_path << "\n";
    return 2;
  }

  ROOT::EnableImplicitMT();

  EventListIO el(list_path);
  ROOT::RDataFrame rdf = el.rdf();

  auto mask_ext = el.mask_for_ext();
  auto mask_mc = el.mask_for_mc_like();
  auto mask_data = el.mask_for_data();

  ROOT::RDF::RNode base =
      rdf.Define("inf_score_0",
                 [](const ROOT::RVec<float>& scores) {
                   if (scores.empty()) return -1.0f;
                   return scores[0];
                 },
                 {"inf_scores"})
          .Define("unit_w", []() { return 1.0; });

  ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
  ROOT::RDF::RNode node_mc = filter_not_by_mask(filter_by_mask(base, mask_mc), mask_ext);
  ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

  std::vector<Entry> entries;
  entries.reserve(include_data ? 3 : 2);

  std::vector<const Entry*> mc;
  std::vector<const Entry*> data;

  ProcessorEntry rec_mc;
  rec_mc.source = Type::kMC;
  entries.emplace_back(make_entry(std::move(node_mc), rec_mc));
  Entry& e_mc = entries.back();
  mc.push_back(&e_mc);

  ProcessorEntry rec_ext;
  rec_ext.source = Type::kExt;
  entries.emplace_back(make_entry(std::move(node_ext), rec_ext));
  Entry& e_ext = entries.back();
  mc.push_back(&e_ext);

  Entry* p_data = nullptr;
  if (include_data) {
    ProcessorEntry rec_data;
    rec_data.source = Type::kData;
    entries.emplace_back(make_entry(std::move(node_data), rec_data));
    p_data = &entries.back();
    data.push_back(p_data);
  }

  if (!extra_sel_expr.empty()) {
    const bool named_column = rdf.HasColumn(extra_sel_expr);

    if (named_column) {
      e_mc.selection.nominal.node =
          e_mc.selection.nominal.node.Filter([](bool pass) { return pass; }, {extra_sel_expr});
      e_ext.selection.nominal.node =
          e_ext.selection.nominal.node.Filter([](bool pass) { return pass; }, {extra_sel_expr});
      if (p_data != nullptr) {
        p_data->selection.nominal.node =
            p_data->selection.nominal.node.Filter([](bool pass) { return pass; }, {extra_sel_expr});
      }
    } else if (is_simple_identifier(extra_sel_expr)) {
      std::cerr << "[plot_model_logit_raw_stack] selection column '" << extra_sel_expr
                << "' is missing; skipping extra selection.\n";
    } else {
      e_mc.selection.nominal.node = e_mc.selection.nominal.node.Filter(extra_sel_expr);
      e_ext.selection.nominal.node = e_ext.selection.nominal.node.Filter(extra_sel_expr);
      if (p_data != nullptr) {
        p_data->selection.nominal.node = p_data->selection.nominal.node.Filter(extra_sel_expr);
      }
    }
  }

  Plotter plotter;
  auto& opt = plotter.options();
  opt.use_log_y = use_logy;
  opt.legend_on_top = true;
  opt.annotate_numbers = true;
  opt.overlay_signal = false;      // avoid internal signal rescaling
  opt.show_ratio = false;          // raw counts vs unweighted MC ratio is usually not meaningful
  opt.show_ratio_band = false;
  opt.normalise_by_bin_width = false;
  opt.y_title = "Entries";
  opt.run_numbers = {"1"};
  opt.image_format = "pdf";

  const double pot_data = el.total_pot_data();
  const double pot_mc = el.total_pot_mc();
  opt.total_protons_on_target = (pot_data > 0.0 ? pot_data : pot_mc);
  opt.beamline = el.beamline_label();

  draw_raw_stack_plot(plotter, mc, data, include_data);

  std::cout << "[plot_model_logit_raw_stack] done\n";
  return 0;
}
