// plot/macro/stack_samples.C
//
// Run with:
//   ./nuxsec macro stack_samples.C
//   ./nuxsec macro stack_samples.C 'stack_samples("reco_neutrino_vertex_sce_z")'
//   ./nuxsec macro stack_samples.C 'stack_samples("reco_neutrino_vertex_sce_z","/path/to/samples.tsv",50,0,1200,"w_nominal","true",false)'
//
// Notes:
//   - This macro loads aggregated samples (samples.tsv -> SampleIO -> original analysis tree)
//   - It runs your analysis column derivations so that "analysis_channels" exists for stacking.
//   - The stack is grouped by "analysis_channels"; expr controls the x-axis variable only.
//   - MC yields are scaled by w_nominal unless an alternative weight is provided.
//   - Output dir/format follow PlotEnv defaults (NUXSEC_PLOT_DIR / NUXSEC_PLOT_FORMAT).

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>

#include "AnalysisConfigService.hh"
#include "ColumnDerivationService.hh"
#include "Plotter.hh"
#include "PlottingHelper.hh"
#include "RDataFrameService.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"


using namespace nu;

int stack_samples_impl(const std::string &expr,
                       const std::string &samples_tsv,
                       int nbins,
                       double xmin,
                       double xmax,
                       const std::string &mc_weight,
                       const std::string &extra_sel,
                       bool use_logy)
{
    ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_samples_tsv() : samples_tsv;
    std::cout << "[stack_samples] samples_tsv=" << list_path << "\n";

    const auto sample_list = read_samples(list_path);

    const auto &analysis = AnalysisConfigService::instance();
    const std::string tree_name = analysis.tree_name();

    std::vector<Entry> entries;
    entries.reserve(sample_list.size());

    std::vector<const Entry *> mc;
    std::vector<const Entry *> data;
    mc.reserve(sample_list.size());
    data.reserve(sample_list.size());

    for (const auto &sl : sample_list)
    {
        SampleIO::Sample sample = SampleIO::read(sl.output_path);

        ROOT::RDataFrame rdf = RDataFrameService::load_sample(sample, tree_name);

        const ProcessorEntry proc_entry = analysis.make_processor(sample);
        const auto &deriver = ColumnDerivationService::instance();
        ROOT::RDF::RNode node = deriver.define(rdf, proc_entry);

        entries.emplace_back(make_entry(std::move(node), proc_entry));
        Entry &entry = entries.back();
        if (!extra_sel.empty())
        {
            entry.selection.nominal.node = entry.selection.nominal.node.Filter(extra_sel);
        }
        entry.pot_nom = pick_pot_nom(sample);
        entry.beamline = SampleIO::beam_mode_name(sample.beam); // "numi" or "bnb"
        entry.period = {};

        const Entry *eptr = &entry;
        if (is_data_origin(sample.origin))
        {
            data.push_back(eptr);
        }
        else
        {
            mc.push_back(eptr);
        }
    }

    const std::string weight = mc_weight.empty() ? "w_nominal" : mc_weight;
    const TH1DModel spec = make_spec(expr, nbins, xmin, xmax, weight);

    Plotter plotter;
    auto &opt = plotter.options();
    opt.use_log_y = use_logy;
    opt.legend_on_top = true;
    opt.annotate_numbers = true;
    opt.show_ratio_band = true;
    opt.x_title = expr;
    opt.y_title = "Events";

    plotter.draw_stack(spec, mc, data);
    return 0;
}


int stack_samples(const std::string &expr = "reco_neutrino_vertex_sce_z",
                  const std::string &samples_tsv = "",
                  int nbins = 50,
                  double xmin = 0.0,
                  double xmax = 2.5,
                  const std::string &mc_weight = "w_nominal",
                  const std::string &extra_sel = "true",
                  bool use_logy = false)
{
    return stack_samples_impl(expr,
                                            samples_tsv,
                                            nbins,
                                            xmin,
                                            xmax,
                                            mc_weight,
                                            extra_sel,
                                            use_logy);
}
