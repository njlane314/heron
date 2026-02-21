/* -- C++ -- */
/**
 *  @file  macro/AnalysisModelExample.C
 *
 *  @brief Example macro showing a user-defined AnalysisModel wired through
 *         ExecutionPolicy and AnalysisContext.
 */

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AnalysisContext.hh"
#include "AnalysisModel.hh"
#include "Dataset.hh"
#include "ExecutionPolicy.hh"

namespace
{

using AnalysisModel = heron::AnalysisModel;
using ExecutionPolicy = nu::ExecutionPolicy;
using Dataset = heron::Dataset;

template <typename TPolicy, typename TServices = std::nullptr_t>
using AnalysisContext = heron::AnalysisContext<TPolicy, TServices>;
class ExampleAnalysisModel : public AnalysisModel
{
  public:
    void configure() override
    {
        m_selection_name = "sel_triggered_muon";
    }

    template <typename TPolicy, typename TServices>
    void configure(const AnalysisContext<TPolicy, TServices> &context)
    {
        const std::string &label = context.metadata().analysis_label;
        if (!label.empty())
        {
            m_selection_name = "sel_" + label;
        }
    }

    void define_channels() override
    {
        channel("inclusive", []() { return 99; }, {"is_data"});
    }

    void define_selections() override
    {
        m_weight = weight("w_nominal", []() { return 1.0; }, {"w_nominal"});
        const auto c_trigger = cut("trigger", []() { return true; }, {"sel_trigger"});
        m_nominal = selection(m_selection_name, c_trigger, m_weight);
    }

    void define_outputs() override
    {
        const auto p_muon = var("muon_p", []() { return 0.0; }, {"mu_p"});

        hist1d("h_muon_p",
               p_muon.name,
               40,
               0.0,
               2.0,
               "Muon momentum;p [GeV];Events",
               m_nominal.name,
               m_weight.name);

        snapshot("events", {"run", "sub", "evt", p_muon.name, m_weight.name}, m_nominal.name);
    }

  private:
    std::string m_selection_name;
    Weight m_weight;
    Selection m_nominal;
};

} // namespace

void AnalysisModelExample(const char *sample_list_path = "scratch/out/out/sample/samples.tsv")
{
    ExecutionPolicy policy = []() { ExecutionPolicy p = ExecutionPolicy::from_env("HERON_PLOT_IMT"); p.nThreads = 4; p.deterministic = true; p.apply("AnalysisModelExample"); return p; }();

    Dataset dataset;
    if (!dataset.load_samples(sample_list_path))
    {
        throw std::runtime_error(std::string("AnalysisModelExample: failed to load sample list: ") + sample_list_path);
    }

    std::cout << "[AnalysisModelExample] loaded_samples=" << dataset.samples().size() << " from=" << sample_list_path << "\n";

    AnalysisContext<ExecutionPolicy> context(policy, "triggered_muon");

    ExampleAnalysisModel model;
    model.initialise(context);

    std::cout << "[AnalysisModelExample] vars=" << model.vars().size()
              << " cuts=" << model.cuts().size()
              << " weights=" << model.weights().size()
              << " selections=" << model.selections().size()
              << " hist1d=" << model.h1().size()
              << " snapshots=" << model.snapshots().size()
              << "\n";
}
