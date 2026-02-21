/* -- C++ -- */
/**
 *  @file  macro/AnalysisModelExample.C
 *
 *  @brief Minimal example macro that loads samples and initialises the
 *         default analysis model with an analysis context.
 */

#include <iostream>
#include <stdexcept>
#include <string>

#include "AnalysisContext.hh"
#include "Dataset.hh"
#include "DefaultAnalysisModel.hh"
#include "ExecutionPolicy.hh"

void AnalysisModelExample(const char *sample_list_path = "scratch/out/out/sample/samples.tsv")
{
    try
    {
        nu::ExecutionPolicy policy{.nThreads = 4, .enableImplicitMT = true, .deterministic = true};
        policy.apply("AnalysisModelExample");

        heron::Dataset dataset;
        if (!dataset.load_samples(sample_list_path))
        {
            throw std::runtime_error(std::string("AnalysisModelExample: failed to load sample list: ") + sample_list_path);
        }

        heron::AnalysisContext<nu::ExecutionPolicy> context(policy, "default");

        DefaultAnalysisModel &model = const_cast<DefaultAnalysisModel &>(DefaultAnalysisModel::instance());
        model.initialise(context);

        std::cout << "[AnalysisModelExample] model=" << model.name()
                  << " tree=" << model.tree_name()
                  << " selections=" << model.selections().size()
                  << " loaded_samples=" << dataset.samples().size()
                  << " from=" << sample_list_path
                  << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "[AnalysisModelExample] " << e.what() << "\n";
    }
}
