/**
 *  @file  lib/NuAna/src/SampleRDF.cxx
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame
 */

#include "NuAna/SampleRDF.h"

#include <utility>

namespace nuana
{

ROOT::RDataFrame SampleRDF::load_sample(const nuio::Sample &sample, const std::string &tree_name)
{
    std::vector<std::string> files = collect_files(sample);
    return ROOT::RDataFrame(tree_name, files);
}

ROOT::RDF::RNode SampleRDF::define_variables(ROOT::RDF::RNode node,
                                             const std::vector<VariableDefinition> &definitions)
{
    ROOT::RDF::RNode updated_node = std::move(node);
    for (const VariableDefinition &definition : definitions)
    {
        updated_node = updated_node.Define(definition.name, definition.expression);
    }

    return updated_node;
}

std::vector<std::string> SampleRDF::collect_files(const nuio::Sample &sample)
{
    std::vector<std::string> files;
    files.reserve(sample.stages.size());
    for (const nuio::SampleStage &stage : sample.stages)
    {
        files.push_back(stage.artio_path);
    }

    return files;
}

} // namespace nuana
