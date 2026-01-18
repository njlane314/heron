/* -- C++ -- */
/**
 *  @file  rdf/src/RDFBuilder.cc
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame.
 */

#include "RDFBuilder.hh"

#include <utility>

namespace nuxsec
{

ROOT::RDataFrame RDFBuilder::load_sample(const Sample &sample, const std::string &tree_name)
{
    std::vector<std::string> files = collect_files(sample);
    return ROOT::RDataFrame(tree_name, files);
}

ROOT::RDF::RNode RDFBuilder::define_variables(ROOT::RDF::RNode node,
                                             const std::vector<ColumnDefinition> &definitions)
{
    ROOT::RDF::RNode updated_node = std::move(node);
    for (const ColumnDefinition &definition : definitions)
    {
        updated_node = updated_node.Define(definition.name, definition.expression);
    }

    return updated_node;
}

std::vector<std::string> RDFBuilder::collect_files(const Sample &sample)
{
    std::vector<std::string> files;
    files.reserve(sample.stages.size());
    for (const SampleStage &stage : sample.stages)
    {
        files.push_back(stage.artio_path);
    }

    return files;
}

} // namespace nuxsec
