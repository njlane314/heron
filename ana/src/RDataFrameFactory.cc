/* -- C++ -- */
/**
 *  @file  ana/src/RDataFrameFactory.cc
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame.
 */

#include "RDataFrameFactory.hh"

#include <utility>

namespace nuxsec
{

ROOT::RDataFrame RDataFrameFactory::load_sample(const SampleIO::Sample &sample, const std::string &tree_name)
{
    std::vector<std::string> files = collect_files(sample);
    return ROOT::RDataFrame(tree_name, files);
}

ROOT::RDF::RNode RDataFrameFactory::define_variables(ROOT::RDF::RNode node,
                                             const std::vector<ColumnDefinition> &definitions)
{
    ROOT::RDF::RNode updated_node = std::move(node);
    for (const ColumnDefinition &definition : definitions)
    {
        updated_node = updated_node.Define(definition.name, definition.expression);
    }

    return updated_node;
}

std::vector<std::string> RDataFrameFactory::collect_files(const SampleIO::Sample &sample)
{
    std::vector<std::string> files;
    files.reserve(sample.inputs.size());
    for (const SampleIO::ProvenanceInput &input : sample.inputs)
    {
        files.push_back(input.art_path);
    }

    return files;
}

} // namespace nuxsec
