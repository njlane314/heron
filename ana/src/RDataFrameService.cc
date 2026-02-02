/* -- C++ -- */
/**
 *  @file  ana/src/RDataFrameService.cc
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame.
 */

#include "RDataFrameService.hh"

#include <utility>

namespace nuxsec
{

ROOT::RDataFrame RDataFrameService::load_sample(const sample::SampleIO::Sample &sample,
                                                const std::string &tree_name)
{
    std::vector<std::string> files = sample::SampleIO::resolve_root_files(sample);
    return ROOT::RDataFrame(tree_name, files);
}

ROOT::RDF::RNode RDataFrameService::define_variables(ROOT::RDF::RNode node,
                                             const std::vector<Column> &definitions)
{
    ROOT::RDF::RNode updated_node = std::move(node);
    for (const Column &definition : definitions)
    {
        updated_node = updated_node.Define(definition.name, definition.expression);
    }

    return updated_node;
}

} // namespace nuxsec
