/* -- C++ -- */
/**
 *  @file  ana/include/RDataFrameService.hh
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame,
 *         covering input configuration and dataframe initialisation.
 */

#ifndef HERON_ANA_RDATA_FRAME_SERVICE_H
#define HERON_ANA_RDATA_FRAME_SERVICE_H

#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "SampleIO.hh"


struct Column
{
    std::string name;
    std::string expression;
    std::string description;
};

class RDataFrameService
{
  public:
    virtual ~RDataFrameService() = default;

    virtual ROOT::RDataFrame load_sample(const SampleIO::Sample &sample,
                                         const std::string &tree_name) const
    {
        std::vector<std::string> files = SampleIO::resolve_root_files(sample);
        return ROOT::RDataFrame(tree_name, files);
    }

    virtual ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                               const std::vector<Column> &definitions) const
    {
        ROOT::RDF::RNode updated_node = std::move(node);
        for (const Column &definition : definitions)
        {
            updated_node = updated_node.Define(definition.name, definition.expression);
        }
        return updated_node;
    }
};


#endif // HERON_ANA_RDATA_FRAME_SERVICE_H
