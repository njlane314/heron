/**
 *  @file  lib/NuAna/include/NuAna/SampleRDF.h
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame
 */

#ifndef NUANA_SAMPLE_RDF_H
#define NUANA_SAMPLE_RDF_H

#include <ROOT/RDataFrame.hxx>

#include <string>
#include <vector>

#include "NuIO/SampleIO.h"

namespace nuana
{

struct VariableDefinition
{
    std::string name;
    std::string expression;
    std::string description;
};

class SampleRDF
{
  public:
    static ROOT::RDataFrame load_sample(const nuio::Sample &sample, const std::string &tree_name);
    static ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                             const std::vector<VariableDefinition> &definitions);
    static std::vector<std::string> collect_files(const nuio::Sample &sample);
};

} // namespace nuana

#endif // NUANA_SAMPLE_RDF_H
