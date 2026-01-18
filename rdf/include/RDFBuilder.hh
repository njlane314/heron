/* -- C++ -- */
/**
 *  @file  rdf/include/RDFBuilder.hh
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame.
 */

#ifndef Nuxsec_RDF_RDF_BUILDER_H_INCLUDED
#define Nuxsec_RDF_RDF_BUILDER_H_INCLUDED

#include <ROOT/RDataFrame.hxx>

#include <string>
#include <vector>

#include "SampleIO.hh"

namespace nuxsec
{

struct ColumnDefinition
{
    std::string name;
    std::string expression;
    std::string description;
};

class RDFBuilder
{
  public:
    static ROOT::RDataFrame load_sample(const Sample &sample, const std::string &tree_name);
    static ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                             const std::vector<ColumnDefinition> &definitions);
    static std::vector<std::string> collect_files(const Sample &sample);
};

} // namespace nuxsec

#endif // Nuxsec_RDF_RDF_BUILDER_H_INCLUDED
