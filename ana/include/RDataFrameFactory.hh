/* -- C++ -- */
/**
 *  @file  ana/include/RDataFrameFactory.hh
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame.
 */

#ifndef NUXSEC_ANA_RDATA_FRAME_FACTORY_H
#define NUXSEC_ANA_RDATA_FRAME_FACTORY_H

#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

#include "SampleIO.hh"

namespace nuxsec
{

struct ColumnDefinition
{
    std::string name;
    std::string expression;
    std::string description;
};

class RDataFrameFactory
{
  public:
    static ROOT::RDataFrame load_sample(const SampleIO::Sample &sample, const std::string &tree_name);
    static ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                             const std::vector<ColumnDefinition> &definitions);
    static std::vector<std::string> collect_files(const SampleIO::Sample &sample);
};

} // namespace nuxsec

#endif // NUXSEC_ANA_RDATA_FRAME_FACTORY_H
