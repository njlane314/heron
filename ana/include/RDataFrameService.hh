/* -- C++ -- */
/**
 *  @file  ana/include/RDataFrameService.hh
 *
 *  @brief Sample loading and variable definitions for ROOT RDataFrame.
 */

#ifndef NUXSEC_ANA_RDATA_FRAME_SERVICE_H
#define NUXSEC_ANA_RDATA_FRAME_SERVICE_H

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

class RDataFrameService
{
  public:
    static ROOT::RDataFrame load_sample(const sample::SampleIO::Sample &sample,
                                        const std::string &tree_name);

    static ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                             const std::vector<ColumnDefinition> &definitions);
};

} // namespace nuxsec

#endif // NUXSEC_ANA_RDATA_FRAME_SERVICE_H
