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
    static ROOT::RDataFrame load_sample(const SampleIO::Sample &sample,
                                        const std::string &tree_name);

    static ROOT::RDF::RNode define_variables(ROOT::RDF::RNode node,
                                             const std::vector<Column> &definitions);
};


#endif // HERON_ANA_RDATA_FRAME_SERVICE_H
