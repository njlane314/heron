/* -- C++ -- */
/**
 *  @file  ana/src/RDataFrameService.cpp
 *
 *  @brief RDataFrame service accessors.
 */

#include "RDataFrameService.hh"

#include "DefaultAnalysisModel.hh"

const RDataFrameService &RDataFrameService::instance()
{
    return DefaultAnalysisModel::instance();
}
