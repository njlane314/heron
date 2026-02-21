/* -- C++ -- */
/**
 *  @file  ana/src/AnalysisConfigService.cpp
 *
 *  @brief Analysis configuration service accessors.
 */

#include "AnalysisConfigService.hh"

#include "DefaultAnalysisModel.hh"

const AnalysisConfigService &AnalysisConfigService::instance()
{
    return DefaultAnalysisModel::instance();
}
