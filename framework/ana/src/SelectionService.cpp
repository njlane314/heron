/* -- C++ -- */
/**
 *  @file  ana/src/SelectionService.cpp
 *
 *  @brief Selection service accessors.
 */

#include "SelectionService.hh"

#include "DefaultAnalysisModel.hh"

const SelectionService &SelectionService::instance()
{
    return DefaultAnalysisModel::instance();
}
