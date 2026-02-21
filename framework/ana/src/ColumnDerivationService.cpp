/* -- C++ -- */
/**
 *  @file  ana/src/ColumnDerivationService.cpp
 *
 *  @brief Column derivation service accessors.
 */

#include "ColumnDerivationService.hh"

#include "DefaultAnalysisModel.hh"

const ColumnDerivationService &ColumnDerivationService::instance()
{
    return DefaultAnalysisModel::instance();
}
