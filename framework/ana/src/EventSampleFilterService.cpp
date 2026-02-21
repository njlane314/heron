/* -- C++ -- */
/**
 *  @file  ana/src/EventSampleFilterService.cpp
 *
 *  @brief Event sample filter service accessors.
 */

#include "EventSampleFilterService.hh"

#include "DefaultAnalysisModel.hh"

const EventSampleFilterService &EventSampleFilterService::instance()
{
    return DefaultAnalysisModel::instance();
}
