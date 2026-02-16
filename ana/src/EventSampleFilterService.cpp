/* -- C++ -- */
/**
 *  @file  ana/src/EventSampleFilterService.cpp
 *
 *  @brief Sample-origin filters for event-level RDataFrame processing.
 */

#include "EventSampleFilterService.hh"



const char *EventSampleFilterService::filter_stage(SampleIO::SampleOrigin origin)
{
    using SampleOrigin = SampleIO::SampleOrigin;

    if (origin == SampleOrigin::kOverlay)
    {
        return "filter_overlay";
    }
    if (origin == SampleOrigin::kStrangeness)
    {
        return "filter_strangeness";
    }
    return nullptr;
}

ROOT::RDF::RNode EventSampleFilterService::apply(ROOT::RDF::RNode node, SampleIO::SampleOrigin origin)
{
    using SampleOrigin = SampleIO::SampleOrigin;

    if (origin == SampleOrigin::kOverlay)
    {
        return node.Filter([](int strange) { return strange == 0; }, {"count_strange"});
    }
    if (origin == SampleOrigin::kStrangeness)
    {
        return node.Filter([](int strange) { return strange > 0; }, {"count_strange"});
    }
    return node;
}
