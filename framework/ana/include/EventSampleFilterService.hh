/* -- C++ -- */
/**
 *  @file  ana/include/EventSampleFilterService.hh
 *
 *  @brief Sample-origin filters for event-level RDataFrame processing.
 */

#ifndef HERON_ANA_EVENT_SAMPLE_FILTER_SERVICE_H
#define HERON_ANA_EVENT_SAMPLE_FILTER_SERVICE_H

#include <ROOT/RDataFrame.hxx>

#include "SampleIO.hh"


class EventSampleFilterService
{
  public:
    virtual ~EventSampleFilterService() = default;

    virtual const char *filter_stage(SampleIO::SampleOrigin origin) const = 0;
    virtual ROOT::RDF::RNode apply(ROOT::RDF::RNode node, SampleIO::SampleOrigin origin) const = 0;
};


#endif // HERON_ANA_EVENT_SAMPLE_FILTER_SERVICE_H
