/* -- C++ -- */
/**
 *  @file  plot/include/PlottingHelper.hh
 *
 *  @brief Plotting helper utilities and plotStackedHist implementation, providing
 *         reusable routines for assembling stacked histograms.
 */

#ifndef HERON_PLOT_PLOTTING_HELPER_H
#define HERON_PLOT_PLOTTING_HELPER_H

#include <string>

#include <ROOT/RDataFrame.hxx>

#include "PlotDescriptors.hh"
#include "SampleIO.hh"


namespace nu
{

std::string env_or(const char *key, const std::string &fallback);
std::string default_samples_tsv();
std::string default_event_list_root();

bool is_data_origin(SampleIO::SampleOrigin o);

double pick_pot_nom(const SampleIO::Sample &s);

Entry make_entry(ROOT::RDF::RNode node, const ProcessorEntry &proc_entry);

TH1DModel make_spec(const std::string &expr,
                    int nbins,
                    double xmin,
                    double xmax,
                    const std::string &weight);

} // namespace nu


#endif // HERON_PLOT_PLOTTING_HELPER_H
