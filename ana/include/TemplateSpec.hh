/* -- C++ -- */
/**
 *  @file  ana/include/TemplateSpec.hh
 *
 *  @brief Template specification helpers for analysis histogram booking.
 */

#ifndef NUXSEC_ANA_TEMPLATE_SPEC_H
#define NUXSEC_ANA_TEMPLATE_SPEC_H

#include <string>
#include <vector>

namespace nuxsec
{

struct TemplateSpec1D
{
    std::string name;
    std::string title;
    std::string selection;
    std::string variable;
    std::string weight;
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 0.0;
};

/// Parse a tab-separated template specification file into 1D template entries.
std::vector<TemplateSpec1D> read_template_spec_1d_tsv(const std::string &path);

} // namespace nuxsec

#endif // NUXSEC_ANA_TEMPLATE_SPEC_H
