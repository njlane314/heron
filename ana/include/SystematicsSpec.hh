/* -- C++ -- */
/**
 *  @file  ana/include/SystematicsSpec.hh
 *
 *  @brief Systematics specification helpers.
 */

#ifndef NUXSEC_ANA_SYSTEMATICS_SPEC_H
#define NUXSEC_ANA_SYSTEMATICS_SPEC_H

#include <string>
#include <vector>

namespace nuxsec
{

struct UnisimSpec
{
    std::string name;
    std::string up_branch;
    std::string dn_branch;
    bool one_sided = false;
    bool log_normal = false;
    bool floatable = false;
};

struct MultisimSpec
{
    std::string name;
    std::string vec_branch;
    std::string cv_branch;
    int max_universes = -1;
    int max_modes = 20;
    double keep_fraction = 0.99;
    bool split_rate_shape = true;
    bool rate_log_normal = true;
};

struct SystematicsConfig
{
    std::vector<UnisimSpec> unisim;
    std::vector<MultisimSpec> multisim;
};

const SystematicsConfig &DefaultSystematics();

} // namespace nuxsec

#endif // NUXSEC_ANA_SYSTEMATICS_SPEC_H
