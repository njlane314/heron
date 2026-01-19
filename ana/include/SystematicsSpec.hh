/* -- C++ -- */
/**
 *  @file  ana/include/SystematicsSpec.hh
 *
 *  @brief Systematics specification helpers.
 */

#ifndef Nuxsec_ANA_SYSTEMATICS_SPEC_H_INCLUDED
#define Nuxsec_ANA_SYSTEMATICS_SPEC_H_INCLUDED

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

#endif // Nuxsec_ANA_SYSTEMATICS_SPEC_H_INCLUDED
