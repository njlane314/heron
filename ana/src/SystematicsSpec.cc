/* -- C++ -- */
/**
 *  @file  ana/src/SystematicsSpec.cc
 *
 *  @brief Default systematics specification list.
 */

#include "SystematicsSpec.hh"

namespace nuxsec
{

const SystematicsConfig &DefaultSystematics()
{
    static const SystematicsConfig cfg{
        {
            {"RPA_CCQE", "knobRPAup", "knobRPAdn", false, false, false},
            {"XSecShape_CCMEC", "knobCCMECup", "knobCCMECdn", false, false, false},
            {"DecayAngMEC", "knobDecayAngMECup", "knobDecayAngMECdn", false, false, false},
            {"Theta_Delta2Npi", "knobThetaDelta2Npiup", "knobThetaDelta2Npidn", false, false, false},
            {"NormCCCOH", "knobNormCCCOHup", "knobNormCCCOHdn", true, true, false},
            {"NormNCCOH", "knobNormNCCOHup", "knobNormNCCOHdn", true, true, false},
        },
        {
            {"ppfx", "weightsPPFX", "ppfx_cv", 600, 30, 0.99, true, true},
            {"genie_all", "weightsGenie", "", -1, 30, 0.99, true, false},
            {"reint", "weightsReint", "", -1, 20, 0.99, true, false},
        }};
    return cfg;
}

} // namespace nuxsec
