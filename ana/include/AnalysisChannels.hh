/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisChannels.hh
 *
 *  @brief Analysis channel definitions used during event categorisation.
 */

#ifndef HERON_ANA_ANALYSIS_CHANNELS_H
#define HERON_ANA_ANALYSIS_CHANNELS_H

#include <cstdlib>

/** \brief Channel definitions for reconstruction-level analysis categorisation. */
class AnalysisChannels
{
  public:
    enum class AnalysisChannel
    {
        Unknown = 0,      ///< Unclassified channel (fallback for unknown categorisation).
        External = 1,     ///< External (non-neutrino) or out-of-volume background.
        OutFV = 2,        ///< Interaction outside the truth fiducial volume.
        MuCC0pi_ge1p = 10,   ///< Muon-neutrino charged-current with 0 pions and ≥1 proton.
        MuCC1pi = 11,        ///< Muon-neutrino charged-current with exactly one charged pion.
        MuCCPi0OrGamma = 12, ///< Muon-neutrino charged-current with π0 or photon activity.
        MuCCNpi = 13,        ///< Muon-neutrino charged-current with more than one pion.
        NC = 14,             ///< Neutral-current interaction in fiducial volume.
        CCS1 = 15,           ///< Charged-current interaction with exactly one strange hadron.
        CCSgt1 = 16,         ///< Charged-current interaction with multiple strange hadrons.
        ECCC = 17,           ///< Electron-neutrino charged-current interaction.
        MuCCOther = 18,      ///< Other muon-neutrino charged-current topologies.
        DataInclusive = 99   ///< Inclusive data channel (non-MC).
    };

    static AnalysisChannel classify_analysis_channel(
        bool in_fiducial,
        int nu_pdg,
        int ccnc,
        int count_strange,
        int n_p,
        int n_pi_minus,
        int n_pi_plus,
        int n_pi0,
        int n_gamma)
    {
        const int npi = n_pi_minus + n_pi_plus;

        if (!in_fiducial)
        {
            if (nu_pdg == 0)
                return AnalysisChannel::External;
            return AnalysisChannel::OutFV;
        }

        if (ccnc == 1)
            return AnalysisChannel::NC;

        if (ccnc == 0 && count_strange > 0)
        {
            if (count_strange == 1)
                return AnalysisChannel::CCS1;
            return AnalysisChannel::CCSgt1;
        }

        if (std::abs(nu_pdg) == 12 && ccnc == 0)
            return AnalysisChannel::ECCC;

        if (std::abs(nu_pdg) == 14 && ccnc == 0)
        {
            if (npi == 0 && n_p > 0)
                return AnalysisChannel::MuCC0pi_ge1p;
            if (npi == 1 && n_pi0 == 0)
                return AnalysisChannel::MuCC1pi;
            if (n_pi0 > 0 || n_gamma >= 2)
                return AnalysisChannel::MuCCPi0OrGamma;
            if (npi > 1)
                return AnalysisChannel::MuCCNpi;
            return AnalysisChannel::MuCCOther;
        }

        return AnalysisChannel::Unknown;
    }

    static int to_int(AnalysisChannel channel)
    {
        return static_cast<int>(channel);
    }
};

#endif // HERON_ANA_ANALYSIS_CHANNELS_H
