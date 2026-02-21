/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisChannels.hh
 *
 *  @brief Analysis channel definitions used during event categorisation.
 */

#ifndef HERON_ANA_ANALYSIS_CHANNELS_H
#define HERON_ANA_ANALYSIS_CHANNELS_H

#include <cmath>
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
        SignalLambdaCCQE = 15,    ///< Signal interaction in CCQE Lambda mode.
        SignalLambdaCCRES = 16,   ///< Signal interaction in CCRES Lambda mode.
        SignalLambdaCCDIS = 17,   ///< Signal interaction in CCDIS Lambda mode.
        SignalLambdaCCOther = 18, ///< Signal interaction in other CC Lambda modes.
        ECCC = 19,                ///< Electron-neutrino charged-current interaction.
        MuCCOther = 20,           ///< Other muon-neutrino charged-current topologies.
        DataInclusive = 99        ///< Inclusive data channel (non-MC).
    };

    static AnalysisChannel classify_lambda_signal_channel(int interaction_type)
    {
        // Support both low-number simb::int_type_ values and Nuance-offset values.
        // Some ntuples persist 0/1/2 (QE/Res/DIS), others store 1001/10xx variants.
        if (interaction_type == 0 || interaction_type == 1001)
            return AnalysisChannel::SignalLambdaCCQE;
        if (interaction_type == 1 || interaction_type == 1073 || interaction_type == 1076)
            return AnalysisChannel::SignalLambdaCCRES;
        if (interaction_type == 2 || interaction_type == 1091)
            return AnalysisChannel::SignalLambdaCCDIS;
        return AnalysisChannel::SignalLambdaCCOther;
    }

    static AnalysisChannel classify_analysis_channel(
        bool in_fiducial,
        int nu_pdg,
        int ccnc,
        int interaction_type,
        int n_p,
        int n_pi_minus,
        int n_pi_plus,
        int n_pi0,
        int n_gamma,
        bool is_nu_mu_cc,
        int lam_pdg,
        float mu_p,
        float p_p,
        float pi_p,
        float lam_decay_sep)
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

        if (is_signal(is_nu_mu_cc, ccnc, in_fiducial, lam_pdg, mu_p, p_p, pi_p, lam_decay_sep))
            return classify_lambda_signal_channel(interaction_type);

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

    static bool is_signal(
        bool is_nu_mu_cc,
        int ccnc,
        bool in_fiducial,
        int lam_pdg,
        float mu_p,
        float p_p,
        float pi_p,
        float lam_decay_sep)
    {
        const float min_mu_p = 0.10f;
        const float min_p_p = 0.30f;
        const float min_pi_p = 0.10f;
        const float min_lam_decay_sep = 0.50f;

        if (!is_nu_mu_cc)
            return false;
        if (ccnc != 0)
            return false;
        if (!in_fiducial)
            return false;
        if (lam_pdg != 3122)
            return false;
        if (!std::isfinite(mu_p) || !std::isfinite(p_p) || !std::isfinite(pi_p) ||
            !std::isfinite(lam_decay_sep))
            return false;
        if (mu_p < min_mu_p || p_p < min_p_p || pi_p < min_pi_p)
            return false;
        return lam_decay_sep >= min_lam_decay_sep;
    }
};

#endif // HERON_ANA_ANALYSIS_CHANNELS_H
