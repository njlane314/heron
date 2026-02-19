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
        SignalLambdaQE = 15,   ///< Signal interaction in Lambda quasi-elastic mode.
        SignalLambdaAssociated = 16, ///< Signal interaction in associated Lambda production mode.
        ECCC = 17,           ///< Electron-neutrino charged-current interaction.
        MuCCOther = 18,      ///< Other muon-neutrino charged-current topologies.
        DataInclusive = 99   ///< Inclusive data channel (non-MC).
    };

    static bool is_lambda_associated_mode(int interaction_mode)
    {
        // Based on simb::int_type_ / Nuance-offset values carried in MCNeutrino::Mode().
        return interaction_mode == 1073 || // kResCCNuKaonPlusLambda0
               interaction_mode == 1076;   // kResCCNuBarKaon0Lambda0
    }

    static AnalysisChannel classify_analysis_channel(
        bool in_fiducial,
        int nu_pdg,
        int ccnc,
        int interaction_mode,
        int n_p,
        int n_pi_minus,
        int n_pi_plus,
        int n_pi0,
        int n_gamma,
        bool is_nu_mu_cc,
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

        if (is_signal(is_nu_mu_cc, ccnc, in_fiducial, mu_p, p_p, pi_p, lam_decay_sep))
        {
            if (is_lambda_associated_mode(interaction_mode))
                return AnalysisChannel::SignalLambdaAssociated;
            return AnalysisChannel::SignalLambdaQE;
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

    static bool is_signal(
        bool is_nu_mu_cc,
        int ccnc,
        bool in_fiducial,
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
        if (!std::isfinite(mu_p) || !std::isfinite(p_p) || !std::isfinite(pi_p) ||
            !std::isfinite(lam_decay_sep))
            return false;
        if (mu_p < min_mu_p || p_p < min_p_p || pi_p < min_pi_p)
            return false;
        return lam_decay_sep >= min_lam_decay_sep;
    }
};

#endif // HERON_ANA_ANALYSIS_CHANNELS_H
