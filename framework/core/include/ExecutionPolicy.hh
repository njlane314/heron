/* -- C++ -- */
/**
 *  @file  core/include/ExecutionPolicy.hh
 *
 *  @brief Lightweight execution policy helper for ROOT macros.
 */

#ifndef HERON_CORE_EXECUTION_POLICY_H
#define HERON_CORE_EXECUTION_POLICY_H

#include <cstdlib>
#include <iostream>
#include <string>

#include <ROOT/RDataFrame.hxx>

#include <TRandom.h>

namespace nu
{
struct ExecutionPolicy
{
    unsigned int nThreads = 0;
    bool enableImplicitMT = false;
    bool deterministic = false;
    bool deterministicMerging = false;

    bool implicit_mt_enabled() const noexcept { return enableImplicitMT; }

    static bool env_enabled(const char *name, bool defaultValue = false)
    {
        const char *v = std::getenv(name);
        if (v == NULL || *v == '\0')
            return defaultValue;

        return std::string(v) != "0";
    }

    static ExecutionPolicy from_env(const char *name = "HERON_PLOT_IMT")
    {
        ExecutionPolicy policy;
        policy.enableImplicitMT = env_enabled(name, false);
        return policy;
    }

    bool apply(const std::string &label = "ExecutionPolicy") const
    {
        if (deterministic && gRandom != NULL)
            gRandom->SetSeed(1);

        if (enableImplicitMT)
        {
            if (ROOT::IsImplicitMTEnabled())
                ROOT::DisableImplicitMT();

            if (nThreads > 0)
                ROOT::EnableImplicitMT(nThreads);
            else
                ROOT::EnableImplicitMT();

            std::cout << "[" << label << "] implicit MT enabled";
            if (nThreads > 0)
                std::cout << " (nThreads=" << nThreads << ")";
            if (deterministic)
                std::cout << ", deterministic seed=1";
            if (deterministicMerging)
                std::cout << ", deterministic merge ordering requested";
            std::cout << "\n";
            return true;
        }

        if (ROOT::IsImplicitMTEnabled())
            ROOT::DisableImplicitMT();

        std::cout << "[" << label << "] implicit MT disabled";
        if (deterministic)
            std::cout << ", deterministic seed=1";
        if (deterministicMerging)
            std::cout << ", deterministic merge ordering requested";
        std::cout << "\n";

        return true;
    }
};

}

#endif
