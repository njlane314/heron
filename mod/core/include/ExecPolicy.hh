/* -- C++ -- */
/**
 *  @file  core/include/ExecPolicy.hh
 *
 *  @brief Lightweight execution policy helper for ROOT macros.
 */

#ifndef HERON_CORE_EXEC_POLICY_H
#define HERON_CORE_EXEC_POLICY_H

#include <iostream>
#include <string>

#include <ROOT/RDataFrame.hxx>

#include <TRandom.h>

namespace nu
{
struct ExecPolicy
{
    unsigned int nThreads = 0;
    bool enableImplicitMT = false;
    bool deterministic = false;

    bool apply(const std::string &label = "ExecPolicy") const
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
            std::cout << "\n";
            return true;
        }

        if (ROOT::IsImplicitMTEnabled())
            ROOT::DisableImplicitMT();

        std::cout << "[" << label << "] implicit MT disabled";
        if (deterministic)
            std::cout << ", deterministic seed=1";
        std::cout << "\n";

        return true;
    }
};
}

#endif
