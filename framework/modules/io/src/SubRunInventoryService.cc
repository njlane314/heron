/* -- C++ -- */
/**
 *  @file  framework/io/src/SubRunInventoryService.cc
 *
 *  @brief Implementation of the SubRun inventory helpers.
 */

#include "SubRunInventoryService.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "ArtFileProvenanceIO.hh"


Summary SubRunInventoryService::scan_subruns(const std::vector<std::string> &files)
{
    Summary out;

    const std::vector<std::string> candidates = {"nuselection/SubRun", "SubRun"};
    std::string tree_path;

    for (const auto &f : files)
    {
        std::unique_ptr<TFile> file(TFile::Open(f.c_str(), "READ"));
        if (!file || file->IsZombie())
        {
            throw std::runtime_error("Failed to open input ROOT file: " + f);
        }

        for (const auto &name : candidates)
        {
            if (dynamic_cast<TTree *>(file->Get(name.c_str())))
            {
                tree_path = name;
                break;
            }
        }
        if (!tree_path.empty())
        {
            break;
        }
    }

    if (tree_path.empty())
    {
        throw std::runtime_error("No input files contained a SubRun tree.");
    }

    TChain chain(tree_path.c_str());
    for (const auto &f : files)
    {
        chain.Add(f.c_str());
    }

    if (!chain.GetBranch("run") || !chain.GetBranch("subRun") || !chain.GetBranch("pot"))
    {
        throw std::runtime_error("SubRun tree missing required branches (run, subRun, pot).");
    }

    Int_t run = 0;
    Int_t subRun = 0;
    Double_t pot = 0.0;

    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("subRun", &subRun);
    chain.SetBranchAddress("pot", &pot);

    const Long64_t n = chain.GetEntries();
    out.n_entries = static_cast<long long>(n);

    std::vector<Subrun> pairs;
    pairs.reserve(static_cast<size_t>(n));

    for (Long64_t i = 0; i < n; ++i)
    {
        chain.GetEntry(i);
        out.pot_sum += static_cast<double>(pot);
        pairs.push_back(Subrun{static_cast<int>(run), static_cast<int>(subRun)});
    }

    std::sort(pairs.begin(), pairs.end(),
              [](const Subrun &a, const Subrun &b)
              {
                  if (a.run != b.run)
                  {
                      return a.run < b.run;
                  }
                  return a.subrun < b.subrun;
              });

    pairs.erase(std::unique(pairs.begin(), pairs.end(),
                            [](const Subrun &a, const Subrun &b)
                            {
                                return a.run == b.run && a.subrun == b.subrun;
                            }),
                pairs.end());

    out.unique_pairs = std::move(pairs);

    return out;
}
