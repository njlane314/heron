/* -- C++ -- */
/**
 *  @file  plot/include/ParticleChannels.hh
 *
 *  @brief Particle-type display properties for particle-level (per-reco-object)
 *         stacked histograms. Categories are based on truth-matched PDG codes.
 */

#ifndef HERON_PLOT_PARTICLE_CHANNELS_H
#define HERON_PLOT_PARTICLE_CHANNELS_H

#include <map>
#include <string>
#include <vector>

#include "RtypesCore.h"
#include "TColor.h"

namespace nu
{

class ParticleChannels
{
  public:
    struct Properties
    {
        int key = 0; // category key (we use the canonical PDG code; 0 means "Other / unmatched")
        std::string plain_name;
        std::string tex_label;
        Color_t fill_colour = kBlack;
        int fill_style = 1001;
    };

    // Map an input PDG code (can be signed) to a plotting category key.
    // Unknown / unmatched values are mapped to 0 ("Other / unmatched").
    static int classify(int pdg)
    {
        const int ap = (pdg < 0) ? -pdg : pdg;
        switch (ap)
        {
        case 13:
            return 13; // mu
        case 11:
            return 11; // e
        case 22:
            return 22; // gamma
        case 2212:
            return 2212; // p
        case 2112:
            return 2112; // n
        case 211:
            return 211; // pi±
        case 111:
            return 111; // pi0
        case 321:
            return 321; // K±
        default:
            return 0;
        }
    }

    static bool matches(int key, int pdg) { return classify(pdg) == key; }

    static const Properties &properties(int key)
    {
        const auto &m = mapping();
        auto it = m.find(key);
        if (it == m.end())
        {
            return m.at(0);
        }
        return it->second;
    }

    static std::string label(int key) { return properties(key).tex_label; }
    static std::string name(int key) { return properties(key).plain_name; }
    static int colour(int key) { return properties(key).fill_colour; }
    static int fill_style(int key) { return properties(key).fill_style; }

    static const std::vector<int> &keys()
    {
        // NOTE: ordering here is only used for "booking"; the final plotted order
        // is yield-sorted in StackedHist.
        static const std::vector<int> v = {13, 2212, 211, 111, 11, 22, 321, 2112, 0};
        return v;
    }

  private:
    static const std::map<int, Properties> &mapping()
    {
        static const std::map<int, Properties> m = {
            {13, {13, "muon", "#mu^{#pm}", kAzure + 2, 1001}},
            {2212, {2212, "proton", "p", kRed - 4, 1001}},
            {211, {211, "pion", "#pi^{#pm}", kOrange + 1, 1001}},
            {111, {111, "pi0", "#pi^{0}", kOrange - 3, 1001}},
            {11, {11, "electron", "e^{#pm}", kMagenta + 1, 1001}},
            {22, {22, "gamma", "#gamma", kYellow - 7, 1001}},
            {321, {321, "kaon", "K^{#pm}", kSpring + 4, 1001}},
            {2112, {2112, "neutron", "n", kCyan + 1, 1001}},
            {0, {0, "other", "Other / unmatched", kGray + 1, 1001}},
        };
        return m;
    }
};

} // namespace nu

#endif // HERON_PLOT_PARTICLE_CHANNELS_H
