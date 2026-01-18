/* -- C++ -- */
/**
 *  @file  plot/include/PlotChannels.hh
 *
 *  @brief Channel display properties for plotting.
 */

#ifndef Nuxsec_PLOT_CHANNELS_H_INCLUDED
#define Nuxsec_PLOT_CHANNELS_H_INCLUDED

#include "RtypesCore.h"
#include "TColor.h"

#include <map>
#include <string>
#include <vector>

namespace nuxsec
{
namespace plot
{

class Channels
{
  public:
    struct Properties
    {
        int key = 0;
        std::string plain_name;
        std::string tex_label;
        Color_t fill_color = kBlack;
        int fill_style = 1001;
    };

    static const Properties &properties(int code)
    {
        const auto &m = mapping();
        auto it = m.find(code);
        if (it == m.end())
        {
            return m.at(99);
        }
        return it->second;
    }

    static std::string label(int code) { return properties(code).tex_label; }
    static std::string name(int code) { return properties(code).plain_name; }
    static int color(int code) { return properties(code).fill_color; }
    static int fill_style(int code) { return properties(code).fill_style; }

    static const std::vector<int> &signal_keys()
    {
        static const std::vector<int> v = {15, 16};
        return v;
    }

    static const std::vector<int> &mc_keys()
    {
        static const std::vector<int> v = {1, 2, 10, 11, 12, 13, 14, 15, 16, 17, 18, 99};
        return v;
    }

  private:
    static const std::map<int, Properties> &mapping()
    {
        static const std::map<int, Properties> m = {
            {0, {0, "data", "Data", kBlack, 1001}},
            {1, {1, "external", "Cosmic", kTeal + 2, 3345}},
            {2, {2, "out_fv", "Out FV", kYellow - 7, 1001}},
            {10, {10, "numu_cc_np0pi", "#nu_{#mu}CC Np0#pi", kRed, 1001}},
            {11, {11, "numu_cc_0pnpi", "#nu_{#mu}CC 0p1#pi^{#pm}", kRed - 7, 1001}},
            {12, {12, "numu_cc_pi0gg", "#nu_{#mu}CC #pi^{0}/#gamma#gamma", kOrange, 1001}},
            {13, {13, "numu_cc_npnpi", "#nu_{#mu}CC multi-#pi^{#pm}", kViolet, 1001}},
            {14, {14, "nc", "#nu_{x}NC", kBlue, 1001}},
            {15, {15, "numu_cc_1s", "#nu_{#mu}CC single-strange", kSpring + 5, 1001}},
            {16, {16, "numu_cc_ms", "#nu_{#mu}CC multi-strange", kGreen + 2, 1001}},
            {17, {17, "nue_cc", "#nu_{e}CC", kMagenta, 1001}},
            {18, {18, "numu_cc_other", "#nu_{#mu}CC Other", kCyan + 2, 1001}},
            {99, {99, "other", "Other", kCyan, 1001}}};
        return m;
    }
};

} // namespace plot
} // namespace nuxsec

#endif // Nuxsec_PLOT_CHANNELS_H_INCLUDED
