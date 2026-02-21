/* -- C++ -- */
/**
 *  @file  framework/plot/include/PlotChannels.hh
 *
 *  @brief Channel display properties for plotting, including colour choices,
 *         labels, and ordering for stacked outputs.
 */

#ifndef HERON_PLOT_CHANNELS_H
#define HERON_PLOT_CHANNELS_H

#include <map>
#include <string>
#include <vector>

#include "RtypesCore.h"
#include "TColor.h"


namespace nu
{

class Channels
{
  public:
    struct Properties
    {
        int key = 0;
        std::string plain_name;
        std::string tex_label;
        Color_t fill_colour = kBlack;
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
    static int colour(int code) { return properties(code).fill_colour; }
    static int fill_style(int code) { return properties(code).fill_style; }

    static const std::vector<int> &signal_keys()
    {
        static const std::vector<int> v = {15, 16, 17, 18};
        return v;
    }

    static const std::vector<int> &mc_keys()
    {
        static const std::vector<int> v = {1, 2, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 99};
        return v;
    }

  private:
    static Color_t colour_from_hex(const char *hex)
    {
        return static_cast<Color_t>(TColor::GetColor(hex));
    }

    static const std::map<int, Properties> &mapping()
    {
        static const std::map<int, Properties> m = {
            {0, {0, "data", "Data", kBlack, 1001}},
            //
            // Palette aligned to the (MicroBooNE-style) scheme in the reference figure:
            //   - signals: greens
            //   - #nu_{#mu} CC: monochrome blues (light #rightarrow dark for subchannels)
            //   - NC: yellow
            //   - #nu_{e} CC: magenta
            //   - out-FV: brown
            //   - cosmic: neutral grey (hatched)
            //   - other: dark red
            //
            {1, {1, "external", "Cosmic", colour_from_hex("#bababa"), 3345}},
            {2, {2, "out_fv", "Out FV", colour_from_hex("#80603e"), 1001}},
            {10, {10, "numu_cc_np0pi", "#nu_{#mu}CC Np0#pi", colour_from_hex("#2f5cf9"), 1001}},
            {11, {11, "numu_cc_0pnpi", "#nu_{#mu}CC 0p1#pi^{#pm}", colour_from_hex("#2347e0"), 1001}},
            {12, {12, "numu_cc_pi0gg", "#nu_{#mu}CC #pi^{0}/#gamma#gamma", colour_from_hex("#1832c7"), 1001}},
            {13, {13, "numu_cc_npnpi", "#nu_{#mu}CC multi-#pi^{#pm}", colour_from_hex("#0c1dae"), 1001}},
            {14, {14, "nc", "#nu_{x}NC", colour_from_hex("#fbcf38"), 1001}},
            {15, {15, "signal_lambda_ccqe", "Signal #Lambda^{0} CCQE (#Lambda^{0} #rightarrow p#pi^{-})", colour_from_hex("#5cfd3f"), 1001}},
            {16, {16, "signal_lambda_ccres", "Signal #Lambda^{0} CCRES (#Lambda^{0} #rightarrow p#pi^{-})", colour_from_hex("#48ca31"), 1001}},
            {17, {17, "signal_lambda_ccdis", "Signal #Lambda^{0} CCDIS (#Lambda^{0} #rightarrow p#pi^{-})", colour_from_hex("#7ae582"), 1001}},
            {18, {18, "signal_lambda_ccother", "Signal #Lambda^{0} CC Other (#Lambda^{0} #rightarrow p#pi^{-})", colour_from_hex("#2dc653"), 1001}},
            {19, {19, "nue_cc", "#nu_{e}CC", colour_from_hex("#c110f9"), 1001}},
            {20, {20, "numu_cc_other", "#nu_{#mu}CC Other", colour_from_hex("#000895"), 1001}},
            {99, {99, "other", "Other", colour_from_hex("#c32910"), 1001}}};
        return m;
    }
};

} // namespace nu


#endif // HERON_PLOT_CHANNELS_H
