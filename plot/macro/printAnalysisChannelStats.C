// plot/macro/printAnalysisChannelStats.C
//
// Run with:
//   ./heron macro printAnalysisChannelStats.C
//   ./heron macro printAnalysisChannelStats.C 'printAnalysisChannelStats("/path/to/event_list.root","true","w_nominal",false)'
//   ./heron macro printAnalysisChannelStats.C 'printAnalysisChannelStats("/path/to/event_list.root","your_selection_here","w_nominal",true)'
//
// Notes:
//   - Input is an event_list_<analysis>.root (same event-list mode as plotStackedHistTrueVertex.C).
//   - Expects input event list to already contain "analysis_channels".
//   - Prints per-channel counts and weighted yields (MC/EXT weighted by mc_weight).
//   - Data (if requested) is printed as counts (no weight).
//   - Channel names/labels follow PlotChannels.hh (nu::Channels).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>
#include <TH1D.h>

#include "EventListIO.hh"
#include "PlotChannels.hh"
#include "PlottingHelper.hh"
#include "SampleCLI.hh"

using namespace nu;

namespace
{
bool looks_like_event_list_root(const std::string &p)
{
    const auto n = p.size();
    if (n < 5 || p.substr(n - 5) != ".root")
        return false;

    std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
    if (!f || f->IsZombie())
        return false;

    const bool has_refs = (f->Get("sample_refs") != nullptr);
    const bool has_events_tree = (f->Get("events") != nullptr);
    const bool has_event_tree_key = (f->Get("event_tree") != nullptr);

    return has_refs && (has_events_tree || has_event_tree_key);
}

struct Totals
{
    long long n = 0;
    double sumw = 0.0;
    double err = 0.0;
};

Totals totals_from_hists(const TH1D &h_counts, const TH1D &h_w)
{
    Totals t;
    t.n = static_cast<long long>(std::llround(h_counts.Integral(1, h_counts.GetNbinsX())));
    double err = 0.0;
    t.sumw = h_w.IntegralAndError(1, h_w.GetNbinsX(), err);
    t.err = err;
    return t;
}


bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto columns = node.GetColumnNames();
    return std::find(columns.begin(), columns.end(), name) != columns.end();
}

bool is_unmapped_channel_code(int code)
{
    const auto &p = Channels::properties(code);
    // Channels::properties returns mapping(99) for unknown codes.
    return (p.key == 99 && code != 99);
}
} // namespace

int printAnalysisChannelStats(const std::string &samples_tsv = "",
                              const std::string &extra_sel = "true",
                              const std::string &mc_weight = "w_nominal",
                              bool include_data = false)
{
    ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;

    std::cout << "[printAnalysisChannelStats] input=" << list_path << "\n";
    std::cout << "[printAnalysisChannelStats] extra_sel=" << (extra_sel.empty() ? std::string("(none)") : extra_sel) << "\n";
    std::cout << "[printAnalysisChannelStats] mc_weight=" << (mc_weight.empty() ? std::string("w_nominal") : mc_weight) << "\n";
    std::cout << "[printAnalysisChannelStats] include_data=" << (include_data ? "true" : "false") << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[printAnalysisChannelStats] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();
    auto mask_data = el.mask_for_data();

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) {
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf;

    ROOT::RDF::RNode node_ext = filter_by_mask(base, mask_ext);
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<size_t>(sid)]);
                                   },
                                   {"sample_id"});

    ROOT::RDF::RNode node_data = filter_by_mask(base, mask_data);

    if (!has_column(base, "analysis_channels"))
    {
        std::cerr << "[printAnalysisChannelStats] missing required column 'analysis_channels' in input event list.\n"
                  << "[printAnalysisChannelStats] this macro expects an event list that already includes derived channels.\n";
        return 1;
    }

    // Apply extra selection.
    if (!extra_sel.empty())
    {
        node_mc = node_mc.Filter(extra_sel);
        node_ext = node_ext.Filter(extra_sel);
        if (include_data)
            node_data = node_data.Filter(extra_sel);
    }

    // Materialise the weight expression as a column (so mc_weight can be an expression).
    const std::string wexpr = mc_weight.empty() ? "w_nominal" : mc_weight;
    node_mc = node_mc.Define("__w__", wexpr);
    node_ext = node_ext.Define("__w__", wexpr);

    // Histogram axis covers Channel-like integer codes (0..99).
    const int nbins = 100;
    const double xmin = -0.5;
    const double xmax = 99.5;

    auto h_mc_n = node_mc.Histo1D({"h_mc_n", "MC;analysis_channels;N", nbins, xmin, xmax}, "analysis_channels");
    auto h_mc_w = node_mc.Histo1D({"h_mc_w", "MC;analysis_channels;sumw", nbins, xmin, xmax}, "analysis_channels", "__w__");

    auto h_ext_n = node_ext.Histo1D({"h_ext_n", "EXT;analysis_channels;N", nbins, xmin, xmax}, "analysis_channels");
    auto h_ext_w = node_ext.Histo1D({"h_ext_w", "EXT;analysis_channels;sumw", nbins, xmin, xmax}, "analysis_channels", "__w__");

    ROOT::RDF::RResultPtr<TH1D> h_data_n;
    if (include_data)
        h_data_n = node_data.Histo1D({"h_data_n", "Data;analysis_channels;N", nbins, xmin, xmax}, "analysis_channels");

    if (include_data)
        ROOT::RDF::RunGraphs({h_mc_n, h_mc_w, h_ext_n, h_ext_w, h_data_n});
    else
        ROOT::RDF::RunGraphs({h_mc_n, h_mc_w, h_ext_n, h_ext_w});

    const Totals t_mc = totals_from_hists(*h_mc_n, *h_mc_w);
    const Totals t_ext = totals_from_hists(*h_ext_n, *h_ext_w);
    const bool show_ext = (t_ext.n != 0) || (std::fabs(t_ext.sumw) > 0.0) || (std::fabs(t_ext.err) > 0.0);

    const double stack_sumw = t_mc.sumw + t_ext.sumw;
    const double stack_err = std::sqrt(t_mc.err * t_mc.err + t_ext.err * t_ext.err);
    const long long stack_n = t_mc.n + t_ext.n;

    long long data_n = 0;
    if (include_data)
        data_n = static_cast<long long>(std::llround(h_data_n->Integral(1, h_data_n->GetNbinsX())));

    std::cout << "\n[printAnalysisChannelStats] totals after selection:\n";
    std::cout << "  MC   : N=" << t_mc.n << "  sumw=" << t_mc.sumw << " +/- " << t_mc.err << "\n";
    if (show_ext)
        std::cout << "  EXT  : N=" << t_ext.n << "  sumw=" << t_ext.sumw << " +/- " << t_ext.err << "\n";
    std::cout << "  STACK: N=" << stack_n << "  sumw=" << stack_sumw << " +/- " << stack_err << "\n";
    if (include_data)
        std::cout << "  DATA : N=" << data_n << "\n";

    // Start from the standard plotting keys (+ 0), then add any additional non-zero codes found.
    std::set<int> codes;
    codes.insert(0);
    for (int k : Channels::mc_keys())
        codes.insert(k);

    for (int code = 0; code <= 99; ++code)
    {
        const int bin = code + 1;
        const double n_any =
            h_mc_n->GetBinContent(bin) + h_ext_n->GetBinContent(bin) + (include_data ? h_data_n->GetBinContent(bin) : 0.0);
        if (n_any != 0.0)
            codes.insert(code);
    }

    std::cout << "\n[printAnalysisChannelStats] per-channel breakdown:\n";

    std::cout << std::left
              << std::setw(6) << "Code"
              << std::setw(22) << "Name"
              << std::setw(30) << "Label"
              << std::right
              << std::setw(12) << "MC N"
              << std::setw(14) << "MC sumw"
              << std::setw(14) << "MC err";
    if (show_ext)
        std::cout << std::setw(12) << "EXT N"
                  << std::setw(14) << "EXT sumw"
                  << std::setw(14) << "EXT err";
    if (include_data)
        std::cout << std::setw(12) << "DATA N";
    std::cout << std::setw(14) << "STACK sumw"
              << std::setw(14) << "STACK err"
              << std::setw(10) << "Frac"
              << std::setw(10) << "Map?"
              << "\n";

    std::cout << std::fixed << std::setprecision(3);

    for (int code : codes)
    {
        const int bin = code + 1;

        const long long n_mc = static_cast<long long>(std::llround(h_mc_n->GetBinContent(bin)));
        const double y_mc = h_mc_w->GetBinContent(bin);
        const double e_mc = h_mc_w->GetBinError(bin);

        const long long n_ext = static_cast<long long>(std::llround(h_ext_n->GetBinContent(bin)));
        const double y_ext = h_ext_w->GetBinContent(bin);
        const double e_ext = h_ext_w->GetBinError(bin);

        long long n_data = 0;
        if (include_data)
            n_data = static_cast<long long>(std::llround(h_data_n->GetBinContent(bin)));

        // Skip codes that are completely empty everywhere.
        if (n_mc == 0 && n_ext == 0 && (!include_data || n_data == 0))
            continue;

        const double y_stack = y_mc + y_ext;
        const double e_stack = std::sqrt(e_mc * e_mc + e_ext * e_ext);
        const double frac = (stack_sumw != 0.0) ? (y_stack / stack_sumw) : 0.0;

        const std::string mapflag = is_unmapped_channel_code(code) ? "no" : "yes";

        std::cout << std::left
                  << std::setw(6) << code
                  << std::setw(22) << Channels::name(code)
                  << std::setw(30) << Channels::label(code)
                  << std::right
                  << std::setw(12) << n_mc
                  << std::setw(14) << y_mc
                  << std::setw(14) << e_mc;
        if (show_ext)
            std::cout << std::setw(12) << n_ext
                      << std::setw(14) << y_ext
                      << std::setw(14) << e_ext;
        if (include_data)
            std::cout << std::setw(12) << n_data;
        std::cout << std::setw(14) << y_stack
                  << std::setw(14) << e_stack
                  << std::setw(10) << frac
                  << std::setw(10) << mapflag
                  << "\n";
    }

    std::cout.flush();
    return 0;
}
