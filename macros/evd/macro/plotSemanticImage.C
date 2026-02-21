/* -- C++ -- */
/// \file evd/macro/plotSemanticImage.C
/// \brief Render semantic images for a specific run/subrun/event.

#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "../../../framework/evd/include/EventDisplay.hh"

using namespace heron::evd;

void plot_semantic_image(const std::string &input_file,
                         int run,
                         int sub,
                         int evt,
                         const std::string &out_dir = "./plots/evd_semantic",
                         const std::string &tree_name = "events")
{
    ROOT::RDataFrame df(tree_name, input_file);

    EventDisplay::BatchOptions opt;
    std::ostringstream selection;
    selection << opt.cols.run << " == " << run
              << " && " << opt.cols.sub << " == " << sub
              << " && " << opt.cols.evt << " == " << evt;
    opt.selection_expr = selection.str();
    opt.n_events = 1;
    opt.out_dir = out_dir;
    opt.mode = EventDisplay::Mode::Semantic;

    EventDisplay::render_from_rdf(df, opt);
}

void plot_semantic_image_by_channel(const std::string &input_file,
                                    int channel,
                                    const std::string &out_dir = "./plots/evd_semantic",
                                    const std::string &tree_name = "events")
{
    ROOT::RDataFrame df(tree_name, input_file);
    auto filtered = df.Filter([channel](int ch) { return ch == channel; },
                              {"analysis_channels"});
    const auto n_rows = static_cast<std::size_t>(filtered.Count().GetValue());
    if (n_rows == 0)
    {
        std::cerr << "[plot_semantic_image_by_channel] No events found for channel "
                  << channel << ".\n";
        return;
    }

    std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<std::size_t> dist(0, n_rows - 1);
    const auto idx = static_cast<ULong64_t>(dist(rng));
    auto picked = filtered.Range(idx, idx + 1);

    EventDisplay::BatchOptions opt;
    opt.n_events = 1;
    opt.out_dir = out_dir;
    opt.mode = EventDisplay::Mode::Semantic;

    EventDisplay::render_from_rdf(picked, opt);
}
