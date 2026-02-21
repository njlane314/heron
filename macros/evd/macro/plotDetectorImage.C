/* -- C++ -- */
/// \file evd/macro/plotDetectorImage.C
/// \brief Render detector images for a specific run/subrun/event.

#include <ROOT/RDataFrame.hxx>

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../../../mod/evd/include/EventDisplay.hh"

using namespace heron::evd;

namespace
{

bool file_exists(const std::string &path)
{
    std::ifstream f(path.c_str());
    return f.good();
}

std::string find_default_event_list_path()
{
    const std::vector<std::string> candidates{
        "/exp/uboone/data/users/nlane/heron/out/event/events.root"};

    for (const auto &path : candidates)
    {
        if (file_exists(path))
            return path;
    }

    return "";
}

} // namespace

void plot_detector_image(const std::string &input_file,
                         int run,
                         int sub,
                         int evt,
                         const std::string &out_dir = "./plots/evd_detector",
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
    opt.mode = EventDisplay::Mode::Detector;

    EventDisplay::render_from_rdf(df, opt);
}

void plot_detector_image_by_channel(const std::string &input_file,
                                    int channel,
                                    const std::string &out_dir = "./plots/evd_detector",
                                    const std::string &tree_name = "events")
{
    ROOT::RDataFrame df(tree_name, input_file);
    auto filtered = df.Filter([channel](int ch) { return ch == channel; },
                              {"analysis_channels"});
    const auto n_rows = static_cast<std::size_t>(filtered.Count().GetValue());
    if (n_rows == 0)
    {
        std::cerr << "[plot_detector_image_by_channel] No events found for channel "
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
    opt.mode = EventDisplay::Mode::Detector;

    EventDisplay::render_from_rdf(picked, opt);
}

void plot_random_detector_image(const std::string &input_file,
                                const std::string &out_dir = "./plots/evd_detector",
                                const std::string &tree_name = "events")
{
    ROOT::RDataFrame df(tree_name, input_file);
    const auto n_rows = static_cast<std::size_t>(df.Count().GetValue());
    if (n_rows == 0)
    {
        std::cerr << "[plot_random_detector_image] No events found in input file '"
                  << input_file << "'.\n";
        return;
    }

    std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<std::size_t> dist(0, n_rows - 1);
    const auto idx = static_cast<ULong64_t>(dist(rng));
    auto picked = df.Range(idx, idx + 1);

    EventDisplay::BatchOptions opt;
    opt.n_events = 1;
    opt.out_dir = out_dir;
    opt.mode = EventDisplay::Mode::Detector;

    EventDisplay::render_from_rdf(picked, opt);
}

/// \brief ROOT entry-point matching this macro's CamelCase filename.
///
/// Call this overload explicitly for default behaviour:
///   heron macro plotDetectorImage.C \
///     'plotDetectorImage("/path/to/event_list.root", 1, 1, 1)'
void plotDetectorImage()
{
    const auto input_file = find_default_event_list_path();
    if (input_file.empty())
    {
        std::cerr
            << "No default event list file found.\n"
            << "Usage: plotDetectorImage(\"/path/to/event_list.root\", run, subrun, event"
            << "[, out_dir[, tree_name]])\n"
            << "       plotDetectorImage(\"/path/to/event_list.root\"[, out_dir[, tree_name]])\n"
            << "       plotDetectorImageByChannel(\"/path/to/event_list.root\", channel"
            << "[, out_dir[, tree_name]])\n";
        return;
    }

    std::cout << "[plotDetectorImage] Using default event list: " << input_file << "\n";
    plot_random_detector_image(input_file);
}

void plotDetectorImage(const std::string &input_file,
                       int run,
                       int sub,
                       int evt,
                       const std::string &out_dir = "./plots/evd_detector",
                       const std::string &tree_name = "events")
{
    plot_detector_image(input_file, run, sub, evt, out_dir, tree_name);
}

/// \brief Entry-point for plotting a random event from an event list file.
void plotDetectorImage(const std::string &input_file,
                       const std::string &out_dir = "./plots/evd_detector",
                       const std::string &tree_name = "events")
{
    plot_random_detector_image(input_file, out_dir, tree_name);
}

/// \brief Alternate entry-point using an analysis channel lookup.
void plotDetectorImageByChannel(const std::string &input_file,
                                int channel,
                                const std::string &out_dir = "./plots/evd_detector",
                                const std::string &tree_name = "events")
{
    plot_detector_image_by_channel(input_file, channel, out_dir, tree_name);
}
