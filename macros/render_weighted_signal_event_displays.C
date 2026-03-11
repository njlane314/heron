// macros/render_weighted_signal_event_displays.C
//
// Render event displays for signal MC events after a score cut.
//
// Important:
//   An event display is per-event. The event weight does not change the image
//   itself; it changes which events are selected for the gallery.
//
// Candidate set:
//   - pure MC only (MC-like with EXT removed)
//   - pass base_sel
//   - pass signal_sel  => signal
//   - pass inf_scores[0] > score_cut
//   - have positive finite weight_col
//
// Selection modes for the gallery:
//   - mode="top"    : choose the highest-weight events deterministically
//   - mode="sample" : weighted sampling without replacement using weight_col
//
// Notes:
//   - The rendered PDF follows the underlying dataframe/tree order of the
//     selected events, not the weight rank. The chosen_events.tsv sidecar gives
//     the chosen rows, weights, and scores explicitly.
//   - By default this targets weighted MC signal. If you also want EXT,
//     add a row-wise EXT flag and keep EXT rows alongside signal MC rows.
//
// Run with:
//   ./heron macro render_weighted_signal_event_displays.C
//   ./heron macro render_weighted_signal_event_displays.C \
//     'render_weighted_signal_event_displays("./scratch/out/event_list.root")'
//
// Examples:
//   ./heron macro render_weighted_signal_event_displays.C \
//     'render_weighted_signal_event_displays("", "sel_muon", "is_signal", "w_nominal", 6.0, 12, "top")'
//
//   ./heron macro render_weighted_signal_event_displays.C \
//     'render_weighted_signal_event_displays("", "sel_muon", "is_signal", "w_nominal", 6.0, 12, "sample", "detector", "plots/sig_score_gt6", "sig_score_gt6.pdf", 12345ULL)'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#if defined(__CLING__)
R__ADD_INCLUDE_PATH(framework/core/include)
R__ADD_INCLUDE_PATH(framework/modules/ana/include)
R__ADD_INCLUDE_PATH(framework/modules/io/include)
R__ADD_INCLUDE_PATH(framework/modules/plot/include)
R__ADD_INCLUDE_PATH(framework/modules/evd/include)
#endif

#include "EventDisplay.hh"
#include "EventListIO.hh"
#include "PlottingHelper.hh"
#include "SelectionService.hh"
#include "include/MacroGuard.hh"
#include "include/MacroIO.hh"

using namespace nu;

namespace
{

bool has_column(ROOT::RDF::RNode node, const std::string &name)
{
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

bool require_columns(ROOT::RDF::RNode node,
                     const std::vector<std::string> &cols,
                     const std::string &label)
{
    bool ok = true;
    for (const auto &c : cols)
    {
        if (!has_column(node, c))
        {
            std::cerr << "[render_weighted_signal_event_displays] missing "
                      << label << " column '" << c << "'.\n";
            ok = false;
        }
    }
    return ok;
}

struct Candidate
{
    ULong64_t entry = 0;
    int sample_id = -1;
    int run = 0;
    int sub = 0;
    int evt = 0;
    double score = 0.0;
    double weight = 0.0;
};

struct EventKey
{
    int sample_id = -1;
    int run = 0;
    int sub = 0;
    int evt = 0;

    bool operator==(const EventKey &other) const
    {
        return std::tie(sample_id, run, sub, evt) == std::tie(other.sample_id, other.run, other.sub, other.evt);
    }
};

struct EventKeyHash
{
    std::size_t operator()(const EventKey &key) const
    {
        std::size_t seed = 0u;
        seed ^= static_cast<std::size_t>(key.sample_id) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.run) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.sub) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.evt) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        return seed;
    }
};

std::vector<std::size_t> pick_top_weight(const std::vector<Candidate> &cand,
                                         std::size_t n_pick)
{
    std::vector<std::size_t> idx(cand.size());
    std::iota(idx.begin(), idx.end(), std::size_t{0});

    std::stable_sort(idx.begin(), idx.end(),
                     [&](std::size_t a, std::size_t b) {
                         return cand[a].weight > cand[b].weight;
                     });

    if (idx.size() > n_pick)
        idx.resize(n_pick);

    return idx;
}

// Efraimidis-Spirakis weighted sampling without replacement.
// For positive weights w_i, draw u_i ~ U(0,1] and keep the largest keys
//
//   key_i = log(u_i) / w_i
//
// so larger weights are more likely to be selected.
std::vector<std::size_t> pick_weighted_sample(const std::vector<Candidate> &cand,
                                              std::size_t n_pick,
                                              unsigned long long seed)
{
    struct KeyedIndex
    {
        double key = -std::numeric_limits<double>::infinity();
        std::size_t idx = 0;
    };

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(std::numeric_limits<double>::min(), 1.0);

    std::vector<KeyedIndex> keyed;
    keyed.reserve(cand.size());

    for (std::size_t i = 0; i < cand.size(); ++i)
    {
        const double w = cand[i].weight;
        if (!(w > 0.0) || !std::isfinite(w))
            continue;

        const double u = uni(rng);
        const double key = std::log(u) / w; // <= 0; values closer to 0 win
        keyed.push_back(KeyedIndex{key, i});
    }

    std::stable_sort(keyed.begin(), keyed.end(),
                     [](const KeyedIndex &a, const KeyedIndex &b) {
                         return a.key > b.key;
                     });

    if (keyed.size() > n_pick)
        keyed.resize(n_pick);

    std::vector<std::size_t> out;
    out.reserve(keyed.size());
    for (const auto &x : keyed)
        out.push_back(x.idx);

    return out;
}

void write_selected_tsv(const std::string &path,
                        const std::vector<Candidate> &cand,
                        const std::vector<std::size_t> &chosen,
                        const std::string &weight_col,
                        const std::string &mode)
{
    std::ofstream ofs(path);
    ofs << "rank\tmode\tentry\tsample_id\trun\tsub\tevt\tinf_score_0\t" << weight_col << "\n";

    for (std::size_t rank = 0; rank < chosen.size(); ++rank)
    {
        const auto &c = cand[chosen[rank]];
        ofs << (rank + 1) << "\t" << mode << "\t"
            << c.entry << "\t"
            << c.sample_id << "\t"
            << c.run << "\t"
            << c.sub << "\t"
            << c.evt << "\t"
            << std::setprecision(17) << c.score << "\t"
            << std::setprecision(17) << c.weight << "\n";
    }
}

} // namespace

int render_weighted_signal_event_displays(
    const std::string &event_list_path = "",
    const std::string &base_sel = "sel_muon",
    const std::string &signal_sel = "is_signal",
    const std::string &weight_col = "w_nominal",
    double score_cut = 6.0,
    unsigned long long n_events = 12,
    const std::string &mode = "top",
    const std::string &display_mode = "detector",
    const std::string &out_dir = "plots/sig_event_displays_score_gt6",
    const std::string &combined_pdf = "sig_event_displays_score_gt6.pdf",
    unsigned long long rng_seed = 12345ULL)
{
    return heron::macro::run_with_guard("render_weighted_signal_event_displays", [&]() -> int {
        const std::string input_path =
            event_list_path.empty() ? default_event_list_root() : event_list_path;

        std::cout << "[render_weighted_signal_event_displays] input="
                  << input_path << "\n";

        if (!looks_like_event_list_root(input_path))
        {
            std::cerr << "[render_weighted_signal_event_displays] input is not an event-list root file: "
                      << input_path << "\n";
            return 1;
        }

        if (n_events == 0)
        {
            std::cerr << "[render_weighted_signal_event_displays] n_events must be > 0\n";
            return 1;
        }

        std::vector<Candidate> cand;
        std::vector<std::size_t> chosen;
        const auto mode_enum = heron::evd::EventDisplay::parse_mode(display_mode);

        {
            ROOT::EnableImplicitMT();

            EventListIO el(input_path);

            ROOT::RDF::RNode base = SelectionService::decorate(el.rdf())
                                    .Define("__entry__",
                                            [](ULong64_t e) { return e; },
                                            {"rdfentry_"})
                                    .Define("inf_score_0",
                                            [](const ROOT::RVec<float> &scores) {
                                                return scores.empty() ? -1.0e9 : static_cast<double>(scores[0]);
                                            },
                                            {"inf_scores"});

        if (!has_column(base, weight_col))
        {
            std::cerr << "[render_weighted_signal_event_displays] missing weight column '"
                      << weight_col << "'.\n";
            return 1;
        }

        if (!require_columns(base,
                             {"sample_id", "run", "sub", "evt"},
                             "event id"))
            return 1;

        if (mode_enum == heron::evd::EventDisplay::Mode::Detector)
        {
            if (!require_columns(base,
                                 {"detector_image_u", "detector_image_v", "detector_image_w"},
                                 "detector image"))
                return 1;
        }
        else
        {
            if (!require_columns(base,
                                 {"semantic_image_u", "semantic_image_v", "semantic_image_w"},
                                 "semantic image"))
                return 1;
        }

        base = base.Define("__is_signal__", signal_sel)
                   .Define("__w__",
                           [](double w) {
                               return std::isfinite(w) ? w : 0.0;
                           },
                           {weight_col});

        auto mask_mc_like = el.mask_for_mc_like();
        auto mask_ext = el.mask_for_ext();

        // Pure MC only: signal = MC-like excluding EXT and requiring signal.
        ROOT::RDF::RNode node = filter_by_sample_mask(base, mask_mc_like, "sample_id");
        node = filter_not_sample_mask(node, mask_ext, "sample_id");

        if (!base_sel.empty())
        {
            if (has_column(node, base_sel))
                node = node.Filter([](bool pass) { return pass; }, {base_sel});
            else
                node = node.Filter(base_sel);
        }

        node = node.Filter([](bool is_sig) { return is_sig; }, {"__is_signal__"})
                   .Filter([score_cut](double s) {
                               return std::isfinite(s) && (s > score_cut);
                           },
                           {"inf_score_0"})
                   .Filter([](double w) {
                               return std::isfinite(w) && (w > 0.0);
                           },
                           {"__w__"});

        auto n_cand_ptr = node.Count();
        auto sumw_ptr = node.Sum<double>("__w__");
        auto entries_ptr = node.Take<ULong64_t>("__entry__");
        auto sample_ids_ptr = node.Take<int>("sample_id");
        auto runs_ptr = node.Take<int>("run");
        auto subs_ptr = node.Take<int>("sub");
        auto evts_ptr = node.Take<int>("evt");
        auto scores_ptr = node.Take<double>("inf_score_0");
        auto weights_ptr = node.Take<double>("__w__");

        const auto n_cand = *n_cand_ptr;
        if (n_cand == 0)
        {
            std::cerr << "[render_weighted_signal_event_displays] no signal MC candidates after"
                      << " base_sel='" << base_sel << "', " << signal_sel << ","
                      << " inf_score_0 > " << score_cut
                      << ", and " << weight_col << " > 0.\n";
            return 1;
        }

        const double sumw = *sumw_ptr;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[render_weighted_signal_event_displays] candidate signal events     = "
                  << n_cand << "\n";
        std::cout << "[render_weighted_signal_event_displays] candidate weighted yield    = "
                  << sumw << "\n";

        const std::vector<ULong64_t> entries = *entries_ptr;
        const std::vector<int> sample_ids = *sample_ids_ptr;
        const std::vector<int> runs = *runs_ptr;
        const std::vector<int> subs = *subs_ptr;
        const std::vector<int> evts = *evts_ptr;
        const std::vector<double> scores = *scores_ptr;
        const std::vector<double> weights = *weights_ptr;

        cand.resize(entries.size());
        for (std::size_t i = 0; i < cand.size(); ++i)
        {
            cand[i].entry = entries[i];
            cand[i].sample_id = sample_ids[i];
            cand[i].run = runs[i];
            cand[i].sub = subs[i];
            cand[i].evt = evts[i];
            cand[i].score = scores[i];
            cand[i].weight = weights[i];
        }

        if (n_events > cand.size())
            n_events = cand.size();

        if (mode == "sample" || mode == "weighted_sample")
        {
            chosen = pick_weighted_sample(cand, static_cast<std::size_t>(n_events), rng_seed);
        }
        else if (mode == "top" || mode == "top_weight")
        {
            chosen = pick_top_weight(cand, static_cast<std::size_t>(n_events));
        }
        else
        {
            std::cerr << "[render_weighted_signal_event_displays] unknown mode='"
                      << mode << "'. Use 'top' or 'sample'.\n";
            return 1;
        }

        if (chosen.empty())
        {
            std::cerr << "[render_weighted_signal_event_displays] no events selected for rendering.\n";
            return 1;
        }

        const double chosen_sumw = [&]() {
            double out = 0.0;
            for (const auto idx : chosen)
                out += cand[idx].weight;
            return out;
        }();

        std::filesystem::create_directories(out_dir);

        const std::string selected_tsv =
            (std::filesystem::path(out_dir) / "chosen_events.tsv").string();
        write_selected_tsv(selected_tsv, cand, chosen, weight_col, mode);

        std::cout << "[render_weighted_signal_event_displays] selected events        = "
                  << chosen.size() << "\n";
        std::cout << "[render_weighted_signal_event_displays] selected weighted yield = "
                  << chosen_sumw;
        if (sumw > 0.0)
            std::cout << " (" << (100.0 * chosen_sumw / sumw) << "% of candidate weight)";
        std::cout << "\n";
        std::cout << "[render_weighted_signal_event_displays] wrote selected event list: "
                  << selected_tsv << "\n";
        std::cout << "[render_weighted_signal_event_displays] note: rendered PDF order follows dataframe order, not chosen rank.\n";

        std::cout << "[render_weighted_signal_event_displays] chosen events:\n";
        for (std::size_t rank = 0; rank < chosen.size(); ++rank)
        {
            const auto &c = cand[chosen[rank]];
            std::cout << "  [" << (rank + 1) << "]"
                      << " run=" << c.run
                      << " sub=" << c.sub
                      << " evt=" << c.evt
                      << " score=" << c.score
                      << " " << weight_col << "=" << c.weight
                      << "\n";
        }

        }

        if (ROOT::IsImplicitMTEnabled())
            ROOT::DisableImplicitMT();

        auto chosen_events = std::make_shared<std::unordered_set<EventKey, EventKeyHash>>();
        for (const auto idx : chosen)
        {
            chosen_events->insert(EventKey{cand[idx].sample_id, cand[idx].run, cand[idx].sub, cand[idx].evt});
        }

        EventListIO el_render(input_path);
        ROOT::RDF::RNode node_chosen =
            SelectionService::decorate(el_render.rdf())
                .Filter([chosen_events](int sample_id, int run, int sub, int evt) {
                            return chosen_events->count(EventKey{sample_id, run, sub, evt}) != 0u;
                        },
                        {"sample_id", "run", "sub", "evt"});

        heron::evd::EventDisplay::BatchOptions opt;
        opt.selection_expr = "";
        opt.n_events = chosen_events->size();
        opt.out_dir = out_dir;
        opt.image_format = "pdf";
        opt.combined_pdf = combined_pdf;
        opt.manifest_path =
            (std::filesystem::path(out_dir) / "render_manifest.json").string();
        opt.mode = mode_enum;
        opt.display.out_dir = out_dir;
        opt.display.use_log_z = true;
        opt.display.show_legend = true;

        heron::evd::EventDisplay::render_from_rdf(node_chosen, opt);

        std::cout << "[render_weighted_signal_event_displays] output dir = "
                  << out_dir << "\n";
        if (!combined_pdf.empty())
        {
            std::cout << "[render_weighted_signal_event_displays] combined pdf = "
                      << (std::filesystem::path(out_dir) / combined_pdf).string() << "\n";
        }
        std::cout << "[render_weighted_signal_event_displays] done\n";
        return 0;
    });
}
