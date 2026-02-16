// plot/macro/plotMuonCandidateMipnessVsLengthByChannel.C
//
// Dedicated macro for post-muon-candidate diagnostics:
//   - x-axis: muon-candidate track length [cm]
//   - y-axis: muon-candidate MIPness (median planes)
//
// It writes three MC-only 2D histograms filtered by analysis channel:
//   - NC (analysis_channels == 14)
//   - nue CC (analysis_channels == 17)
//   - inclusive numu CC (analysis_channels in {10,11,12,13,15,16,18})
//
// Run with:
//   ./heron macro plotMuonCandidateMipnessVsLengthByChannel.C
//   ./heron macro plotMuonCandidateMipnessVsLengthByChannel.C \
//       'plotMuonCandidateMipnessVsLengthByChannel("./scratch/out/event_list_myana.root")'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "EventListIO.hh"
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

bool implicit_mt_enabled()
{
    const char *env = std::getenv("HERON_PLOT_IMT");
    return env != nullptr && std::string(env) != "0";
}

std::string getenv_or(const char *key, const std::string &fallback)
{
    const char *v = std::getenv(key);
    if (!v || std::string(v).empty())
        return fallback;
    return std::string(v);
}

std::string sanitize_for_filename(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (const char c : s)
    {
        const bool ok =
            (c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') ||
            (c == '_');
        out.push_back(ok ? c : '_');
    }
    return out;
}

bool has_column(const std::vector<std::string> &cols, const std::string &name)
{
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

std::string channel_is_one_of_expr(const std::vector<int> &channels)
{
    if (channels.empty())
        return "false";

    std::string expr = "(";
    for (std::size_t i = 0; i < channels.size(); ++i)
    {
        if (i > 0)
            expr += " || ";
        expr += "analysis_channels == " + std::to_string(channels[i]);
    }
    expr += ")";
    return expr;
}

std::string format_decimal(double value, int precision = 1)
{
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(precision);
    os << value;
    return os.str();
}

double weighted_sum(const ROOT::RDF::RNode &node,
                    const std::string &weight_col,
                    const std::string &selection)
{
    ROOT::RDF::RNode n = node;
    if (!selection.empty())
        n = n.Filter(selection);
    return n.Sum<double>(weight_col).GetValue();
}

ROOT::RDF::RNode define_mu_cand_columns(ROOT::RDF::RNode node,
                                        float mu_score_cut,
                                        float mu_len_cut_cm,
                                        float mu_dist_cut_cm,
                                        unsigned mu_req_gen)
{
    using ROOT::RVec;

    node = node.Define(
        "mu_cand_idx",
        [mu_score_cut, mu_len_cut_cm, mu_dist_cut_cm, mu_req_gen](const RVec<float> &scores,
                                                                  const RVec<float> &lengths,
                                                                  const RVec<float> &distances,
                                                                  const RVec<unsigned> &generations) -> int {
            const std::size_t n = std::min(std::min(scores.size(), lengths.size()),
                                           std::min(distances.size(), generations.size()));
            int best = -1;
            float best_len = -std::numeric_limits<float>::infinity();

            for (std::size_t i = 0; i < n; ++i)
            {
                const float s = scores[i];
                const float L = lengths[i];
                const float d = distances[i];

                if (!(std::isfinite(s) && std::isfinite(L) && std::isfinite(d)))
                    continue;

                const bool ok = (s > mu_score_cut) &&
                                (L > mu_len_cut_cm) &&
                                (d < mu_dist_cut_cm) &&
                                (generations[i] == mu_req_gen);

                if (!ok)
                    continue;

                if (L > best_len)
                {
                    best_len = L;
                    best = static_cast<int>(i);
                }
            }
            return best;
        },
        {"track_shower_scores", "track_length", "track_distance_to_vertex", "pfp_generations"});

    node = node.Define(
        "mu_cand_length",
        [](int idx, const RVec<float> &v) -> float {
            if (idx < 0 || static_cast<std::size_t>(idx) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(idx)];
        },
        {"mu_cand_idx", "track_length"});

    return node;
}

ROOT::RDF::RNode define_mipness_columns(ROOT::RDF::RNode node,
                                        double Mref_mev_per_cm,
                                        double tail_weight_w)
{
    using ROOT::RVec;

    auto D_from_stats = [Mref_mev_per_cm, tail_weight_w](float T70, float Q50, float Q90) -> double {
        if (!(std::isfinite(T70) && std::isfinite(Q50) && std::isfinite(Q90)))
            return std::numeric_limits<double>::quiet_NaN();
        if (!(T70 > 0.0f && Q50 > 0.0f && Q90 > 0.0f))
            return std::numeric_limits<double>::quiet_NaN();
        return std::log(static_cast<double>(T70) / Mref_mev_per_cm) +
               tail_weight_w * std::log(static_cast<double>(Q90) / static_cast<double>(Q50));
    };

    auto mipness_med_at = [D_from_stats](int idx,
                                         const RVec<float> &T70u, const RVec<float> &Q50u, const RVec<float> &Q90u,
                                         const RVec<float> &T70v, const RVec<float> &Q50v, const RVec<float> &Q90v,
                                         const RVec<float> &T70y, const RVec<float> &Q50y, const RVec<float> &Q90y) -> float {
        auto get = [](int i, const RVec<float> &v) -> float {
            if (i < 0 || static_cast<std::size_t>(i) >= v.size())
                return std::numeric_limits<float>::quiet_NaN();
            return v[static_cast<std::size_t>(i)];
        };

        std::vector<double> Ds;
        Ds.reserve(3);

        const double Du = D_from_stats(get(idx, T70u), get(idx, Q50u), get(idx, Q90u));
        const double Dv = D_from_stats(get(idx, T70v), get(idx, Q50v), get(idx, Q90v));
        const double Dy = D_from_stats(get(idx, T70y), get(idx, Q50y), get(idx, Q90y));

        if (std::isfinite(Du))
            Ds.push_back(Du);
        if (std::isfinite(Dv))
            Ds.push_back(Dv);
        if (std::isfinite(Dy))
            Ds.push_back(Dy);

        if (Ds.empty())
            return std::numeric_limits<float>::quiet_NaN();

        std::sort(Ds.begin(), Ds.end());
        const double Dmed = Ds[Ds.size() / 2];
        return static_cast<float>(std::exp(-Dmed));
    };

    node = node.Define(
        "mu_cand_mipness_med",
        mipness_med_at,
        {"mu_cand_idx",
         "track_dedx_T70_u", "track_dedx_Q50_u", "track_dedx_Q90_u",
         "track_dedx_T70_v", "track_dedx_Q50_v", "track_dedx_Q90_v",
         "track_dedx_T70_y", "track_dedx_Q50_y", "track_dedx_Q90_y"});

    return node;
}

void save_2d(const ROOT::RDF::RNode &node,
             const std::string &x,
             const std::string &y,
             const std::string &w,
             const std::string &sel,
             const std::string &tag,
             int nbx, double xmin, double xmax,
             int nby, double ymin, double ymax,
             double x_threshold,
             double y_threshold)
{
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(255);

    const std::string out_dir = getenv_or("HERON_PLOT_DIR", "./scratch/plots");
    const std::string out_fmt = getenv_or("HERON_PLOT_FORMAT", "pdf");
    gSystem->mkdir(out_dir.c_str(), /*recursive=*/true);

    ROOT::RDF::RNode n = node;
    if (!sel.empty())
        n = n.Filter(sel);

    const std::string hname = "h2_" + sanitize_for_filename(tag);
    auto h2 = n.Histo2D({hname.c_str(), "", nbx, xmin, xmax, nby, ymin, ymax}, x, y, w);

    const double sum_w = weighted_sum(n, w, "");
    const double sum_x = weighted_sum(n, w, "(" + x + " > " + format_decimal(x_threshold, 3) + ")");
    const double sum_y = weighted_sum(n, w, "(" + y + " > " + format_decimal(y_threshold, 3) + ")");
    const double sum_xy = weighted_sum(n, w,
                                       "(" + x + " > " + format_decimal(x_threshold, 3) + ") && (" +
                                           y + " > " + format_decimal(y_threshold, 3) + ")");

    const double frac_x = (sum_w > 0.0) ? (100.0 * sum_x / sum_w) : 0.0;
    const double frac_y = (sum_w > 0.0) ? (100.0 * sum_y / sum_w) : 0.0;
    const double frac_xy = (sum_w > 0.0) ? (100.0 * sum_xy / sum_w) : 0.0;

    TCanvas c(("c2_" + sanitize_for_filename(tag)).c_str(), "", 760, 700);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.16);
    c.SetBottomMargin(0.11);
    c.SetTopMargin(0.08);
    c.SetTicks(1, 1);
    c.SetLogz(true);

    h2->SetTitle(";Muon-candidate length [cm];Minimally Ionising Calorimetry");
    h2->GetXaxis()->SetTitleOffset(1.1);
    h2->GetYaxis()->SetTitleOffset(1.3);
    h2->GetZaxis()->SetTitle("Events");
    h2->GetZaxis()->SetTitleOffset(1.1);
    h2->SetMinimum(0.5);
    h2->Draw("COLZ");

    const std::string summary =
        "In-range #Sigmaw: " + format_decimal(sum_w, 1) +
        "; L > " + format_decimal(x_threshold, 1) + " cm: " + format_decimal(frac_x, 1) + "%" +
        "; MIPness > " + format_decimal(y_threshold, 2) + ": " + format_decimal(frac_y, 1) + "%" +
        "; both: " + format_decimal(frac_xy, 1) + "%";

    std::cout << "[plotMuonCandidateMipnessVsLengthByChannel] " << tag << ": " << summary << "\n";

    const std::string out_path = out_dir + "/" + sanitize_for_filename(tag) + "." + out_fmt;
    c.SaveAs(out_path.c_str());
    std::cout << "[plotMuonCandidateMipnessVsLengthByChannel] wrote " << out_path << "\n";
}

} // namespace

int plotMuonCandidateMipnessVsLengthByChannel(const std::string &samples_tsv = "",
                                              const std::string &extra_sel = "true",
                                              double Mref_mev_per_cm = 2.10,
                                              double mip_tail_weight = 1.0,
                                              float mu_score_cut = 0.5f,
                                              float mu_len_cut_cm = 10.0f,
                                              float mu_dist_cut_cm = 4.0f,
                                              unsigned mu_req_gen = 2u,
                                              double mipness_threshold = 0.50)
{
    if (implicit_mt_enabled())
        ROOT::EnableImplicitMT();

    const std::string list_path = samples_tsv.empty() ? default_event_list_root() : samples_tsv;
    std::cout << "[plotMuonCandidateMipnessVsLengthByChannel] input=" << list_path << "\n";

    if (!looks_like_event_list_root(list_path))
    {
        std::cerr << "[plotMuonCandidateMipnessVsLengthByChannel] input is not an event list root file: " << list_path << "\n";
        return 1;
    }

    EventListIO el(list_path);
    ROOT::RDataFrame rdf = el.rdf();

    const std::vector<std::string> cols = rdf.GetColumnNames();
    const std::vector<std::string> required = {
        "sel_triggered_muon",
        "analysis_channels",
        "w_nominal",
        "track_shower_scores",
        "track_length",
        "track_distance_to_vertex",
        "pfp_generations",
        "track_dedx_T70_u",
        "track_dedx_T70_v",
        "track_dedx_T70_y",
        "track_dedx_Q50_u",
        "track_dedx_Q50_v",
        "track_dedx_Q50_y",
        "track_dedx_Q90_u",
        "track_dedx_Q90_v",
        "track_dedx_Q90_y"};

    std::vector<std::string> missing;
    for (const auto &c : required)
        if (!has_column(cols, c))
            missing.push_back(c);

    if (!missing.empty())
    {
        std::cerr << "[plotMuonCandidateMipnessVsLengthByChannel] missing required event-list columns:\n";
        for (const auto &c : missing)
            std::cerr << "  - " << c << "\n";
        return 2;
    }

    auto mask_ext = el.mask_for_ext();
    auto mask_mc = el.mask_for_mc_like();

    auto filter_by_mask = [](ROOT::RDF::RNode n, std::shared_ptr<const std::vector<char>> mask) {
        return n.Filter(
            [mask](int sid) {
                return sid >= 0
                       && sid < static_cast<int>(mask->size())
                       && (*mask)[static_cast<std::size_t>(sid)];
            },
            {"sample_id"});
    };

    ROOT::RDF::RNode base = rdf;
    ROOT::RDF::RNode node_mc = filter_by_mask(base, mask_mc)
                                   .Filter([mask_ext](int sid) {
                                       return !(sid >= 0
                                                && sid < static_cast<int>(mask_ext->size())
                                                && (*mask_ext)[static_cast<std::size_t>(sid)]);
                                   },
                                   {"sample_id"});

    node_mc = define_mu_cand_columns(node_mc, mu_score_cut, mu_len_cut_cm, mu_dist_cut_cm, mu_req_gen);
    node_mc = define_mipness_columns(node_mc, Mref_mev_per_cm, mip_tail_weight);

    const std::string sel_post_mu = "sel_triggered_muon && (mu_cand_idx >= 0)";
    const std::string sel_finite_mip = "(mu_cand_mipness_med == mu_cand_mipness_med)";
    const std::string sel_base = "(" + sel_post_mu + ") && (" + sel_finite_mip + ") && (" + extra_sel + ")";

    const std::string sel_ch_nc = channel_is_one_of_expr({14});
    const std::string sel_ch_nue_cc = channel_is_one_of_expr({17});
    const std::string sel_ch_numu_cc_inclusive = channel_is_one_of_expr({10, 11, 12, 13, 15, 16, 18});

    save_2d(node_mc,
            "mu_cand_length",
            "mu_cand_mipness_med",
            "w_nominal",
            sel_base + " && " + sel_ch_nc,
            "mipness2d_post_mucand_length_vs_mipnessmed_ch_nc",
            56, 0.0, 700.0,
            56, 0.0, 1.4,
            static_cast<double>(mu_len_cut_cm), mipness_threshold);

    save_2d(node_mc,
            "mu_cand_length",
            "mu_cand_mipness_med",
            "w_nominal",
            sel_base + " && " + sel_ch_nue_cc,
            "mipness2d_post_mucand_length_vs_mipnessmed_ch_nuecc",
            56, 0.0, 700.0,
            56, 0.0, 1.4,
            static_cast<double>(mu_len_cut_cm), mipness_threshold);

    save_2d(node_mc,
            "mu_cand_length",
            "mu_cand_mipness_med",
            "w_nominal",
            sel_base + " && " + sel_ch_numu_cc_inclusive,
            "mipness2d_post_mucand_length_vs_mipnessmed_ch_numucc_inclusive",
            56, 0.0, 700.0,
            56, 0.0, 1.4,
            static_cast<double>(mu_len_cut_cm), mipness_threshold);

    return 0;
}
