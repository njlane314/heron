/* -- C++ -- */
/**
 *  @file  ana/src/SystematicsBuilder.cc
 *
 *  @brief Implementation of systematic template building.
 */

#include "SystematicsBuilder.hh"

#include <TFile.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <ROOT/RDFHelpers.hxx>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "ColumnDerivationService.hh"
#include "RDataFrameFactory.hh"
#include "SampleIO.hh"
#include "TemplateIO.hh"

namespace nuxsec
{

namespace
{

inline float ushort_to_weight(unsigned short w)
{
    return static_cast<float>(w) * 0.001f;
}

inline double dot(const std::vector<double> &a, const std::vector<double> &b)
{
    double s = 0.0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i)
    {
        s += a[i] * b[i];
    }
    return s;
}

inline double rms(const std::vector<double> &x)
{
    if (x.empty())
    {
        return 0.0;
    }
    double s2 = 0.0;
    for (double v : x)
    {
        s2 += v * v;
    }
    return std::sqrt(s2 / static_cast<double>(x.size()));
}

} // namespace

bool SystematicsBuilder::is_varied_sample_kind(SampleIO::SampleKind k, const Options &opt)
{
    switch (k)
    {
    case SampleIO::SampleKind::kOverlay:
        return opt.include_overlay;
    case SampleIO::SampleKind::kDirt:
        return opt.include_dirt;
    case SampleIO::SampleKind::kStrangeness:
        return opt.include_strangeness;
    default:
        return false;
    }
}

std::unique_ptr<TH1D> SystematicsBuilder::read_nominal_hist(const std::string &root_path,
                                                            const std::string &sample_name,
                                                            const std::string &hist_name)
{
    std::unique_ptr<TFile> f(TFile::Open(root_path.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
        throw std::runtime_error("read_nominal_hist: cannot open " + root_path);
    }

    const std::string key = sample_name + "/hists/" + hist_name;
    TH1 *obj = dynamic_cast<TH1 *>(f->Get(key.c_str()));
    if (!obj)
    {
        throw std::runtime_error("read_nominal_hist: missing TH1 at " + key);
    }

    auto *h = dynamic_cast<TH1D *>(obj);
    if (!h)
    {
        throw std::runtime_error("read_nominal_hist: expected TH1D at " + key);
    }

    std::unique_ptr<TH1D> out(dynamic_cast<TH1D *>(h->Clone(("nom_" + hist_name).c_str())));
    out->SetDirectory(nullptr);
    return out;
}

void SystematicsBuilder::write_one_syst(const std::string &root_path,
                                        const std::string &sample_name,
                                        const std::string &syst_name,
                                        const std::string &variation,
                                        const std::vector<std::pair<std::string, const TH1 *>> &hists)
{
    nuxsec::TemplateIO::write_syst_histograms(root_path, sample_name, syst_name, variation, hists);
}

void SystematicsBuilder::clamp_non_negative(TH1D &h)
{
    const int nb = h.GetNbinsX();
    for (int i = 1; i <= nb; ++i)
    {
        if (h.GetBinContent(i) < 0.0)
        {
            h.SetBinContent(i, 0.0);
        }
    }
}

int SystematicsBuilder::detect_universe_count(const SampleIO::Sample &sample,
                                              const std::string &tree_name,
                                              const std::string &vec_branch)
{
    ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, tree_name);

    auto v = rdf.Range(0, 1).Take<std::vector<unsigned short>>(vec_branch).GetValue();
    if (v.empty())
    {
        return 0;
    }
    return static_cast<int>(v[0].size());
}

void SystematicsBuilder::build_unisim(const SampleIO::Sample &sample,
                                      const std::string &tree_name,
                                      const std::vector<TemplateSpec1D> &specs,
                                      const std::string &template_root_path,
                                      const UnisimSpec &uspec)
{
    ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, tree_name);

    nuxsec::ProcessorEntry proc = {};
    switch (sample.kind)
    {
    case SampleIO::SampleKind::kOverlay:
    case SampleIO::SampleKind::kDirt:
    case SampleIO::SampleKind::kStrangeness:
        proc.source = SourceKind::kMC;
        proc.pot_nom = sample.db_tortgt_pot_sum;
        proc.pot_eqv = sample.subrun_pot_sum;
        break;
    default:
        return;
    }

    ROOT::RDF::RNode node = nuxsec::ColumnDerivationService::instance().define(rdf, proc);

    const std::string wup = "__w_unisim_up_" + uspec.name;
    const std::string wdn = "__w_unisim_dn_" + uspec.name;

    node = node.Define(
        wup,
        [](float w, double k) { return static_cast<float>(static_cast<double>(w) * k); },
        {"w_template", uspec.up_branch});

    if (uspec.one_sided)
    {
        node = node.Define(wdn, [](float w) { return w; }, {"w_template"});
    }
    else
    {
        node = node.Define(
            wdn,
            [](float w, double k) { return static_cast<float>(static_cast<double>(w) * k); },
            {"w_template", uspec.dn_branch});
    }

    std::vector<ROOT::RDF::RResultPtr<TH1D>> hup;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hdn;
    std::vector<ROOT::RDF::RResultHandle> handles;
    hup.reserve(specs.size());
    hdn.reserve(specs.size());
    handles.reserve(2 * specs.size());

    for (const auto &hs : specs)
    {
        ROOT::RDF::RNode f = node;
        if (!hs.selection.empty())
        {
            f = f.Filter(hs.selection, hs.name);
        }

        ROOT::RDF::TH1DModel m_up((hs.name + "_up").c_str(), hs.title.c_str(),
                                 hs.nbins, hs.xmin, hs.xmax);
        ROOT::RDF::TH1DModel m_dn((hs.name + "_dn").c_str(), hs.title.c_str(),
                                 hs.nbins, hs.xmin, hs.xmax);

        auto hu = f.Histo1D(m_up, hs.variable, wup);
        auto hd = f.Histo1D(m_dn, hs.variable, wdn);

        hup.push_back(hu);
        hdn.push_back(hd);
        handles.emplace_back(hu);
        handles.emplace_back(hd);
    }

    ROOT::RDF::RunGraphs(handles);

    std::vector<std::pair<std::string, const TH1 *>> pos;
    std::vector<std::pair<std::string, const TH1 *>> neg;
    pos.reserve(specs.size());
    neg.reserve(specs.size());

    for (size_t i = 0; i < specs.size(); ++i)
    {
        const TH1D &hu = hup[i].GetValue();
        const TH1D &hd = hdn[i].GetValue();
        pos.emplace_back(specs[i].name, &hu);
        neg.emplace_back(specs[i].name, &hd);
    }

    write_one_syst(template_root_path, sample.sample_name, uspec.name, "pos", pos);
    write_one_syst(template_root_path, sample.sample_name, uspec.name, "neg", neg);

    nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, uspec.name, "type", "unisim");
    nuxsec::TemplateIO::write_syst_flag_meta(template_root_path,
                                                 uspec.name,
                                                 "log_normal",
                                                 uspec.log_normal ? "1" : "0");
    nuxsec::TemplateIO::write_syst_flag_meta(template_root_path,
                                                 uspec.name,
                                                 "floatable",
                                                 uspec.floatable ? "1" : "0");
}

void SystematicsBuilder::build_multisim_joint_eigenmodes(const std::vector<SampleIO::Sample> &samples,
                                                         const std::vector<std::string> &sample_names,
                                                         const std::string &tree_name,
                                                         const std::vector<TemplateSpec1D> &specs,
                                                         const std::string &template_root_path,
                                                         const MultisimSpec &mspec,
                                                         const Options &opt)
{
    if (samples.empty() || specs.empty())
    {
        return;
    }

    int U = mspec.max_universes;
    if (U <= 0)
    {
        U = detect_universe_count(samples.front(), tree_name, mspec.vec_branch);
    }
    if (U <= 0)
    {
        throw std::runtime_error("Multisim " + mspec.name + ": could not determine universe count");
    }

    struct Block
    {
        std::string sample;
        std::string hist;
        int nbins;
        int offset;
    };

    std::vector<Block> blocks;
    blocks.reserve(sample_names.size() * specs.size());

    int L = 0;
    for (const auto &sname : sample_names)
    {
        for (const auto &hs : specs)
        {
            blocks.push_back({sname, hs.name, hs.nbins, L});
            L += hs.nbins;
        }
    }

    std::vector<double> T0(static_cast<size_t>(L), 0.0);
    for (const auto &b : blocks)
    {
        std::unique_ptr<TH1D> h = read_nominal_hist(template_root_path, b.sample, b.hist);
        for (int i = 0; i < b.nbins; ++i)
        {
            T0[static_cast<size_t>(b.offset + i)] = h->GetBinContent(i + 1);
        }
    }

    const double T0dotT0 = dot(T0, T0);
    if (T0dotT0 <= 0.0)
    {
        throw std::runtime_error("Multisim " + mspec.name + ": nominal stack has zero norm");
    }

    struct SampleBook
    {
        std::string name;
        std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h;
        std::vector<ROOT::RDF::RResultHandle> handles;
    };

    std::vector<SampleBook> books;
    books.reserve(samples.size());

    for (size_t sidx = 0; sidx < samples.size(); ++sidx)
    {
        const SampleIO::Sample &sample = samples[sidx];

        ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, tree_name);

        nuxsec::ProcessorEntry proc = {};
        proc.source = SourceKind::kMC;
        proc.pot_nom = sample.db_tortgt_pot_sum;
        proc.pot_eqv = sample.subrun_pot_sum;

        ROOT::RDF::RNode node = nuxsec::ColumnDerivationService::instance().define(rdf, proc);

        SampleBook sb;
        sb.name = sample.sample_name;
        sb.h.resize(static_cast<size_t>(U));
        sb.handles.reserve(static_cast<size_t>(U) * specs.size());

        for (int u = 0; u < U; ++u)
        {
            const std::string wcol = "__w_" + mspec.name + "_u" + std::to_string(u);

            if (mspec.cv_branch.empty())
            {
                node = node.Define(
                    wcol,
                    [u](float w, const std::vector<unsigned short> &wv) {
                        const float uni = (u < static_cast<int>(wv.size()))
                                              ? ushort_to_weight(wv[static_cast<size_t>(u)])
                                              : 1.0f;
                        return w * uni;
                    },
                    {"w_template", mspec.vec_branch});
            }
            else
            {
                node = node.Define(
                    wcol,
                    [u](float w, const std::vector<unsigned short> &wv, float cv) {
                        const float uni = (u < static_cast<int>(wv.size()))
                                              ? ushort_to_weight(wv[static_cast<size_t>(u)])
                                              : 1.0f;
                        const float denom = (cv > 0.0f) ? cv : 1.0f;
                        return w * (uni / denom);
                    },
                    {"w_template", mspec.vec_branch, mspec.cv_branch});
            }

            sb.h[static_cast<size_t>(u)].reserve(specs.size());

            for (size_t hidx = 0; hidx < specs.size(); ++hidx)
            {
                const auto &hs = specs[hidx];

                ROOT::RDF::RNode f = node;
                if (!hs.selection.empty())
                {
                    f = f.Filter(hs.selection, hs.name);
                }

                const std::string model_name = hs.name + "_u" + std::to_string(u);
                ROOT::RDF::TH1DModel m(model_name.c_str(), hs.title.c_str(),
                                       hs.nbins, hs.xmin, hs.xmax);

                auto hist = f.Histo1D(m, hs.variable, wcol);
                sb.h[static_cast<size_t>(u)].push_back(hist);
                sb.handles.emplace_back(hist);
            }
        }

        ROOT::RDF::RunGraphs(sb.handles);
        books.push_back(std::move(sb));
    }

    TMatrixDSym Vshape(L);
    std::vector<double> thetas;
    thetas.reserve(static_cast<size_t>(U));

    for (int u = 0; u < U; ++u)
    {
        std::vector<double> Tu(static_cast<size_t>(L), 0.0);

        int cursor = 0;
        for (size_t sidx = 0; sidx < books.size(); ++sidx)
        {
            for (size_t hidx = 0; hidx < specs.size(); ++hidx)
            {
                const TH1D &h = books[sidx].h[static_cast<size_t>(u)][hidx].GetValue();
                const int nb = h.GetNbinsX();
                for (int i = 0; i < nb; ++i)
                {
                    Tu[static_cast<size_t>(cursor + i)] = h.GetBinContent(i + 1);
                }
                cursor += nb;
            }
        }

        double alpha = 0.0;
        if (mspec.split_rate_shape)
        {
            std::vector<double> d(static_cast<size_t>(L), 0.0);
            for (int i = 0; i < L; ++i)
            {
                d[static_cast<size_t>(i)] = Tu[static_cast<size_t>(i)] - T0[static_cast<size_t>(i)];
            }
            alpha = dot(T0, d) / T0dotT0;

            const double s = 1.0 + alpha;
            if (s > 0.0)
            {
                thetas.push_back(std::log(s));
            }
            else
            {
                thetas.push_back(0.0);
            }

            std::vector<double> R(static_cast<size_t>(L), 0.0);
            for (int i = 0; i < L; ++i)
            {
                R[static_cast<size_t>(i)] = Tu[static_cast<size_t>(i)] - s * T0[static_cast<size_t>(i)];
            }

            for (int i = 0; i < L; ++i)
            {
                const double Ri = R[static_cast<size_t>(i)];
                for (int j = 0; j <= i; ++j)
                {
                    Vshape(i, j) += Ri * R[static_cast<size_t>(j)];
                }
            }
        }
        else
        {
            std::vector<double> diff(static_cast<size_t>(L), 0.0);
            for (int i = 0; i < L; ++i)
            {
                diff[static_cast<size_t>(i)] = T0[static_cast<size_t>(i)] - Tu[static_cast<size_t>(i)];
            }

            for (int i = 0; i < L; ++i)
            {
                const double di = diff[static_cast<size_t>(i)];
                for (int j = 0; j <= i; ++j)
                {
                    Vshape(i, j) += di * diff[static_cast<size_t>(j)];
                }
            }
        }
    }

    const double invU = 1.0 / static_cast<double>(U);
    Vshape *= invU;

    if (mspec.split_rate_shape && mspec.rate_log_normal)
    {
        const double sigma_theta = rms(thetas);

        const std::string rate_name = mspec.name + "_rate";
        nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, rate_name, "type", "rate");
        nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, rate_name, "log_normal", "1");
        nuxsec::TemplateIO::write_syst_flag_meta(template_root_path,
                                                     rate_name,
                                                     "sigma_theta",
                                                     std::to_string(sigma_theta));

        const double up_scale = std::exp(+sigma_theta);
        const double dn_scale = std::exp(-sigma_theta);

        for (const auto &b : blocks)
        {
            std::unique_ptr<TH1D> hnom = read_nominal_hist(template_root_path, b.sample, b.hist);

            std::unique_ptr<TH1D> hup(dynamic_cast<TH1D *>(hnom->Clone((b.hist + "_rate_up").c_str())));
            std::unique_ptr<TH1D> hdn(dynamic_cast<TH1D *>(hnom->Clone((b.hist + "_rate_dn").c_str())));
            hup->SetDirectory(nullptr);
            hdn->SetDirectory(nullptr);

            hup->Scale(up_scale);
            hdn->Scale(dn_scale);

            write_one_syst(template_root_path, b.sample, rate_name, "pos", {{b.hist, hup.get()}});
            write_one_syst(template_root_path, b.sample, rate_name, "neg", {{b.hist, hdn.get()}});
        }
    }

    TMatrixDSymEigen eig(Vshape);
    const TVectorD evals = eig.GetEigenValues();
    const TMatrixD evecs = eig.GetEigenVectors();

    std::vector<int> idx(static_cast<size_t>(L));
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return evals[a] > evals[b]; });

    double total_var = 0.0;
    for (int i = 0; i < L; ++i)
    {
        if (evals[i] > 0.0)
        {
            total_var += evals[i];
        }
    }
    if (total_var <= 0.0)
    {
        return;
    }

    int nmodes = 0;
    double cum = 0.0;
    for (int k = 0; k < L; ++k)
    {
        const int i = idx[static_cast<size_t>(k)];
        const double lam = evals[i];
        if (lam <= 0.0)
        {
            continue;
        }
        cum += lam;
        ++nmodes;
        if (nmodes >= mspec.max_modes)
        {
            break;
        }
        if ((cum / total_var) >= mspec.keep_fraction)
        {
            break;
        }
    }

    nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, mspec.name, "type", "multisim_eigen");
    nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, mspec.name, "nuniv", std::to_string(U));
    nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, mspec.name, "nmodes", std::to_string(nmodes));

    for (int m = 0; m < nmodes; ++m)
    {
        const int col = idx[static_cast<size_t>(m)];
        const double lam = std::max(0.0, evals[col]);
        const double sdev = std::sqrt(lam);

        const std::string mode_name =
            mspec.name + "_mode" + (m < 10 ? "0" : "") + std::to_string(m);

        nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, mode_name, "type", "eigenmode");
        nuxsec::TemplateIO::write_syst_flag_meta(template_root_path, mode_name, "parent", mspec.name);

        for (const auto &b : blocks)
        {
            std::unique_ptr<TH1D> hnom = read_nominal_hist(template_root_path, b.sample, b.hist);

            std::unique_ptr<TH1D> hup(dynamic_cast<TH1D *>(hnom->Clone((b.hist + "_up").c_str())));
            std::unique_ptr<TH1D> hdn(dynamic_cast<TH1D *>(hnom->Clone((b.hist + "_dn").c_str())));
            hup->SetDirectory(nullptr);
            hdn->SetDirectory(nullptr);

            for (int i = 0; i < b.nbins; ++i)
            {
                const int gi = b.offset + i;
                const double delta = sdev * evecs(gi, col);

                const double nom = hnom->GetBinContent(i + 1);
                const double upv = nom + delta;
                const double dnv = nom - delta;

                hup->SetBinContent(i + 1, upv);
                hdn->SetBinContent(i + 1, dnv);
                hup->SetBinError(i + 1, hnom->GetBinError(i + 1));
                hdn->SetBinError(i + 1, hnom->GetBinError(i + 1));
            }

            if (opt.clamp_negative_bins)
            {
                clamp_non_negative(*hup);
                clamp_non_negative(*hdn);
            }

            write_one_syst(template_root_path, b.sample, mode_name, "pos", {{b.hist, hup.get()}});
            write_one_syst(template_root_path, b.sample, mode_name, "neg", {{b.hist, hdn.get()}});
        }
    }
}

void SystematicsBuilder::build_all(const std::vector<SampleListEntry> &entries,
                                   const std::string &tree_name,
                                   const std::vector<TemplateSpec1D> &fit_specs,
                                   const std::string &template_root_path,
                                   const SystematicsConfig &cfg,
                                   const Options &opt)
{
    if (opt.nthreads > 1)
    {
        ROOT::EnableImplicitMT(opt.nthreads);
    }

    std::vector<SampleIO::Sample> mc_samples;
    std::vector<std::string> mc_names;

    for (const auto &e : entries)
    {
        SampleIO::Sample s = SampleIO::read(e.output_path);
        if (!is_varied_sample_kind(s.kind, opt))
        {
            continue;
        }
        mc_names.push_back(s.sample_name);
        mc_samples.push_back(std::move(s));
    }

    for (const auto &us : cfg.unisim)
    {
        for (const auto &s : mc_samples)
        {
            build_unisim(s, tree_name, fit_specs, template_root_path, us);
        }
    }

    for (const auto &ms : cfg.multisim)
    {
        build_multisim_joint_eigenmodes(mc_samples, mc_names, tree_name, fit_specs, template_root_path, ms, opt);
    }
}

} // namespace nuxsec
