/* -- C++ -- */
/**
 *  @file  io/src/TemplateIO.cc
 *
 *  @brief Implementation for the NuXSec template container and ROOT IO helpers.
 */

#include "TemplateIO.hh"

#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TParameter.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace
{

TDirectory *get_or_make_dir(TDirectory *parent, const std::string &name)
{
    if (!parent)
    {
        return nullptr;
    }
    if (name.empty())
    {
        return parent;
    }
    if (auto *existing = parent->GetDirectory(name.c_str()))
    {
        return existing;
    }
    return parent->mkdir(name.c_str());
}

}

namespace nuxsec
{

const double TemplateIO::kIgnoreExposureScale = 1.0e3;

TemplateIO::TemplateIO()
    : TNamed()
{
    SetName("Default");
    m_mutable = false;
    m_nx = 1;
    m_ny = 1;
    m_true_bin_count = 1;
    m_bin_fraction.assign(1, 0.0);
    m_bin_stat.assign(1, 0.0);
}

TemplateIO::TemplateIO(const char *name, int n_x, double x_min, double x_max,
                       int n_y, double y_min, double y_max)
    : TNamed()
{
    SetName(name);
    m_mutable = true;

    m_nx = (n_x > 0 ? n_x : 1);
    m_ny = (n_y > 1 ? n_y : 1);
    m_true_bin_count = m_nx * m_ny;

    m_x_min = x_min;
    m_x_max = x_max;
    m_y_min = y_min;
    m_y_max = y_max;

    m_total_yield = 0.0;
    m_bin_fraction.assign(m_true_bin_count, 0.0);
    m_bin_stat.assign(m_true_bin_count, 0.0);
}

TemplateIO::TemplateIO(const TemplateIO &rhs)
    : TNamed(rhs)
{
    *this = rhs;
}

TemplateIO &TemplateIO::operator=(const TemplateIO &rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    TNamed::operator=(rhs);

    m_mutable = false;

    m_nx = rhs.m_nx;
    m_ny = rhs.m_ny;
    m_true_bin_count = rhs.m_true_bin_count;
    m_x_min = rhs.m_x_min;
    m_x_max = rhs.m_x_max;
    m_y_min = rhs.m_y_min;
    m_y_max = rhs.m_y_max;

    m_total_yield = rhs.m_total_yield;
    m_bin_fraction = rhs.m_bin_fraction;
    m_bin_stat = rhs.m_bin_stat;

    m_syst_names = rhs.m_syst_names;
    m_syst = rhs.m_syst;
    m_float_flag = rhs.m_float_flag;
    m_log_normal_flag = rhs.m_log_normal_flag;

    m_syst_index_outer = rhs.m_syst_index_outer;
    m_syst_index_inner = rhs.m_syst_index_inner;

    m_linearised = rhs.m_linearised;
    m_lin_bin_fraction = rhs.m_lin_bin_fraction;
    m_lin_bin_stat = rhs.m_lin_bin_stat;
    m_efficiency_sums = rhs.m_efficiency_sums;

    return *this;
}

TemplateIO::~TemplateIO() = default;

double TemplateIO::BinFraction(int ix, int iy) const
{
    if (ix < 0 || ix >= m_nx)
    {
        return 0.0;
    }
    if (m_ny > 1)
    {
        if (iy < 0 || iy >= m_ny)
        {
            return 0.0;
        }
        return m_bin_fraction[Index(ix, iy)];
    }
    return m_bin_fraction[ix];
}

double TemplateIO::BinYield(int ix, int iy) const
{
    return BinFraction(ix, iy) * m_total_yield;
}

double TemplateIO::BinStatError(int ix, int iy) const
{
    if (ix < 0 || ix >= m_nx)
    {
        return 0.0;
    }
    if (m_ny > 1)
    {
        if (iy < 0 || iy >= m_ny)
        {
            return 0.0;
        }
        return m_bin_stat[Index(ix, iy)];
    }
    return m_bin_stat[ix];
}

void TemplateIO::SetTotalYield(double yield)
{
    if (m_mutable)
    {
        m_total_yield = yield;
    }
}

void TemplateIO::SetBinFraction(double fraction, int ix, int iy)
{
    if (!m_mutable)
    {
        return;
    }
    if (ix < 0 || ix >= m_nx)
    {
        return;
    }
    if (m_ny > 1)
    {
        if (iy < 0 || iy >= m_ny)
        {
            return;
        }
        m_bin_fraction[Index(ix, iy)] = fraction;
        return;
    }
    m_bin_fraction[ix] = fraction;
}

void TemplateIO::SetBinStatError(double error, int ix, int iy)
{
    if (!m_mutable)
    {
        return;
    }
    if (ix < 0 || ix >= m_nx)
    {
        return;
    }
    if (m_ny > 1)
    {
        if (iy < 0 || iy >= m_ny)
        {
            return;
        }
        m_bin_stat[Index(ix, iy)] = error;
        return;
    }
    m_bin_stat[ix] = error;
}

bool TemplateIO::CheckStats(const TH1 *hist, int verbose) const
{
    if (!hist)
    {
        return false;
    }
    if (!hist->GetSumw2N() && verbose > 0)
    {
        ::Warning("TemplateIO::CheckStats",
                  "Input hist %s has no Sumw2; stat errors may be wrong.", hist->GetName());
    }
    return true;
}

bool TemplateIO::FillFromTH1(const TH1 *hist, int rebin)
{
    if (!hist)
    {
        return false;
    }
    if (m_ny > 1)
    {
        return false;
    }

    std::unique_ptr<TH1> tmp(static_cast<TH1 *>(hist->Clone("nuxsec_tmp_1d")));
    if (!tmp)
    {
        return false;
    }

    if (tmp->GetNbinsX() != m_nx)
    {
        if (rebin > 0 && (tmp->GetNbinsX() / rebin) == m_nx)
        {
            tmp->Rebin(rebin);
        }
        else
        {
            return false;
        }
    }

    if (!CheckStats(tmp.get(), 1))
    {
        return false;
    }

    const double sum = tmp->Integral(1, m_nx);
    m_total_yield = sum;

    if (sum <= 0.0)
    {
        std::fill(m_bin_fraction.begin(), m_bin_fraction.end(), 0.0);
        std::fill(m_bin_stat.begin(), m_bin_stat.end(), 0.0);
        return true;
    }

    for (int i = 0; i < m_nx; ++i)
    {
        const double content = tmp->GetBinContent(i + 1);
        const double error = tmp->GetBinError(i + 1);

        if (content < 0.0)
        {
            return false;
        }

        const double frac = content / sum;
        const double ef = error / sum;

        m_bin_fraction[i] = frac;
        m_bin_stat[i] = ef;
    }

    return true;
}

bool TemplateIO::FillFromTH2(const TH2 *hist, int rebin_x, int rebin_y)
{
    if (!hist)
    {
        return false;
    }
    if (m_ny < 2)
    {
        return false;
    }

    std::unique_ptr<TH2> tmp(static_cast<TH2 *>(hist->Clone("nuxsec_tmp_2d")));
    if (!tmp)
    {
        return false;
    }

    if (tmp->GetNbinsX() != m_nx || tmp->GetNbinsY() != m_ny)
    {
        if (rebin_x > 0 && rebin_y > 0 && (tmp->GetNbinsX() / rebin_x) == m_nx
            && (tmp->GetNbinsY() / rebin_y) == m_ny)
        {
            tmp->Rebin2D(rebin_x, rebin_y);
        }
        else
        {
            return false;
        }
    }

    if (!CheckStats(tmp.get(), 1))
    {
        return false;
    }

    const double sum = tmp->Integral(1, m_nx, 1, m_ny);
    m_total_yield = sum;

    if (sum <= 0.0)
    {
        std::fill(m_bin_fraction.begin(), m_bin_fraction.end(), 0.0);
        std::fill(m_bin_stat.begin(), m_bin_stat.end(), 0.0);
        return true;
    }

    for (int iy = 0; iy < m_ny; ++iy)
    {
        for (int ix = 0; ix < m_nx; ++ix)
        {
            const int ibin = Index(ix, iy);
            const double content = tmp->GetBinContent(ix + 1, iy + 1);
            const double error = tmp->GetBinError(ix + 1, iy + 1);

            if (content < 0.0)
            {
                return false;
            }

            const double frac = content / sum;
            const double ef = error / sum;

            m_bin_fraction[ibin] = frac;
            m_bin_stat[ibin] = ef;
        }
    }

    return true;
}

bool TemplateIO::HasSystematic(const std::string &name) const
{
    return (SystIndex(name) >= 0);
}

int TemplateIO::SystIndex(const std::string &name) const
{
    if (name.empty())
    {
        return -1;
    }

    for (size_t i = 0; i < m_syst_names.size(); ++i)
    {
        if (m_syst_names[i] == name)
        {
            return static_cast<int>(i);
        }
    }

    return -1;
}

void TemplateIO::AddSystematicFrac(const char *name, const TH1 *pos_frac, const TH1 *neg_frac)
{
    if (!m_mutable)
    {
        return;
    }
    if (!pos_frac || !neg_frac)
    {
        return;
    }
    if (!name || !name[0])
    {
        return;
    }
    if (pos_frac->GetNbinsX() != m_nx || neg_frac->GetNbinsX() != m_nx)
    {
        return;
    }
    if (m_ny > 1)
    {
        return;
    }

    const std::string n(name);
    if (HasSystematic(n))
    {
        return;
    }

    m_syst_names.push_back(n);
    m_syst.push_back(std::vector<NuXSecSystBin>(m_true_bin_count));
    m_float_flag.push_back(0);
    m_log_normal_flag.push_back(0);

    for (int i = 0; i < m_nx; ++i)
    {
        NuXSecSystBin &sys_bin = m_syst.back()[i];
        const double pos = pos_frac->GetBinContent(i + 1);
        const double neg = neg_frac->GetBinContent(i + 1);

        sys_bin.sigma_pos = pos;
        sys_bin.sigma_neg = neg;
        sys_bin.asym = (pos != neg);
    }
}

void TemplateIO::AddSystematicFrac2D(const char *name, const TH2 *pos_frac, const TH2 *neg_frac)
{
    if (!m_mutable)
    {
        return;
    }
    if (!pos_frac || !neg_frac)
    {
        return;
    }
    if (!name || !name[0])
    {
        return;
    }
    if (pos_frac->GetNbinsX() != m_nx || pos_frac->GetNbinsY() != m_ny)
    {
        return;
    }
    if (neg_frac->GetNbinsX() != m_nx || neg_frac->GetNbinsY() != m_ny)
    {
        return;
    }

    const std::string n(name);
    if (HasSystematic(n))
    {
        return;
    }

    m_syst_names.push_back(n);
    m_syst.push_back(std::vector<NuXSecSystBin>(m_true_bin_count));
    m_float_flag.push_back(0);
    m_log_normal_flag.push_back(0);

    for (int iy = 0; iy < m_ny; ++iy)
    {
        for (int ix = 0; ix < m_nx; ++ix)
        {
            NuXSecSystBin &sys_bin = m_syst.back()[Index(ix, iy)];
            const double pos = pos_frac->GetBinContent(ix + 1, iy + 1);
            const double neg = neg_frac->GetBinContent(ix + 1, iy + 1);

            sys_bin.sigma_pos = pos;
            sys_bin.sigma_neg = neg;
            sys_bin.asym = (pos != neg);
        }
    }
}

bool TemplateIO::GetFloatFlag(const std::string &name) const
{
    const int index = SystIndex(name);
    if (index < 0)
    {
        return false;
    }
    return (m_float_flag[index] != 0);
}

bool TemplateIO::GetLogNormalFlag(const std::string &name) const
{
    const int index = SystIndex(name);
    if (index < 0)
    {
        return false;
    }
    return (m_log_normal_flag[index] != 0);
}

void TemplateIO::SetFloatFlag(const std::string &name, bool on)
{
    if (!m_mutable)
    {
        return;
    }
    const int index = SystIndex(name);
    if (index < 0)
    {
        return;
    }
    m_float_flag[index] = (on ? 1 : 0);
}

void TemplateIO::SetLogNormalFlag(const std::string &name, bool on)
{
    if (!m_mutable)
    {
        return;
    }
    const int index = SystIndex(name);
    if (index < 0)
    {
        return;
    }
    m_log_normal_flag[index] = (on ? 1 : 0);
}

void TemplateIO::Linearise(const std::vector<std::string> &global_syst_order)
{
    if (!m_mutable)
    {
        return;
    }
    if (global_syst_order.size() != m_syst_names.size())
    {
        return;
    }

    if (!m_syst_index_outer.empty() || !m_syst_index_inner.empty())
    {
        return;
    }

    const int n_syst = NSyst();
    m_syst_index_outer.assign(n_syst, -1);
    m_syst_index_inner.assign(n_syst, -1);

    std::vector<int> indices(n_syst, -1);
    for (int i = 0; i < n_syst; ++i)
    {
        const std::string &name = m_syst_names[i];
        for (int j = 0; j < n_syst; ++j)
        {
            if (name == global_syst_order[j])
            {
                indices[i] = j;
                break;
            }
        }
        if (indices[i] < 0)
        {
            return;
        }
    }

    std::vector<int> sorted_indices = indices;
    std::sort(sorted_indices.begin(), sorted_indices.end());
    if (sorted_indices.front() != 0 || sorted_indices.back() != n_syst - 1)
    {
        return;
    }

    for (int i = 0; i < n_syst; ++i)
    {
        const int outer = indices[i];
        m_syst_index_outer[i] = outer;
        m_syst_index_inner[outer] = i;
    }

    m_lin_bin_fraction = m_bin_fraction;
    m_lin_bin_stat = m_bin_stat;
    m_efficiency_sums.assign(m_true_bin_count, 0.0);

    for (int is = 0; is < n_syst; ++is)
    {
        const int index = m_syst_index_inner[is];
        if (index < 0)
        {
            return;
        }

        const std::vector<NuXSecSystBin> &syst = m_syst[index];
        for (int ib = 0; ib < m_true_bin_count; ++ib)
        {
            const NuXSecSystBin &sys_bin = syst[ib];
            const double delta_eff = 0.5 * (std::abs(sys_bin.sigma_pos) + std::abs(sys_bin.sigma_neg));
            m_efficiency_sums[ib] += delta_eff;
        }
    }

    m_linearised = true;
}

double TemplateIO::BinFractionVaried(int ix, int iy, const double *pulls) const
{
    if (!pulls)
    {
        return 0.0;
    }
    if (ix < 0 || ix >= m_nx)
    {
        return 0.0;
    }
    if (m_ny > 1)
    {
        if (iy < 0 || iy >= m_ny)
        {
            return 0.0;
        }
    }

    const int ibin = Index(ix, (m_ny > 1) ? iy : 0);

    if (!m_linearised || m_syst_index_outer.empty() || m_syst_index_inner.empty())
    {
        double out = m_bin_fraction[ibin];
        for (size_t is = 0; is < m_syst.size(); ++is)
        {
            const NuXSecSystBin &sys_bin = m_syst[is][ibin];
            const double delta = sys_bin.asym ? AsymDelta(pulls[is], sys_bin.sigma_pos, sys_bin.sigma_neg)
                                              : (pulls[is] * sys_bin.sigma_pos);
            out += (m_bin_fraction[ibin] * delta);
        }
        return out;
    }

    double sum = m_efficiency_sums[ibin];
    if (sum <= 0.0)
    {
        return m_lin_bin_fraction[ibin];
    }

    double var = 0.0;
    for (size_t is = 0; is < m_syst.size(); ++is)
    {
        const int outer = m_syst_index_outer[is];
        const NuXSecSystBin &sys_bin = m_syst[is][ibin];
        const double delta = sys_bin.asym ? AsymDelta(pulls[outer], sys_bin.sigma_pos, sys_bin.sigma_neg)
                                          : (pulls[outer] * sys_bin.sigma_pos);
        var += delta;
    }

    return m_lin_bin_fraction[ibin] * (1.0 + var);
}

void TemplateIO::PrepareExclusionSums()
{
    for (size_t is = 0; is < m_syst.size(); ++is)
    {
        for (int ib = 0; ib < m_true_bin_count; ++ib)
        {
            NuXSecSystBin &sys_bin = m_syst[is][ib];
            sys_bin.exclusion_sum = 0.0;
        }
    }

    for (int ib = 0; ib < m_true_bin_count; ++ib)
    {
        double sum = 0.0;
        for (size_t is = 0; is < m_syst.size(); ++is)
        {
            NuXSecSystBin &sys_bin = m_syst[is][ib];
            const double pos = std::abs(sys_bin.sigma_pos);
            const double neg = std::abs(sys_bin.sigma_neg);
            sum += (pos > neg) ? pos : neg;
        }
        for (size_t is = 0; is < m_syst.size(); ++is)
        {
            NuXSecSystBin &sys_bin = m_syst[is][ib];
            sys_bin.exclusion_sum = sum - std::max(std::abs(sys_bin.sigma_pos), std::abs(sys_bin.sigma_neg));
        }
    }
}

TH1 *TemplateIO::MakeTH1(const std::string &title) const
{
    if (m_ny > 1)
    {
        return nullptr;
    }

    const std::string hname = std::string(GetName()) + "_th1";
    auto *hist = new TH1D(hname.c_str(), title.c_str(), m_nx, m_x_min, m_x_max);
    hist->Sumw2();

    for (int i = 0; i < m_nx; ++i)
    {
        const double content = m_total_yield * m_bin_fraction[i];
        const double error = m_total_yield * m_bin_stat[i];
        hist->SetBinContent(i + 1, content);
        hist->SetBinError(i + 1, error);
    }

    return hist;
}

TH2 *TemplateIO::MakeTH2(const std::string &title) const
{
    if (m_ny < 2)
    {
        return nullptr;
    }

    const std::string hname = std::string(GetName()) + "_th2";
    auto *hist = new TH2D(hname.c_str(), title.c_str(), m_nx, m_x_min, m_x_max,
                          m_ny, m_y_min, m_y_max);
    hist->Sumw2();

    for (int iy = 0; iy < m_ny; ++iy)
    {
        for (int ix = 0; ix < m_nx; ++ix)
        {
            const int ibin = Index(ix, iy);
            const double content = m_total_yield * m_bin_fraction[ibin];
            const double error = m_total_yield * m_bin_stat[ibin];
            hist->SetBinContent(ix + 1, iy + 1, content);
            hist->SetBinError(ix + 1, iy + 1, error);
        }
    }

    return hist;
}

void TemplateIO::write_histograms(const std::string &root_path,
                                  const std::string &sample_name,
                                  const std::vector<std::pair<std::string, const TH1 *>> &hists,
                                  const TemplateWriteOptions &opt)
{
    if (hists.empty())
    {
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateIO::write_histograms: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *hdir = get_or_make_dir(sdir, "hists");
    if (!hdir)
    {
        throw std::runtime_error("TemplateIO::write_histograms: failed to create directories");
    }

    hdir->cd();
    for (const auto &entry : hists)
    {
        const std::string &hname = entry.first;
        const TH1 *hptr = entry.second;
        if (!hptr)
        {
            continue;
        }

        std::unique_ptr<TH1> clone(dynamic_cast<TH1 *>(hptr->Clone(hname.c_str())));
        if (!clone)
        {
            throw std::runtime_error("TemplateIO::write_histograms: failed to clone TH1 " + hname);
        }

        clone->SetDirectory(hdir);
        const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
        clone->Write(hname.c_str(), wopt);
    }

    fout->Write();
    fout->Close();
}

void TemplateIO::write_string_meta(const std::string &root_path,
                                   const std::string &sample_name,
                                   const std::string &key,
                                   const std::string &value,
                                   const TemplateWriteOptions &opt)
{
    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateIO::write_string_meta: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *mdir = get_or_make_dir(sdir, "meta");
    if (!mdir)
    {
        throw std::runtime_error("TemplateIO::write_string_meta: failed to create meta directory");
    }

    mdir->cd();
    TNamed named(key.c_str(), value.c_str());
    const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
    named.Write(key.c_str(), wopt);

    fout->Write();
    fout->Close();
}

void TemplateIO::write_double_meta(const std::string &root_path,
                                   const std::string &sample_name,
                                   const std::string &key,
                                   double value,
                                   const TemplateWriteOptions &opt)
{
    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateIO::write_double_meta: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *mdir = get_or_make_dir(sdir, "meta");
    if (!mdir)
    {
        throw std::runtime_error("TemplateIO::write_double_meta: failed to create meta directory");
    }

    mdir->cd();
    TParameter<double> par(key.c_str(), value);
    const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
    par.Write(key.c_str(), wopt);

    fout->Write();
    fout->Close();
}

void TemplateIO::write_syst_histograms(const std::string &root_path,
                                       const std::string &sample_name,
                                       const std::string &syst_name,
                                       const std::string &variation,
                                       const std::vector<std::pair<std::string, const TH1 *>> &hists,
                                       const TemplateWriteOptions &opt)
{
    if (hists.empty())
    {
        return;
    }
    if (syst_name.empty() || variation.empty())
    {
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateIO::write_syst_histograms: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *sysd = get_or_make_dir(sdir, "systs");
    TDirectory *thiss = get_or_make_dir(sysd, syst_name);
    TDirectory *vard = get_or_make_dir(thiss, variation);
    TDirectory *hdir = get_or_make_dir(vard, "hists");
    if (!hdir)
    {
        throw std::runtime_error("TemplateIO::write_syst_histograms: failed to create directories");
    }

    hdir->cd();
    for (const auto &entry : hists)
    {
        const std::string &hname = entry.first;
        const TH1 *hptr = entry.second;
        if (!hptr)
        {
            continue;
        }

        std::unique_ptr<TH1> clone(dynamic_cast<TH1 *>(hptr->Clone(hname.c_str())));
        if (!clone)
        {
            throw std::runtime_error("TemplateIO::write_syst_histograms: failed to clone TH1 " + hname);
        }

        clone->SetDirectory(hdir);
        const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
        clone->Write(hname.c_str(), wopt);
    }

    fout->Write();
    fout->Close();
}

void TemplateIO::write_syst_flag_meta(const std::string &root_path,
                                      const std::string &syst_name,
                                      const std::string &key,
                                      const std::string &value,
                                      const TemplateWriteOptions &opt)
{
    if (syst_name.empty())
    {
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateIO::write_syst_flag_meta: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *gdir = get_or_make_dir(top, "__global__");
    TDirectory *mdir = get_or_make_dir(gdir, "meta");
    TDirectory *sdir = get_or_make_dir(mdir, "systs");
    TDirectory *thiss = get_or_make_dir(sdir, syst_name);
    if (!thiss)
    {
        throw std::runtime_error("TemplateIO::write_syst_flag_meta: failed to create syst meta directories");
    }

    thiss->cd();
    TNamed named(key.c_str(), value.c_str());
    const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
    named.Write(key.c_str(), wopt);

    fout->Write();
    fout->Close();
}

} // namespace nuxsec
