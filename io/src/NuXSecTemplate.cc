/* -- C++ -- */
/**
 *  @file  io/src/NuXSecTemplate.cc
 *
 *  @brief Implementation for the NuXSec template container.
 */

#include "NuXSecTemplate.hh"

#include <TError.h>
#include <TGenericClassInfo.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace nuxsec
{

ClassImp(NuXSecTemplate);

const double NuXSecTemplate::kIgnoreExposureScale = 1.0e3;

NuXSecTemplate::NuXSecTemplate()
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

NuXSecTemplate::NuXSecTemplate(const char *name, int n_x, double x_min, double x_max,
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

NuXSecTemplate::NuXSecTemplate(const NuXSecTemplate &rhs)
    : TNamed(rhs)
{
    *this = rhs;
}

NuXSecTemplate &NuXSecTemplate::operator=(const NuXSecTemplate &rhs)
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

NuXSecTemplate::~NuXSecTemplate() = default;

double NuXSecTemplate::BinFraction(int ix, int iy) const
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

double NuXSecTemplate::BinYield(int ix, int iy) const
{
    return BinFraction(ix, iy) * m_total_yield;
}

double NuXSecTemplate::BinStatError(int ix, int iy) const
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

void NuXSecTemplate::SetTotalYield(double yield)
{
    if (m_mutable)
    {
        m_total_yield = yield;
    }
}

void NuXSecTemplate::SetBinFraction(double fraction, int ix, int iy)
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

void NuXSecTemplate::SetBinStatError(double error, int ix, int iy)
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

bool NuXSecTemplate::CheckStats(const TH1 *hist, int verbose) const
{
    if (!hist)
    {
        return false;
    }
    if (!hist->GetSumw2N() && verbose > 0)
    {
        ::Warning("NuXSecTemplate::CheckStats",
                  "Input hist %s has no Sumw2; stat errors may be wrong.", hist->GetName());
    }
    return true;
}

bool NuXSecTemplate::FillFromTH1(const TH1 *hist, int rebin)
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
        m_bin_fraction[i] = content / sum;
        m_bin_stat[i] = error;
    }
    return true;
}

bool NuXSecTemplate::FillFromTH2(const TH2 *hist, int rebin_x, int rebin_y)
{
    if (!hist)
    {
        return false;
    }
    if (m_ny <= 1)
    {
        return false;
    }

    std::unique_ptr<TH2> tmp(static_cast<TH2 *>(hist->Clone("nuxsec_tmp_2d")));
    if (!tmp)
    {
        return false;
    }

    if (tmp->GetNbinsX() != m_nx)
    {
        if (rebin_x > 0 && (tmp->GetNbinsX() / rebin_x) == m_nx)
        {
            tmp->RebinX(rebin_x);
        }
        else
        {
            return false;
        }
    }
    if (tmp->GetNbinsY() != m_ny)
    {
        if (rebin_y > 0 && (tmp->GetNbinsY() / rebin_y) == m_ny)
        {
            tmp->RebinY(rebin_y);
        }
        else
        {
            return false;
        }
    }

    const double sum = tmp->Integral(1, m_nx, 1, m_ny);
    m_total_yield = sum;

    if (sum <= 0.0)
    {
        std::fill(m_bin_fraction.begin(), m_bin_fraction.end(), 0.0);
        std::fill(m_bin_stat.begin(), m_bin_stat.end(), 0.0);
        return true;
    }

    for (int ix = 0; ix < m_nx; ++ix)
    {
        for (int iy = 0; iy < m_ny; ++iy)
        {
            const double content = tmp->GetBinContent(ix + 1, iy + 1);
            const double error = tmp->GetBinError(ix + 1, iy + 1);
            if (content < 0.0)
            {
                return false;
            }
            const int k = Index(ix, iy);
            m_bin_fraction[k] = content / sum;
            m_bin_stat[k] = error;
        }
    }
    return true;
}

bool NuXSecTemplate::HasSystematic(const std::string &name) const
{
    return (SystIndex(name) >= 0);
}

int NuXSecTemplate::SystIndex(const std::string &name) const
{
    for (int i = 0; i < static_cast<int>(m_syst_names.size()); ++i)
    {
        if (m_syst_names[i] == name)
        {
            return i;
        }
    }
    return -1;
}

void NuXSecTemplate::AddSystematicFrac(const char *name, const TH1 *pos_frac, const TH1 *neg_frac)
{
    if (!name || !pos_frac || !neg_frac)
    {
        return;
    }
    if (m_ny > 1)
    {
        return;
    }
    if (HasSystematic(name))
    {
        return;
    }
    if (pos_frac->GetNbinsX() != m_nx || neg_frac->GetNbinsX() != m_nx)
    {
        return;
    }

    const int s = static_cast<int>(m_syst_names.size());
    m_syst_names.push_back(name);
    m_syst.emplace_back(m_true_bin_count);
    m_float_flag.push_back(0);
    m_log_normal_flag.push_back(0);

    for (int i = 0; i < m_nx; ++i)
    {
        NuXSecSystBin syst_bin;
        syst_bin.sigma_pos = pos_frac->GetBinContent(i + 1);
        syst_bin.sigma_neg = neg_frac->GetBinContent(i + 1);
        syst_bin.asym = (syst_bin.sigma_pos != syst_bin.sigma_neg);
        m_syst[s][i] = syst_bin;
    }
}

void NuXSecTemplate::AddSystematicFrac2D(const char *name, const TH2 *pos_frac, const TH2 *neg_frac)
{
    if (!name || !pos_frac || !neg_frac)
    {
        return;
    }
    if (m_ny <= 1)
    {
        return;
    }
    if (HasSystematic(name))
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

    const int s = static_cast<int>(m_syst_names.size());
    m_syst_names.push_back(name);
    m_syst.emplace_back(m_true_bin_count);
    m_float_flag.push_back(0);
    m_log_normal_flag.push_back(0);

    for (int ix = 0; ix < m_nx; ++ix)
    {
        for (int iy = 0; iy < m_ny; ++iy)
        {
            NuXSecSystBin syst_bin;
            syst_bin.sigma_pos = pos_frac->GetBinContent(ix + 1, iy + 1);
            syst_bin.sigma_neg = neg_frac->GetBinContent(ix + 1, iy + 1);
            syst_bin.asym = (syst_bin.sigma_pos != syst_bin.sigma_neg);
            m_syst[s][Index(ix, iy)] = syst_bin;
        }
    }
}

bool NuXSecTemplate::GetFloatFlag(const std::string &name) const
{
    const int index = SystIndex(name);
    if (index < 0)
    {
        return false;
    }
    return (m_float_flag[index] & 0x1);
}

bool NuXSecTemplate::GetLogNormalFlag(const std::string &name) const
{
    const int index = SystIndex(name);
    if (index < 0)
    {
        return false;
    }
    return (m_log_normal_flag[index] & 0x1);
}

void NuXSecTemplate::SetFloatFlag(const std::string &name, bool on)
{
    const int index = SystIndex(name);
    if (index < 0)
    {
        return;
    }
    if (on)
    {
        m_float_flag[index] |= 0x1;
    }
    else
    {
        m_float_flag[index] &= ~0x1;
    }
}

void NuXSecTemplate::SetLogNormalFlag(const std::string &name, bool on)
{
    const int index = SystIndex(name);
    if (index < 0)
    {
        return;
    }
    if (on)
    {
        m_log_normal_flag[index] |= 0x1;
    }
    else
    {
        m_log_normal_flag[index] &= ~0x1;
    }
}

void NuXSecTemplate::Linearise(const std::vector<std::string> &global_syst_order)
{
    m_lin_bin_fraction = m_bin_fraction;
    m_lin_bin_stat = m_bin_stat;

    m_syst_index_outer.assign(m_syst_names.size(), -1);
    m_syst_index_inner.assign(global_syst_order.size(), -1);

    for (int i = 0; i < static_cast<int>(m_syst_names.size()); ++i)
    {
        for (int g = 0; g < static_cast<int>(global_syst_order.size()); ++g)
        {
            if (m_syst_names[i] == global_syst_order[g])
            {
                m_syst_index_outer[i] = g;
                m_syst_index_inner[g] = i;
                break;
            }
        }
    }

    m_efficiency_sums.assign(m_true_bin_count, 0.0);
    m_linearised = true;
}

double NuXSecTemplate::BinFractionVaried(int ix, int iy, const double *pulls) const
{
    if (!m_linearised)
    {
        return 0.0;
    }
    if (ix < 0 || ix >= m_nx)
    {
        return 0.0;
    }
    if (m_ny > 1 && (iy < 0 || iy >= m_ny))
    {
        return 0.0;
    }

    const int k = (m_ny > 1) ? Index(ix, iy) : ix;
    const double nominal = m_lin_bin_fraction[k];
    if (nominal <= 0.0)
    {
        return 0.0;
    }

    double delta = 0.0;
    for (int s = 0; s < static_cast<int>(m_syst.size()); ++s)
    {
        const int g = m_syst_index_outer[s];
        if (g < 0)
        {
            continue;
        }
        const double pull = pulls[g];

        const NuXSecSystBin &syst_bin = m_syst[s][k];
        double local_delta = AsymDelta(pull, syst_bin.sigma_pos, syst_bin.sigma_neg);

        if (m_log_normal_flag[s] & 0x1)
        {
            local_delta = std::exp(local_delta) - 1.0;
        }

        delta += local_delta;
    }

    const double varied = nominal * (1.0 + delta);
    return (varied > 0.0) ? varied : 0.0;
}

void NuXSecTemplate::PrepareExclusionSums()
{
}

TH1 *NuXSecTemplate::MakeTH1(const std::string &title) const
{
    if (m_ny > 1)
    {
        return nullptr;
    }

    const std::string name = std::string(GetName()) + (title.empty() ? "" : ("-" + title));
    auto *hist = new TH1D(name.c_str(), name.c_str(), m_nx, m_x_min, m_x_max);
    hist->Sumw2();

    for (int i = 0; i < m_nx; ++i)
    {
        hist->SetBinContent(i + 1, BinYield(i));
        hist->SetBinError(i + 1, BinStatError(i));
    }
    return hist;
}

TH2 *NuXSecTemplate::MakeTH2(const std::string &title) const
{
    if (m_ny <= 1)
    {
        return nullptr;
    }

    const std::string name = title.empty() ? std::string(GetName()) : title;
    auto *hist = new TH2D(name.c_str(), name.c_str(), m_nx, m_x_min, m_x_max, m_ny, m_y_min, m_y_max);
    hist->Sumw2();

    for (int ix = 0; ix < m_nx; ++ix)
    {
        for (int iy = 0; iy < m_ny; ++iy)
        {
            hist->SetBinContent(ix + 1, iy + 1, BinYield(ix, iy));
            hist->SetBinError(ix + 1, iy + 1, BinStatError(ix, iy));
        }
    }
    return hist;
}

} // namespace nuxsec
