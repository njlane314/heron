/* -- C++ -- */
/**
 *  @file  io/include/NuXSecTemplate.hh
 *
 *  @brief Collie-style template container for NuXSec fits.
 */

#ifndef NUXSEC_IO_NUXSEC_TEMPLATE_H
#define NUXSEC_IO_NUXSEC_TEMPLATE_H

#include <TNamed.h>

#include <string>
#include <vector>

class TH1;
class TH2;

namespace nuxsec
{

struct NuXSecSystBin
{
    double sigma_pos = 0.0;
    double sigma_neg = 0.0;
    bool asym = false;

    double base_delta_efficiency = 0.0;
    double exclusion_sum = 0.0;

};

/** \brief CollieDistribution-like template representation for NuXSec. */
class NuXSecTemplate : public TNamed
{
  public:
    static const double kIgnoreExposureScale;

    NuXSecTemplate();
    NuXSecTemplate(const char *name, int n_x, double x_min, double x_max,
                   int n_y = 1, double y_min = 0, double y_max = 1);

    NuXSecTemplate(const NuXSecTemplate &rhs);
    NuXSecTemplate &operator=(const NuXSecTemplate &rhs);
    ~NuXSecTemplate() override;

    int NX() const { return m_nx; }
    int NY() const { return m_ny; }
    int NBins() const { return m_true_bin_count; }

    double XMin() const { return m_x_min; }
    double XMax() const { return m_x_max; }
    double YMin() const { return m_y_min; }
    double YMax() const { return m_y_max; }

    double BinFraction(int ix, int iy = -1) const;
    double BinYield(int ix, int iy = -1) const;
    double BinStatError(int ix, int iy = -1) const;

    double TotalYield() const { return m_total_yield; }
    void SetTotalYield(double yield);

    void SetBinFraction(double fraction, int ix, int iy = -1);
    void SetBinStatError(double error, int ix, int iy = -1);

    bool FillFromTH1(const TH1 *hist, int rebin = 0);
    bool FillFromTH2(const TH2 *hist, int rebin_x = 0, int rebin_y = 0);

    int NSyst() const { return static_cast<int>(m_syst_names.size()); }
    const std::string &SystName(int index) const { return m_syst_names.at(index); }

    bool HasSystematic(const std::string &name) const;
    int SystIndex(const std::string &name) const;

    void AddSystematicFrac(const char *name, const TH1 *pos_frac, const TH1 *neg_frac);
    void AddSystematicFrac2D(const char *name, const TH2 *pos_frac, const TH2 *neg_frac);

    bool GetFloatFlag(const std::string &name) const;
    bool GetLogNormalFlag(const std::string &name) const;
    void SetFloatFlag(const std::string &name, bool on);
    void SetLogNormalFlag(const std::string &name, bool on);

    void Linearise(const std::vector<std::string> &global_syst_order);
    bool IsLinearised() const { return m_linearised; }

    double BinFractionVaried(int ix, int iy, const double *pulls) const;

    void PrepareExclusionSums();

    TH1 *MakeTH1(const std::string &title = "") const;
    TH2 *MakeTH2(const std::string &title = "") const;

  private:
    int Index(int ix, int iy) const
    {
        return (m_ny > 1) ? (ix + iy * m_nx) : ix;
    }

    static inline double AsymDelta(double pull, double sigma_pos, double sigma_neg)
    {
        return (pull >= 0.0) ? (pull * sigma_pos) : (pull * sigma_neg);
    }

    bool CheckStats(const TH1 *hist, int verbose = 1) const;

  private:
    bool m_mutable = false;

    int m_nx = 1;
    int m_ny = 1;
    int m_true_bin_count = 1;
    double m_x_min = -1.0;
    double m_x_max = -1.0;
    double m_y_min = -1.0;
    double m_y_max = -1.0;

    double m_total_yield = 0.0;
    std::vector<double> m_bin_fraction;
    std::vector<double> m_bin_stat;

    std::vector<std::string> m_syst_names;
    std::vector<std::vector<NuXSecSystBin>> m_syst;
    std::vector<unsigned char> m_float_flag;
    std::vector<unsigned char> m_log_normal_flag;

    std::vector<int> m_syst_index_outer;
    std::vector<int> m_syst_index_inner;

    bool m_linearised = false;
    std::vector<double> m_lin_bin_fraction;
    std::vector<double> m_lin_bin_stat;
    std::vector<double> m_efficiency_sums;

};

} // namespace nuxsec

#endif // NUXSEC_IO_NUXSEC_TEMPLATE_H
