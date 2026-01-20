/* -- C++ -- */
#ifndef NUXSEC_ANA_SYSTEMATICSTOOLS_H
#define NUXSEC_ANA_SYSTEMATICSTOOLS_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "TH1.h"
#include "TDecompSVD.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace nuxsec {

// -------------------------
// Slice quality (Eq. 10.1)
// -------------------------
struct SliceQuality {
  double purity = 0.0;       // P_s
  double completeness = 0.0; // C_s
};

inline SliceQuality ComputeSliceQuality(unsigned n_hits_slice_and_nu,
                                        unsigned n_hits_slice,
                                        unsigned n_hits_nu) {
  SliceQuality out;
  out.purity = (n_hits_slice > 0u) ? (double(n_hits_slice_and_nu) / double(n_hits_slice)) : 0.0;
  out.completeness = (n_hits_nu > 0u) ? (double(n_hits_slice_and_nu) / double(n_hits_nu)) : 0.0;
  return out;
}

// -------------------------
// Histogram/vector helpers
// -------------------------
inline std::vector<double> ToVector(const TH1& h) {
  const int nb = h.GetNbinsX();
  std::vector<double> v(nb, 0.0);
  for (int i = 1; i <= nb; ++i) v[i - 1] = h.GetBinContent(i);
  return v;
}

inline void AssertSameSize(const std::vector<double>& a,
                           const std::vector<double>& b,
                           const std::string& what) {
  if (a.size() != b.size()) {
    throw std::runtime_error(what + ": size mismatch " + std::to_string(a.size()) +
                             " vs " + std::to_string(b.size()));
  }
}

inline double Dot(const std::vector<double>& a, const std::vector<double>& b) {
  AssertSameSize(a, b, "Dot");
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// -----------------------------------------
// Remove global rate shift (Eq. 9.11)
// alpha = (T0^T (Tu - T0)) / (T0^T T0)
// s     = 1 + alpha
// R     = Tu - s T0
// -----------------------------------------
struct RateResidual {
  double alpha = 0.0;
  double scale = 1.0;
  std::vector<double> residual;
};

inline RateResidual RemoveGlobalRate(const std::vector<double>& nominal,
                                     const std::vector<double>& varied) {
  AssertSameSize(nominal, varied, "RemoveGlobalRate");

  const double denom = Dot(nominal, nominal);

  RateResidual out;
  if (denom <= 0.0) {
    // Nominal is empty; can't project meaningfully.
    out.alpha = 0.0;
    out.scale = 1.0;
    out.residual.resize(nominal.size());
    for (std::size_t i = 0; i < nominal.size(); ++i) out.residual[i] = varied[i] - nominal[i];
    return out;
  }

  std::vector<double> delta(nominal.size(), 0.0);
  for (std::size_t i = 0; i < nominal.size(); ++i) delta[i] = varied[i] - nominal[i];

  out.alpha = Dot(nominal, delta) / denom;
  out.scale = 1.0 + out.alpha;

  out.residual.resize(nominal.size());
  for (std::size_t i = 0; i < nominal.size(); ++i) out.residual[i] = varied[i] - out.scale * nominal[i];
  return out;
}

// -------------------------
// Envelope + coverage
// -------------------------
struct Envelope {
  std::vector<double> lo;
  std::vector<double> hi;
};

inline Envelope MakeEnvelope(const std::vector<std::vector<double>>& vecs) {
  if (vecs.empty()) throw std::runtime_error("MakeEnvelope: empty input");

  const std::size_t nb = vecs.front().size();
  for (const auto& v : vecs) {
    if (v.size() != nb) throw std::runtime_error("MakeEnvelope: inconsistent vector sizes");
  }

  Envelope env;
  env.lo.assign(nb, +std::numeric_limits<double>::infinity());
  env.hi.assign(nb, -std::numeric_limits<double>::infinity());

  for (const auto& v : vecs) {
    for (std::size_t i = 0; i < nb; ++i) {
      env.lo[i] = std::min(env.lo[i], v[i]);
      env.hi[i] = std::max(env.hi[i], v[i]);
    }
  }
  return env;
}

struct CoverageResult {
  double coverage = 0.0;         // weighted fraction inside
  double uncovered_weight = 0.0; // weighted fraction outside
  std::vector<int> uncovered_bins; // 1-based bin indices outside (after tolerance)
};

inline CoverageResult Coverage(const std::vector<double>& x,
                               const Envelope& env,
                               const std::optional<std::vector<double>>& weights = std::nullopt,
                               double tol = 0.0) {
  AssertSameSize(x, env.lo, "Coverage(x, env.lo)");
  AssertSameSize(x, env.hi, "Coverage(x, env.hi)");
  if (weights) AssertSameSize(x, *weights, "Coverage(weights)");

  double wsum = 0.0;
  double win = 0.0;

  CoverageResult out;

  for (std::size_t i = 0; i < x.size(); ++i) {
    const double w = weights ? (*weights)[i] : 1.0;
    wsum += w;

    const bool inside = (x[i] >= env.lo[i] - tol) && (x[i] <= env.hi[i] + tol);
    if (inside) {
      win += w;
    } else {
      out.uncovered_bins.push_back(int(i) + 1);
    }
  }

  if (wsum > 0.0) {
    out.coverage = win / wsum;
    out.uncovered_weight = 1.0 - out.coverage;
  } else {
    out.coverage = 0.0;
    out.uncovered_weight = 0.0;
  }
  return out;
}

// -------------------------------------------------
// Covariance from an ensemble of residual vectors
// V = (1/U) sum_u r_u r_u^T   (Eq. 9.12 form)
// (Use r_u = (T0 - Tu) for CV-anchored multisim,
//  Use r_u = R^(u) for detvar shape-only covariance)
// -------------------------------------------------
inline TMatrixD CovarianceFromResiduals(const std::vector<std::vector<double>>& residuals) {
  if (residuals.empty()) throw std::runtime_error("CovarianceFromResiduals: empty input");

  const std::size_t nb = residuals.front().size();
  for (const auto& r : residuals) {
    if (r.size() != nb) throw std::runtime_error("CovarianceFromResiduals: inconsistent vector sizes");
  }

  TMatrixD V(int(nb), int(nb));
  V.Zero();

  for (const auto& r : residuals) {
    for (std::size_t i = 0; i < nb; ++i) {
      for (std::size_t j = 0; j < nb; ++j) {
        V(int(i), int(j)) += r[i] * r[j];
      }
    }
  }

  V *= (1.0 / double(residuals.size()));
  return V;
}

// -------------------------
// SVD pseudo-inverse + chi2
// -------------------------
inline TMatrixD PseudoInverseSVD(const TMatrixD& m, double rcond = 1e-12) {
  if (m.GetNrows() != m.GetNcols()) throw std::runtime_error("PseudoInverseSVD: matrix not square");

  TDecompSVD svd(m);
  const Bool_t ok = svd.Decompose();
  if (!ok) throw std::runtime_error("PseudoInverseSVD: SVD decomposition failed");

  const TVectorD s = svd.GetSig();
  TMatrixD U = svd.GetU();
  TMatrixD V = svd.GetV();

  const int n = s.GetNrows();
  double smax = 0.0;
  for (int i = 0; i < n; ++i) smax = std::max(smax, s[i]);
  const double thresh = rcond * smax;

  TMatrixD Sinv(n, n);
  Sinv.Zero();
  for (int i = 0; i < n; ++i) {
    if (s[i] > thresh) Sinv(i, i) = 1.0 / s[i];
  }

  // m = U S V^T  =>  m^+ = V S^+ U^T
  TMatrixD UT(TMatrixD::kTransposed, U);
  TMatrixD out = V * Sinv * UT;
  return out;
}

inline double Chi2(const std::vector<double>& r, const TMatrixD& Vinv) {
  const int n = Vinv.GetNrows();
  if (Vinv.GetNcols() != n) throw std::runtime_error("Chi2: Vinv not square");
  if (int(r.size()) != n) throw std::runtime_error("Chi2: size mismatch");

  TVectorD rv(n);
  for (int i = 0; i < n; ++i) rv[i] = r[std::size_t(i)];

  TVectorD tmp = Vinv * rv;
  return rv * tmp;
}

} // namespace nuxsec

#endif
