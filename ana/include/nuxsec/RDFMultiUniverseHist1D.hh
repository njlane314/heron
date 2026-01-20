/* -- C++ -- */
#ifndef NUXSEC_ANA_RDFMULTIUNIVERSEHIST1D_H
#define NUXSEC_ANA_RDFMULTIUNIVERSEHIST1D_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"

namespace nuxsec {

/// Configuration for building a 1D multisim histogram set.
struct Hist1DConfig {
  std::string name;
  std::string title;
  int nbins = 20;
  double xmin = 0.0;
  double xmax = 1.0;
};

struct MultiUniverseHist1D {
  std::unique_ptr<TH1D> nominal;                // T^(0)
  std::vector<std::unique_ptr<TH1D>> universes; // T^(s,u), size U
};

// Weight vector column must be ROOT::RVec<float> (common for multisim branches).
inline MultiUniverseHist1D BuildMultiUniverseHist1D_FloatWeights(ROOT::RDF::RNode df,
                                                                 const Hist1DConfig& cfg,
                                                                 const std::string& zCol,
                                                                 const std::string& wCvCol,
                                                                 const std::string& wSysVecCol,
                                                                 std::size_t nUniverses) {
  if (nUniverses == 0) throw std::runtime_error("BuildMultiUniverseHist1D: nUniverses=0");

  const unsigned nSlots = df.GetNSlots();

  // Slot-local histograms (thread-safe fill).
  std::vector<std::unique_ptr<TH1D>> hNomSlot(nSlots);
  std::vector<std::vector<std::unique_ptr<TH1D>>> hUniSlot(nSlots);

  for (unsigned s = 0; s < nSlots; ++s) {
    hNomSlot[s] = std::make_unique<TH1D>(
        (cfg.name + "_nom_slot" + std::to_string(s)).c_str(),
        cfg.title.c_str(),
        cfg.nbins, cfg.xmin, cfg.xmax);
    hNomSlot[s]->SetDirectory(nullptr);

    hUniSlot[s].resize(nUniverses);
    for (std::size_t u = 0; u < nUniverses; ++u) {
      hUniSlot[s][u] = std::make_unique<TH1D>(
          (cfg.name + "_u" + std::to_string(u) + "_slot" + std::to_string(s)).c_str(),
          cfg.title.c_str(),
          cfg.nbins, cfg.xmin, cfg.xmax);
      hUniSlot[s][u]->SetDirectory(nullptr);
    }
  }

  // Fill.
  df.ForeachSlot(
      [&](unsigned slot, double z, double wcv, const ROOT::RVec<float>& wsys) {
        hNomSlot[slot]->Fill(z, wcv);

        const std::size_t n = std::min<std::size_t>(wsys.size(), nUniverses);
        for (std::size_t u = 0; u < n; ++u) {
          hUniSlot[slot][u]->Fill(z, wcv * double(wsys[u]));
        }
      },
      {zCol, wCvCol, wSysVecCol});

  // Merge slots.
  MultiUniverseHist1D out;
  out.nominal = std::make_unique<TH1D>(
      (cfg.name + "_nom").c_str(), cfg.title.c_str(), cfg.nbins, cfg.xmin, cfg.xmax);
  out.nominal->SetDirectory(nullptr);

  out.universes.resize(nUniverses);
  for (std::size_t u = 0; u < nUniverses; ++u) {
    out.universes[u] = std::make_unique<TH1D>(
        (cfg.name + "_u" + std::to_string(u)).c_str(), cfg.title.c_str(), cfg.nbins, cfg.xmin,
        cfg.xmax);
    out.universes[u]->SetDirectory(nullptr);
  }

  for (unsigned s = 0; s < nSlots; ++s) {
    out.nominal->Add(hNomSlot[s].get());
    for (std::size_t u = 0; u < nUniverses; ++u) out.universes[u]->Add(hUniSlot[s][u].get());
  }

  return out;
}

} // namespace nuxsec

#endif
