/* -- C++ -- */
/**
 *  @file  ana/include/SystematicsBuilder.hh
 *
 *  @brief Build systematic templates from unisim and multisim variations.
 */

#ifndef Nuxsec_ANA_SYSTEMATICS_BUILDER_H_INCLUDED
#define Nuxsec_ANA_SYSTEMATICS_BUILDER_H_INCLUDED

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>

#include "SampleTypes.hh"
#include "SystematicsSpec.hh"
#include "TemplateSpec.hh"

namespace nuxsec
{

struct SampleListEntry
{
    std::string sample_name;
    std::string sample_kind;
    std::string beam_mode;
    std::string output_path;
};

class SystematicsBuilder final
{
  public:
    struct Options
    {
        int nthreads = 1;
        bool include_overlay = true;
        bool include_dirt = true;
        bool include_strangeness = true;
        bool clamp_negative_bins = true;
    };

    static void BuildAll(const std::vector<SampleListEntry> &entries,
                         const std::string &tree_name,
                         const std::vector<TemplateSpec1D> &fit_specs,
                         const std::string &template_root_path,
                         const SystematicsConfig &cfg,
                         const Options &opt);

  private:
    static bool IsVariedSampleKind(SampleKind k, const Options &opt);

    static std::unique_ptr<TH1D> ReadNominalHist(const std::string &root_path,
                                                 const std::string &sample_name,
                                                 const std::string &hist_name);

    static void WriteOneSyst(const std::string &root_path,
                             const std::string &sample_name,
                             const std::string &syst_name,
                             const std::string &variation,
                             const std::vector<std::pair<std::string, const TH1 *>> &hists);

    static void BuildUnisim(const Sample &sample,
                            const std::string &tree_name,
                            const std::vector<TemplateSpec1D> &specs,
                            const std::string &template_root_path,
                            const UnisimSpec &spec);

    static void BuildMultisimJointEigenmodes(const std::vector<Sample> &samples,
                                             const std::vector<std::string> &sample_names,
                                             const std::string &tree_name,
                                             const std::vector<TemplateSpec1D> &specs,
                                             const std::string &template_root_path,
                                             const MultisimSpec &mspec,
                                             const Options &opt);

    static int DetectUniverseCount(const Sample &sample,
                                   const std::string &tree_name,
                                   const std::string &vec_branch);

    static void ClampNonNegative(TH1D &h);
};

} // namespace nuxsec

#endif // Nuxsec_ANA_SYSTEMATICS_BUILDER_H_INCLUDED
