/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisProcessor.hh
 *
 *  @brief Variable definitions for analysis RDataFrame processing
 */

#ifndef Nuxsec_ANA_ANALYSIS_PROCESSOR_H_INCLUDED
#define Nuxsec_ANA_ANALYSIS_PROCESSOR_H_INCLUDED

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace nuxsec
{

enum class Source {
    Data,
    Ext,
    MC,
    Unknown
};

enum class Slice {
    None,
    BeamInclusive,
    StrangenessInclusive
};

enum class Channel : std::uint8_t {
    OutFV = 1,
    External = 2,
    MuCC0pi_ge1p = 10,
    MuCC1pi = 11,
    MuCCPi0OrGamma = 12,
    MuCCNpi = 13,
    NC = 14,
    CCS1 = 15,
    CCSgt1 = 16,
    ECCC = 17,
    MuCCOther = 18,
    DataInclusive = 20,
    Unknown = 99
};

inline const char *ChannelLabel(Channel c)
{
    switch (c) {
    case Channel::OutFV:
        return "Out-FV";
    case Channel::External:
        return "External";
    case Channel::MuCC0pi_ge1p:
        return "CC0#pi, #geq1p";
    case Channel::MuCC1pi:
        return "CC1#pi^{#pm}";
    case Channel::MuCCPi0OrGamma:
        return "CC #pi^{0} / #gamma-rich";
    case Channel::MuCCNpi:
        return "CC N#pi^{#pm} (N>1)";
    case Channel::NC:
        return "NC (in-FV)";
    case Channel::CCS1:
        return "CC + 1 strange";
    case Channel::CCSgt1:
        return "CC + >1 strange";
    case Channel::ECCC:
        return "#nu_{e} CC (no strange)";
    case Channel::MuCCOther:
        return "CC other";
    case Channel::DataInclusive:
        return "Data (incl.)";
    default:
        return "Other";
    }
}

namespace sample
{

enum class origin {
    data,
    beam,
    strangeness,
    ext,
    dirt,
    unknown
};

inline origin origin_from(std::string_view s)
{
    if (s == "data")
        return origin::data;
    if (s == "beam")
        return origin::beam;
    if (s == "strangeness")
        return origin::strangeness;
    if (s == "ext" || s == "external")
        return origin::ext;
    if (s == "dirt")
        return origin::dirt;
    if (s == "mc")
        return origin::beam;
    return origin::unknown;
}

inline Source to_source(origin o)
{
    switch (o) {
    case origin::data:
        return Source::Data;
    case origin::ext:
        return Source::Ext;
    case origin::unknown:
        return Source::Unknown;
    default:
        return Source::MC;
    }
}

inline Slice to_slice(origin o)
{
    switch (o) {
    case origin::beam:
        return Slice::BeamInclusive;
    case origin::strangeness:
        return Slice::StrangenessInclusive;
    default:
        return Slice::None;
    }
}

inline origin from_source_slice(Source src, Slice sl)
{
    if (src == Source::Data)
        return origin::data;
    if (src == Source::Ext)
        return origin::ext;
    if (src == Source::Unknown)
        return origin::unknown;
    switch (sl) {
    case Slice::StrangenessInclusive:
        return origin::strangeness;
    case Slice::BeamInclusive:
        return origin::beam;
    default:
        return origin::beam;
    }
}

} // namespace sample

struct Frame
{
    std::shared_ptr<ROOT::RDataFrame> df;
    mutable std::optional<ROOT::RDF::RNode> node;

    Frame() = default;
    Frame(std::shared_ptr<ROOT::RDataFrame> df_in, ROOT::RDF::RNode node_in)
        : df(std::move(df_in))
        , node(std::move(node_in))
    {
    }

    auto report() const
    {
        if (!node)
            throw std::runtime_error("Frame::report: node is not initialised");
        return node->Report();
    }

    ROOT::RDF::RNode rnode() const
    {
        if (!node)
            throw std::runtime_error("Frame::rnode: node is not initialised");
        return *node;
    }
};

struct Entry
{
    std::string beamline;
    std::string period;
    Source source = Source::Unknown;
    Slice slice = Slice::BeamInclusive;
    sample::origin kind = sample::origin::unknown;
    std::vector<std::string> files;
    std::string file;

    double pot_nom = 0.0;
    double pot_eqv = 0.0;
    double trig_nom = 0.0;
    double trig_eqv = 0.0;

    Frame nominal;
    std::unordered_map<std::string, Frame> detvars;

    ROOT::RDF::RNode rnode() const { return nominal.rnode(); }
    const Frame *detvar(const std::string &tag) const
    {
        auto it = detvars.find(tag);
        return it == detvars.end() ? nullptr : &it->second;
    }
    std::vector<std::string> variation_tags() const
    {
        std::vector<std::string> out;
        out.reserve(detvars.size());
        for (auto const &kv : detvars)
            out.push_back(kv.first);
        return out;
    }
};

struct ProcessorEntry
{
    Source source = Source::Unknown;
    double pot_nom = 0.0;
    double pot_eqv = 0.0;
    double trig_nom = 0.0;
    double trig_eqv = 0.0;
};

/** \brief Apply analysis variable definitions to an RDataFrame. */
class AnalysisProcessor
{
  public:
    ROOT::RDF::RNode Run(ROOT::RDF::RNode node, const ProcessorEntry &rec) const;
    static const AnalysisProcessor &Processor();

  private:
    static constexpr double kRecognisedPurityMin;
    static constexpr double kRecognisedCompletenessMin;

    static constexpr float kTrainingFraction;
    static constexpr bool kTrainingIncludeExt;

    static bool IsInTruthVolume(float x, float y, float z) noexcept;
    static bool IsInRecoVolume(float x, float y, float z) noexcept;
};

} // namespace nuxsec

#endif // Nuxsec_ANA_ANALYSIS_PROCESSOR_H_INCLUDED
