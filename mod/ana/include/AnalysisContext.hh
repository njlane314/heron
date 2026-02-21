/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisContext.hh
 *
 *  @brief Lightweight analysis macro context that bundles analysis policy and
 *         configured services for one macro invocation.
 */

#ifndef HERON_ANA_ANALYSIS_CONTEXT_H
#define HERON_ANA_ANALYSIS_CONTEXT_H

#include <memory>
#include <utility>

namespace heron
{

/**
 *  \brief Shared macro context initialised once at macro entry.
 *
 *  This class is intentionally small and ownership-focused so macros can do
 *  this at the top level:
 *
 *  \code
 *  AnalysisContext ctx(pol, MakeServices("services.fcl"));
 *  \endcode
 *
 *  Example `services.fcl` content (adapt to your local service types):
 *
 *  \code
 *  services: {
 *    analysis_config: {
 *      name: "numi-analysis"
 *      tree_name: "events"
 *    }
 *
 *    selection: {
 *      channel: "inclusive_mu_cc"
 *      min_muon_track_length_cm: 30.0
 *    }
 *
 *    columns: {
 *      enable_reco_vertex_distance: true
 *      enable_muon_candidate_score: true
 *    }
 *  }
 *  \endcode
 */
template <typename TPolicy, typename TServices>
class AnalysisContext
{
  public:
    using policy_type = TPolicy;
    using services_type = TServices;

    AnalysisContext(const policy_type &policy, services_type services)
        : m_policy(policy)
        , m_services(std::make_shared<services_type>(std::move(services)))
    {
    }

    const policy_type &policy() const noexcept { return m_policy; }
    const services_type &services() const noexcept { return *m_services; }
    const services_type *operator->() const noexcept { return m_services.get(); }

  private:
    const policy_type &m_policy;
    std::shared_ptr<services_type> m_services;
};

} // namespace heron

#endif // HERON_ANA_ANALYSIS_CONTEXT_H
