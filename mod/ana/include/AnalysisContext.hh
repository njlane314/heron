/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisContext.hh
 *
 *  @brief Analysis macro context that centralises execution policy, services,
 *         provenance metadata and output destinations.
 */

#ifndef HERON_ANA_ANALYSIS_CONTEXT_H
#define HERON_ANA_ANALYSIS_CONTEXT_H

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

namespace heron
{

/**
 *  \brief Shared macro context initialised once at macro entry.
 *
 *  AnaContext is the macro "environment": analyses read policy/services/
 *  metadata/outputs from this object, rather than steering thread policy,
 *  process metadata, or output paths directly.
 */
template <typename TPolicy, typename TServices = std::nullptr_t>
class AnalysisContext
{
  public:
    using policy_type = TPolicy;
    using services_type = TServices;

    struct Metadata
    {
        std::string schema_version;
        std::string config_tag;
        std::string git_hash;
        std::string analysis_label;
    };

    struct Outputs
    {
        std::string artefact_dir;
        std::string plot_dir;
        std::string table_dir;
    };

    AnalysisContext(const policy_type &policy, services_type services)
        : m_policy(policy)
        , m_services(std::make_shared<services_type>(std::move(services)))
        , m_metadata(default_metadata("analysis"))
        , m_outputs(default_outputs())
    {
    }

    AnalysisContext(const policy_type &policy,
                    services_type services,
                    Metadata metadata,
                    Outputs outputs)
        : m_policy(policy)
        , m_services(std::make_shared<services_type>(std::move(services)))
        , m_metadata(std::move(metadata))
        , m_outputs(std::move(outputs))
    {
    }

    explicit AnalysisContext(const policy_type &policy,
                             const std::string &analysis_label = "analysis")
        : m_policy(policy)
        , m_services(std::make_shared<services_type>(services_type()))
        , m_metadata(default_metadata(analysis_label))
        , m_outputs(default_outputs())
    {
    }

    const policy_type &policy() const noexcept { return m_policy; }
    const services_type &services() const noexcept { return *m_services; }
    const services_type *operator->() const noexcept { return m_services.get(); }

    const Metadata &metadata() const noexcept { return m_metadata; }
    const Outputs &outputs() const noexcept { return m_outputs; }

    bool implicit_mt_enabled() const noexcept
    {
        return m_policy.implicit_mt_enabled();
    }

    bool apply_runtime(const std::string &label = "AnalysisContext") const
    {
        return m_policy.apply(label);
    }

  private:
    static std::string getenv_or_default(const char *name, const char *fallback)
    {
        const char *v = std::getenv(name);
        if (v == NULL || *v == '\0')
            return fallback;

        return std::string(v);
    }

    static Metadata default_metadata(const std::string &analysis_label)
    {
        Metadata m;
        m.schema_version = getenv_or_default("HERON_SCHEMA_VERSION", "unknown");
        m.config_tag = getenv_or_default("HERON_CONFIG_TAG", "default");
        m.git_hash = getenv_or_default("HERON_GIT_HASH", "unknown");
        m.analysis_label = analysis_label;
        return m;
    }

    static Outputs default_outputs()
    {
        Outputs o;
        o.artefact_dir = getenv_or_default("HERON_ARTIFACT_DIR", "./scratch/artifacts");
        o.plot_dir = getenv_or_default("HERON_PLOT_DIR", "./scratch/plots");
        o.table_dir = getenv_or_default("HERON_TABLE_DIR", "./scratch/tables");
        return o;
    }

    policy_type m_policy;
    std::shared_ptr<services_type> m_services;
    Metadata m_metadata;
    Outputs m_outputs;
};

} // namespace heron


#endif // HERON_ANA_ANALYSIS_CONTEXT_H
