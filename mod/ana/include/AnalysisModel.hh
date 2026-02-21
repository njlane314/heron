/* -- C++ -- */
/**
 *  @file  ana/include/AnalysisModel.hh
 *
 *  @brief Declarative base model for analysis variables, selections,
 *         histograms, and snapshots.
 */

#ifndef HERON_ANA_ANALYSIS_MODEL_H
#define HERON_ANA_ANALYSIS_MODEL_H

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "AnalysisContext.hh"

namespace heron
{

/**
 *  \brief Named variable declaration with optional dependency names.
 */
template <typename TValue>
struct Var
{
    std::string name;
    std::function<TValue()> expression;
    std::vector<std::string> dependencies;
};

/**
 *  \brief Named cut predicate with optional dependency names.
 */
struct Cut
{
    std::string name;
    std::function<bool()> predicate;
    std::vector<std::string> dependencies;
};

/**
 *  \brief Named weight expression with optional dependency names.
 */
struct Weight
{
    std::string name;
    std::function<double()> expression;
    std::vector<std::string> dependencies;

    bool enabled = false;
};

/**
 *  \brief Named selection composed from one cut and an optional weight.
 */
struct Selection
{
    std::string name;
    Cut cut;
    Weight weight;
};

/**
 *  \brief Descriptor for a one-dimensional histogram declaration.
 */
struct Hist1DSpec
{
    std::string name;
    std::string title;
    std::string variable;
    int bins = 0;
    double x_min = 0.0;
    double x_max = 0.0;
    std::string selection;
    std::string weight;
};

/**
 *  \brief Descriptor for a snapshot declaration.
 */
struct SnapshotSpec
{
    std::string name;
    std::vector<std::string> columns;
    std::string selection;
};

/**
 *  \brief Base class for declarative analysis models.
 *
 *  Typical usage in a concrete model:
 *
 *  \code
 *  class MuonSelectionModel : public heron::AnalysisModel
 *  {
 *    public:
 *      void define() override
 *      {
 *          const auto p_muon_p = var("muon_p", []() { return 0.0; }, {"reco_muon_p"});
 *          const auto c_fid = cut("fiducial", []() { return true; }, {"reco_vertex_x", "reco_vertex_y", "reco_vertex_z"});
 *          const auto w_cv = weight("cv", []() { return 1.0; }, {"event_weight_cv"});
 *          const auto s_nominal = selection("nominal", c_fid, w_cv);
 *
 *          hist1d("h_muon_p", p_muon_p.name, 40, 0.0, 2.0, "Muon momentum;p [GeV];Events", s_nominal.name, w_cv.name);
 *          snapshot("selected_events", {"run", "subrun", "event", p_muon_p.name}, s_nominal.name);
 *      }
 *  };
 *  \endcode
 */
class AnalysisModel
{
  public:
    virtual ~AnalysisModel() = default;

    /// Declare model variables, cuts, histograms, and snapshots.
    virtual void define() = 0;

    /// Optional hook for external configuration/services.
    virtual void configure() {}

    /// Optional context-aware configuration hook.
    template <typename TPolicy, typename TServices>
    void configure(const AnalysisContext<TPolicy, TServices> &context)
    {
        (void) context;
        configure();
    }

    /// Reset and build the model using default configuration.
    bool initialise()
    {
        clear();
        configure();
        define();
        return true;
    }

    /// Reset and build the model using an analysis context.
    template <typename TPolicy, typename TServices>
    bool initialise(const AnalysisContext<TPolicy, TServices> &context)
    {
        clear();
        configure(context);
        define();
        return true;
    }

    const std::vector<Var<double> > &vars() const noexcept { return m_vars; }
    const std::vector<Cut> &cuts() const noexcept { return m_cuts; }
    const std::vector<Weight> &weights() const noexcept { return m_weights; }
    const std::vector<Selection> &selections() const noexcept { return m_selections; }
    const std::vector<Hist1DSpec> &h1() const noexcept { return m_h1; }
    const std::vector<SnapshotSpec> &snapshots() const noexcept { return m_snapshots; }

  protected:
    template <typename TExpression>
    Var<double> var(std::string name, TExpression expression, std::vector<std::string> dependencies = {})
    {
        Var<double> v;
        v.name = std::move(name);
        v.expression = std::function<double()>(std::move(expression));
        v.dependencies = std::move(dependencies);
        m_vars.push_back(v);
        return v;
    }

    template <typename TPredicate>
    Cut cut(std::string name, TPredicate predicate, std::vector<std::string> dependencies = {})
    {
        Cut c;
        c.name = std::move(name);
        c.predicate = std::function<bool()>(std::move(predicate));
        c.dependencies = std::move(dependencies);
        m_cuts.push_back(c);
        return c;
    }

    template <typename TExpression>
    Weight weight(std::string name, TExpression expression, std::vector<std::string> dependencies = {})
    {
        Weight w;
        w.name = std::move(name);
        w.expression = std::function<double()>(std::move(expression));
        w.dependencies = std::move(dependencies);
        w.enabled = true;
        m_weights.push_back(w);
        return w;
    }

    Selection selection(std::string name, Cut c, Weight w = {});
    Hist1DSpec hist1d(std::string name,
                      std::string variable,
                      int bins,
                      double x_min,
                      double x_max,
                      std::string title = "",
                      std::string selection = "",
                      std::string weight = "");
    SnapshotSpec snapshot(std::string name, std::vector<std::string> columns, std::string selection = "");
    void clear() noexcept;

  private:
    std::vector<Var<double> > m_vars;
    std::vector<Cut> m_cuts;
    std::vector<Weight> m_weights;
    std::vector<Selection> m_selections;
    std::vector<Hist1DSpec> m_h1;
    std::vector<SnapshotSpec> m_snapshots;
};

} // namespace heron

#endif // HERON_ANA_ANALYSIS_MODEL_H
