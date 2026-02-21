/* -- C++ -- */
/**
 *  @file  ana/src/AnalysisModel.cpp
 *
 *  @brief Declarative base model helpers.
 */

#include "AnalysisModel.hh"

#include <utility>

namespace heron
{

Selection AnalysisModel::selection(std::string name, Cut c, Weight w)
{
    Selection s;
    s.name = std::move(name);
    s.cut = std::move(c);
    s.weight = std::move(w);
    m_selections.push_back(s);
    return s;
}

Hist1DSpec AnalysisModel::hist1d(std::string name,
                            std::string variable,
                            int bins,
                            double x_min,
                            double x_max,
                            std::string title,
                            std::string selection,
                            std::string weight)
{
    Hist1DSpec h;
    h.name = std::move(name);
    h.title = std::move(title);
    h.variable = std::move(variable);
    h.bins = bins;
    h.x_min = x_min;
    h.x_max = x_max;
    h.selection = std::move(selection);
    h.weight = std::move(weight);
    m_h1.push_back(h);
    return h;
}

SnapshotSpec AnalysisModel::snapshot(std::string name, std::vector<std::string> columns, std::string selection)
{
    SnapshotSpec s;
    s.name = std::move(name);
    s.columns = std::move(columns);
    s.selection = std::move(selection);
    m_snapshots.push_back(s);
    return s;
}


void AnalysisModel::clear() noexcept
{
    m_vars.clear();
    m_cuts.clear();
    m_weights.clear();
    m_selections.clear();
    m_h1.clear();
    m_snapshots.clear();
}

} // namespace heron
