/* -- C++ -- */
/// \file macros/include/MacroColumns.hh
/// \brief Header-only column checks shared by ROOT macros.

#ifndef heron_macro_MACROCOLUMNS_H
#define heron_macro_MACROCOLUMNS_H

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace heron {
namespace macro {

inline bool has_column(const std::vector<std::string> &available_columns, const std::string &column_name)
{
  for (std::vector<std::string>::const_iterator it = available_columns.begin(); it != available_columns.end(); ++it) {
    if (*it == column_name) return true;
  }

  return false;
}

inline std::vector<std::string> missing_required_columns(const std::vector<std::string> &available_columns,
                                                         const std::vector<std::string> &required_columns)
{
  std::vector<std::string> missing_columns;
  for (std::vector<std::string>::const_iterator it = required_columns.begin(); it != required_columns.end(); ++it) {
    if (!has_column(available_columns, *it)) missing_columns.push_back(*it);
  }

  return missing_columns;
}

inline std::vector<std::string> missing_required_columns(const std::unordered_set<std::string> &available_columns,
                                                         const std::vector<std::string> &required_columns)
{
  std::vector<std::string> missing_columns;
  for (std::vector<std::string>::const_iterator it = required_columns.begin(); it != required_columns.end(); ++it) {
    if (available_columns.find(*it) == available_columns.end()) missing_columns.push_back(*it);
  }

  return missing_columns;
}

inline void print_missing_columns(const std::vector<std::string> &missing_columns,
                                  std::ostream &stream = std::cerr)
{
  if (missing_columns.empty()) return;

  stream << "missing required columns:";
  for (std::vector<std::string>::const_iterator it = missing_columns.begin(); it != missing_columns.end(); ++it) {
    stream << " " << *it;
  }
  stream << "\n";
}

inline int require_columns(const std::unordered_set<std::string> &available_columns,
                           const std::vector<std::string> &required_columns,
                           const std::string &macro_name,
                           const std::string &label,
                           std::ostream &stream = std::cerr)
{
  const std::vector<std::string> missing_columns = missing_required_columns(available_columns, required_columns);
  if (missing_columns.empty()) return 0;

  stream << "[" << macro_name << "] missing required columns for " << label << ":\n";
  for (std::vector<std::string>::const_iterator it = missing_columns.begin(); it != missing_columns.end(); ++it) {
    stream << "  - " << *it << "\n";
  }

  return 1;
}

} // namespace macro
} // namespace heron

#endif
