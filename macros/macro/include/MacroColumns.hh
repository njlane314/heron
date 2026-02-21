/* -- C++ -- */
/// \file macros/macro/include/MacroColumns.hh
/// \brief Header-only helpers for tree/RDataFrame column validation.

#ifndef heron_macro_MACROCOLUMNS_H
#define heron_macro_MACROCOLUMNS_H

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

namespace heron {
namespace macro {

inline bool has_column(const std::vector<std::string> &columns, const std::string &column)
{
  return std::find(columns.begin(), columns.end(), column) != columns.end();
}

inline std::vector<std::string>
missing_required_columns(const std::vector<std::string> &columns,
                         const std::initializer_list<std::string> &required_columns)
{
  std::vector<std::string> missing_columns;
  for (const std::string &required_column : required_columns) {
    if (!has_column(columns, required_column)) {
      missing_columns.push_back(required_column);
    }
  }

  return missing_columns;
}

inline void print_missing_columns(const std::string &macro_name,
                                  const std::vector<std::string> &missing_columns)
{
  if (missing_columns.empty()) {
    return;
  }

  std::cerr << "[" << macro_name << "] missing required columns:";
  for (const std::string &column_name : missing_columns) {
    std::cerr << " " << column_name;
  }
  std::cerr << "\n";
}

} // namespace macro
} // namespace heron

#endif
