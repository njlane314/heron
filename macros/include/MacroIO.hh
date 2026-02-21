/* -- C++ -- */
/// \file macros/include/MacroIO.hh
/// \brief Header-only I/O helpers shared by ROOT macros.

#ifndef heron_macro_MACROIO_H
#define heron_macro_MACROIO_H

#include <algorithm>
#include <cctype>
#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"

namespace heron {
namespace macro {

inline bool looks_like_root_file(const std::string &path)
{
  if (path.size() < 5) return false;

  std::string suffix = path.substr(path.size() - 5);
  std::transform(suffix.begin(), suffix.end(), suffix.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  return suffix == ".root";
}

inline std::unique_ptr<TFile> open_root_file_read(const std::string &path)
{
  TFile *p_file = TFile::Open(path.c_str(), "READ");
  if (p_file == NULL) return std::unique_ptr<TFile>();

  if (p_file->IsZombie()) {
    p_file->Close();
    delete p_file;
    return std::unique_ptr<TFile>();
  }

  return std::unique_ptr<TFile>(p_file);
}

inline bool validate_root_input_path(const std::string &path, std::ostream &stream = std::cerr)
{
  if (path.empty()) {
    stream << "input path is empty\n";
    return false;
  }

  if (!looks_like_root_file(path)) {
    stream << "input path does not look like a ROOT file: " << path << "\n";
    return false;
  }

  return true;
}

inline bool looks_like_event_list_root(const std::string &path)
{
  std::unique_ptr<TFile> p_file = open_root_file_read(path);
  if (!p_file) return false;

  const bool has_sample_refs = (p_file->Get("sample_refs") != NULL);
  const bool has_events_tree = (p_file->Get("events") != NULL);
  const bool has_event_tree_key = (p_file->Get("event_tree") != NULL);

  return has_sample_refs && (has_events_tree || has_event_tree_key);
}

} // namespace macro
} // namespace heron

#endif
