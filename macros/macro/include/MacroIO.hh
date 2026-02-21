/* -- C++ -- */
/// \file macros/macro/include/MacroIO.hh
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

} // namespace macro
} // namespace heron

#endif
