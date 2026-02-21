/* -- C++ -- */
/// \file macros/macro/include/MacroIO.hh
/// \brief Header-only helpers for macro-level input/output checks.

#ifndef heron_macro_MACROIO_H
#define heron_macro_MACROIO_H

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>

#include "TFile.h"

namespace heron {
namespace macro {

inline bool looks_like_root_file(const std::string &path)
{
  if (path.size() < 5) {
    return false;
  }

  std::string extension = path.substr(path.size() - 5);
  std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });

  return extension == ".root";
}

inline std::unique_ptr<TFile> open_root_file_read(const std::string &path)
{
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    return std::unique_ptr<TFile>();
  }

  return file;
}

inline bool validate_root_input_path(const std::string &path)
{
  if (path.empty()) {
    return false;
  }

  return looks_like_root_file(path);
}

} // namespace macro
} // namespace heron

#endif
