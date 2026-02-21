/* -- C++ -- */
/// \file macros/include/MacroEnv.hh
/// \brief Header-only environment and path helpers shared by ROOT macros.

#ifndef heron_macro_MACROENV_H
#define heron_macro_MACROENV_H

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace heron {
namespace macro {

inline std::string getenv_or(const char *key, const std::string &fallback)
{
  const char *value = std::getenv(key);
  if (!value) return fallback;

  const std::string text(value);
  if (text.empty()) return fallback;

  return text;
}

inline bool env_truthy(const char *value)
{
  if (!value) return false;

  const std::string text(value);
  return (text == "1" || text == "true" || text == "TRUE" || text == "yes" || text == "YES" ||
          text == "on" || text == "ON");
}

inline std::string plot_out_dir()
{
  return getenv_or("HERON_PLOT_DIR", "./scratch/plots");
}

inline std::string plot_out_fmt()
{
  return getenv_or("HERON_PLOT_FORMAT", "pdf");
}

inline bool file_exists(const std::string &path)
{
  std::ifstream file(path.c_str());
  return file.good();
}

inline std::string find_first_existing_path(const std::vector<std::string> &candidates)
{
  for (std::vector<std::string>::const_iterator it = candidates.begin(); it != candidates.end(); ++it) {
    if (file_exists(*it)) return *it;
  }

  return "";
}

inline std::string find_default_event_list_path()
{
  const std::vector<std::string> candidates{
      "/exp/uboone/data/users/nlane/heron/out/event/events.root"};

  return find_first_existing_path(candidates);
}

} // namespace macro
} // namespace heron

#endif
