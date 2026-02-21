/* -- C++ -- */
/// \file macros/include/MacroGuard.hh
/// \brief Shared helpers to provide consistent exception guards for ROOT macros.

#ifndef heron_macro_MACROGUARD_H
#define heron_macro_MACROGUARD_H

#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

namespace heron {
namespace macro {

inline void log_guard_exception(const std::string &macro_name, const std::exception &ex)
{
  std::cerr << "[" << macro_name << "] exception: " << ex.what() << "\n";
}

inline void log_guard_unknown_exception(const std::string &macro_name)
{
  std::cerr << "[" << macro_name << "] unknown exception\n";
}

template <typename Callable>
inline decltype(auto) run_with_guard(const std::string &macro_name, Callable &&callable)
{
  try {
    return std::invoke(std::forward<Callable>(callable));
  } catch (const std::exception &ex) {
    log_guard_exception(macro_name, ex);
    throw;
  } catch (...) {
    log_guard_unknown_exception(macro_name);
    throw;
  }
}

template <typename Callable>
inline bool run_with_guard_no_throw(const std::string &macro_name, Callable &&callable)
{
  try {
    std::invoke(std::forward<Callable>(callable));
    return true;
  } catch (const std::exception &ex) {
    log_guard_exception(macro_name, ex);
    return false;
  } catch (...) {
    log_guard_unknown_exception(macro_name);
    return false;
  }
}

} // namespace macro
} // namespace heron

#endif
