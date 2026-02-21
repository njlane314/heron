/* -- C++ -- */
/// \file macros/macro/MacroGuard.hh
/// \brief Shared helper to provide a consistent exception guard for ROOT macros.

#ifndef heron_macro_MACROGUARD_H
#define heron_macro_MACROGUARD_H

#include <exception>
#include <iostream>
#include <string>
#include <utility>

namespace heron {
namespace macro {

template <typename Callable>
inline decltype(auto) run_with_guard(const std::string &macro_name, Callable &&callable)
{
  try {
    return std::forward<Callable>(callable)();
  } catch (const std::exception &ex) {
    std::cerr << "[" << macro_name << "] exception: " << ex.what() << "\n";
    throw;
  } catch (...) {
    std::cerr << "[" << macro_name << "] unknown exception\n";
    throw;
  }
}

} // namespace macro
} // namespace heron

#endif
