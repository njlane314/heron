/* -- C++ -- */
/// \file macros/macro/include/MacroStyle.hh
/// \brief Header-only helper for common plotting style defaults.

#ifndef heron_macro_MACROSTYLE_H
#define heron_macro_MACROSTYLE_H

#include "TStyle.h"

namespace heron {
namespace macro {

inline void apply_default_plot_style(TStyle *p_style)
{
  if (p_style == NULL) {
    return;
  }

  p_style->SetOptStat(0);
  p_style->SetTitleBorderSize(0);
  p_style->SetLegendBorderSize(0);
  p_style->SetPadTickX(1);
  p_style->SetPadTickY(1);
}

} // namespace macro
} // namespace heron

#endif
