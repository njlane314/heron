#include "TBranch.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"

#include <cstdio>
#include <cctype>
#include <algorithm>
#include <map>
#include <set>
#include "include/MacroGuard.hh"

namespace {

struct MultisimSummary {
  std::set<int> universe_ids;
  int histogram_count;
  bool has_1d;
  bool has_2d;
  TString example_1d;
  TString example_2d;

  MultisimSummary() : histogram_count(0), has_1d(false), has_2d(false) {}
};

void print_indent(const int level) {
  for (int i = 0; i < level; ++i) {
    std::printf("  ");
  }
}

TString make_relative_path(const TString& object_path, const TString& file_path) {
  TString rel = object_path;
  if (rel.BeginsWith(file_path)) {
    rel.Remove(0, file_path.Length());
    if (rel.BeginsWith("/")) {
      rel.Remove(0, 1);
    }
  }
  return rel;
}

TString format_histogram_summary_line(const TString& rel_path, TObject* object) {
  if (object == NULL) {
    return "";
  }

  // TH2 inherits from TH1, so check TH2 first.
  if (object->InheritsFrom("TH2")) {
    TH2* h2 = dynamic_cast<TH2*>(object);
    if (h2 == NULL) {
      return "";
    }

    const TAxis* x = h2->GetXaxis();
    const TAxis* y = h2->GetYaxis();

    TString line;
    line.Form("%s [%s] Title=\"%s\" X=\"%s\" bins=%d range=[%.6g,%.6g] Y=\"%s\" bins=%d range=[%.6g,%.6g]",
              rel_path.Data(),
              h2->ClassName(),
              h2->GetTitle(),
              x ? x->GetTitle() : "",
              x ? x->GetNbins() : 0,
              x ? x->GetXmin() : 0.0,
              x ? x->GetXmax() : 0.0,
              y ? y->GetTitle() : "",
              y ? y->GetNbins() : 0,
              y ? y->GetXmin() : 0.0,
              y ? y->GetXmax() : 0.0);
    return line;
  }

  if (object->InheritsFrom("TH1")) {
    TH1* h1 = dynamic_cast<TH1*>(object);
    if (h1 == NULL) {
      return "";
    }

    const TAxis* x = h1->GetXaxis();

    TString line;
    line.Form("%s [%s] Title=\"%s\" X=\"%s\" bins=%d range=[%.6g,%.6g]",
              rel_path.Data(),
              h1->ClassName(),
              h1->GetTitle(),
              x ? x->GetTitle() : "",
              x ? x->GetNbins() : 0,
              x ? x->GetXmin() : 0.0,
              x ? x->GetXmax() : 0.0);
    return line;
  }

  return "";
}

void print_branch_list(const TObjArray* branches, const int indent_level) {
  if (branches == NULL) {
    return;
  }

  for (int i = 0; i < branches->GetEntries(); ++i) {
    TBranch* branch = dynamic_cast<TBranch*>(branches->At(i));
    if (branch == NULL) {
      continue;
    }

    print_indent(indent_level);
    std::printf("- %s\n", branch->GetName());
    print_branch_list(branch->GetListOfBranches(), indent_level + 1);
  }
}

void print_tree(const TString& object_path, TTree* tree) {
  if (tree == NULL) {
    return;
  }

  std::printf("\nTree: %s\n", object_path.Data());
  std::printf("  Entries: %lld\n", static_cast<Long64_t>(tree->GetEntries()));
  std::printf("  Branches:\n");
  print_branch_list(tree->GetListOfBranches(), 2);
}

int parse_universe_id(const TString& histogram_name) {
  const Ssiz_t uni_position = histogram_name.Index("_Uni_");
  if (uni_position < 0) {
    return -1;
  }

  const Ssiz_t id_start = uni_position + 5;
  Ssiz_t id_end = id_start;
  while (id_end < histogram_name.Length() && std::isdigit(histogram_name[id_end])) {
    ++id_end;
  }

  if (id_end == id_start) {
    return -1;
  }

  return TString(histogram_name(id_start, id_end - id_start)).Atoi();
}

TString normalise_multisim_type(const TString& histogram_name) {
  const Ssiz_t uni_position = histogram_name.Index("_Uni_");
  if (uni_position < 0) {
    return "";
  }

  TString prefix = histogram_name(0, uni_position);
  const Ssiz_t id_start = uni_position + 5;
  Ssiz_t suffix_start = id_start;
  while (suffix_start < histogram_name.Length() && std::isdigit(histogram_name[suffix_start])) {
    ++suffix_start;
  }

  TString suffix = histogram_name(suffix_start, histogram_name.Length() - suffix_start);
  if (suffix.EndsWith("_2D")) {
    suffix.Remove(suffix.Length() - 3);
  }

  return prefix + suffix;
}

void scan_directory(TDirectory* directory,
                    const TString& directory_path,
                    const TString& file_path,
                    int& tree_count,
                    int& object_count,
                    std::map<TString, int>& class_counts,
                    std::set<TString>& directory_paths,
                    std::map<TString, MultisimSummary>& multisim_summaries,
                    std::set<TString>& non_multisim_hist_summaries) {
  if (directory == NULL) {
    return;
  }

  directory_paths.insert(directory_path);

  TIter next_key(directory->GetListOfKeys());
  TKey* key = NULL;
  while ((key = dynamic_cast<TKey*>(next_key())) != NULL) {
    TClass* key_class = TClass::GetClass(key->GetClassName());
    const bool is_tree_key = key_class != NULL && key_class->InheritsFrom(TTree::Class());
    const bool is_directory_key = key_class != NULL && key_class->InheritsFrom(TDirectory::Class());

    TObject* object = key->ReadObj();
    if (object == NULL) {
      continue;
    }

    TString object_path = directory_path;
    if (!object_path.EndsWith("/")) {
      object_path += "/";
    }
    object_path += key->GetName();

    ++object_count;
    class_counts[key->GetClassName()] += 1;

    const bool is_multisim_path = object_path.Contains("/Multisims/");
    const bool is_hist = object->InheritsFrom("TH1") || object->InheritsFrom("TH2");
    if (is_hist) {
      const TString rel_path = make_relative_path(object_path, file_path);
      const TString line = format_histogram_summary_line(rel_path, object);
      if (!line.IsNull()) {
        if (!is_multisim_path) {
          // This is the list you can paste back for us to identify the CV / parent / detsmear / other plots.
          non_multisim_hist_summaries.insert(line);
        }
      }
    }

    if (object_path.Contains("/Multisims/")) {
      const bool is_hist1d = object->InheritsFrom("TH1") && !object->InheritsFrom("TH2");
      const bool is_hist2d = object->InheritsFrom("TH2");

      if (is_hist1d || is_hist2d) {
        const TString multisim_type = normalise_multisim_type(key->GetName());
        if (!multisim_type.IsNull()) {
          MultisimSummary& summary = multisim_summaries[multisim_type];
          summary.histogram_count += 1;
          summary.has_1d = summary.has_1d || is_hist1d;
          summary.has_2d = summary.has_2d || is_hist2d;

          const int universe_id = parse_universe_id(key->GetName());
          if (universe_id >= 0) {
            summary.universe_ids.insert(universe_id);
          }

          // Grab one fully-qualified example (incl. axes/binning) for each multisim type.
          const TString rel_path = make_relative_path(object_path, file_path);
          const TString example_line = format_histogram_summary_line(rel_path, object);
          if (!example_line.IsNull()) {
            if (is_hist1d && summary.example_1d.IsNull()) {
              summary.example_1d = example_line;
            }
            if (is_hist2d && summary.example_2d.IsNull()) {
              summary.example_2d = example_line;
            }
          }
        }
      }
    }

    if (is_tree_key || object->InheritsFrom(TTree::Class())) {
      TTree* tree = dynamic_cast<TTree*>(object);
      print_tree(object_path, tree);
      ++tree_count;
      continue;
    }

    if (is_directory_key || object->InheritsFrom(TDirectory::Class())) {
      TDirectory* sub_directory = dynamic_cast<TDirectory*>(object);
      scan_directory(sub_directory,
                     object_path,
                     file_path,
                     tree_count,
                     object_count,
                     class_counts,
                     directory_paths,
                     multisim_summaries,
                     non_multisim_hist_summaries);
    }
  }
}

void print_directory_summary(const std::set<TString>& directory_paths) {
  std::printf("Directories:\n");
  for (std::set<TString>::const_iterator it = directory_paths.begin(); it != directory_paths.end(); ++it) {
    std::printf("  - %s\n", it->Data());
  }
}

void print_multisim_summary(const std::map<TString, MultisimSummary>& multisim_summaries) {
  std::printf("Multisim types:\n");
  if (multisim_summaries.empty()) {
    std::printf("  - none found\n");
    return;
  }

  for (std::map<TString, MultisimSummary>::const_iterator it = multisim_summaries.begin();
       it != multisim_summaries.end();
       ++it) {
    const MultisimSummary& summary = it->second;
    std::printf("  - %s (universes=%lu, histograms=%d, 1D=%s, 2D=%s)\n",
                it->first.Data(),
                static_cast<unsigned long>(summary.universe_ids.size()),
                summary.histogram_count,
                summary.has_1d ? "yes" : "no",
                summary.has_2d ? "yes" : "no");
    if (!summary.example_1d.IsNull()) {
      std::printf("      example 1D: %s\n", summary.example_1d.Data());
    }
    if (!summary.example_2d.IsNull()) {
      std::printf("      example 2D: %s\n", summary.example_2d.Data());
    }
  }
}

void print_non_multisim_histogram_summary(const std::set<TString>& hist_summaries) {
  std::printf("Non-multisim histograms (excluding /Multisims/*):\n");
  if (hist_summaries.empty()) {
    std::printf("  - none found\n");
    return;
  }

  std::printf("  - count: %lu\n", static_cast<unsigned long>(hist_summaries.size()));
  for (std::set<TString>::const_iterator it = hist_summaries.begin(); it != hist_summaries.end(); ++it) {
    std::printf("  - %s\n", it->Data());
  }
}

void print_class_summary(const std::map<TString, int>& class_counts) {
  std::printf("Object types:\n");
  for (std::map<TString, int>::const_iterator it = class_counts.begin(); it != class_counts.end(); ++it) {
    std::printf("  - %s: %d\n", it->first.Data(), it->second);
  }
}

void print_file_tree_summary(const char* file_path) {
  if (file_path == NULL || file_path[0] == '\0') {
    std::printf("[printFluxInputTreesAndBranches] empty file path provided\n");
    return;
  }

  TFile input_file(file_path, "READ");
  if (input_file.IsZombie()) {
    std::printf("[printFluxInputTreesAndBranches] failed to open file: %s\n", file_path);
    return;
  }

  std::printf("\n============================================================\n");
  std::printf("File: %s\n", file_path);
  std::printf("============================================================\n");

  int tree_count = 0;
  int object_count = 0;
  std::map<TString, int> class_counts;
  std::set<TString> directory_paths;
  std::map<TString, MultisimSummary> multisim_summaries;
  std::set<TString> non_multisim_hist_summaries;
  scan_directory(&input_file,
                 input_file.GetName(),
                 input_file.GetName(),
                 tree_count,
                 object_count,
                 class_counts,
                 directory_paths,
                 multisim_summaries,
                 non_multisim_hist_summaries);

  std::printf("\n");
  print_directory_summary(directory_paths);
  std::printf("\n");
  print_multisim_summary(multisim_summaries);
  std::printf("\n");
  print_non_multisim_histogram_summary(non_multisim_hist_summaries);
  std::printf("\n");
  print_class_summary(class_counts);
  std::printf("\nSummary: objects=%d trees=%d\n", object_count, tree_count);

  if (tree_count == 0) {
    std::printf("No TTrees were found in this file.\n");
  }
}

} // namespace

void printFluxInputTreesAndBranches(
  const char* fhc_file = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root",
  const char* rhc_file = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root"
) {
  heron::macro::run_with_guard("printFluxInputTreesAndBranches", [&]() {
  print_file_tree_summary(fhc_file);
  print_file_tree_summary(rhc_file);

  });
}
