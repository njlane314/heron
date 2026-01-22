# Scope
This file captures the coding conventions and structural patterns used in the Nuxsec C++ codebase.

## Directory layout
- Modules live in top-level folders (io/, sample/, pot/, rdf/, sel/, ana/, syst/, stat/) with include/, src/,
  and macro/ subfolders.
- Apps/ contains small executables (apps/src).
- lib/ is reserved for build outputs (shared libraries).

## File layout & naming
- Use .hh for headers and .cc for implementation files.
- Use the nuxsec namespace for library code.
- Class and type names use descriptive PascalCase that keeps entities explicit
  (e.g., ProvenanceInput, LogicalSample, Dataset, ChannelDef, Selection, HistDef).
- Global/static variables and macros are defined in implementation files when appropriate.

## Header conventions
- Start headers with /* -- C++ -- */.
- Use concise include guards in the form #ifndef NUXSEC_<MODULE>_<NAME>_H / #define NUXSEC_<MODULE>_<NAME>_H.
- Forward declare classes when possible, include system/ROOT headers after that.
- using namespace std; appears in headers and is acceptable in this codebase.
- Group and order includes with blank lines between groups:
  - Headers (.hh): C/C++ standard library headers, then external dependencies (e.g., ROOT, sqlite3), then project headers.
  - Sources (.cc): corresponding header first, then C/C++ standard library headers, then external dependencies, then project headers.
  - Keep includes alphabetised within each group.

## Documentation & comments
- Use Doxygen-style comments:
  - /** \brief … */ for classes.
  - /// for method/field comments and short descriptions.
  - Multi-line class/method docs often include \ingroup tags and detailed text blocks.

## Formatting & style
- Brace style: opening brace on the same line as the declaration.
- Indentation mixes two spaces and tabs; preserve existing indentation style in the file being edited.
- Inline methods are commonly defined in headers.
- Prefer explicit NULL checks and early returns in implementation code.
- Use British-English spelling in code, comments, and documentation (e.g., normalisation).

## Naming conventions
- Member variables frequently use prefixes:
  - p_ for pointers (e.g., p_histograms, p_container).
  - m_ for member data (e.g., m_hist, m_name).
- Function and method names use lowercase with underscores between words (e.g., add_histogram, store_histograms, get_unique_id).

## Constants & macros
- Macro guards and simple numeric constants are defined in .cc files.
- Inline macro helpers like MIN, MAX, SQR, ABS follow the existing pattern.

## API patterns
- Prefer clear, imperative method names and explicit parameters.
- Methods often return bool for success/failure and use early-exit patterns.
- Use standard containers (std::vector, std::map) with iterators, consistent with existing code.

## Error handling
- Use simple checks and return/assert rather than exceptions.
- Keep behavior consistent with existing patterns (e.g., return false on invalid inputs).

## External dependencies
- ROOT types (e.g., TFile, TH1, TH2) are used directly in headers/implementation.
- Include ROOT headers in headers where needed.
