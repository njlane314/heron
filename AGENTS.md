# Repository Formatting Rules

All source files in this repository must follow the formatting shown in `nuaggregator.cxx` and the example below.

## Required style

- Every `.cxx`, `.cc`, `.cpp`, `.c`, and `.h/.hpp` file must start with a Doxygen file header of the form:
  ```
  /**
   *  @file  <relative path>
   *
   *  @brief <short description>
   */
  ```
- Use Allman braces (opening brace on its own line).
- Indent with 4 spaces.
- Keep namespace braces on their own lines.
- Align pointer and reference symbols to the left (e.g., `const std::string &name`).
- Prefer one statement per line and avoid compacting blocks into single-line statements.

## Suggested clang-format

If you use `clang-format`, apply the following style:

```
{BasedOnStyle: LLVM, IndentWidth: 4, BreakBeforeBraces: Allman, AllowShortFunctionsOnASingleLine: None, ColumnLimit: 0, PointerAlignment: Left, ReferenceAlignment: Left, NamespaceIndentation: None}
```
