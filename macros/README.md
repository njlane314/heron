# Macro standardisation guide

This directory uses thin ROOT macros as an interface layer and keeps reusable behaviour in compiled backend services.

## Macro template (required for new macros)

Use this layout for new plotting macros:

1. Resolve default input path (`default_event_list_root`).
2. Validate event-list input (`looks_like_event_list_root`).
3. Build event-list RDF via `EventListIO`.
4. Build sample-split nodes with:
   - `filter_by_sample_mask`
   - `filter_not_sample_mask`
5. Keep macro-level logic focused on plot configuration and output naming only.

```cpp
int plotMySummary(const std::string &event_list_path = "")
{
    const std::string input = event_list_path.empty() ? default_event_list_root() : event_list_path;
    if (!looks_like_event_list_root(input))
    {
        std::cerr << "[plotMySummary] invalid event list root input: " << input << "\n";
        return 1;
    }

    EventListIO el(input);
    ROOT::RDataFrame rdf = el.rdf();

    auto mask_mc = el.mask_for_mc_like();
    auto mask_ext = el.mask_for_ext();

    ROOT::RDF::RNode node_mc = filter_by_sample_mask(rdf, mask_mc);
    node_mc = filter_not_sample_mask(node_mc, mask_ext);

    // Plot-specific logic only.
    return 0;
}
```

## Style rule for macro interfaces

- Keep macro function signatures small and explicit.
- Prefer one macro entry function per file.
- Avoid duplicating generic helpers that already exist in `framework/plot` and `framework/io`.
- Prefer imperative method names and early returns for validation failures.

## Standardised manifest entrypoints

A logical macro command should be added to `macros/macro/library/manifest.tsv` so users can run a stable alias even if a file moves.

## Legacy migration policy

Migration should proceed in batches:

1. Migrate macros that duplicate input validation and mask filtering first.
2. Migrate one macro per category (cut-flow, efficiency, stacked histogram) before broader rollout.
3. Keep output names and selection semantics unchanged while migrating.
4. Mark migrated macros by using backend helpers directly in their implementation.
