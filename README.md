# Nuxsec

## Planned features

### Shared libraries

- `libArtIO.so`
  - ArtIO.root data model + ROOT dictionaries + I/O helpers (e.g. ArtProvenance, ArtProvenanceIO).
- `libNuIO.so`
  - NuIO bank schemas + dictionaries + manifest/loader helpers (friend-tree assembly).
- `libNuSel.so`
  - Selection definitions evaluated on NuIO.Evt variables; produces/consumes NuIO.Sel.
- `libNuWgt.so`
  - Weight-pack extraction/standardisation; produces/consumes NuIO.Wgt.
- `libNuAna.so`
  - Histogram definitions and filling utilities used by calculators (nominal and systematics share identical channel/bin definitions).
- `libNuSyst.so`
  - Universe loops, error bands, covariance builders.
- `libNuPlot.so`
  - Plotting style and rendering helpers.

### Typical workflow

1. Aggregate LArSoft outputs.
   - Run `artIOaggregator` to create `ArtIO.root`.
2. Build the event store.
   - Run `nuIOcondenser.exe` to assemble `NuIO.*.root` banks.
3. Run nominal analysis.
   - Run `nuAnaCalc.exe` to write/update `NuResults.root` (`Nominal/`).
4. Run systematics later (optional).
   - If not already present, generate `NuIO.Wgt` via `nuIOcondenser.exe`.
   - Run `nuSystCalc.exe` to append `Syst/` products into the existing `NuResults.root`.
5. Plot.
   - Run `nuPlot.exe` from `NuResults.root`.
