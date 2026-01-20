# Nuxsec

ROOT-based utilities for a neutrino cross-section analysis pipeline built around explicit analysis entities
(Fragment → Sample → Dataset/RDF → Channel/Category → Selection → Template/Plot).

## Core concepts (pipeline map)

- **LArSoft job outputs** (partitions / shards)
- **Art provenance aggregation** (stage-level manifests)
- **Logical samples** (collection of shards + merged view)
- **Exposure / POT accounting** (per shard, per logical sample, target POT scaling)
- **RDataFrame construction** (source trees, friend trees, derived columns, weights)
- **Channels / categories** (analysis topology, truth categories, control regions)
- **Selections and template specs** (cutflow, histograms, response matrices, xsec outputs)
- **Template production and plotting** (binned templates, stacked plots)

## Repository structure

This is a COLLIE-like module layout. Each module is built as its own shared library.

```
nuxsec/
  io/      # LArSoft output discovery, file manifests, provenance extraction
  ana/     # analysis-level definitions and ROOT::RDataFrame sources + derived columns
  plot/    # plotting utilities for stacked histograms and diagnostic outputs
  apps/    # unified CLI (aggregators, template makers)
  scripts/ # environment helpers
```

## Macros and executables

### What the macros are used for

Macro folders (for example, `io/macro/` and `ana/macro/`) are reserved for interactive ROOT
scripts that perform post-processing and visualisation of outputs produced by the compiled
applications. Empty `.keep` placeholders are committed so the macro directories stay tracked
in the repository. These macros are intended for quick inspection and plotting workflows, such
as reading output ROOT files, dumping diagnostic values, and generating plots for systematic
variations or fit-response checks.

### How they differ from applications

Applications are compiled executables built by the `Makefile` that perform the core Nuxsec
workflows: aggregating inputs, building RDF datasets, and producing the ROOT outputs used by
downstream analysis. Macros, in contrast, are not standalone executables—they are run
interactively inside ROOT (for example, `root -l` followed by `.L myMacro.C`) to inspect or
visualise those outputs. In short, applications do the heavy calculations, while macros focus
on analysis and plotting of the results.

## Requirements

- C++17 compiler (e.g. `g++`)
- ROOT (for `root-config` and runtime I/O)
- sqlite3 development headers/libs

## Build

```bash
make
```

This produces:

- `build/lib/libNuXsecIO.so`
- `build/lib/libNuXsecSample.so`
- `build/lib/libNuXsecAna.so`
- `build/lib/libNuXsecPlot.so`
- `build/bin/nuxsec`
- `./nuxsec` (wrapper script that runs `build/bin/nuxsec`)

## Analysis processing

The `libNuXsecAna` library provides `nuxsec::ColumnDerivationService` and RDF construction helpers
for defining analysis-level columns (weights, fiducial checks, channel classifications) on `ROOT::RDF::RNode`.

The `libNuXsecPlot` library and `nuxsec` CLI build binned template histograms from
aggregated samples and a compiled analysis definition, serving as the inputs to plotting and downstream
cross-section fits.

## Runtime environment

```bash
source .setup.sh
```

This adds `build/bin` to your `PATH` and `build/lib` to your `LD_LIBRARY_PATH` so the
executables can locate the shared libraries at runtime.

If you prefer to run from the repository root without modifying your `PATH`, use the
wrapper script:

```bash
./nuxsec <command> [args...]
```

## Prepare file lists

File lists are newline-delimited paths to ROOT files (blank lines and `#` comments are ignored).

```bash
cat > data.list <<'LIST'
# stage input files
/path/to/input1.root
/path/to/input2.root
LIST
```

## Example command-line workflow (sealed analysis, no runtime selections/tree/templates)

Assume you run from the repo top-level (`nuxsec/`) and you already have per-stage filelists
produced by your partitioning step.

### 1) Stage → Art provenance ROOT (per partition/stage)

Each stage filelist is a plain text file containing the input art/selection ROOT files for that stage.

```bash
# Data (kind can be omitted if your first file is nuselection_data.root and you rely on auto-detect)
nuxsec artio "data_bnb_run1a:inputs/filelists/data_bnb_run1a.txt:Data:BNB"

# EXT
nuxsec artio "ext_bnb_run1a:inputs/filelists/ext_bnb_run1a.txt:EXT:BNB"

# MC overlay
nuxsec artio "overlay_bnb_run1a:inputs/filelists/overlay_bnb_run1a.txt:Overlay:BNB"

# Dirt
nuxsec artio "dirt_bnb_run1a:inputs/filelists/dirt_bnb_run1a.txt:Dirt:BNB"
```

Outputs (by code convention):

```
build/out/art/art_prov_data_bnb_run1a.root
build/out/art/art_prov_ext_bnb_run1a.root
build/out/art/art_prov_overlay_bnb_run1a.root
build/out/art/art_prov_dirt_bnb_run1a.root
```

Repeat `nuxsec artio ...` for each partition/stage (run1b, run1c, etc.).

### 2) Art provenance ROOT → Sample ROOT (group stages into samples)

Create per-sample filelists that contain the art provenance ROOT outputs from step (1):

```bash
ls build/out/art/art_prov_data_bnb_run1*.root    > inputs/samples/data_bnb_run1.txt
ls build/out/art/art_prov_ext_bnb_run1*.root     > inputs/samples/ext_bnb_run1.txt
ls build/out/art/art_prov_overlay_bnb_run1*.root > inputs/samples/overlay_bnb_run1.txt
ls build/out/art/art_prov_dirt_bnb_run1*.root    > inputs/samples/dirt_bnb_run1.txt
```

Then aggregate each sample:

```bash
nuxsec sample "data_bnb_run1:inputs/samples/data_bnb_run1.txt"
nuxsec sample "ext_bnb_run1:inputs/samples/ext_bnb_run1.txt"
nuxsec sample "overlay_bnb_run1:inputs/samples/overlay_bnb_run1.txt"
nuxsec sample "dirt_bnb_run1:inputs/samples/dirt_bnb_run1.txt"
```

Outputs:

```
build/out/sample/sample_root_<sample>.root
build/out/sample/samples.tsv
```

The TSV is the only input list you pass downstream.

### 3) Samples → Analysis object (templates ROOT), sealed analysis

The compiled analysis definition in this repository is `nuxsec_default_v1` with tree
name `MyTree`. The template maker writes templates plus metadata (analysis name, tree name,
compiled template definitions, sample bookkeeping) into a single ROOT artifact.

```bash
nuxsec template-make build/out/sample/samples.tsv \
  build/out/ana/analysis_default.root 8
```

### 4) Systematics and plots

The CLI now includes a macro-driven plotting stage. The default plot macro generates a POT
timeline plot and writes images under `build/out/plot/`.

```bash
# Plot POT timeline (writes build/out/plot/pot_timeline.png)
nuxsec plot

# Override the output stem if you want a different filename or directory.
nuxsec plot build/out/plot/pot_timeline_run1
```

Fixed outputs (by code convention):

```
build/out/plot/pot_timeline.png
```

Systematics still use macros or extensions; there is no `systs` or `list-analyses` command
in this repository yet.

Shell completion for these commands is available in `scripts/nuxsec-completion.bash` (source it
in your shell profile or session).

### What the user “touches” in this UX

- Filelists for stages (already produced by your partitioning step).
- Per-sample grouping filelists (or a small helper script).
- `samples.tsv` is produced automatically and becomes the stable handoff.
- After that, users only run:
  - `nuxsec template-make …` (required)
  - `nuxsec plot …` (optional, macro-driven plotting outputs).
