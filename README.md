# heron — HERON — Histogram and Event Relay for Orchestrated Normalisation

Documentation: [https://heron.github.io/heron/](https://njlane314.github.io/heron/)

ROOT-based utilities for a neutrino cross-section analysis pipeline. The codebase formalises the analysis
entities (aggregation → sample → dataset/RDF → channel/category → selection → template/plot) and provides
compiled drivers to move between them.

The core design goal is to keep analysis structure explicit: inputs are aggregated into samples with
recorded provenance, and samples feed into a compiled event-level analysis. This makes it easier to audit
which inputs contributed to each step and to reproduce analysis or training artefacts later.

## Architecture 

Each top-level module builds a shared library and exposes typed services or analysis entities. The pipeline
also supports producing a CNN training snapshot: sample-level ROOT outputs can be materialised, split, and
saved as an offline training set before event-level aggregation or plotting.

```
framework/io/    LArSoft output discovery, file manifests, provenance, and run databases
framework/ana/   analysis configuration, selection logic, and ROOT::RDataFrame column derivations
framework/plot/  stacked-histogram and channel plotting helpers
framework/apps/  CLI entrypoints that orchestrate the pipeline
framework/evd/   event display utilities
```


### Runtime 

- `$HERON_OUTPUT_DIR/art/` stores provenance ROOT outputs from `heron art` (the variable must be set).
- `scratch/out/<set>/sample/` stores per-sample ROOT outputs and `samples.tsv` produced by `heron sample`.
- `scratch/out/<set>/event/` stores event-level ROOT outputs produced by `heron event`.

The `<set>` segment defaults to `out` and is controlled by `HERON_SET` or `heron --set`.

## Requirements

- C++17 compiler (e.g. `g++`)
- ROOT (for `root-config` and runtime I/O)
- sqlite3 development headers/libs

## Build

```bash
make all -j12
```

This produces shared libraries and drivers:

- `build/lib/libheronIO.so`
- `build/lib/libheronAna.so`
- `build/lib/libheronPlot.so`
- `build/bin/heron`
- `build/bin/heronArtFileIOdriver`
- `build/bin/heronSampleIOdriver`
- `build/bin/heronEventIOdriver`
- `./heron` (wrapper script that runs `build/bin/heron`)

## CLI Overview

```bash
heron -h
```

```text
Usage: heron <command> [args]

Commands:
  art         Aggregate art provenance for an input
  sample      Aggregate Sample ROOT files from art provenance
  event       Build event-level output from aggregated samples
  paths       Print resolved workspace paths
  env         Print environment exports for a workspace

Global options:
  -S, --set   Workspace selector (default: out)

Run 'heron <command> --help' for command-specific usage.
```

## Run Environment

```bash
source .container.sh
source .setup.sh
```

This adds `build/bin` to your `PATH` and `build/lib` to your `LD_LIBRARY_PATH` so the
executables can locate the shared libraries at runtime.

If you prefer to run from the repository root without modifying your `PATH`, use the
wrapper script:

```bash
./heron <command> [args...]
```

### Environment Variables

- `HERON_SET` selects the active workspace (default: `out`).
- `HERON_OUT_BASE` overrides the base output directory; if unset, `HERON_OUTPUT_DIR` is used before falling back to `<repo>/scratch/out`.
- Temporary snapshot staging is written to `/exp/uboone/data/users/$USER/staging`; `USER` must be set.
- `HERON_OUTPUT_DIR` is required by `heron art`; outputs are written to `$HERON_OUTPUT_DIR/art`.
- `HERON_SAMPLE_DIR` and `HERON_EVENT_DIR` override per-stage output directories for `sample` and `event`.
- `HERON_REPO_ROOT` can be set to override the repo discovery used by the CLI.
- `HERON_TREE_NAME` selects the input tree name for the event builder (default: `Events`).

## Input Files

File lists are newline-delimited paths to ROOT files (blank lines and `#` comments are ignored):

```bash
cat > data.list <<'LIST'
# input file list
/path/to/input1.root
/path/to/input2.root
LIST
```

## Minimal Workflow

Assume you run from the repo root and already have per-input filelists from your partitioning step.
Choose a workspace either by exporting `HERON_SET` or using `heron --set` in each command.

1) **Input → art provenance ROOT (per partition/input)**

```bash
heron --set template art "sample_a:inputs/filelists/sample_a.txt:Data:Beam"
heron --set template art "sample_b:inputs/filelists/sample_b.txt:EXT:Beam"
heron --set template art "sample_c:inputs/filelists/sample_c.txt:Overlay:Beam"
heron --set template art "sample_d:inputs/filelists/sample_d.txt:Dirt:Beam"
```

Outputs (by code convention):

```
scratch/out/template/art/art_prov_sample_a.root
scratch/out/template/art/art_prov_sample_b.root
scratch/out/template/art/art_prov_sample_c.root
scratch/out/template/art/art_prov_sample_d.root
```

2) **Art provenance ROOT → sample ROOT (group inputs into samples)**

Create per-sample filelists containing the art provenance outputs from step (1):

```bash
mkdir -p scratch/out/lists
ls scratch/out/template/art/art_prov_sample_a*.root > scratch/out/lists/sample_a.txt
ls scratch/out/template/art/art_prov_sample_b*.root > scratch/out/lists/sample_b.txt
ls scratch/out/template/art/art_prov_sample_c*.root > scratch/out/lists/sample_c.txt
ls scratch/out/template/art/art_prov_sample_d*.root > scratch/out/lists/sample_d.txt
```

Then aggregate each sample:

```bash
heron --set template sample "sample_a:scratch/out/lists/sample_a.txt"
heron --set template sample "sample_b:scratch/out/lists/sample_b.txt"
heron --set template sample "sample_c:scratch/out/lists/sample_c.txt"
heron --set template sample "sample_d:scratch/out/lists/sample_d.txt"
```

Outputs:

```
scratch/out/template/sample/sample_root_<sample>.root
scratch/out/template/sample/samples.tsv
```

Use the resulting `samples.tsv` downstream for event-level aggregation and plotting.

**Training vs template workspaces**

Keep train/template outputs separated by selecting the workspace instead of moving files:

```bash
heron --set template sample "sample_a:scratch/out/lists/sample_a.txt"
heron --set template event scratch/out/template/event/events.root

heron --set train sample "sample_a:scratch/out/lists/sample_a.txt"
heron --set train macro plotTrainingQA.C
```

Use `heron paths` to print resolved locations or `eval "$(heron env train)"` to switch a shell.

3) **Samples → event-level output (compiled analysis)**

The compiled analysis definition in this repository is `heron_default` with tree
name `Events` by default. Override the input tree name by exporting `HERON_TREE_NAME`
before running the CLI. The event builder writes a single ROOT file containing the
event-level tree plus metadata for the aggregated samples.
If you provide only the output path, the CLI uses the active workspace's
`samples.tsv` automatically. Provide a selection expression as the final
argument to filter events before writing the output.

```bash
heron --set template event scratch/out/template/event/events.root
```

To override the event output schema, pass a columns TSV as the final positional
argument. The TSV expects `type` and `name` columns (see `framework/apps/config/event_columns.tsv`).
If you only want to provide columns, pass `true` as the selection placeholder.

```bash
heron --set template event scratch/out/template/event/events.root true framework/apps/config/event_columns.tsv
```

Selection strings can reference selection columns derived by the SelectionService.
Examples:

```bash
heron --set template event scratch/out/template/event/events.root sel_muon
heron --set template event scratch/out/template/event/events.root sel_inclusive_mu_cc
heron --set template event scratch/out/template/event/events.root sel_triggered_muon
heron --set template event scratch/out/template/event/events.root sel_triggered_slice
heron --set template event scratch/out/template/event/events.root sel_reco_fv
```

4) **Plotting via macros**

Macro support and the in-repository plotting macro collection are temporarily removed from the CLI path.
For standalone ROOT usage, see `framework/macro/AnalysisModelExample.C`, which demonstrates loading
a `samples.tsv` list via `SampleIO`, calling `ExecutionPolicy` first, constructing an
`AnalysisContext`, and then initialising a user-defined `AnalysisModel` subclass.
