# nuxsec

ROOT-based utilities for a neutrino cross-section analysis pipeline. The codebase formalises the analysis
entities (aggregation → sample → dataset/RDF → channel/category → selection → template/plot) and provides
compiled drivers to move between them.

The core design goal is to keep analysis structure explicit: inputs are aggregated into samples with
recorded provenance, samples feed into a compiled event-level analysis, and plots/macros operate on the
resulting outputs. This makes it easier to audit which inputs contributed to each step and to reproduce
analysis or training artefacts later.

## Architecture 

Each top-level module builds a shared library and exposes typed services or analysis entities. The pipeline
also supports producing a CNN training snapshot: sample-level ROOT outputs can be materialised, split, and
saved as an offline training set before event-level aggregation or plotting.

```
io/    LArSoft output discovery, file manifests, provenance, and run databases
ana/   analysis configuration, selection logic, and ROOT::RDataFrame column derivations
plot/  stacked-histogram and channel plotting helpers
apps/  CLI entrypoints that orchestrate the pipeline
```

### Runtime 

- `scratch/out/<set>/art/` stores provenance ROOT outputs from `nuxsec art`.
- `scratch/out/<set>/sample/` stores per-sample ROOT outputs and `samples.tsv` produced by `nuxsec sample`.
- `scratch/out/<set>/event/` stores event-level ROOT outputs produced by `nuxsec event`.
- `scratch/plot/<set>/` stores plot outputs produced by `nuxsec macro` (configurable via `NUXSEC_PLOT_DIR`).

The `<set>` segment defaults to `template` and is controlled by `NUXSEC_SET` or `nuxsec --set`.

## Requirements

- C++17 compiler (e.g. `g++`)
- ROOT (for `root-config` and runtime I/O)
- sqlite3 development headers/libs

## Build

```bash
make all -j12
```

This produces shared libraries and drivers:

- `build/lib/libNuxsecIO.so`
- `build/lib/libNuxsecAna.so`
- `build/lib/libNuxsecPlot.so`
- `build/bin/nuxsec`
- `build/bin/nuxsecArtFileIOdriver`
- `build/bin/nuxsecSampleIOdriver`
- `build/bin/nuxsecEventIOdriver`
- `./nuxsec` (wrapper script that runs `build/bin/nuxsec`)

## CLI Overview

```bash
nuxsec -h
```

```text
Usage: nuxsec <command> [args]

Commands:
  art         Aggregate art provenance for an input
  sample      Aggregate Sample ROOT files from art provenance
  event       Build event-level output from aggregated samples
  macro       Run plot macros
  paths       Print resolved workspace paths
  env         Print environment exports for a workspace

Global options:
  -S, --set   Workspace selector (default: template)

Run 'nuxsec <command> --help' for command-specific usage.
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
./nuxsec <command> [args...]
```

### Environment Variables

- `NUXSEC_SET` selects the active workspace (default: `template`).
- `NUXSEC_OUT_BASE` overrides the base output directory (default: `<repo>/scratch/out`).
- `NUXSEC_PLOT_BASE` overrides the plot base directory (default: `<repo>/scratch/plot`).
- `NUXSEC_ART_DIR`, `NUXSEC_SAMPLE_DIR`, and `NUXSEC_EVENT_DIR` override per-stage output directories.
- `NUXSEC_PLOT_DIR` and `NUXSEC_PLOT_FORMAT` control plot output location and file extension.
- `NUXSEC_REPO_ROOT` can be set to override the repo discovery used by the CLI.
- `NUXSEC_TREE_NAME` selects the input tree name for the event builder (default: `Events`).

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
Choose a workspace either by exporting `NUXSEC_SET` or using `nuxsec --set` in each command.

1) **Input → art provenance ROOT (per partition/input)**

```bash
nuxsec --set template art "sample_a:inputs/filelists/sample_a.txt:Data:Beam"
nuxsec --set template art "sample_b:inputs/filelists/sample_b.txt:EXT:Beam"
nuxsec --set template art "sample_c:inputs/filelists/sample_c.txt:Overlay:Beam"
nuxsec --set template art "sample_d:inputs/filelists/sample_d.txt:Dirt:Beam"
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
nuxsec --set template sample "sample_a:scratch/out/lists/sample_a.txt"
nuxsec --set template sample "sample_b:scratch/out/lists/sample_b.txt"
nuxsec --set template sample "sample_c:scratch/out/lists/sample_c.txt"
nuxsec --set template sample "sample_d:scratch/out/lists/sample_d.txt"
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
nuxsec --set template sample "sample_a:scratch/out/lists/sample_a.txt"
nuxsec --set template event scratch/out/template/event/events.root

nuxsec --set train sample "sample_a:scratch/out/lists/sample_a.txt"
nuxsec --set train macro plotTrainingQA.C
```

Use `nuxsec paths` to print resolved locations or `eval "$(nuxsec env train)"` to switch a shell.

3) **Samples → event-level output (compiled analysis)**

The compiled analysis definition in this repository is `nuxsec_default` with tree
name `Events` by default. Override the input tree name by exporting `NUXSEC_TREE_NAME`
before running the CLI. The event builder writes a single ROOT file containing the
event-level tree plus metadata for the aggregated samples.
If you provide only the output path, the CLI uses the active workspace's
`samples.tsv` automatically. Provide a selection expression as the final
argument to filter events before writing the output.

```bash
nuxsec --set template event scratch/out/template/event/events.root
```

Example with an explicit selection:

```bash
nuxsec --set template event scratch/out/template/event/events.root "sel_triggered_muon"
```

Muon-neutrino selection example:

```bash
nuxsec --set template event scratch/out/template/event/events.root "sel_inclusive_mu_cc"
```

4) **Plotting via macros**

Plotting is macro-driven. Use the `nuxsec macro` helper to run a plot macro
(and optionally a specific function inside it).

```bash
nuxsec --set template macro plotPotSimple.C
```

Shell completion for these commands is available in `scripts/nuxsec-completion.bash` (source it
in your shell profile or session).
