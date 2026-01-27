# nuxsec

ROOT-based utilities for a neutrino cross-section analysis pipeline. The codebase formalises the analysis
entities (aggregation → sample → dataset/RDF → channel/category → selection → template/plot) and provides
compiled drivers to move between them.

The core design goal is to keep analysis structure explicit: inputs are aggregated into samples with
recorded provenance, samples feed into a compiled event-level analysis, and plots/macros operate on the
resulting outputs. This makes it easier to audit which inputs contributed to each step and to reproduce
analysis or training artefacts later.

## Architecture at a glance

Each top-level module builds a shared library and exposes typed services or analysis entities. The pipeline
also supports producing a CNN training snapshot: sample-level ROOT outputs can be materialised, split, and
saved as an offline training set before event-level aggregation or plotting.

```
io/    LArSoft output discovery, file manifests, provenance, and run databases
ana/   analysis configuration, selection logic, and ROOT::RDataFrame column derivations
plot/  stacked-histogram and channel plotting helpers
apps/  CLI entrypoints that orchestrate the pipeline
```

### Key library surfaces (headers)

**IO services**
- `ArtFileProvenanceIO` for art provenance ingestion and ROOT output management.【F:io/include/ArtFileProvenanceIO.hh†L67-L172】
- `SampleIO` for sample aggregation metadata (origin, beam mode, tree names) used by both analysis and
  CNN-training sample materialisation workflows.【F:io/include/SampleIO.hh†L20-L117】
- `EventIO` for event-level ROOT I/O, including open modes and metadata handling.【F:io/include/EventIO.hh†L41-L170】
- `RunDatabaseService` and `SubRunInventoryService` for run/subrun tracking and SQLite-backed state.【F:io/include/RunDatabaseService.hh†L37-L171】【F:io/include/SubRunInventoryService.hh†L23-L146】
- `NormalisationService` for exposure and POT normalisation handling.【F:io/include/NormalisationService.hh†L20-L109】

**Analysis services**
- `AnalysisConfigService` for analysis configuration and on-disk metadata.【F:ana/include/AnalysisConfigService.hh†L21-L126】
- `Selection` for declarative selection/plot definitions, presets, and categorisation logic.【F:ana/include/Selection.hh†L33-L195】
- `RDataFrameService` for RDF sources and derived column orchestration.【F:ana/include/RDataFrameService.hh†L28-L91】
- `ColumnDerivationService` for channel-aware column additions and derived variables.【F:ana/include/ColumnDerivationService.hh†L17-L110】

**Plotting utilities**
- `Plotter` and `StackedHist` for channel-aware stacked histograms and outputs.【F:plot/include/Plotter.hh†L21-L94】【F:plot/include/StackedHist.hh†L30-L88】
- `Channels` and `PlotDescriptors` for plot metadata, cut direction, and visual configuration.【F:plot/include/PlotChannels.hh†L23-L78】【F:plot/include/PlotDescriptors.hh†L20-L100】

### Runtime artefacts (by convention)

- `build/out/art/` stores provenance ROOT outputs from `nuxsec art`.
- `build/out/sample/` stores per-sample ROOT outputs and `samples.tsv` produced by `nuxsec sample`.
- `build/out/event/` stores event-level ROOT outputs produced by `nuxsec event`.
- `build/plot/` stores plot outputs produced by `nuxsec macro` (configurable via `NUXSEC_PLOT_DIR`).

## Requirements

- C++17 compiler (e.g. `g++`)
- ROOT (for `root-config` and runtime I/O)
- sqlite3 development headers/libs

## Build

```bash
source .container.sh
source .setup.sh
```

```bash
make
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

## CLI overview

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

Run 'nuxsec <command> --help' for command-specific usage.
```

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

### Environment knobs

- `NUXSEC_REPO_ROOT` can be set to override the repo discovery used by the CLI.
- `NUXSEC_TREE_NAME` selects the input tree name for the event builder (default: `Events`).
- `NUXSEC_PLOT_DIR` and `NUXSEC_PLOT_FORMAT` control plot output location and file extension.

## Input file lists

File lists are newline-delimited paths to ROOT files (blank lines and `#` comments are ignored):

```bash
cat > data.list <<'LIST'
# input file list
/path/to/input1.root
/path/to/input2.root
LIST
```

## Minimal workflow (analysis + CNN training snapshots)

Assume you run from the repo root and already have per-input filelists from your partitioning step.

1) **Input → art provenance ROOT (per partition/input)**

```bash
nuxsec art "sample_a:inputs/filelists/sample_a.txt:Data:Beam"
nuxsec art "sample_b:inputs/filelists/sample_b.txt:EXT:Beam"
nuxsec art "sample_c:inputs/filelists/sample_c.txt:Overlay:Beam"
nuxsec art "sample_d:inputs/filelists/sample_d.txt:Dirt:Beam"
```

Outputs (by code convention):

```
build/out/art/art_prov_sample_a.root
build/out/art/art_prov_sample_b.root
build/out/art/art_prov_sample_c.root
build/out/art/art_prov_sample_d.root
```

2) **Art provenance ROOT → sample ROOT (group inputs into samples)**

Create per-sample filelists containing the art provenance outputs from step (1):

```bash
mkdir -p build/out/lists
ls build/out/art/art_prov_sample_a*.root > build/out/lists/sample_a.txt
ls build/out/art/art_prov_sample_b*.root > build/out/lists/sample_b.txt
ls build/out/art/art_prov_sample_c*.root > build/out/lists/sample_c.txt
ls build/out/art/art_prov_sample_d*.root > build/out/lists/sample_d.txt
```

Then aggregate each sample:

```bash
nuxsec sample "sample_a:build/out/lists/sample_a.txt"
nuxsec sample "sample_b:build/out/lists/sample_b.txt"
nuxsec sample "sample_c:build/out/lists/sample_c.txt"
nuxsec sample "sample_d:build/out/lists/sample_d.txt"
```

Outputs:

```
build/out/sample/sample_root_<sample>.root
build/out/sample/samples.tsv
```

The `build/out/sample/` artefacts can feed event-level aggregation or be snapshotted for
CNN training (for example, copying the per-sample ROOT outputs into a dedicated training
dataset directory alongside a curated `samples.tsv`).

**Training vs template sample sets (recommended handoff)**

Maintain two disjoint sample aggregations and stage them into separate directories so the
training set never overlaps the template/plotting set.

```bash
mkdir -p build/out/lists_train build/out/lists_template

# Train lists (only training partitions)
for origin in sample_a sample_b sample_c sample_d; do
  ls build/out/art/art_prov_${origin}_train*.root > build/out/lists_train/${origin}_train.txt
done

# Template lists (disjoint from training)
for origin in sample_a sample_b sample_c sample_d; do
  ls build/out/art/art_prov_${origin}_template*.root > build/out/lists_template/${origin}_template.txt
done
```

Aggregate each set and archive the outputs:

```bash
# Training samples/TSV
for origin in sample_a sample_b sample_c sample_d; do
  nuxsec sample "${origin}_train:build/out/lists_train/${origin}_train.txt"
done
mkdir -p build/out/sample_train
mv build/out/sample/sample_root_* build/out/sample/samples.tsv build/out/sample_train/

# Template samples/TSV
for origin in sample_a sample_b sample_c sample_d; do
  nuxsec sample "${origin}_template:build/out/lists_template/${origin}_template.txt"
done
mkdir -p build/out/sample_template
mv build/out/sample/sample_root_* build/out/sample/samples.tsv build/out/sample_template/
```

Use the training snapshot for CNN workflows, and pass
`build/out/sample_template/samples.tsv` downstream for event-level aggregation and plotting.

3) **Samples → event-level output (compiled analysis)**

The compiled analysis definition in this repository is `nuxsec_default` with tree
name `Events` by default. Override the input tree name by exporting `NUXSEC_TREE_NAME`
before running the CLI. The event builder writes a single ROOT file containing the
event-level tree plus metadata for the aggregated samples.

```bash
nuxsec event build/out/sample_template/samples.tsv build/out/event/events.root
```

4) **Plotting via macros**

Plotting is macro-driven. Use the `nuxsec macro` helper to run a plot macro
(and optionally a specific function inside it).

```bash
nuxsec macro plotPotSimple.C
```

Shell completion for these commands is available in `scripts/nuxsec-completion.bash` (source it
in your shell profile or session).
