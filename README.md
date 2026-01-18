# Nuxsec

ROOT-based utilities for a neutrino cross-section analysis pipeline built around explicit analysis entities
(Fragment → Sample → Dataset → Channel/Category → Selection → Product).

## Core concepts (pipeline map)

- **LArSoft job outputs** (partitions / shards)
- **Logical samples** (collection of shards + merged view)
- **Exposure / POT accounting** (per shard, per logical sample, target POT scaling)
- **RDataFrame construction** (source trees, friend trees, derived columns, weights)
- **Channels / categories** (analysis topology, truth categories, control regions)
- **Selections and products** (cutflow, histograms, response matrices, xsec outputs)

## Repository structure

This is a COLLIE-like module layout. Each module is built as its own shared library.

```
nuxsec/
  io/      # LArSoft output discovery, file manifests, provenance extraction
  ana/     # analysis-level definitions and ROOT::RDataFrame sources + derived columns
  apps/    # small CLIs (aggregators, RDF builders)
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

- `build/lib/libNuxsecIO.so`
- `build/lib/libNuxsecSample.so`
- `build/lib/libNuxsecAna.so`
- `build/bin/nuxsecArtIOaggregator`
- `build/bin/nuxsecSampleIOaggregator`
- `build/bin/nuxsecSampleRDFbuilder`

## Analysis processing

The `libNuXsecAna` library provides `nuxsec::AnalysisRdfDefinitions` and RDF construction helpers
for defining analysis-level columns (weights, fiducial checks, channel classifications) on `ROOT::RDF::RNode`
instances used by `nuxsecSampleRDFbuilder`.

## Runtime environment

```bash
source scripts/setup.sh
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

## Run the aggregators

```bash
build/bin/nuxsecArtIOaggregator my_stage:data.list
# writes ./ArtFileProvenance_my_stage.root
```

```bash
build/bin/nuxsecSampleIOaggregator my_sample:data.list
# writes ./SampleRootIO_my_sample.root
# updates ./SampleRootIO_samples.tsv
```
