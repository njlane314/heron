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
instances used by `sampleRDFbuilder`.

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
