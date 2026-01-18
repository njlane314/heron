# Nuxsec

ROOT-based provenance aggregation utilities for NuIO samples.

## Requirements

- C++17 compiler (e.g. `g++`)
- ROOT (for `root-config` and runtime I/O)
- sqlite3 development headers/libs

## Build the shared libraries

```bash
make
```

This produces:

- `build/lib/libNuIO.so`
- `build/lib/libNuAna.so`

## Analysis processing

The NuAna library provides `nuana::NuAnalysisProcessor` for defining analysis-level
columns (weights, fiducial checks, channel classifications) on `ROOT::RDF::RNode`
instances used by `sampleRDFbuilder`.

## Runtime environment

```bash
source scripts/setup.sh
```

## Build the command-line tools

```bash
g++ -std=c++17 -O2 -Wall -Wextra \
  $(root-config --cflags) \
  -I./lib/NuIO/include \
  -L./build/lib -lNuIO \
  -lsqlite3 $(root-config --libs) \
  -o artIOaggregator \
  bin/artIOaggregator/artIOaggregator.cxx

g++ -std=c++17 -O2 -Wall -Wextra \
  $(root-config --cflags) \
  -I./lib/NuIO/include \
  -L./build/lib -lNuIO \
  -lsqlite3 $(root-config --libs) \
  -o sampleIOaggregator \
  bin/sampleIOaggregator/sampleIOaggregator.cxx
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
./artIOaggregator my_stage:data.list
# writes ./ArtIO_my_stage.root
```

```bash
./sampleIOaggregator my_sample:data.list
# writes ./SampleIO_my_sample.root
# updates ./SampleIO_samples.tsv
```
