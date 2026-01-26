#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR="build/out/sample-lists"
ART_DIR="build/out/art"

mkdir -p "${OUTPUT_DIR}"

ls "${ART_DIR}"/art_prov_beam_s*.root > "${OUTPUT_DIR}/numi_fhc_run1_sample_beam.list"
ls "${ART_DIR}"/art_prov_strangeness_run1_fhc.root > "${OUTPUT_DIR}/numi_fhc_run1_sample_strangeness.list"
