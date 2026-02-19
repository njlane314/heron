#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR="${HERON_OUTPUT_DIR}/sample-lists"
ART_DIR="${HERON_OUTPUT_DIR}/art"

if [[ ! -d "${ART_DIR}" && -n "${HERON_SET:-}" && -d "${HERON_OUTPUT_DIR}/${HERON_SET}/art" ]]
then
    ART_DIR="${HERON_OUTPUT_DIR}/${HERON_SET}/art"
fi

mkdir -p "${OUTPUT_DIR}"
shopt -s nullglob

write_list() {
    local pattern="$1"
    local list="${OUTPUT_DIR}/$2"
    local matches=()

    matches=("${ART_DIR}"/${pattern})

    if [[ ${#matches[@]} -eq 0 ]]
    then
        echo "Warning: no art files matched ${pattern} in ${ART_DIR}" >&2
        return 0
    fi

    printf "%s\n" "${matches[@]}" > "${list}"
}

if [[ ! -d "${ART_DIR}" ]]
then
    echo "Warning: missing art directory: ${ART_DIR}" >&2
    exit 0
fi

write_list "art_prov_beam_s*.root" "numi_fhc_run1_sample_beam.list"
write_list "art_prov_strangeness_run1_fhc.root" "numi_fhc_run1_sample_strangeness.list"

shopt -u nullglob
