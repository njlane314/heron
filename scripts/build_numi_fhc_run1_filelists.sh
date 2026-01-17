#!/usr/bin/env bash
set -euo pipefail

USER_NAME="${USER_NAME:-nlane}"
RELEASE="${RELEASE:-v08_00_00_82}"
NAME="${NAME:-numi_fhc_run1}"

BASE="/pnfs/uboone/scratch/users/${USER_NAME}/ntuples/${RELEASE}/${NAME}"
OUTPUT_DIR="inputs/numi_fhc_run1/filelists"

mkdir -p "${OUTPUT_DIR}"

write_list() {
    local stage="$1"
    local dir="$2"
    local list="${OUTPUT_DIR}/${stage}.list"

    if [[ ! -d "${dir}" ]]
    then
        echo "Missing output directory: ${dir}" >&2
        return 1
    fi

    find "${dir}" -maxdepth 1 -type f -name "*.root" | sort > "${list}"

    if [[ ! -s "${list}" ]]
    then
        echo "Warning: no ROOT files found in ${dir}" >&2
    fi
}

write_list "beam_s0" "${BASE}/beam/s0/out"
write_list "beam_s1" "${BASE}/beam/s1/out"
write_list "beam_s2" "${BASE}/beam/s2/out"
write_list "beam_s3" "${BASE}/beam/s3/out"
write_list "dirt_s0" "${BASE}/dirt/s0/out"
write_list "dirt_s1" "${BASE}/dirt/s1/out"
write_list "strangeness_run1_fhc" "${BASE}/strangeness/run1_fhc/out"
write_list "ext_p0" "${BASE}/ext/p0/out"
write_list "ext_p1" "${BASE}/ext/p1/out"
write_list "ext_p2" "${BASE}/ext/p2/out"
