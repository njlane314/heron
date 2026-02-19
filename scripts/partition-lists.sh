#!/usr/bin/env bash
set -euo pipefail

USER_NAME="${USER_NAME:-nlane}"
RELEASE="${RELEASE:-v08_00_00_82}"
NAME="${NAME:-numi_fhc_run1}"

BASE="/pnfs/uboone/scratch/users/${USER_NAME}/ntuples/${RELEASE}/${NAME}"
OUTPUT_DIR="${HERON_OUTPUT_DIR}/inputs"
FILE_PATTERN="${FILE_PATTERN:-nu_selection_data.root,nu_selection.root}"

mkdir -p "${OUTPUT_DIR}"

write_list() {
    local stage="$1"
    local dir="$2"
    local list="${OUTPUT_DIR}/${stage}.list"
    local scratch_list
    local pattern
    local find_args=()
    local patterns=()
    scratch_list="$(mktemp -p "${OUTPUT_DIR}" "scratch_list.XXXXXX")"

    if [[ ! -d "${dir}" ]]
    then
        echo "Warning: missing output directory: ${dir}" >&2
        rm -f "${scratch_list}"
        return 0
    fi

    IFS=',' read -r -a patterns <<< "${FILE_PATTERN}"
    for pattern in "${patterns[@]}"
    do
        if [[ -z "${pattern}" ]]
        then
            continue
        fi

        if [[ ${#find_args[@]} -gt 0 ]]
        then
            find_args+=(-o)
        fi
        find_args+=(-name "${pattern}")
    done

    if [[ ${#find_args[@]} -eq 0 ]]
    then
        echo "Warning: no file patterns supplied" >&2
        rm -f "${scratch_list}"
        return 0
    fi

    find "${dir}" -type f \( "${find_args[@]}" \) | sort > "${scratch_list}"

    if [[ ! -s "${scratch_list}" ]]
    then
        echo "Warning: no ROOT files found in ${dir}" >&2
        rm -f "${scratch_list}"
        return 0
    fi

    mv "${scratch_list}" "${list}"
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
