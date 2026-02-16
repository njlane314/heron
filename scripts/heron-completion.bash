#!/usr/bin/env bash

_heron()
{
  local cur prev
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  local commands="art sample event macro paths env help -h --help"

  _heron_find_root()
  {
    local dir
    local base
    local parent
    local exe
    local exe_dir

    if [[ -n "${HERON_REPO_ROOT:-}" && ( -d "${HERON_REPO_ROOT}/plot/macro" || -d "${HERON_REPO_ROOT}/evd/macro" || -d "${HERON_REPO_ROOT}/standalone/macro" || -d "${HERON_REPO_ROOT}/io/macro" ) ]]; then
      printf "%s" "${HERON_REPO_ROOT}"
      return 0
    fi

    if [[ -n "${HERON_ROOT:-}" && ( -d "${HERON_ROOT}/plot/macro" || -d "${HERON_ROOT}/evd/macro" || -d "${HERON_ROOT}/standalone/macro" || -d "${HERON_ROOT}/io/macro" ) ]]; then
      printf "%s" "${HERON_ROOT}"
      return 0
    fi

    exe="$(command -v heron 2>/dev/null || true)"
    if [[ -n "${exe}" ]]; then
      exe_dir="$(dirname "$(readlink -f "${exe}" 2>/dev/null || printf "%s" "${exe}")")"
      dir="${exe_dir}"
      while [[ -n "${dir}" && "${dir}" != "/" ]]; do
        if [[ -d "${dir}/plot/macro" || -d "${dir}/evd/macro" || -d "${dir}/standalone/macro" || -d "${dir}/io/macro" ]]; then
          printf "%s" "${dir}"
          return 0
        fi
        dir="$(dirname "${dir}")"
      done
    fi

    dir="${PWD}"
    while [[ "${dir}" != "/" ]]; do
      if [[ -d "${dir}/plot/macro" || -d "${dir}/evd/macro" || -d "${dir}/standalone/macro" || -d "${dir}/io/macro" ]]; then
        printf "%s" "${dir}"
        return 0
      fi
      base="$(basename "${dir}")"
      if [[ "${base}" == "plot" && -d "${dir}/macro" ]]; then
        printf "%s" "$(dirname "${dir}")"
        return 0
      fi
      if [[ "${base}" == "evd" && -d "${dir}/macro" ]]; then
        printf "%s" "$(dirname "${dir}")"
        return 0
      fi
      if [[ "${base}" == "io" && -d "${dir}/macro" ]]; then
        printf "%s" "$(dirname "${dir}")"
        return 0
      fi
      if [[ "${base}" == "macro" ]]; then
        parent="$(dirname "${dir}")"
        if [[ "$(basename "${parent}")" == "plot" || "$(basename "${parent}")" == "evd" || "$(basename "${parent}")" == "io" ]]; then
          printf "%s" "$(dirname "${parent}")"
          return 0
        fi
        if [[ "$(basename "${parent}")" == "standalone" ]]; then
          printf "%s" "$(dirname "${parent}")"
          return 0
        fi
      fi
      dir="$(dirname "${dir}")"
    done
    return 1
  }

  _heron_list_macros()
  {
    local repo_root
    local macro_dir
    local macro
    local evd_dir

    repo_root="$(_heron_find_root 2>/dev/null || true)"
    if [[ -n "${repo_root}" ]]; then
      macro_dir="${repo_root}/plot/macro"
      if [[ -d "${macro_dir}" ]]; then
        for macro in "${macro_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      macro_dir="${repo_root}/standalone/macro"
      if [[ -d "${macro_dir}" ]]; then
        for macro in "${macro_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      evd_dir="${repo_root}/evd/macro"
      if [[ -d "${evd_dir}" ]]; then
        for macro in "${evd_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      macro_dir="${repo_root}/io/macro"
      if [[ -d "${macro_dir}" ]]; then
        for macro in "${macro_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      return 0
    fi

    heron macro list 2>/dev/null | awk '/\.C$/ {print $1}'
  }

  if [[ ${COMP_CWORD} -le 1 ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${prev}" == "help" || "${prev}" == "-h" || "${prev}" == "--help" ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${COMP_WORDS[1]}" == "macro" ]]; then
    local macros
    macros="$(_heron_list_macros)"
    if [[ ${COMP_CWORD} -eq 2 ]]; then
      COMPREPLY=( $(compgen -W "${macros}" -- "${cur}") )
      return 0
    fi
  fi

  COMPREPLY=()
  return 0
}

complete -F _heron heron
