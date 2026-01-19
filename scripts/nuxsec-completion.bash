#!/usr/bin/env bash

_nuxsec()
{
  local cur prev
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  local commands="a art artio artio-aggregate s samp sample sample-aggregate t tpl template template-make macro m help -h --help"

  _nuxsec_find_root()
  {
    local dir
    local base
    local parent
    local exe
    local exe_dir

    if [[ -n "${NUXSEC_REPO_ROOT:-}" && -d "${NUXSEC_REPO_ROOT}/plot/macro" ]]; then
      printf "%s" "${NUXSEC_REPO_ROOT}"
      return 0
    fi

    if [[ -n "${NUXSEC_ROOT:-}" && -d "${NUXSEC_ROOT}/plot/macro" ]]; then
      printf "%s" "${NUXSEC_ROOT}"
      return 0
    fi

    exe="$(command -v nuxsec 2>/dev/null || true)"
    if [[ -n "${exe}" ]]; then
      exe_dir="$(dirname "$(readlink -f "${exe}" 2>/dev/null || printf "%s" "${exe}")")"
      dir="${exe_dir}"
      while [[ -n "${dir}" && "${dir}" != "/" ]]; do
        if [[ -d "${dir}/plot/macro" ]]; then
          printf "%s" "${dir}"
          return 0
        fi
        dir="$(dirname "${dir}")"
      done
    fi

    dir="${PWD}"
    while [[ "${dir}" != "/" ]]; do
      if [[ -d "${dir}/plot/macro" ]]; then
        printf "%s" "${dir}"
        return 0
      fi
      base="$(basename "${dir}")"
      if [[ "${base}" == "plot" && -d "${dir}/macro" ]]; then
        printf "%s" "$(dirname "${dir}")"
        return 0
      fi
      if [[ "${base}" == "macro" ]]; then
        parent="$(dirname "${dir}")"
        if [[ "$(basename "${parent}")" == "plot" ]]; then
          printf "%s" "$(dirname "${parent}")"
          return 0
        fi
      fi
      dir="$(dirname "${dir}")"
    done
    return 1
  }

  _nuxsec_list_macros()
  {
    nuxsec macro list 2>/dev/null | awk '/\.C$/ {print $1}'
  }

  if [[ ${COMP_CWORD} -le 1 ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${prev}" == "help" || "${prev}" == "-h" || "${prev}" == "--help" ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${COMP_WORDS[1]}" == "macro" || "${COMP_WORDS[1]}" == "m" ]]; then
    local macros
    macros="$(_nuxsec_list_macros)"
    if [[ ${COMP_CWORD} -eq 2 ]]; then
      COMPREPLY=( $(compgen -W "${macros}" -- "${cur}") )
      return 0
    fi
  fi

  COMPREPLY=()
  return 0
}

complete -F _nuxsec nuxsec
