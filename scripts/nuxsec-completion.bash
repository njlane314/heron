#!/usr/bin/env bash

_nuxsec()
{
  local cur prev
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  local commands="a art artio artio-aggregate s samp sample sample-aggregate t tpl template template-make macro m plot help -h --help"

  _nuxsec_find_root()
  {
    local dir="${PWD}"
    while [[ "${dir}" != "/" ]]; do
      if [[ -d "${dir}/plot/macro" ]]; then
        printf "%s" "${dir}"
        return 0
      fi
      dir="$(dirname "${dir}")"
    done
    return 1
  }

  _nuxsec_list_macros()
  {
    local root
    root="$(_nuxsec_find_root)" || return 0
    local macro_dir="${root}/plot/macro"
    local macro
    shopt -s nullglob
    for macro in "${macro_dir}"/*.C; do
      printf "%s " "$(basename "${macro}")"
    done
    shopt -u nullglob
  }

  if [[ ${COMP_CWORD} -le 1 ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${prev}" == "help" || "${prev}" == "-h" || "${prev}" == "--help" ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${COMP_WORDS[1]}" == "macro" || "${COMP_WORDS[1]}" == "m" || "${COMP_WORDS[1]}" == "plot" ]]; then
    local macro_commands="run list"
    local macros
    macros="$(_nuxsec_list_macros)"
    if [[ ${COMP_CWORD} -eq 2 ]]; then
      COMPREPLY=( $(compgen -W "${macro_commands} ${macros}" -- "${cur}") )
      return 0
    fi
    if [[ "${prev}" == "run" ]]; then
      COMPREPLY=( $(compgen -W "${macros}" -- "${cur}") )
      return 0
    fi
    COMPREPLY=( $(compgen -W "${macros}" -- "${cur}") )
    return 0
  fi

  COMPREPLY=()
  return 0
}

complete -F _nuxsec nuxsec
