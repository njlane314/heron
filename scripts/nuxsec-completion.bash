#!/usr/bin/env bash

_nuxsec()
{
  local cur prev
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  local commands="a art artio artio-aggregate s samp sample sample-aggregate t tpl template template-make plot help -h --help"

  if [[ ${COMP_CWORD} -le 1 ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${prev}" == "help" || "${prev}" == "-h" || "${prev}" == "--help" ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  COMPREPLY=()
  return 0
}

complete -F _nuxsec nuxsec
