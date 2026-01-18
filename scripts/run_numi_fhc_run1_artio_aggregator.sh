#!/usr/bin/env bash
while read -r spec; do
  [[ -z "${spec}" || "${spec}" == \#* ]] && continue
  build/bin/nuxsec artio-aggregate "${spec}"
done < scripts/numi_fhc_run1_artio.args
