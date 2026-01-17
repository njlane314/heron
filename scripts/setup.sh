#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

export LD_LIBRARY_PATH="${REPO_ROOT}/build/lib:${LD_LIBRARY_PATH:-}"

echo "LD_LIBRARY_PATH updated to include ${REPO_ROOT}/build/lib"
