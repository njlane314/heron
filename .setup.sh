source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup sam_web_client

setup root v6_28_12 -q e20:p3915:prof
setup cmake v3_20_0 -q Linux64bit+3.10-2.17
setup nlohmann_json v3_11_2
setup libtorch v1_13_1b -q e20:prof
setup eigen v3_4_0

HERON_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export HERON_ROOT

source "${HERON_ROOT}/scripts/heron-completion.bash"

export PATH="${HERON_ROOT}/build/bin:${PATH}"
export LD_LIBRARY_PATH="${HERON_ROOT}/build/lib:${HERON_ROOT}/build/framework:${LD_LIBRARY_PATH}"
export ROOT_INCLUDE_PATH="${HERON_ROOT}/build/framework:${ROOT_INCLUDE_PATH}"

if [ -d "/exp/uboone/data/users/${USER}/cache" ]; then
  export HERON_TMPDIR="/exp/uboone/data/users/${USER}/cache/heron_tmp"
  mkdir -p "${HERON_TMPDIR}"
fi
