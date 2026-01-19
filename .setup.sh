source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup sam_web_client

setup root v6_28_12 -q e20:p3915:prof
setup cmake v3_20_0 -q Linux64bit+3.10-2.17
setup nlohmann_json v3_11_2
setup tbb v2021_9_0 -q e20
setup libtorch v1_13_1b -q e20:prof
setup eigen v3_4_0

NUXSEC_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NUXSEC_ROOT

source "${NUXSEC_ROOT}/scripts/nuxsec-completion.bash"

export PATH="${NUXSEC_ROOT}/build/bin:${PATH}"
export LD_LIBRARY_PATH="${NUXSEC_ROOT}/build/lib:${NUXSEC_ROOT}/build/framework:${LD_LIBRARY_PATH}"
export ROOT_INCLUDE_PATH="${NUXSEC_ROOT}/build/framework:${ROOT_INCLUDE_PATH}"
