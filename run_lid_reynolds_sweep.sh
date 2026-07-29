#!/usr/bin/env bash
set -euo pipefail

# Run the 300x300 lid cavity at a list of Reynolds numbers.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NX=${NX:-300}
NY=${NY:-300}
IMAGES=${IMAGES:-1}
WALL_VELOCITY=${WALL_VELOCITY:-0.1}
STEP_LIMIT=${STEP_LIMIT:-1600000}
REYNOLDS_NUMBERS=${REYNOLDS_NUMBERS:-"100 500 1000 2500 5000"}
RESULTS_DIR=${RESULTS_DIR:-"${ROOT}/runs/lid_reynolds_$(date +%Y%m%d_%H%M%S)"}

mkdir -p "${RESULTS_DIR}"
for reynolds in ${REYNOLDS_NUMBERS}; do
  viscosity=$(awk -v u="${WALL_VELOCITY}" -v n="${NX}" -v re="${reynolds}" \
    'BEGIN {printf "%.16g", u*n/re}')
  omega=$(awk -v nu="${viscosity}" \
    'BEGIN {printf "%.16g", 1.0/(0.5+3.0*nu)}')
  NX=${NX} NY=${NY} IMAGES=${IMAGES} LBM_STEP_LIMIT=${STEP_LIMIT} \
    RUN_DIR="${RESULTS_DIR}/Re_${reynolds}" \
    "${ROOT}/run_local.sh" lid "${omega}" "${WALL_VELOCITY}"
done
echo "Finished. Results are in ${RESULTS_DIR}."
