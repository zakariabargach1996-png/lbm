#!/bin/bash
set -euo pipefail

IMAGES=${IMAGES:-4}
WALL_VELOCITY=${WALL_VELOCITY:-0.1}
STEP_LIMIT=${STEP_LIMIT:-1600000}
REYNOLDS_NUMBERS=${REYNOLDS_NUMBERS:-"100 500 1000 2500 5000"}

if (( IMAGES < 1 || IMAGES > 4 )); then
  echo "IMAGES must be between 1 and 4." >&2
  exit 2
fi
if (( STEP_LIMIT < 1 )); then
  echo "STEP_LIMIT must be a positive integer." >&2
  exit 2
fi

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RESULTS_DIR=${RESULTS_DIR:-"${ROOT}/lid_reynolds_results_300x300"}
BUILD_DIR="${RESULTS_DIR}/build"
mkdir -p "${BUILD_DIR}/mod"

set +u
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
set -u

ifx -O3 -fpp -coarray=shared \
  -DLBM_NX=300 -DLBM_NY=300 \
  -module "${BUILD_DIR}/mod" \
  "${ROOT}/src/lbm_prams.f90" \
  "${ROOT}/src/domain_decomposition.f90" \
  "${ROOT}/src/lbm_solver.f90" \
  "${ROOT}/src/lbm_output.f90" \
  "${ROOT}/app/main.f90" \
  -o "${BUILD_DIR}/milestone4_lid"

export FOR_COARRAY_NUM_IMAGES=${IMAGES}
export OMP_NUM_THREADS=1
export LBM_STEP_LIMIT=${STEP_LIMIT}

for reynolds in ${REYNOLDS_NUMBERS}; do
  run_dir="${RESULTS_DIR}/Re_${reynolds}"
  mkdir -p "${run_dir}/figures"

  viscosity=$(awk -v u="${WALL_VELOCITY}" -v n=300 -v re="${reynolds}" \
    'BEGIN {printf "%.16g", u*n/re}')
  omega=$(awk -v nu="${viscosity}" \
    'BEGIN {printf "%.16g", 1.0/(0.5+3.0*nu)}')

  {
    echo "Requested Reynolds number: ${reynolds}"
    echo "Grid: 300 x 300"
    echo "Wall velocity: ${WALL_VELOCITY}"
    echo "Viscosity: ${viscosity}"
    echo "Omega: ${omega}"
    echo "Step limit: ${STEP_LIMIT}"
  } | tee "${run_dir}/parameters.txt"

  (
    cd "${run_dir}"
    nice -n 5 "${BUILD_DIR}/milestone4_lid" \
      lid "${omega}" "${WALL_VELOCITY}" 2>&1 | tee simulation.log
  )

  MPLCONFIGDIR=/tmp/matplotlib-lbm conda run -n meep \
    python "${ROOT}/plots.py" --data-dir "${run_dir}" \
    --output-dir "${run_dir}/figures"

  conda run -n meep python "${ROOT}/test/validate_outputs.py" \
    --data-dir "${run_dir}"
done

echo "Finished lid sweep in ${RESULTS_DIR}"
