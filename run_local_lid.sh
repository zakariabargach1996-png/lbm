#!/bin/bash
set -euo pipefail

# Conservative defaults for a laptop/desktop run. At omega=1 and U=1/18,
# nu=1/6 and the 300-wide cavity has Re=100.
IMAGES=${IMAGES:-4}
OMEGA=${OMEGA:-1.0}
WALL_VELOCITY=${WALL_VELOCITY:-0.0555555555555556}

if (( IMAGES < 1 || IMAGES > 4 )); then
  echo "IMAGES must be between 1 and 4 for the conservative local run." >&2
  exit 2
fi

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUN_DIR=${RUN_DIR:-"${ROOT}/local_lid_output_$(date +%Y%m%d_%H%M%S)"}
mkdir -p "${RUN_DIR}/mod" "${RUN_DIR}/figures"

# oneAPI's setup script probes optional components and is not nounset-clean.
set +u
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
set -u

NU=$(awk -v omega="${OMEGA}" 'BEGIN {print (1.0/omega-0.5)/3.0}')
RE=$(awk -v u="${WALL_VELOCITY}" -v nu="${NU}" 'BEGIN {print u*300.0/nu}')

echo "Local moving-lid run"
echo "  domain:       300 x 300"
echo "  images:       ${IMAGES}"
echo "  omega:        ${OMEGA}"
echo "  viscosity:    ${NU}"
echo "  wall speed:   ${WALL_VELOCITY}"
echo "  Reynolds:     ${RE}"
echo "  output:       ${RUN_DIR}"
echo "  memory use is expected to stay well below 1 GiB"

ifx -O3 -fpp -coarray=shared \
  -DLBM_NX=300 -DLBM_NY=300 \
  -module "${RUN_DIR}/mod" \
  "${ROOT}/src/lbm_prams.f90" \
  "${ROOT}/src/domain_decomposition.f90" \
  "${ROOT}/src/lbm_solver.f90" \
  "${ROOT}/src/lbm_output.f90" \
  "${ROOT}/app/main.f90" \
  -o "${RUN_DIR}/milestone4_local"

export FOR_COARRAY_NUM_IMAGES=${IMAGES}
export OMP_NUM_THREADS=1

cd "${RUN_DIR}"
# Lower priority keeps the desktop responsive while the CPU-bound solver runs.
nice -n 5 ./milestone4_local lid "${OMEGA}" "${WALL_VELOCITY}" \
  2>&1 | tee simulation.log

MPLCONFIGDIR=/tmp/matplotlib-lbm conda run -n meep \
  python "${ROOT}/plots.py" --data-dir "${RUN_DIR}" \
  --output-dir "${RUN_DIR}/figures"

conda run -n meep python "${ROOT}/test/validate_outputs.py" \
  --data-dir "${RUN_DIR}"

echo "Finished. Open ${RUN_DIR}/figures/moving_lid_evolution.png"
