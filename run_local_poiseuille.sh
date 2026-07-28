#!/bin/bash
set -euo pipefail

# Small defaults keep the steady-state run practical on a laptop/desktop.
# Override NX, NY, IMAGES, or OMEGA in the environment when needed.
NX=${NX:-64}
NY=${NY:-64}
IMAGES=${IMAGES:-4}
OMEGA=${OMEGA:-1.0}

if (( NX < 4 || NY < 4 )); then
  echo "NX and NY must both be at least 4." >&2
  exit 2
fi

if (( IMAGES < 1 || IMAGES > 4 )); then
  echo "IMAGES must be between 1 and 4 for the conservative local run." >&2
  exit 2
fi

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUN_DIR=${RUN_DIR:-"${ROOT}/local_poiseuille_output_$(date +%Y%m%d_%H%M%S)"}
mkdir -p "${RUN_DIR}/mod" "${RUN_DIR}/figures"

# oneAPI's setup script probes optional components and is not nounset-clean.
set +u
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
set -u

NU=$(awk -v omega="${OMEGA}" 'BEGIN {print (1.0/omega-0.5)/3.0}')

echo "Local Poiseuille run"
echo "  domain:       ${NX} x ${NY}"
echo "  images:       ${IMAGES} (automatic 2D decomposition)"
echo "  omega:        ${OMEGA}"
echo "  viscosity:    ${NU}"
echo "  target umax:  0.05"
echo "  output:       ${RUN_DIR}"

ifx -O3 -fpp -coarray=shared \
  -DLBM_NX="${NX}" -DLBM_NY="${NY}" \
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
nice -n 5 ./milestone4_local poiseuille "${OMEGA}" \
  2>&1 | tee simulation.log

MPLCONFIGDIR=/tmp/matplotlib-lbm conda run -n meep \
  python "${ROOT}/plots.py" --data-dir "${RUN_DIR}" \
  --output-dir "${RUN_DIR}/figures"

conda run -n meep python "${ROOT}/test/validate_outputs.py" \
  --data-dir "${RUN_DIR}"

echo "Finished. Open ${RUN_DIR}/figures/poiseuille_evolution.png"
