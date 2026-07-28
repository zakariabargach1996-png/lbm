#!/bin/bash
set -euo pipefail

# The report sweep uses the requested 300x300 lattice.
NX=${NX:-300}
NY=${NY:-300}
IMAGES=${IMAGES:-4}
OMEGAS=${OMEGAS:-"0.6 0.8 1.0 1.2 1.4 1.6"}

if (( NX < 16 || NY < 16 )); then
  echo "NX and NY must both be at least 16." >&2
  exit 2
fi
if (( IMAGES < 1 || IMAGES > 4 )); then
  echo "IMAGES must be between 1 and 4 for the local run." >&2
  exit 2
fi

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RESULTS_DIR=${RESULTS_DIR:-"${ROOT}/shear_viscosity_sweep_$(date +%Y%m%d_%H%M%S)"}
BUILD_DIR="${RESULTS_DIR}/build"
mkdir -p "${BUILD_DIR}/mod" "${RESULTS_DIR}/figures"

# oneAPI's setup script probes optional components and is not nounset-clean.
set +u
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
set -u

ifx -O3 -fpp -coarray=shared \
  -DLBM_NX="${NX}" -DLBM_NY="${NY}" \
  -module "${BUILD_DIR}/mod" \
  "${ROOT}/src/lbm_prams.f90" \
  "${ROOT}/src/domain_decomposition.f90" \
  "${ROOT}/src/lbm_solver.f90" \
  "${ROOT}/src/lbm_output.f90" \
  "${ROOT}/app/main.f90" \
  -o "${BUILD_DIR}/milestone4_shear"

export FOR_COARRAY_NUM_IMAGES=${IMAGES}
export OMP_NUM_THREADS=1

{
  echo "grid=${NX}x${NY}"
  echo "images=${IMAGES}"
  echo "omegas=${OMEGAS}"
} > "${RESULTS_DIR}/sweep_parameters.txt"

for omega in ${OMEGAS}; do
  if ! awk -v value="${omega}" 'BEGIN { exit !(value > 0.0 && value < 2.0) }'; then
    echo "Invalid omega '${omega}'; every omega must satisfy 0 < omega < 2." >&2
    exit 2
  fi

  omega_label=$(printf '%0.3f' "${omega}")
  run_dir="${RESULTS_DIR}/omega_${omega_label}"
  mkdir -p "${run_dir}"
  viscosity=$(awk -v value="${omega}" \
    'BEGIN { printf "%.16g", (1.0/value - 0.5)/3.0 }')

  echo "Running omega=${omega_label}, applied viscosity=${viscosity}"
  {
    echo "omega=${omega}"
    echo "applied_viscosity=${viscosity}"
    echo "grid=${NX}x${NY}"
    echo "images=${IMAGES}"
  } > "${run_dir}/parameters.txt"

  (
    cd "${run_dir}"
    nice -n 5 "${BUILD_DIR}/milestone4_shear" shear "${omega}" \
      2>&1 | tee simulation.log
  )

  conda run -n meep python "${ROOT}/test/validate_outputs.py" \
    --data-dir "${run_dir}"
done

MPLCONFIGDIR=/tmp/matplotlib-lbm conda run -n meep \
  python "${ROOT}/shear_sweep_plots.py" \
  --data-dir "${RESULTS_DIR}" \
  --output-dir "${RESULTS_DIR}/figures"

echo "Finished shear-wave sweep in ${RESULTS_DIR}"
echo "Figures:"
echo "  ${RESULTS_DIR}/figures/applied_vs_measured_viscosity.png"
echo "  ${RESULTS_DIR}/figures/applied_vs_measured_viscosity_by_timestep.png"
echo "  ${RESULTS_DIR}/figures/measured_vs_theory_decay.png"
