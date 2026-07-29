#!/usr/bin/env bash
set -euo pipefail

# Run the shear-wave case for several relaxation parameters.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NX=${NX:-300}
NY=${NY:-300}
IMAGES=${IMAGES:-1}
OMEGAS=${OMEGAS:-"0.6 0.8 1.0 1.2 1.4 1.6"}
RESULTS_DIR=${RESULTS_DIR:-"${ROOT}/runs/shear_sweep_$(date +%Y%m%d_%H%M%S)"}
PYTHON=${PYTHON:-python3}
PLOT=${PLOT:-0}

mkdir -p "${RESULTS_DIR}"
for omega in ${OMEGAS}; do
  if ! awk -v value="${omega}" 'BEGIN { exit !(value > 0.0 && value < 2.0) }'; then
    echo "Invalid omega '${omega}'; expected 0 < omega < 2." >&2
    exit 2
  fi
  label=$(printf '%0.3f' "${omega}")
  NX=${NX} NY=${NY} IMAGES=${IMAGES} \
    RUN_DIR="${RESULTS_DIR}/omega_${label}" \
    "${ROOT}/run_local.sh" shear "${omega}"
done

if (( PLOT )); then
  "${PYTHON}" "${ROOT}/shear_sweep_plots.py" \
    --data-dir "${RESULTS_DIR}" --output-dir "${RESULTS_DIR}/figures"
fi
echo "Finished. Results are in ${RESULTS_DIR}."
