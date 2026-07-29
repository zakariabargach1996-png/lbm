#!/usr/bin/env bash
set -euo pipefail

# Run the 300x300 moving-lid case through the portable local driver.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export NX=${NX:-300}
export NY=${NY:-300}
export IMAGES=${IMAGES:-1}
export RUN_DIR=${RUN_DIR:-"${ROOT}/runs/local_lid_$(date +%Y%m%d_%H%M%S)"}

exec "${ROOT}/run_local.sh" lid \
  "${OMEGA:-1.0}" "${WALL_VELOCITY:-0.0555555555555556}" "${STEPS:-}"
