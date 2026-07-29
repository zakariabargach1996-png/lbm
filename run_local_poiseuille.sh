#!/usr/bin/env bash
set -euo pipefail

# Run the Poiseuille case through the portable local driver.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export NX=${NX:-64}
export NY=${NY:-64}
export IMAGES=${IMAGES:-1}
export RUN_DIR=${RUN_DIR:-"${ROOT}/runs/local_poiseuille_$(date +%Y%m%d_%H%M%S)"}

exec "${ROOT}/run_local.sh" poiseuille \
  "${OMEGA:-1.0}" "${WALL_VELOCITY:-0.05}" "${STEPS:-}"
