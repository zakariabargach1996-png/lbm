#!/bin/bash
set -euo pipefail

# Submit repeated fixed-image placement benchmarks for any solver case.
CASE=${CASE:-lid}
IMAGES=${IMAGES:-384}
NODE_COUNTS=${NODE_COUNTS:-"16 24 48"}
MAX_TASKS_PER_NODE=${MAX_TASKS_PER_NODE:-32}
SIZES=${SIZES:-"10000x10000"}
STEPS=${STEPS:-5000}
REPEATS=${REPEATS:-3}
PARTITION=${PARTITION:-cpu_il}
TIME_LIMIT=${TIME_LIMIT:-02:00:00}
MEM_PER_NODE=${MEM_PER_NODE:-8G}
EXCLUSIVE=${EXCLUSIVE:-1}
DRY_RUN=${DRY_RUN:-0}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if [[ ! "${CASE}" =~ ^(shear|density|couette|poiseuille|lid)$ ]]; then
  echo "CASE must be shear, density, couette, poiseuille, or lid." >&2
  exit 2
fi
for variable in IMAGES MAX_TASKS_PER_NODE REPEATS STEPS; do
  value=${!variable}
  if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${variable} must be a positive integer." >&2
    exit 2
  fi
done
if [[ "${DRY_RUN}" != 0 && "${DRY_RUN}" != 1 ]] ||
   [[ "${EXCLUSIVE}" != 0 && "${EXCLUSIVE}" != 1 ]]; then
  echo "DRY_RUN and EXCLUSIVE must each be either 0 or 1." >&2
  exit 2
fi

submitted=0
for size in ${SIZES}; do
  if ! [[ "${size}" =~ ^[1-9][0-9]*x[1-9][0-9]*$ ]]; then
    echo "Invalid size '${size}'; expected NXxNY." >&2
    exit 2
  fi
  NX=${size%x*}
  NY=${size#*x}
  for nodes in ${NODE_COUNTS}; do
    if ! [[ "${nodes}" =~ ^[1-9][0-9]*$ ]]; then
      echo "Invalid node count '${nodes}'; expected a positive integer." >&2
      exit 2
    fi
    if (( IMAGES % nodes != 0 )); then
      echo "Invalid placement: ${IMAGES} images cannot be divided evenly over ${nodes} nodes." >&2
      exit 2
    fi
    tasks_per_node=$((IMAGES / nodes))
    if (( tasks_per_node > MAX_TASKS_PER_NODE )); then
      echo "Invalid placement: ${tasks_per_node} images/node exceeds MAX_TASKS_PER_NODE=${MAX_TASKS_PER_NODE}." >&2
      exit 2
    fi

    # No decomposition can give a non-empty block to more than NX*NY images.
    if (( IMAGES > NX * NY )); then
      echo "Skipping ${size}: ${IMAGES} images exceed ${NX}x${NY} cells."
      continue
    fi

    for ((repeat = 1; repeat <= REPEATS; repeat++)); do
      command=(
        sbatch
        --partition="${PARTITION}"
        --time="${TIME_LIMIT}"
        --nodes="${nodes}"
        --ntasks="${IMAGES}"
        --ntasks-per-node="${tasks_per_node}"
        --cpus-per-task=1
        --mem="${MEM_PER_NODE}"
        --job-name="${CASE}_${NX}x${NY}_i${IMAGES}_n${nodes}"
        --export="ALL,CASE=${CASE},NIMAGES=${IMAGES},NODES=${nodes},TASKS_PER_NODE=${tasks_per_node},NX=${NX},NY=${NY},STEPS=${STEPS},REPEAT=${repeat}"
      )
      if (( EXCLUSIVE )); then
        command+=(--exclusive)
      fi
      command+=("${ROOT}/lid_performance2.slurm")
      printf 'Submitting %-12s images=%-4d nodes=%-3d images/node=%-2d repeat=%d/%d\n' \
        "${size}" "${IMAGES}" "${nodes}" "${tasks_per_node}" "${repeat}" "${REPEATS}"
      if (( DRY_RUN )); then
        printf '  '
        printf '%q ' "${command[@]}"
        printf '\n'
      else
        "${command[@]}"
      fi
      submitted=$((submitted + 1))
    done
  done
done

if (( DRY_RUN )); then
  echo "${submitted} scaling job(s) planned."
else
  echo "${submitted} scaling job(s) submitted."
fi
