#!/usr/bin/env bash
#SBATCH --job-name=lbm_perf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=lbm_perf_%j.out
#SBATCH --error=lbm_perf_%j.err

set -euo pipefail

# Run a fixed-step, output-free benchmark and record its Slurm placement.

CASE=${CASE:-lid}
NX=${NX:-300}
NY=${NY:-300}
STEPS=${STEPS:-5000}
NIMAGES=${NIMAGES:-${SLURM_NTASKS:-1}}
NODES=${NODES:-${SLURM_JOB_NUM_NODES:-1}}
TASKS_PER_NODE=${TASKS_PER_NODE:-$((NIMAGES / NODES))}
REPEAT=${REPEAT:-1}
ROOT=${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
RUN_DIR=${RUN_DIR:-"${ROOT}/cluster_runs/${CASE}_${NX}x${NY}_n${NODES}_tpn${TASKS_PER_NODE}_r${REPEAT}_${SLURM_JOB_ID:-manual}"}

if [[ ! "${CASE}" =~ ^(shear|density|couette|poiseuille|lid)$ ]]; then
  echo "CASE must name one solver case for a performance job." >&2
  exit 2
fi
if (( NIMAGES != NODES * TASKS_PER_NODE )); then
  echo "NIMAGES must equal NODES * TASKS_PER_NODE." >&2
  exit 2
fi
if [[ -n "${SLURM_NTASKS:-}" ]] && (( NIMAGES != SLURM_NTASKS )); then
  echo "NIMAGES must equal the allocated task count (${SLURM_NTASKS})." >&2
  exit 2
fi

export CASE NX NY STEPS NIMAGES RUN_DIR
export PERFORMANCE_ONLY=1
export COARRAY_MODE=${COARRAY_MODE:-distributed}
export I_MPI_PIN=${I_MPI_PIN:-1}
export I_MPI_PIN_DOMAIN=${I_MPI_PIN_DOMAIN:-core}
export I_MPI_PERHOST=${I_MPI_PERHOST:-${TASKS_PER_NODE}}

bash "${ROOT}/bwunicluster.slurm"

result_dir="${RUN_DIR}/${CASE}"
performance_file="${result_dir}/performance.csv"
if [[ ! -f "${performance_file}" ]]; then
  echo "Expected performance file was not produced: ${performance_file}" >&2
  exit 1
fi

{
  echo "job_id,nodes,tasks_per_node,images,nx,ny,steps,repeat,nodelist"
  echo "${SLURM_JOB_ID:-manual},${NODES},${TASKS_PER_NODE},${NIMAGES},${NX},${NY},${STEPS},${REPEAT},${SLURM_JOB_NODELIST:-local}"
} > "${RUN_DIR}/allocation.csv"

{
  IFS= read -r performance_header
  IFS= read -r performance_row
} < "${performance_file}"
printf '%s,nodes,tasks_per_node,repeat,job_id\n' "${performance_header}" \
  > "${RUN_DIR}/placement.csv"
printf '%s,%d,%d,%d,%s\n' "${performance_row}" "${NODES}" \
  "${TASKS_PER_NODE}" "${REPEAT}" "${SLURM_JOB_ID:-manual}" \
  >> "${RUN_DIR}/placement.csv"

echo "Performance data: ${RUN_DIR}/placement.csv"
