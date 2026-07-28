#!/bin/bash
# Submit one coarray image per SLURM task.  For example, this uses four nodes
# and 32 images on each node:
# sbatch --nodes=4 --ntasks=128 --ntasks-per-node=32 \
#   --export=ALL,NIMAGES=128,NODES=4,TASKS_PER_NODE=32,NX=1200,NY=1200 \
#   lid_performance.slurm
#SBATCH --job-name=lbm_lid_perf
#SBATCH --partition=cpu_il
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=00:30:00
#SBATCH --mem=8gb
#SBATCH --output=lid_perf_%j.out
#SBATCH --error=lid_perf_%j.err

set -euo pipefail

NIMAGES=${NIMAGES:-${SLURM_NTASKS}}
NX=${NX:-300}
NY=${NY:-300}
STEPS=${STEPS:-5000}
OMEGA=${OMEGA:-1.0}
WALL_VELOCITY=${WALL_VELOCITY:-0.0555555555555556}
NODES=${NODES:-${SLURM_JOB_NUM_NODES}}
TASKS_PER_NODE=${TASKS_PER_NODE:-1}
REPEAT=${REPEAT:-1}

if (( NIMAGES != SLURM_NTASKS )); then
  echo "NIMAGES must equal the requested --ntasks (${SLURM_NTASKS})." >&2
  exit 2
fi
if (( NODES != SLURM_JOB_NUM_NODES )); then
  echo "NODES must equal the allocated node count (${SLURM_JOB_NUM_NODES})." >&2
  exit 2
fi
if (( NIMAGES != NODES * TASKS_PER_NODE )); then
  echo "NIMAGES must equal NODES * TASKS_PER_NODE." >&2
  exit 2
fi

module purge
module load compiler/intel/2023.1.0_llvm
module load mpi/impi/2021.11

ROOT=${SLURM_SUBMIT_DIR}
RUN_DIR=${ROOT}/cluster_runs/lid_${NX}x${NY}_n${NODES}_tpn${TASKS_PER_NODE}_r${REPEAT}_${SLURM_JOB_ID}
mkdir -p "${RUN_DIR}/mod"

# Distributed coarrays are required once images span more than one node.
ifx -O3 -fpp -coarray=distributed -DLBM_PERFORMANCE_ONLY \
  -DLBM_NX="${NX}" -DLBM_NY="${NY}" -module "${RUN_DIR}/mod" \
  "${ROOT}/src/lbm_prams.f90" \
  "${ROOT}/src/domain_decomposition.f90" \
  "${ROOT}/src/lbm_solver.f90" \
  "${ROOT}/src/lbm_output.f90" \
  "${ROOT}/app/main.f90" \
  -o "${RUN_DIR}/milestone4_cluster"

export OMP_NUM_THREADS=1
export FOR_COARRAY_NUM_IMAGES=${NIMAGES}
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_PERHOST=${TASKS_PER_NODE}
# Intel's coarray runtime starts mpiexec.hydra itself.  These settings make
# Hydra use the SLURM allocation and place the requested number of images on
# each node.
export I_MPI_HYDRA_BOOTSTRAP=slurm
cd "${RUN_DIR}"

cat > allocation.csv <<EOF
job_id,nodes,tasks_per_node,images,nx,ny,steps,repeat,nodelist
${SLURM_JOB_ID},${NODES},${TASKS_PER_NODE},${NIMAGES},${NX},${NY},${STEPS},${REPEAT},${SLURM_JOB_NODELIST}
EOF

echo "Job ID:          ${SLURM_JOB_ID}"
echo "Nodes:           ${NODES} (${SLURM_JOB_NODELIST})"
echo "Tasks per node:  ${TASKS_PER_NODE}"
echo "Coarray images:  ${NIMAGES}"
echo "Grid and steps:  ${NX} x ${NY}, ${STEPS}"
echo "Repeat:          ${REPEAT}"

# Do not wrap this in mpirun: Intel's distributed-coarray runtime launches the
# images according to FOR_COARRAY_NUM_IMAGES.  An outer mpirun would instead
# start N independent coarray programs.
./milestone4_cluster lid "${OMEGA}" "${WALL_VELOCITY}" "${STEPS}" \
  2>&1 | tee simulation.log

{
  IFS= read -r performance_header
  IFS= read -r performance_row
} < performance.csv
printf '%s,nodes,tasks_per_node,repeat,job_id\n' "${performance_header}" > placement.csv
printf '%s,%d,%d,%d,%s\n' "${performance_row}" "${NODES}" \
  "${TASKS_PER_NODE}" "${REPEAT}" "${SLURM_JOB_ID}" >> placement.csv

echo "Performance row: ${RUN_DIR}/performance.csv"
echo "Placement row:   ${RUN_DIR}/placement.csv"
