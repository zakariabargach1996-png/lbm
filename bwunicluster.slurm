#!/bin/bash
#SBATCH --job-name=lbm_coarray
#SBATCH --partition=cpu_il
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00
#SBATCH --mem=8gb
#SBATCH --output=lbm_%j.out
#SBATCH --error=lbm_%j.err

set -euo pipefail

# Values may be overridden when submitting, for example:
# sbatch --export=ALL,CASE=lid,OMEGA=1.5,WALL_VELOCITY=0.05 bwunicluster.slurm
CASE=${CASE:-poiseuille}
if [[ "${CASE}" == "lid" ]]; then
  # Keep the assigned cavity size fixed while omega and wall speed vary.
  NX=300
  NY=300
else
  NX=${NX:-512}
  NY=${NY:-512}
fi
OMEGA=${OMEGA:-1.4285714285714286}
WALL_VELOCITY=${WALL_VELOCITY:-0.05}
STEPS=${STEPS:-0}
OMEGAS=${OMEGAS:-"0.6 0.8 1.0 1.2 1.4 1.6 1.7"}

module purge
module load compiler/intel
module load mpi/impi

cd "${SLURM_SUBMIT_DIR}"

echo "Job ID:       ${SLURM_JOB_ID}"
echo "Node:         ${SLURM_JOB_NODELIST}"
echo "Case:         ${CASE}"
echo "Grid:         ${NX} x ${NY}"
echo "Omega:        ${OMEGA}"
echo "Wall speed:   ${WALL_VELOCITY}"
echo "Coarray images: ${SLURM_CPUS_PER_TASK}"
echo "Compiler:     $(command -v ifx)"

# Build directly so the compute node does not need FPM installed. The source
# order matters because domain_decomposition uses the lbm_params module.
ifx -O3 -fpp -coarray=shared \
  -DLBM_NX="${NX}" -DLBM_NY="${NY}" \
  src/lbm_prams.f90 \
  src/domain_decomposition.f90 \
  src/lbm_solver.f90 \
  src/lbm_output.f90 \
  app/main.f90 \
  -o milestone4_cluster

# The shared coarray runtime creates this many images on the allocated node.
export FOR_COARRAY_NUM_IMAGES=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=1

if [[ "${CASE}" == "viscosity-sweep" ]]; then
  for sweep_omega in ${OMEGAS}; do
    echo "Running shear wave with omega=${sweep_omega}"
    ./milestone4_cluster shear "${sweep_omega}" "${WALL_VELOCITY}"
  done
else
  if (( STEPS > 0 )); then
    ./milestone4_cluster "${CASE}" "${OMEGA}" "${WALL_VELOCITY}" "${STEPS}"
  else
    ./milestone4_cluster "${CASE}" "${OMEGA}" "${WALL_VELOCITY}"
  fi
fi

echo "Finished at $(date)"
