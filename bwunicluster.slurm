#!/usr/bin/env bash
#SBATCH --job-name=lbm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=lbm_%j.out
#SBATCH --error=lbm_%j.err

set -euo pipefail

# Build and run one or all cases inside a Slurm allocation.

CASE=${CASE:-poiseuille}
NX=${NX:-128}
NY=${NY:-128}
OMEGA=${OMEGA:-1.4285714285714286}
WALL_VELOCITY=${WALL_VELOCITY:-0.05}
STEPS=${STEPS:-}
PERFORMANCE_ONLY=${PERFORMANCE_ONLY:-0}
NIMAGES=${NIMAGES:-${SLURM_NTASKS:-1}}
COARRAY_MODE=${COARRAY_MODE:-distributed}
FC=${FC:-ifx}
ROOT=${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
RUN_DIR=${RUN_DIR:-"${ROOT}/cluster_runs/${CASE}_${NX}x${NY}_${SLURM_JOB_ID:-manual}"}

if [[ ! "${CASE}" =~ ^(shear|density|couette|poiseuille|lid|all)$ ]]; then
  echo "CASE must be shear, density, couette, poiseuille, lid, or all." >&2
  exit 2
fi
for item in NX NY NIMAGES; do
  value=${!item}
  if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${item} must be a positive integer." >&2
    exit 2
  fi
done
if [[ -n "${STEPS}" ]] && ! [[ "${STEPS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "STEPS must be a positive integer when set." >&2
  exit 2
fi
if [[ "${PERFORMANCE_ONLY}" != 0 && "${PERFORMANCE_ONLY}" != 1 ]]; then
  echo "PERFORMANCE_ONLY must be 0 or 1." >&2
  exit 2
fi
if (( PERFORMANCE_ONLY )) && [[ -z "${STEPS}" ]]; then
  echo "PERFORMANCE_ONLY=1 requires STEPS." >&2
  exit 2
fi

# Sites can provide their own environment setup, for example:
# MODULES="compiler/intel mpi/impi" or ENV_SETUP=/path/to/setup.sh.
if [[ -n "${MODULES:-}" ]]; then
  module purge
  for package in ${MODULES}; do
    module load "${package}"
  done
fi
if [[ -n "${ENV_SETUP:-}" ]]; then
  set +u
  source "${ENV_SETUP}"
  set -u
fi
if ! command -v "${FC}" >/dev/null 2>&1; then
  echo "Compiler '${FC}' is unavailable; set MODULES, ENV_SETUP, or FC." >&2
  exit 127
fi

mkdir -p "${RUN_DIR}/build/mod"
compiler_name=$(basename "${FC}")
case "${compiler_name}" in
  ifx|ifort)
    compile_command=(
      "${FC}" -O3 -fpp "-coarray=${COARRAY_MODE}"
      "-DLBM_NX=${NX}" "-DLBM_NY=${NY}"
      -module "${RUN_DIR}/build/mod"
    )
    launcher=direct
    ;;
  caf)
    compile_command=(
      "${FC}" -O3 -cpp "-DLBM_NX=${NX}" "-DLBM_NY=${NY}"
      -J "${RUN_DIR}/build/mod"
    )
    launcher=cafrun
    ;;
  gfortran)
    if (( NIMAGES != 1 )); then
      echo "Plain gfortran supports only NIMAGES=1; use caf or Intel Fortran." >&2
      exit 2
    fi
    compile_command=(
      "${FC}" -O3 -cpp -fcoarray=single
      "-DLBM_NX=${NX}" "-DLBM_NY=${NY}"
      -J "${RUN_DIR}/build/mod"
    )
    launcher=direct
    ;;
  *)
    echo "FC must name ifx, ifort, caf, or gfortran." >&2
    exit 2
    ;;
esac
if (( PERFORMANCE_ONLY )); then
  compile_command+=(-DLBM_PERFORMANCE_ONLY)
fi
compile_command+=(
  "${ROOT}/src/lbm_prams.f90"
  "${ROOT}/src/domain_decomposition.f90"
  "${ROOT}/src/lbm_solver.f90"
  "${ROOT}/src/lbm_output.f90"
  "${ROOT}/app/main.f90"
  -o "${RUN_DIR}/build/lbm_solver"
)
(
  cd "${RUN_DIR}/build"
  "${compile_command[@]}"
)

export FOR_COARRAY_NUM_IMAGES=${NIMAGES}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export I_MPI_HYDRA_BOOTSTRAP=${I_MPI_HYDRA_BOOTSTRAP:-slurm}

run_one() {
  local selected_case=$1
  local output_dir="${RUN_DIR}/${selected_case}"
  local command=(
    "${RUN_DIR}/build/lbm_solver"
    "${selected_case}" "${OMEGA}" "${WALL_VELOCITY}"
  )
  [[ -n "${STEPS}" ]] && command+=("${STEPS}")
  mkdir -p "${output_dir}"
  echo "Running ${selected_case}: ${NX}x${NY}, ${NIMAGES} image(s)."
  (
    cd "${output_dir}"
    if [[ "${launcher}" == "cafrun" ]]; then
      cafrun -n "${NIMAGES}" "${command[@]}" 2>&1 | tee simulation.log
    else
      "${command[@]}" 2>&1 | tee simulation.log
    fi
  )
}

if [[ "${CASE}" == "all" ]]; then
  for selected_case in shear density couette poiseuille lid; do
    run_one "${selected_case}"
  done
else
  run_one "${CASE}"
fi

echo "Finished. Results are in ${RUN_DIR}."
