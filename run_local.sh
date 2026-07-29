#!/usr/bin/env bash
set -euo pipefail

# Build and run any solver case with Intel Fortran, OpenCoarrays, or GNU Fortran.

usage() {
  echo "Usage: $0 CASE [OMEGA] [WALL_VELOCITY] [MAX_STEPS]"
  echo "CASE: shear, density, couette, poiseuille, lid, or all"
  echo "Environment: NX NY IMAGES RUN_DIR FC BUILD_TYPE PERFORMANCE_ONLY"
}

case_name=${1:-}
omega=${2:-1.4285714285714286}
wall_velocity=${3:-0.05}
max_steps=${4:-}

if [[ -z "${case_name}" || "${case_name}" == "-h" || "${case_name}" == "--help" ]]; then
  usage
  [[ -n "${case_name}" ]] && exit 0
  exit 2
fi
if [[ ! "${case_name}" =~ ^(shear|density|couette|poiseuille|lid|all)$ ]]; then
  echo "Unknown case: ${case_name}" >&2
  usage >&2
  exit 2
fi

NX=${NX:-128}
NY=${NY:-128}
IMAGES=${IMAGES:-1}
BUILD_TYPE=${BUILD_TYPE:-release}
PERFORMANCE_ONLY=${PERFORMANCE_ONLY:-0}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUN_DIR=${RUN_DIR:-"${ROOT}/runs/local_${case_name}_$(date +%Y%m%d_%H%M%S)"}

for item in NX NY IMAGES; do
  value=${!item}
  if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${item} must be a positive integer." >&2
    exit 2
  fi
done
if (( NX < 4 || NY < 4 )); then
  echo "NX and NY must both be at least 4." >&2
  exit 2
fi
if [[ "${PERFORMANCE_ONLY}" != 0 && "${PERFORMANCE_ONLY}" != 1 ]]; then
  echo "PERFORMANCE_ONLY must be 0 or 1." >&2
  exit 2
fi
if [[ -n "${max_steps}" ]] && ! [[ "${max_steps}" =~ ^[1-9][0-9]*$ ]]; then
  echo "MAX_STEPS must be a positive integer." >&2
  exit 2
fi
if (( PERFORMANCE_ONLY )) && [[ -z "${max_steps}" ]]; then
  echo "PERFORMANCE_ONLY=1 requires MAX_STEPS." >&2
  exit 2
fi

compiler=${FC:-}
if [[ -z "${compiler}" ]]; then
  if command -v ifx >/dev/null 2>&1; then
    compiler=ifx
  elif command -v caf >/dev/null 2>&1 && command -v cafrun >/dev/null 2>&1; then
    compiler=caf
  elif command -v gfortran >/dev/null 2>&1; then
    compiler=gfortran
  else
    echo "No supported Fortran compiler found (ifx, caf, or gfortran)." >&2
    exit 127
  fi
fi
if ! command -v "${compiler}" >/dev/null 2>&1; then
  echo "Compiler not found: ${compiler}" >&2
  exit 127
fi

compiler_name=$(basename "${compiler}")
case "${compiler_name}" in
  ifx)
    preprocess_flag=-fpp
    coarray_flag=-coarray=shared
    launcher=direct
    ;;
  caf)
    preprocess_flag=-cpp
    coarray_flag=
    launcher=cafrun
    ;;
  gfortran)
    preprocess_flag=-cpp
    coarray_flag=-fcoarray=single
    launcher=direct
    if (( IMAGES != 1 )); then
      echo "Plain gfortran supports one image here; install OpenCoarrays or use IMAGES=1." >&2
      exit 2
    fi
    ;;
  *)
    echo "Unsupported FC=${compiler}; use ifx, caf, or gfortran." >&2
    exit 2
    ;;
esac

case "${BUILD_TYPE}" in
  release) optimization=-O3 ;;
  debug) optimization=-O0 ;;
  *) echo "BUILD_TYPE must be release or debug." >&2; exit 2 ;;
esac

mkdir -p "${RUN_DIR}/build/mod"
binary="${RUN_DIR}/build/lbm_solver"
compile_command=(
  "${compiler}" "${optimization}" "${preprocess_flag}"
  "-DLBM_NX=${NX}" "-DLBM_NY=${NY}"
)
[[ -n "${coarray_flag}" ]] && compile_command+=("${coarray_flag}")
if [[ "${compiler_name}" == "ifx" ]]; then
  compile_command+=(-module "${RUN_DIR}/build/mod")
else
  compile_command+=(-J "${RUN_DIR}/build/mod")
fi
if (( PERFORMANCE_ONLY )); then
  compile_command+=(-DLBM_PERFORMANCE_ONLY)
fi
compile_command+=(
  "${ROOT}/src/lbm_prams.f90"
  "${ROOT}/src/domain_decomposition.f90"
  "${ROOT}/src/lbm_solver.f90"
  "${ROOT}/src/lbm_output.f90"
  "${ROOT}/app/main.f90"
  -o "${binary}"
)

echo "Building ${NX}x${NY} solver with ${compiler_name}."
(
  cd "${RUN_DIR}/build"
  "${compile_command[@]}"
)

run_one() {
  local selected_case=$1
  local output_dir="${RUN_DIR}/${selected_case}"
  local command=("${binary}" "${selected_case}" "${omega}" "${wall_velocity}")
  [[ -n "${max_steps}" ]] && command+=("${max_steps}")
  mkdir -p "${output_dir}"

  echo "Running ${selected_case}: images=${IMAGES}, output=${output_dir}"
  (
    cd "${output_dir}"
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
    if [[ "${launcher}" == "cafrun" ]]; then
      cafrun -n "${IMAGES}" "${command[@]}" 2>&1 | tee simulation.log
    else
      export FOR_COARRAY_NUM_IMAGES=${IMAGES}
      "${command[@]}" 2>&1 | tee simulation.log
    fi
  )
}

if [[ "${case_name}" == "all" ]]; then
  for selected_case in shear density couette poiseuille lid; do
    run_one "${selected_case}"
  done
else
  run_one "${case_name}"
fi

echo "Finished. Results are in ${RUN_DIR}."
