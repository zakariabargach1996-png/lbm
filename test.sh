#!/usr/bin/env bash
set -euo pipefail

# Build unit tests and optionally run a short simulation of every case.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
FC=${FC:-}
if [[ -z "${FC}" ]]; then
  if command -v ifx >/dev/null 2>&1; then
    FC=ifx
  elif command -v gfortran >/dev/null 2>&1; then
    FC=gfortran
  else
    echo "No supported Fortran compiler found." >&2
    exit 127
  fi
fi
case "$(basename "${FC}")" in
  ifx) compiler_flags=(-fpp -coarray=shared); module_flag=-module ;;
  gfortran) compiler_flags=(-cpp -fcoarray=single); module_flag=-J ;;
  *) echo "FC must name ifx or gfortran for unit tests." >&2; exit 2 ;;
esac

build_dir=$(mktemp -d "${TMPDIR:-/tmp}/lbm-test.XXXXXXXX")
trap 'rm -rf "${build_dir}"' EXIT
(
  cd "${build_dir}"
  "${FC}" "${compiler_flags[@]}" "${module_flag}" "${build_dir}" \
    "${ROOT}/src/lbm_prams.f90" "${ROOT}/test/check.f90" \
    -o "${build_dir}/check"
)
FOR_COARRAY_NUM_IMAGES=1 "${build_dir}/check"

if [[ "${SMOKE:-0}" == 1 ]]; then
  NX=${NX:-16} NY=${NY:-16} IMAGES=${IMAGES:-1} \
    RUN_DIR=${RUN_DIR:-"${ROOT}/runs/smoke_test"} \
    "${ROOT}/run_local.sh" all 1.0 0.02 "${STEPS:-20}"
fi
