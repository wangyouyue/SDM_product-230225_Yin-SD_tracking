#!/usr/bin/env bash
# Source this helper before running Python jobs on SQUID.
# It suppresses module pager output and prefers the user's sdm_env conda
# environment for NetCDF-heavy GMD2026 analysis.

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

export PAGER=cat
export LESS="${LESS:-FRX}"
export LMOD_PAGER=cat
export MODULES_PAGER=cat

module_log="${MODULE_QUIET_LOG:-$(mktemp "${TMPDIR:-/tmp}/gmd2026_module.XXXXXX")}"
cleanup_module_log=0
if [ -z "${MODULE_QUIET_LOG:-}" ]; then
  cleanup_module_log=1
fi
: > "$module_log"

if command -v module >/dev/null 2>&1; then
  if ! module --silent purge >>"$module_log" 2>&1; then
    if ! module -s purge >>"$module_log" 2>&1; then
      module purge >>"$module_log" 2>&1 || true
    fi
  fi

  if ! module --silent load BasePy/2026 >>"$module_log" 2>&1; then
    if ! module -s load BasePy/2026 >>"$module_log" 2>&1; then
      if ! module load BasePy/2026 >>"$module_log" 2>&1; then
        echo "Failed to load BasePy/2026; see ${module_log}" >&2
        sed -n '1,160p' "$module_log" >&2 || true
        return 1 2>/dev/null || exit 1
      fi
    fi
  fi
fi

if [ "$cleanup_module_log" = "1" ]; then
  rm -f "$module_log"
fi

preferred_conda_env="${GMD2026_CONDA_ENV:-sdm_env}"
if [ -z "${GMD2026_PYTHON:-}" ] && [ -n "${preferred_conda_env}" ]; then
  conda_setup_candidates=()
  if command -v conda >/dev/null 2>&1; then
    conda_base="$(conda info --base 2>/dev/null || true)"
    if [ -n "${conda_base}" ]; then
      conda_setup_candidates+=("${conda_base}/etc/profile.d/conda.sh")
    fi
  fi
  conda_setup_candidates+=(
    "${HOME}/miniconda3/etc/profile.d/conda.sh"
    "${HOME}/anaconda3/etc/profile.d/conda.sh"
  )
  for conda_setup in "${conda_setup_candidates[@]}"; do
    if [ -f "${conda_setup}" ]; then
      # shellcheck disable=SC1090
      source "${conda_setup}"
      break
    fi
  done
  if command -v conda >/dev/null 2>&1; then
    if conda activate "${preferred_conda_env}" >/dev/null 2>&1; then
      export GMD2026_PYTHON="$(command -v python)"
      export GMD2026_PYTHON_ENV="conda:${preferred_conda_env}"
    fi
  fi
fi

default_gmd2026_venv="${GMD2026_NETCDF4_VENV:-${HOME}/venvs/gmd2026-netcdf4}"
if [ -z "${GMD2026_PYTHON:-}" ] && [ -x "${default_gmd2026_venv}/bin/python" ]; then
  export GMD2026_PYTHON="${default_gmd2026_venv}/bin/python"
  export GMD2026_PYTHON_ENV="venv:${default_gmd2026_venv}"
fi

if [ -z "${GMD2026_PYTHON:-}" ]; then
  if command -v python3 >/dev/null 2>&1; then
    export GMD2026_PYTHON="$(command -v python3)"
    export GMD2026_PYTHON_ENV="${GMD2026_PYTHON_ENV:-python3}"
  else
    echo "Failed to locate python3 or ${default_gmd2026_venv}/bin/python" >&2
    return 1 2>/dev/null || exit 1
  fi
fi
