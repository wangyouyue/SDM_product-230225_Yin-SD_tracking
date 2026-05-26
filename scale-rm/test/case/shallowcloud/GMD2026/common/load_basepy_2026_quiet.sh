#!/usr/bin/env bash
# Source this helper before running Python jobs on SQUID.
# It suppresses module pager output and prefers a user venv with netCDF4.

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

default_gmd2026_venv="${GMD2026_NETCDF4_VENV:-${HOME}/venvs/gmd2026-netcdf4}"
if [ -z "${GMD2026_PYTHON:-}" ] && [ -x "${default_gmd2026_venv}/bin/python" ]; then
  export GMD2026_PYTHON="${default_gmd2026_venv}/bin/python"
fi

if [ -z "${GMD2026_PYTHON:-}" ]; then
  if command -v python3 >/dev/null 2>&1; then
    export GMD2026_PYTHON="$(command -v python3)"
  else
    echo "Failed to locate python3 or ${default_gmd2026_venv}/bin/python" >&2
    return 1 2>/dev/null || exit 1
  fi
fi
