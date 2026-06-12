#!/usr/bin/env bash
set -euo pipefail

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

export PAGER=cat
export LESS="${LESS:-FRX}"
export LMOD_PAGER=cat
export MODULES_PAGER=cat
export SCALE_SYS="${SCALE_SYS:-Linux64-intel-impi}"

if command -v module >/dev/null 2>&1; then
  module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1 >/dev/null 2>&1 || true
fi

root_dir=$(cd "$(dirname "$0")/.." && pwd); status=0
while IFS= read -r mf; do
  d=$(dirname "$mf"); echo "== ${d#$root_dir/} =="
  top=$(cd "$d/../../../../../../.." && pwd); scale=$(cd "$d/../../../../../.." && pwd)
  echo "TOPDIR=$top"; echo "SCALE_RM_DIR=$scale"
  [ -f "$top/Mkinclude" ] || { echo "ERROR: missing $top/Mkinclude" >&2; status=1; }
  [ -d "$scale/src" ] || { echo "ERROR: missing $scale/src" >&2; status=1; }
  [ -d "$top/bin" ] || echo "INFO: $top/bin is not present before build; this is expected before running build_squid.sh" >&2
  (cd "$d" && make -n allclean >/dev/null) || { echo "ERROR: make -n allclean failed in $d" >&2; status=1; }
done < <(find "$root_dir" -mindepth 3 -maxdepth 3 -name Makefile | sort)
exit "$status"
