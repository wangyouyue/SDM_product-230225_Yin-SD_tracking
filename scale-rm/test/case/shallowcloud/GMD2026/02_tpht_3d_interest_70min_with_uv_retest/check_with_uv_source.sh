#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
scale_rm_root=$(cd ../../../../.. && pwd)
mkinit="${scale_rm_root}/src/preprocess/mod_mkinit.f90"

fail=0

check_file_contains() {
  local file="$1"
  local pattern="$2"
  local label="$3"
  if ! grep -Fq "$pattern" "$file"; then
    printf 'FAIL: %s not found in %s\n' "$label" "$file" >&2
    fail=1
  fi
}

check_file_absent() {
  local file="$1"
  local pattern="$2"
  local label="$3"
  if grep -Fq "$pattern" "$file"; then
    printf 'FAIL: stale no_uv %s remains in %s\n' "$label" "$file" >&2
    fail=1
  fi
}

check_file_contains "$mkinit" 'velx(k,i,j) =  3.0_RP + 4.3 * GRID_CZ(k)*1.E-3_RP' 'DYCOMS RF02 initial velx with UV'
check_file_contains "$mkinit" 'vely(k,i,j) = -9.0_RP + 5.6 * GRID_CZ(k)*1.E-3_RP' 'DYCOMS RF02 initial vely with UV'
check_file_absent "$mkinit" 'velx(k,i,j) = 0.0_RP ! 3.0_RP + 4.3 * GRID_CZ(k)*1.E-3_RP' 'initial velx'
check_file_absent "$mkinit" 'vely(k,i,j) = 0.0_RP !-9.0_RP + 5.6 * GRID_CZ(k)*1.E-3_RP' 'initial vely'

for mod_user in fw_discovery/code/mod_user.f90 bw_reconstruction/code/mod_user.f90; do
  check_file_contains "$mod_user" 'U_GEOS(k) =  3.0_RP + 4.3 * CZ(k)*1.E-3_RP' 'DYCOMS geostrophic U_GEOS with UV'
  check_file_contains "$mod_user" 'V_GEOS(k) = -9.0_RP + 5.6 * CZ(k)*1.E-3_RP' 'DYCOMS geostrophic V_GEOS with UV'
  check_file_absent "$mod_user" 'U_GEOS(k) = 0.0_RP ! 3.0_RP + 4.3 * CZ(k)*1.E-3_RP' 'U_GEOS'
  check_file_absent "$mod_user" 'V_GEOS(k) = 0.0_RP ! -9.0_RP + 5.6 * CZ(k)*1.E-3_RP' 'V_GEOS'
done

if [ "$fail" != "0" ]; then
  cat >&2 <<'EOF'

The with-UV TPHT retest requires both:
  1. scale-rm/src/preprocess/mod_mkinit.f90 uses the standard DYCOMS RF02 initial wind.
  2. fw_discovery/code/mod_user.f90 and bw_reconstruction/code/mod_user.f90 use the standard DYCOMS geostrophic wind.

Run ./apply_with_uv_source_patch.sh before rebuilding, then rerun this check.
EOF
  exit 1
fi

echo "PASS: with-UV DYCOMS initial wind and geostrophic forcing are configured."
