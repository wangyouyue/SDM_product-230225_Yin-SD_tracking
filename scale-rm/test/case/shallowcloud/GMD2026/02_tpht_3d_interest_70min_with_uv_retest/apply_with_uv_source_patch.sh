#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkinit="$(cd ../../../../.. && pwd)/src/preprocess/mod_mkinit.f90"

if grep -Fq 'velx(k,i,j) =  3.0_RP + 4.3 * GRID_CZ(k)*1.E-3_RP' "$mkinit" \
  && grep -Fq 'vely(k,i,j) = -9.0_RP + 5.6 * GRID_CZ(k)*1.E-3_RP' "$mkinit"; then
  echo "mod_mkinit.f90 already uses standard DYCOMS RF02 initial UV."
  exit 0
fi

perl -0pi -e 's/velx\(k,i,j\) = 0\.0_RP ! 3\.0_RP \+ 4\.3 \* GRID_CZ\(k\)\*1\.E-3_RP/velx(k,i,j) =  3.0_RP + 4.3 * GRID_CZ(k)*1.E-3_RP/' "$mkinit"
perl -0pi -e 's/vely\(k,i,j\) = 0\.0_RP !-9\.0_RP \+ 5\.6 \* GRID_CZ\(k\)\*1\.E-3_RP/vely(k,i,j) = -9.0_RP + 5.6 * GRID_CZ(k)*1.E-3_RP/' "$mkinit"

bash ./check_with_uv_source.sh
