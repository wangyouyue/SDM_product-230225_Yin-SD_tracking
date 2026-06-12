#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
source ../common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" prepare_hook_smoke_cases.py

(
  cd fw_all_hook_2min
  bash ../../common/build_squid.sh
)

topdir=$(cd ../../../../../.. && pwd)
if [ -x "${topdir}/bin/scale-rm" ] && [ -x "${topdir}/bin/scale-rm_init" ]; then
  exe_dir="${topdir}/bin"
elif [ -x "fw_all_hook_2min/scale-rm" ] && [ -x "fw_all_hook_2min/scale-rm_init" ]; then
  exe_dir="$(pwd)/fw_all_hook_2min"
else
  echo "Could not find built scale-rm and scale-rm_init executables." >&2
  echo "Checked ${topdir}/bin and $(pwd)/fw_all_hook_2min." >&2
  exit 127
fi

staging_dir="$(pwd)/.hook_executables"
rm -rf "$staging_dir"
mkdir -p "$staging_dir"
for exe in scale-rm scale-rm_init; do
  cp -f "${exe_dir}/${exe}" "$staging_dir/"
  chmod +x "${staging_dir}/${exe}"
done
exe_dir="$staging_dir"

for case_dir in \
  shared_init_hook_2min \
  nt_nolog_hook_2min \
  fw_all_hook_2min \
  fw_all_unhooked_2min \
  bw_selected_hook_2min \
  tpht_ids_hook_2min \
  tpht_bw_hook_2min \
  tpht_bw_missing_ids_failfast; do
  (
    cd "$case_dir"
    for exe in scale-rm scale-rm_init; do
      rm -f "./${exe}"
      cp -f "${exe_dir}/${exe}" .
      chmod +x "./${exe}"
    done
  )
done
