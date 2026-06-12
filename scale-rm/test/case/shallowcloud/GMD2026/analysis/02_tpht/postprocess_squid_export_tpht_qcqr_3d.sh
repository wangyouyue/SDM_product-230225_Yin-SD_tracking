#!/usr/bin/env bash
# Export GMD2026 TPHT BW history.pe* QC/QR fields for 3-D rendering.
set -euo pipefail

if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck source=/dev/null
  source /etc/profile.d/modules.sh
fi

root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/02_tpht) root="$(cd "$root/../.." && pwd)" ;;
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac

cd "$root"

env_helper="${GMD2026_ENV_HELPER:-${root}/common/load_basepy_2026_quiet.sh}"
if [ ! -f "$env_helper" ]; then
  echo "GMD2026 environment helper is missing: ${env_helper}" >&2
  exit 2
fi
# shellcheck source=/dev/null
source "$env_helper"

gmd_outdir="${GMD2026_GMD_OUTDIR:-${root}/analysis_outputs_for_GMD}"
case_dir="${GMD2026_TPHT_BW_CASE_DIR:-${root}/02_tpht_3d_interest_70min/bw_reconstruction}"
export_outdir="${GMD2026_TPHT_3D_EXPORT_DIR:-${gmd_outdir}/video_sources/tpht_3d}"
mkdir -p "${gmd_outdir}/logs" "$export_outdir"

echo "GMD2026 TPHT 3-D QC/QR export"
echo "root=${root}"
echo "case_dir=${case_dir}"
echo "export_outdir=${export_outdir}"
echo "python=${GMD2026_PYTHON}"
date "+start_time=%Y-%m-%dT%H:%M:%S%z"

"${GMD2026_PYTHON}" analysis/02_tpht/export_02_tpht_qcqr_3d_for_render_squid.py \
  --root "$root" \
  --case-dir "$case_dir" \
  --outdir "$export_outdir" \
  ${GMD2026_TPHT_3D_EXPORT_ARGS:-}

date "+end_time=%Y-%m-%dT%H:%M:%S%z"
