#!/usr/bin/env bash
# Render a GMD2026 TPHT 3-D diagnostic video from exported QC/QR fields.
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

analysis_outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
gmd_outdir="${GMD2026_GMD_OUTDIR:-${root}/analysis_outputs_for_GMD}"
export_file="${GMD2026_TPHT_3D_EXPORT_FILE:-${gmd_outdir}/video_sources/tpht_3d/tpht_bw_qc_qr_3d_fullgrid.nc}"
mkdir -p "${gmd_outdir}/logs" "${gmd_outdir}/videos/tpht_3d"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${gmd_outdir}/logs/matplotlib}"
mkdir -p "$MPLCONFIGDIR"

echo "GMD2026 TPHT 3-D predecessor video"
echo "root=${root}"
echo "analysis_outdir=${analysis_outdir}"
echo "gmd_outdir=${gmd_outdir}"
echo "export_file=${export_file}"
echo "python=${GMD2026_PYTHON}"
date "+start_time=%Y-%m-%dT%H:%M:%S%z"

"${GMD2026_PYTHON}" analysis/02_tpht/render_02_tpht_predecessor_3d_video.py \
  --root "$root" \
  --analysis-outdir "$analysis_outdir" \
  --gmd-outdir "$gmd_outdir" \
  --export-file "$export_file" \
  ${GMD2026_TPHT_3D_RENDER_ARGS:-}

date "+end_time=%Y-%m-%dT%H:%M:%S%z"
