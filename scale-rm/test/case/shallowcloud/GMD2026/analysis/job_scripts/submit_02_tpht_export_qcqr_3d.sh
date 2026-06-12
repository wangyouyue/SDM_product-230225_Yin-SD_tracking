#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=8
#PBS -l elapstim_req=12:00:00
#PBS -N gmd02_3dexp
set -euo pipefail
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac
analysis_outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
gmd_outdir="${GMD2026_GMD_OUTDIR:-${root}/analysis_outputs_for_GMD}"
export GMD2026_ROOT="$root" GMD2026_OUTDIR="$analysis_outdir" GMD2026_GMD_OUTDIR="$gmd_outdir" PYTHONUNBUFFERED=1
export GMD2026_ANALYSIS_JOB_NAME="02_tpht_export_qcqr_3d"
export GMD2026_ANALYSIS_COMMAND='bash analysis/02_tpht/postprocess_squid_export_tpht_qcqr_3d.sh'
exec "${root}/analysis/job_scripts/_run_analysis_job.sh"
