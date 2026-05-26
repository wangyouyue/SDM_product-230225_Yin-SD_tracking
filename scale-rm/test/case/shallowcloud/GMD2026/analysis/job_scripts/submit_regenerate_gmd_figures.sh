#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=1
#PBS -l elapstim_req=02:00:00
#PBS -N gmd_gmd_figs
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
export GMD2026_ANALYSIS_JOB_NAME="regenerate_gmd_figures"
export GMD2026_ANALYSIS_COMMAND='bash analysis/regenerate_gmd_figures.sh --root "${GMD2026_ROOT}" --input "${GMD2026_OUTDIR}" --outdir "${GMD2026_GMD_OUTDIR}" ${GMD2026_REGENERATE_ARGS:-}'
exec "${root}/analysis/job_scripts/_run_analysis_job.sh"
