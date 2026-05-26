#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=1
#PBS -l elapstim_req=04:00:00
#PBS -N gmd02_hist
set -euo pipefail
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac
outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
export GMD2026_ROOT="$root" GMD2026_OUTDIR="$outdir" PYTHONUNBUFFERED=1
export GMD2026_ANALYSIS_JOB_NAME="02_tpht_history_sanity"
export GMD2026_ANALYSIS_COMMAND='bash analysis/02_tpht/run_02_tpht_history_sanity.sh --root "${GMD2026_ROOT}" --outdir "${GMD2026_OUTDIR}" ${GMD2026_HISTORY_SANITY_ARGS:-}'
exec "${root}/analysis/job_scripts/_run_analysis_job.sh"
