#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=4
#PBS -l elapstim_req=03:00:00
#PBS -N gmd_analysis_p
set -euo pipefail
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac
outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
workers="${GMD2026_WORKERS:-4}"
export GMD2026_ROOT="$root" GMD2026_OUTDIR="$outdir" PYTHONUNBUFFERED=1
export GMD2026_ANALYSIS_JOB_NAME="parallel_analysis"
export GMD2026_ANALYSIS_COMMAND='bash analysis/run_all_analysis.sh --root "${GMD2026_ROOT}" --outdir "${GMD2026_OUTDIR}" --workers '"${workers}"' ${GMD2026_ANALYSIS_ARGS:-}'
exec "${root}/analysis/job_scripts/_run_analysis_job.sh"
