#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=8
#PBS -l elapstim_req=12:00:00
#PBS -N gmd02_lpar
set -euo pipefail
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac
outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
export GMD2026_ROOT="$root" GMD2026_OUTDIR="$outdir" PYTHONUNBUFFERED=1
export TPHT_LIGHT_WORKERS="${TPHT_LIGHT_WORKERS:-8}"
export GMD2026_ANALYSIS_JOB_NAME="02_tpht_light_parallel"
export GMD2026_ANALYSIS_COMMAND='bash analysis/02_tpht/run_02_tpht_light_parallel.sh --root "${GMD2026_ROOT}" --outdir "${GMD2026_OUTDIR}" ${GMD2026_ANALYSIS_ARGS:-}'
exec "${root}/analysis/job_scripts/_run_analysis_job.sh"
