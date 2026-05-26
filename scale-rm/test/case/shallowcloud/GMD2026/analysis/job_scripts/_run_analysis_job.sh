#!/usr/bin/env bash
# Shared execution wrapper for GMD2026 analysis qsub jobs.
set -euo pipefail

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

export PYTHONUNBUFFERED=1

job_name="${GMD2026_ANALYSIS_JOB_NAME:-analysis_job}"
root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac
outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
command="${GMD2026_ANALYSIS_COMMAND:-}"

if [ -z "$command" ]; then
  echo "GMD2026_ANALYSIS_COMMAND is not set" >&2
  exit 2
fi

cd "$root"
mkdir -p "${outdir}/logs"

exec >"${outdir}/logs/${job_name}.out" 2>"${outdir}/logs/${job_name}.err"

echo "GMD2026 analysis job: ${job_name}"
echo "root=${root}"
echo "outdir=${outdir}"
echo "command=${command}"
date "+start_time=%Y-%m-%dT%H:%M:%S%z"

env_helper="${GMD2026_ENV_HELPER:-${root}/common/load_basepy_2026_quiet.sh}"
if [ ! -f "${env_helper}" ]; then
  echo "GMD2026 environment helper is missing: ${env_helper}" >&2
  exit 2
fi
source "${env_helper}"

start_epoch=$(date +%s)
set +e
eval "$command"
status=$?
set -e
end_epoch=$(date +%s)

cat >"${outdir}/logs/${job_name}.analysis_job_metrics.json" <<EOF
{
  "job_name": "${job_name}",
  "start_epoch": ${start_epoch},
  "end_epoch": ${end_epoch},
  "wallclock_s": $((end_epoch - start_epoch)),
  "exit_status": ${status}
}
EOF

date "+end_time=%Y-%m-%dT%H:%M:%S%z"
exit "$status"
