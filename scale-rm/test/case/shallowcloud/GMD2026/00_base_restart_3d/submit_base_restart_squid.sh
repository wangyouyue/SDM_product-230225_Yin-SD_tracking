#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ../common/submit_utils.sh
job=$(submit_job_in_dir base_restart_3d squid_run.sh)
failure_report=$(submit_after_in_dir "$job" ../common serial_collect_failures.sh)
echo "Base restart job: $job"
echo "Failure-summary job: $failure_report"
