#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ../common/submit_utils.sh
fw_job=$(submit_job_in_dir fw_discovery squid_run.sh)
merge_job=$(submit_after_in_dir "$fw_job" postprocess serial_merge_ids.sh)
bw_job=$(submit_after_in_dir "$merge_job" bw_reconstruction squid_run.sh)
check_job=$(submit_after_in_dir "$bw_job" postprocess serial_check_tpht.sh)
failure_report=$(submit_after_in_dir "$check_job" ../common serial_collect_failures.sh)
echo "FW job: $fw_job"
echo "Merge job: $merge_job"
echo "BW job: $bw_job"
echo "Check job: $check_job"
echo "Failure-summary job: $failure_report"
