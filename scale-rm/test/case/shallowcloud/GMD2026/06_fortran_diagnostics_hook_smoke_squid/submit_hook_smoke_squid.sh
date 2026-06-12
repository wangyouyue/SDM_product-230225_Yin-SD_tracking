#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
source ../common/submit_utils.sh

require_case_executables() {
  local case_dir="$1"
  for exe in scale-rm scale-rm_init; do
    if [ ! -x "${case_dir}/${exe}" ]; then
      echo "Missing ${case_dir}/${exe}." >&2
      echo "Run 'bash build_hook_smoke_squid.sh' before submitting, or use PREPARE_CASES=1 bash submit_hook_smoke_squid.sh." >&2
      exit 127
    fi
  done
}

case_dirs=(
  shared_init_hook_2min
  nt_nolog_hook_2min
  fw_all_hook_2min
  fw_all_unhooked_2min
  bw_selected_hook_2min
  tpht_ids_hook_2min
  tpht_bw_missing_ids_failfast
  tpht_bw_hook_2min
)

if [ "${PREPARE_CASES:-0}" = "1" ]; then
  bash build_hook_smoke_squid.sh
fi

for case_dir in "${case_dirs[@]}"; do
  require_case_executables "$case_dir"
done

jobs=()
last_job=""
for case_dir in "${case_dirs[@]:0:7}"; do
  if [ -z "$last_job" ]; then
    job=$(submit_job_in_dir "$case_dir" squid_run.sh)
  else
    job=$(submit_after_in_dir "$last_job" "$case_dir" squid_run.sh)
  fi
  jobs+=("$job")
  last_job="$job"
done

merge_job=$(submit_after_in_dir "$last_job" postprocess serial_merge_hook_ids.sh)
tpht_bw_job=$(submit_after_in_dir "$merge_job" tpht_bw_hook_2min squid_run.sh)
jobs+=("$merge_job" "$tpht_bw_job")
check_job=$(submit_after_in_dir "$tpht_bw_job" postprocess serial_check_hook.sh)
failure_report=$(submit_after_in_dir "$check_job" ../common serial_collect_failures.sh)

printf "Hook smoke chained jobs: %s\n" "${jobs[*]}"
echo "TPHT merge job: $merge_job"
echo "TPHT BW job: $tpht_bw_job"
echo "Hook diagnostic checker job: $check_job"
echo "Failure-summary job: $failure_report"
