#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ../common/submit_utils.sh
cases=()
for c in ref_full_stratified ref_full_random stratified_f*_s?? random_f*_s??; do
  [ -d "$c" ] && cases+=("$c")
done
mapfile -t jobs < <(submit_chain_in_dirs squid_run.sh "${cases[@]}")
last_job="${jobs[$((${#jobs[@]} - 1))]}"
post=$(submit_after_in_dir "$last_job" postprocess parallel_evaluate.sh)
failure_report=$(submit_after_in_dir "$post" ../common serial_collect_failures.sh)
printf "Model jobs: %s\n" "${jobs[*]}"
echo "Postprocess job: $post"
echo "Failure-summary job: $failure_report"
