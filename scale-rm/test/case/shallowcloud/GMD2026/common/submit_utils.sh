#!/usr/bin/env bash
set -euo pipefail

extract_qsub_jobid() {
  awk '
    {
      for (i = 1; i <= NF; i++) {
        token = $i
        gsub(/[,;]$/, "", token)
        if (token ~ /^[0-9]+([.][A-Za-z0-9_.-]+)?$/) {
          print token
          exit
        }
      }
    }
  '
}

submit_job() {
  local script="$1"
  local output jobid status
  set +e
  output=$(qsub "$script" 2>&1)
  status=$?
  set -e
  if [ "$status" -ne 0 ]; then
    printf "%s\n" "$output" >&2
    return "$status"
  fi
  jobid=$(printf "%s\n" "$output" | extract_qsub_jobid)
  if [ -z "$jobid" ]; then
    printf "Could not parse qsub job id from output:\n%s\n" "$output" >&2
    return 1
  fi
  printf "%s\n" "$jobid"
}

submit_after() {
  local dep="$1"
  local script="$2"
  local output jobid status
  set +e
  output=$(qsub --after "$dep" "$script" 2>&1)
  status=$?
  set -e
  jobid=$(printf "%s\n" "$output" | extract_qsub_jobid)
  if [ "$status" -eq 0 ] && [ -n "$jobid" ]; then
    printf "%s\n" "$jobid"
    return 0
  fi

  # SQUID/NQSV installations differ in where they expect --after relative to
  # the jobscript.  Retry with the script first before reporting failure.
  local first_output="$output"
  set +e
  output=$(qsub "$script" --after "$dep" 2>&1)
  status=$?
  set -e
  jobid=$(printf "%s\n" "$output" | extract_qsub_jobid)
  if [ "$status" -ne 0 ] || [ -z "$jobid" ]; then
    printf "Could not submit dependency job. First qsub output:\n%s\nSecond qsub output:\n%s\n" \
      "$first_output" "$output" >&2
    [ "$status" -ne 0 ] && return "$status"
    return 1
  fi
  printf "%s\n" "$jobid"
}

submit_job_in_dir() {
  local dir="$1"
  local script="${2:-squid_run.sh}"
  (cd "$dir" && submit_job "$script")
}

submit_after_in_dir() {
  local dep="$1"
  local dir="$2"
  local script="${3:-squid_run.sh}"
  (cd "$dir" && submit_after "$dep" "$script")
}

submit_chain_in_dirs() {
  local script="${1:-squid_run.sh}"
  shift
  local last_job=""
  local job=""
  for dir in "$@"; do
    if [ -z "$last_job" ]; then
      job=$(submit_job_in_dir "$dir" "$script")
    else
      job=$(submit_after_in_dir "$last_job" "$dir" "$script")
    fi
    printf "%s\n" "$job"
    last_job="$job"
  done
}
