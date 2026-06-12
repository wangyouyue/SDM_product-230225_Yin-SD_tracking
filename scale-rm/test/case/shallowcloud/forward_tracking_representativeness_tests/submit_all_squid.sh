#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

submitted=0
for job in */squid_run.sh; do
  [ -f "${job}" ] || continue
  case_dir="${job%/squid_run.sh}"
  (
    cd "${case_dir}"
    qsub squid_run.sh
  )
  submitted=$((submitted + 1))
done

echo "submitted_squid_jobs=${submitted}"
