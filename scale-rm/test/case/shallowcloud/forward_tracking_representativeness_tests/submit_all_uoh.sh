#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

submitted=0
for job in */UoH_run.pbs; do
  [ -f "${job}" ] || continue
  case_dir="${job%/UoH_run.pbs}"
  (
    cd "${case_dir}"
    qsub UoH_run.pbs
  )
  submitted=$((submitted + 1))
done

echo "submitted_uoh_jobs=${submitted}"
