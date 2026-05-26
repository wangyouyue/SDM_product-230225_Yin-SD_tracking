#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=1
#PBS -l elapstim_req=00:10:00
#PBS -N gmd_collect_failures
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PBS_O_WORKDIR:-$script_dir}"
gmd_root="$(cd .. && pwd)"
cd "$gmd_root"
mkdir -p job_status
source common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" common/collect_job_failures.py > job_status/collect_job_failures.log 2>&1
