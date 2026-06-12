#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=4
#PBS -l elapstim_req=00:30:00
#PBS -T intmpi
#PBS -N tpht_bw_missing_ids_failfast
set -euo pipefail

cd "${PBS_O_WORKDIR:-$(pwd)}"
export GMD_CASE_GROUP="06_fortran_diagnostics_hook_smoke_squid"
export GMD_CASE_NAME="tpht_bw_missing_ids_failfast"
export GMD_RUN_INIT=0
export GMD_RUN_MAIN=1
export GMD_USE_SHARED_INIT=1
export GMD_SHARED_INIT_DIR="../shared_init_hook_2min"
export GMD_EXPECT_FAILURE=1
source ../../common/squid_model_job.sh
