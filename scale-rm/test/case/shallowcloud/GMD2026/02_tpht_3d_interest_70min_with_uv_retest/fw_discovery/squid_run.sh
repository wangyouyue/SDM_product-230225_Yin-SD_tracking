#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=64
#PBS -l elapstim_req=12:00:00
#PBS -T intmpi
#PBS -N tpht_uv_fw
set -euo pipefail

cd "${PBS_O_WORKDIR:-$(pwd)}"
export GMD_CASE_GROUP="02_tpht_3d_interest_70min_with_uv_retest"
export GMD_CASE_NAME="fw_discovery"
export GMD_RUN_INIT=1
export GMD_RUN_MAIN=1
export GMD_USE_SHARED_INIT=0
export GMD_SHARED_INIT_DIR=""
export GMD_EXPECT_FAILURE=0
source ../../common/squid_model_job.sh
