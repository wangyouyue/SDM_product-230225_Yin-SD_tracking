#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=64
#PBS -l elapstim_req=12:00:00
#PBS -T intmpi
#PBS -N nt_coallog
set -euo pipefail

cd "${PBS_O_WORKDIR:-$(pwd)}"
export GMD_CASE_GROUP="01_bench_3d_samp_30min"
export GMD_CASE_NAME="nt_coallog"
export GMD_RUN_INIT=0
export GMD_RUN_MAIN=1
export GMD_USE_SHARED_INIT=1
export GMD_SHARED_INIT_DIR="../nt_nolog"
export GMD_EXPECT_FAILURE=0
source ../../common/squid_model_job.sh
