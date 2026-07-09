#!/bin/bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -b 1
#PBS -l cpunum_job=8
#PBS -l elapstim_req=04:00:00
#PBS -T intmpi
#PBS -N CTPHTv12io

set -euo pipefail

cd "$PBS_O_WORKDIR"

source /etc/profile.d/modules.sh
module load BaseCPU/2024
module load inteloneAPI/2023.2
module load hdf5/1.12.3.mpi
module load netcdf-c/4.9.2
module load netcdf-fortran/4.6.1

export OMP_NUM_THREADS=1
export I_MPI_PIN=0
export COLD_TPHT_PYTHON=python3.11
export COLD_TPHT_PROD_RANKS="${COLD_TPHT_PROD_RANKS:-4}"
export COLD_TPHT_PROD_DURATION_SEC="${COLD_TPHT_PROD_DURATION_SEC:-30.0}"
export COLD_TPHT_VALIDATION_STAMP="prod_io_${PBS_JOBID:-manual}_$(date +%Y%m%dT%H%M%S)"

bash ./run_v12_production_io_threshold.sh
