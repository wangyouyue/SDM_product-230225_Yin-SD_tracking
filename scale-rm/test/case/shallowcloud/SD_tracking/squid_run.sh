#!/bin/csh
#PBS -q SQUID
#PBS --group=hp250136 
#PBS -m b
#PBS -b 10 
#PBS -l cpunum_job=40
#PBS -l elapstim_req=10:00:00
#PBS -T intmpi
#PBS -N fw_tracking

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}

#-----------program execution------------

module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1

# run
mpiexec -n 400 ./scale-rm_init init.conf || exit
mpiexec -n 400 ./scale-rm run.conf || exit
