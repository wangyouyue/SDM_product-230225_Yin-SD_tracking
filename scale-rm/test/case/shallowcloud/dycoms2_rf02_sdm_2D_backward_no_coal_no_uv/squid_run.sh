#!/bin/csh
#PBS -q SQUID
#PBS --group=hp250136 
#PBS -m b
#PBS -b 1 
#PBS -l cpunum_job=40
#PBS -l elapstim_req=40:00:00
#PBS -T intmpi
#PBS -N bw_tracking

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}

#-----------program execution------------

module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1

# run
mpiexec -n 40 ./scale-rm_init init.conf || exit
mpiexec -n 40 ./scale-rm run.conf || exit
