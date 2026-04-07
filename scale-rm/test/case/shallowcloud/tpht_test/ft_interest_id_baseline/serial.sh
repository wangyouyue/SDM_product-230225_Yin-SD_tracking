#!/bin/bash
#------- qsub option -----------
#PBS -q SQUID
#PBS --group=hp250136
#PBS -l elapstim_req=24:00:00

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}

#------- Program execution -----------
module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1

# run
time ./particle_tracer_opt
