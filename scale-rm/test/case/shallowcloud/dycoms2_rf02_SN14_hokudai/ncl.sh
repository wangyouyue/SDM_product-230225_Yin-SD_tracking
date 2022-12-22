#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=ha
#PJM -L node=10
#PJM --mpi proc=400
#PJM -L elapse=10:00:00
#PJM -g r22775
#PJM -j
#------- Program execution -------#
module load intel impi hdf5 netcdf netcdf-fortran
ulimit -u 2048
export OPENBLAS_NUM_THREADS=1
source ~/.bashrc
conda activate ncl_stable
# run
rm -rf SN14_history_*.ncl
MPI=3600
core_num=400
let p=MPI/core_num
i=0
while(( $i<$(($core_num-1)) ))
do
    cp ./SN14_history.ncl ./SN14_history_$(($i*$p)).ncl
    sed -i "s/0,MPI/$(($i*$p)),$((($i+1)*$p))/" ./SN14_history_$(($i*$p)).ncl
    ncl SN14_history_$(($i*$p)).ncl &
    let "i++"
done
cp ./SN14_history.ncl ./SN14_history_$(($i*$p)).ncl
sed -i "s/0,MPI/$(($i*$p)),$((($i+1)*$p))/" ./SN14_history_$(($i*$p)).ncl
ncl SN14_history_$(($i*$p)).ncl
wait
