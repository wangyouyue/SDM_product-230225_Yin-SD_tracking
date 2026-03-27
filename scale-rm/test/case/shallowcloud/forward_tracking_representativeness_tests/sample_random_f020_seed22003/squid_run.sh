#!/bin/csh
#PBS -q SQUID
#PBS --group=hp250136 
#PBS -m b
#PBS -b 1 
#PBS -l cpunum_job=40
#PBS -l elapstim_req=40:00:00
#PBS -T intmpi
#PBS -N sample_random_f020_seed22003
#PBS -M yinchongzhi@gmail.com

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}

#-----------program execution------------

module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1

# Run-time diagnostics for performance comparison
set RUN_START_EPOCH=`date +%s`
set INIT_TIMED = 0
if ( -x /usr/bin/time ) then
  /usr/bin/time -v -o time_scale_rm_init.log mpiexec -n 40 ./scale-rm_init init.conf
  if ( $status == 0 ) then
    set INIT_TIMED = 1
  else
    echo "[perf] warning: /usr/bin/time -v is unavailable, fallback to direct run."
  endif
endif
if ( $INIT_TIMED == 0 ) then
  mpiexec -n 40 ./scale-rm_init init.conf
endif
if ( $status != 0 ) exit $status

set MAIN_TIMED = 0
if ( -x /usr/bin/time ) then
  /usr/bin/time -v -o time_scale_rm_main.log mpiexec -n 40 ./scale-rm run.conf
  if ( $status == 0 ) then
    set MAIN_TIMED = 1
  else
    echo "[perf] warning: /usr/bin/time -v is unavailable, fallback to direct run."
  endif
endif
if ( $MAIN_TIMED == 0 ) then
  mpiexec -n 40 ./scale-rm run.conf
endif
if ( $status != 0 ) exit $status

set RUN_END_EPOCH=`date +%s`
@ RUN_ELAPSED_SEC = $RUN_END_EPOCH - $RUN_START_EPOCH
echo "[perf] total_elapsed_sec=${RUN_ELAPSED_SEC}"

# Try scheduler-level resource statistics if PBS tools are available
if ( $?PBS_JOBID ) then
  echo "=== Scheduler resource usage (PBS) ==="
  if ( `which qstat >& /dev/null; echo $status` == 0 ) then
    qstat -fx "${PBS_JOBID}" | egrep "resources_used|Resource_List"
  else if ( `which tracejob >& /dev/null; echo $status` == 0 ) then
    tracejob -n 1 "${PBS_JOBID}" | egrep -i "resources_used|mem|vmem|walltime|cput"
  else
    echo "qstat/tracejob is not available in this environment."
  endif
endif
