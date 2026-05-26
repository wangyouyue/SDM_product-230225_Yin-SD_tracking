#!/usr/bin/env bash
set -euo pipefail
rm -f LOG.pe* *.o* *.e* job_metrics.json time_scale_rm_init.log time_scale_rm_main.log
rm -f SD_selected_NetCDF_* SD_all_NetCDF_* SD_coal_output_NetCDF_* history*.nc history.pe*
rm -f monitor.pe* boundary.pe* refstate.pe* random_number_restart* superdroplet_restart*
rm -f dycom_restart_rf02* init_dycom_rf02*
rm -rf fw_output bw_output fw_tracking restart_output logs
