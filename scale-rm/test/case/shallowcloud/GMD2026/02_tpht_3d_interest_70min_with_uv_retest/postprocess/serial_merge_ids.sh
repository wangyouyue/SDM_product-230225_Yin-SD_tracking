#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=8
#PBS -l elapstim_req=12:00:00
#PBS -N tpht_uv_merge
set -euo pipefail

cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs

if [ -n "${GMD_PYTHON_VENV:-}" ]; then
  source "${GMD_PYTHON_VENV}/bin/activate"
elif [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
  module purge
  module load BasePy/2026
fi

python3 merge_tracking_interest_ids.py \
  --input-glob "../fw_discovery/fw_tracking/tracking_interest_ids.pe*.ids" \
  --output "../fw_discovery/fw_tracking/tracking_interest_ids_merged.ids" \
  --rank-bucket-basename "../fw_discovery/fw_tracking/tracking_interest_ids_dedup" \
  --workers "${TPHT_MERGE_WORKERS:-8}" \
  > logs/merge_tracking_interest_ids.log 2>&1
