#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=1
#PBS -l elapstim_req=00:30:00
#PBS -N gmd_submit_analysis
set -euo pipefail

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

export PYTHONUNBUFFERED=1
root="${GMD2026_ROOT:-${PBS_O_WORKDIR:-$(pwd)}}"
case "$root" in
  */analysis/job_scripts) root="$(cd "$root/../.." && pwd)" ;;
  */analysis) root="$(cd "$root/.." && pwd)" ;;
esac
outdir="${GMD2026_OUTDIR:-${root}/analysis_outputs}"
export GMD2026_ROOT="$root" GMD2026_OUTDIR="$outdir"

cd "$root"
mkdir -p "${outdir}/logs"
log_file="${outdir}/logs/submit_analysis_then_plots.log"
exec > >(tee -a "$log_file") 2>&1
start_epoch=$(date +%s)

extract_jobid() {
  awk '{for (i=1; i<=NF; i++) {gsub(/[,;]$/, "", $i); if ($i ~ /^[0-9]+([.][A-Za-z0-9_.-]+)?$/) {print $i; exit}}}'
}

submit_job() {
  local script="$1" output jobid
  output=$(qsub "$script")
  jobid=$(printf "%s\n" "$output" | extract_jobid)
  if [ -z "$jobid" ]; then
    echo "Could not parse job id from qsub output for ${script}: ${output}" >&2
    exit 1
  fi
  printf "%s\n" "$jobid"
}

submit_after_colon() {
  local deps="$1" script="$2" output status jobid
  set +e
  output=$(qsub --after "$deps" "$script" 2>&1)
  status=$?
  set -e
  jobid=$(printf "%s\n" "$output" | extract_jobid)
  if [ "$status" -eq 0 ] && [ -n "$jobid" ]; then
    printf "%s\n" "$jobid"
    return 0
  fi
  echo "qsub --after colon dependency failed; falling back to polling aggregator." >&2
  echo "$output" >&2
  return 1
}

submit_polling_aggregator() {
  local deps="$1" script_dir="$2" wait_script jobid output
  wait_script="${outdir}/logs/wait_then_plots_$(date +%Y%m%dT%H%M%S).sh"
  cat >"$wait_script" <<EOF
#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=1
#PBS -l elapstim_req=04:00:00
#PBS -N gmd_wait_plots
set -euo pipefail
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
export GMD2026_ROOT="${root}"
export GMD2026_OUTDIR="${outdir}"
cd "${root}"
for jobid in ${deps//:/ }; do
  while qstat "\$jobid" >/dev/null 2>&1; do
    sleep 60
  done
done
bash "${script_dir}/submit_all_plots.sh"
EOF
  output=$(qsub "$wait_script")
  jobid=$(printf "%s\n" "$output" | extract_jobid)
  if [ -z "$jobid" ]; then
    echo "Could not parse aggregator job id from qsub output: ${output}" >&2
    exit 1
  fi
  printf "%s\n" "$jobid"
}

script_dir="${root}/analysis/job_scripts"
job00=$(submit_job "${script_dir}/submit_00_restart.sh")
job01=$(submit_job "${script_dir}/submit_01_benchmark.sh")
job02l=$(submit_job "${script_dir}/submit_02_tpht_summary.sh")
job02h=$(submit_job "${script_dir}/submit_02_tpht_heavy.sh")
job03=$(submit_job "${script_dir}/submit_03_sampling.sh")
job04=$(submit_job "${script_dir}/submit_04_sdnc_scaling.sh")
job05=$(submit_job "${script_dir}/submit_05_outint_io.sh")

deps="${job00}:${job01}:${job02l}:${job02h}:${job03}:${job04}:${job05}"
if plot_job=$(submit_after_colon "$deps" "${script_dir}/submit_all_plots.sh"); then
  dependency_mode="qsub --after colon"
else
  plot_job=$(submit_polling_aggregator "$deps" "$script_dir")
  dependency_mode="polling aggregator"
fi

end_epoch=$(date +%s)
cat >"${outdir}/logs/submit_analysis_then_plots.analysis_job_metrics.json" <<EOF
{
  "job_name": "submit_analysis_then_plots",
  "start_epoch": ${start_epoch},
  "end_epoch": ${end_epoch},
  "wallclock_s": $((end_epoch - start_epoch)),
  "exit_status": 0
}
EOF

cat <<EOF
Submitted GMD2026 analysis jobs:
  00_restart: ${job00}
  01_benchmark: ${job01}
  02_tpht_light: ${job02l}
  02_tpht_heavy: ${job02h}
  03_sampling: ${job03}
  04_sdnc_scaling: ${job04}
  05_outint_io: ${job05}
  plots: ${plot_job}
Dependency mode: ${dependency_mode}
EOF
