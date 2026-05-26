#!/usr/bin/env bash
# Common SQUID model runner for all GMD2026 model cases.
# Source this file from a case-level squid_run.sh.

set -uo pipefail

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
  export PAGER=cat
  export LESS="${LESS:-FRX}"
  export LMOD_PAGER=cat
  export MODULES_PAGER=cat
  module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1
else
  echo "WARNING: /etc/profile.d/modules.sh not found; assuming a local test environment." >&2
fi

cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs
start_epoch=$(date +%s)
status=0

if [ "${CLEAN_OUTPUTS:-0}" = "1" ]; then
  ./clean_outputs.sh
  mkdir -p logs
fi

repo_root=$(cd ../../../../../../.. && pwd)
for exe in scale-rm scale-rm_init; do
  if [ "$status" != "0" ]; then
    break
  fi
  if [ ! -x "./${exe}" ]; then
    if [ ! -x "${repo_root}/bin/${exe}" ]; then
      echo "Missing executable: ${repo_root}/bin/${exe}" >&2
      status=127
      break
    fi
    cp -f "${repo_root}/bin/${exe}" .
    chmod +x "./${exe}"
  fi
done

makefile_ranks=$(
  awk '$1 == "TPROC" && $2 == "=" {print $3; exit}' Makefile 2>/dev/null
)
mpi_ranks="${MPI_RANKS:-${makefile_ranks:-1}}"
omp_threads="${OMP_NUM_THREADS:-1}"

run_with_time() {
  local log_file="$1"
  shift
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -v -o "$log_file" "$@"
  else
    "$@"
  fi
}

extract_namelist_value() {
  local key="$1"
  local file="${2:-run.conf}"
  [ -f "$file" ] || return 0
  awk -F= -v key="$key" '
    {
      lhs = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", lhs)
      if (lhs == key) {
        value = $2
        sub(/!.*/, "", value)
        sub(/,.*/, "", value)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
        gsub(/^"/, "", value)
        gsub(/"$/, "", value)
        print value
        exit
      }
    }
  ' "$file"
}

ensure_parent_dir() {
  local basename_value="$1"
  local parent_dir
  [ -n "$basename_value" ] || return 0
  parent_dir=$(dirname "$basename_value")
  if [ -n "$parent_dir" ] && [ "$parent_dir" != "." ]; then
    mkdir -p "$parent_dir"
  fi
}

prepare_configured_output_dirs() {
  local config_file key value
  for config_file in init.conf run.conf; do
    [ -f "$config_file" ] || continue
    for key in \
      ATMOS_RESTART_OUT_BASENAME \
      LAND_RESTART_OUT_BASENAME \
      OCEAN_RESTART_OUT_BASENAME \
      URBAN_RESTART_OUT_BASENAME \
      ATMOS_BOUNDARY_OUT_BASENAME \
      RANDOM_OUT_BASENAME \
      SD_OUT_BASENAME \
      HISTORY_DEFAULT_BASENAME \
      tracking_id_output_basename; do
      value=$(extract_namelist_value "$key" "$config_file")
      ensure_parent_dir "$value"
    done
  done
}

prepare_tpht_id_output() {
  local output_base
  output_base=$(extract_namelist_value tracking_id_output_basename run.conf)
  if [ -n "$output_base" ]; then
    ensure_parent_dir "$output_base"
    rm -f "${output_base}".pe*.ids "${output_base}"_merged.ids "${output_base}"_dedup.pe*.ids
  fi
}

validate_tracking_id_input() {
  local input_base rank rank_file missing
  input_base=$(extract_namelist_value tracking_id_input_basename run.conf)
  [ -n "$input_base" ] || return 0

  if [[ ! "$mpi_ranks" =~ ^[0-9]+$ ]] || [ "$mpi_ranks" -le 0 ]; then
    echo "Invalid MPI_RANKS/TPROC value for TPHT ID preflight: ${mpi_ranks}" >&2
    return 127
  fi

  missing=0
  for ((rank = 0; rank < mpi_ranks; rank++)); do
    printf -v rank_file "%s.pe%06d.ids" "$input_base" "$rank"
    if [ ! -f "$rank_file" ]; then
      echo "Missing TPHT BW ID handoff file: ${rank_file}" >&2
      missing=1
    fi
  done

  if [ "$missing" != "0" ]; then
    echo "Refusing to start BW reconstruction because tracking_id_input_basename is incomplete: ${input_base}" >&2
    return 127
  fi
}

copy_shared_init() {
  local shared_init_dir="${GMD_SHARED_INIT_DIR:-../shared_init_hook_2min}"
  shopt -s nullglob
  local init_files=("${shared_init_dir}"/init_dycom_rf02_*.nc)
  local random_files=("${shared_init_dir}"/random_number_init.pe*)
  shopt -u nullglob

  if [ "${#init_files[@]}" -eq 0 ] || [ "${#random_files[@]}" -eq 0 ]; then
    echo "Missing shared init files in ${shared_init_dir}. Submit the configured shared-init provider case first." >&2
    return 127
  fi
  cp -f "${init_files[@]}" .
  cp -f "${random_files[@]}" .
}

if [ "$status" = "0" ] && [ "${GMD_USE_SHARED_INIT:-0}" = "1" ]; then
  copy_shared_init || status=$?
fi

if [ "$status" = "0" ]; then
  prepare_configured_output_dirs || status=$?
fi

if [ "$status" = "0" ] && [ "${GMD_RUN_INIT:-1}" = "1" ]; then
  run_with_time time_scale_rm_init.log mpiexec -n "$mpi_ranks" ./scale-rm_init init.conf || status=$?
fi

if [ "$status" = "0" ] && [ "${GMD_RUN_MAIN:-1}" = "1" ]; then
  prepare_tpht_id_output
  validate_tracking_id_input || status=$?
fi

if [ "$status" = "0" ] && [ "${GMD_RUN_MAIN:-1}" = "1" ]; then
  run_with_time time_scale_rm_main.log mpiexec -n "$mpi_ranks" ./scale-rm run.conf || status=$?
fi

end_epoch=$(date +%s)
wallclock_s=$((end_epoch - start_epoch))
node_count="${NQSV_NODE_NUM:-${PBS_NUM_NODES:-0}}"
pbs_jobid="${PBS_JOBID:-unknown}"
case_name="${GMD_CASE_NAME:-$(basename "$(pwd)")}"
case_group="${GMD_CASE_GROUP:-$(basename "$(dirname "$(pwd)")")}"
gmd_root="${GMD_ROOT:-$(cd ../.. && pwd)}"
status_dir="${GMD_STATUS_DIR:-${gmd_root}/job_status}"
mkdir -p "$status_dir"

status_label="PASS"
if [ "${GMD_EXPECT_FAILURE:-0}" = "1" ] && [ "$status" != "0" ]; then
  status_label="EXPECTED_FAILURE"
elif [ "${GMD_EXPECT_FAILURE:-0}" = "1" ] && [ "$status" = "0" ]; then
  status_label="UNEXPECTED_PASS"
elif [ "$status" != "0" ]; then
  status_label="FAIL"
fi

cat > job_metrics.json <<JSON
{"case_name":"${case_name}","pbs_jobid":"${pbs_jobid}","start_epoch":${start_epoch},"end_epoch":${end_epoch},"wallclock_s":${wallclock_s},"mpi_ranks":${mpi_ranks},"omp_threads":${omp_threads},"node_count":${node_count},"exit_status":${status}}
JSON

cat > "${status_dir}/${case_group}__${case_name}.json" <<JSON
{"case_group":"${case_group}","case_name":"${case_name}","case_dir":"$(pwd)","pbs_jobid":"${pbs_jobid}","start_epoch":${start_epoch},"end_epoch":${end_epoch},"wallclock_s":${wallclock_s},"mpi_ranks":${mpi_ranks},"omp_threads":${omp_threads},"node_count":${node_count},"exit_status":${status},"expected_failure":${GMD_EXPECT_FAILURE:-0},"status":"${status_label}","job_metrics":"$(pwd)/job_metrics.json"}
JSON

events_file="${status_dir}/job_events.tsv"
append_job_event() {
  if [ ! -f "$events_file" ]; then
    printf "case_group\tcase_name\tstatus\texit_status\texpected_failure\tpbs_jobid\twallclock_s\tcase_dir\n" > "$events_file"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$case_group" "$case_name" "$status_label" "$status" "${GMD_EXPECT_FAILURE:-0}" \
    "$pbs_jobid" "$wallclock_s" "$(pwd)" >> "$events_file"
}

if command -v flock >/dev/null 2>&1; then
  (
    flock 9
    append_job_event
  ) 9>"${status_dir}/job_events.lock"
else
  append_job_event
fi

if [ "${GMD_EXPECT_FAILURE:-0}" = "1" ]; then
  [ "$status" = "0" ] && exit 1 || exit 0
fi
exit "$status"
