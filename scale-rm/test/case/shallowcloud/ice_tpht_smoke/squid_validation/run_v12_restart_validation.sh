#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SMOKE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO_ROOT=$(cd "$SMOKE_DIR/../../../../.." && pwd)
COLD_CASE="$SMOKE_DIR/fw_event_type_smoke"
STAMP=${COLD_TPHT_VALIDATION_STAMP:-restart_$(date +%Y%m%dT%H%M%S)}
RESULT_DIR="$SCRIPT_DIR/results/$STAMP"
LOG="$RESULT_DIR/run.log"
REPORT="$RESULT_DIR/validation_report.md"
PYTHON_BIN=${COLD_TPHT_PYTHON:-python3.11}
FAILURES=0

mkdir -p "$RESULT_DIR"
exec > >(tee "$LOG") 2>&1

record() {
  printf '%s\n' "$*" >> "$REPORT"
}

run_cmd() {
  local label=$1
  shift
  echo
  echo "### $label"
  echo "CMD: $*"
  "$@"
  local status=$?
  echo "EXIT: $status"
  record "| $label | $status |"
  if [ "$status" -ne 0 ]; then
    FAILURES=$((FAILURES + 1))
  fi
  return "$status"
}

clean_cold_case() {
  (
    cd "$COLD_CASE" || exit 1
    ./clean.sh
    mkdir -p fw_tracking fw_output fw_event_output
    mkdir -p restart_mid_output restart_resume_output restart_cont_output restart_legacy_output
  )
}

archive_case_outputs() {
  local case_name=$1
  local case_dir=$2
  local dest="$RESULT_DIR/$case_name"
  mkdir -p "$dest"
  (
    cd "$case_dir" || exit 1
    cp -p LOG.pe* "$dest"/ 2>/dev/null || true
    cp -p SD_*NetCDF* "$dest"/ 2>/dev/null || true
    cp -p SD_*history* "$dest"/ 2>/dev/null || true
    cp -p restart_*_output/SD_*NetCDF* "$dest"/ 2>/dev/null || true
  )
}

run_cold_case() {
  local case_name=$1
  local ranks=$2
  local init_conf=$3
  local run_conf=$4
  clean_cold_case
  (
    set -e
    cd "$COLD_CASE" || exit 1
    mpiexec -n "$ranks" ./scale-rm_init "$init_conf"
    mpiexec -n "$ranks" ./scale-rm "$run_conf"
  )
  local status=$?
  archive_case_outputs "$case_name" "$COLD_CASE"
  return "$status"
}

ensure_cold_binaries() {
  (
    set -e
    cd "$COLD_CASE" || exit 1
    if [ ! -x ./scale-rm ] || [ ! -x ./scale-rm_init ]; then
      env SCALE_ENABLE_SDM=T make -j4
    fi
  )
}

save_current_ordinary_outputs() {
  local dest=$1
  mkdir -p "$dest"
  cp -p "$COLD_CASE"/SD_all_NetCDF_*.pe* "$dest"/ 2>/dev/null || true
  local count
  count=$(find "$dest" -maxdepth 1 -name 'SD_all_NetCDF_*.pe*' -type f | wc -l | tr -d ' ')
  echo "ordinary_output_files_saved=$count"
  if [ "$count" -eq 0 ]; then
    echo "ERROR: no ordinary SD_all_NetCDF files were saved to $dest" >&2
    return 1
  fi
  return 0
}

record "# Cold TPHT v1.2 SQUID Restart/Spatial Restart Validation"
record ""
record "- timestamp: $STAMP"
record "- workdir: $(pwd)"
record "- git_commit: $(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
record "- python: $($PYTHON_BIN --version 2>&1)"
record ""
record "## Command Status"
record ""
record "| step | exit_code |"
record "|---|---:|"

run_cmd "generate validation configs" "$PYTHON_BIN" "$SCRIPT_DIR/make_v12_validation_configs.py"
run_cmd "ensure cold smoke binaries" ensure_cold_binaries

run_cmd "restart continuous smoke" run_cold_case restart_continuous 1 init.conf run_restart_continuous.conf
run_cmd "archive restart continuous ordinary output" save_current_ordinary_outputs "$RESULT_DIR/restart_continuous_exact"
run_cmd "restart part1 smoke" run_cold_case restart_part1 1 init.conf run_restart_part1.conf
run_cmd "restart part2 smoke" bash -lc "cd '$COLD_CASE' && mpiexec -n 1 ./scale-rm run_restart_part2.conf"
archive_case_outputs "restart_part2" "$COLD_CASE"
run_cmd "restart exact validator" "$PYTHON_BIN" "$SMOKE_DIR/compare_cold_interval_fields.py" --left-glob "$RESULT_DIR/restart_continuous_exact/SD_all_NetCDF_*.pe*" --right-glob "$COLD_CASE/SD_all_NetCDF_*.pe*"

run_cmd "spatial restart continuous smoke" run_cold_case spatial_restart_continuous 1 init.conf run_v12_spatial_restart_continuous.conf
run_cmd "archive spatial restart continuous ordinary output" save_current_ordinary_outputs "$RESULT_DIR/spatial_restart_continuous_exact"
run_cmd "spatial restart part1 smoke" run_cold_case spatial_restart_part1 1 init.conf run_v12_spatial_restart_part1.conf
run_cmd "spatial restart part2 smoke" bash -lc "cd '$COLD_CASE' && mpiexec -n 1 ./scale-rm run_v12_spatial_restart_part2.conf"
archive_case_outputs "spatial_restart_part2" "$COLD_CASE"
run_cmd "spatial restart exact validator" "$PYTHON_BIN" "$SMOKE_DIR/compare_cold_interval_fields.py" --left-glob "$RESULT_DIR/spatial_restart_continuous_exact/SD_all_NetCDF_*.pe*" --right-glob "$COLD_CASE/SD_all_NetCDF_*.pe*"
run_cmd "spatial restart nonzero validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_spatial_visit_flag.py --glob 'SD_all_NetCDF_*.pe*' --expect-any"

record ""
record "## Final Status"
record ""
record "- required_failures: $FAILURES"

echo
echo "VALIDATION_REPORT=$REPORT"
echo "REQUIRED_FAILURES=$FAILURES"

if [ "$FAILURES" -ne 0 ]; then
  exit 1
fi
exit 0
