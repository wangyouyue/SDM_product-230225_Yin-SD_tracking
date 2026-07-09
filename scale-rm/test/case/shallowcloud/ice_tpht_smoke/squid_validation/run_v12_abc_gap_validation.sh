#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SMOKE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO_ROOT=$(cd "$SMOKE_DIR/../../../../.." && pwd)
FW_CASE="$SMOKE_DIR/fw_event_type_smoke"
BW_CASE="$SMOKE_DIR/bw_event_type_smoke"
STAMP=${COLD_TPHT_VALIDATION_STAMP:-abc_$(date +%Y%m%dT%H%M%S)}
RESULT_DIR="$SCRIPT_DIR/results/$STAMP"
LOG="$RESULT_DIR/run.log"
REPORT="$RESULT_DIR/validation_report.md"
PYTHON_BIN=${COLD_TPHT_PYTHON:-python3.11}
FAILURES=0
MPI_LAUNCHER_TEXT=${COLD_TPHT_MPI_LAUNCHER:-"mpiexec --oversubscribe"}
read -r -a MPI_LAUNCHER <<< "$MPI_LAUNCHER_TEXT"

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

clean_fw_case() {
  (
    cd "$FW_CASE" || exit 1
    ./clean.sh
    rm -f SD_all_NetCDF_* SD_all_history.pe* SD_selected_history.pe* SD_lifecycle_NetCDF_*
    mkdir -p fw_tracking fw_output fw_event_output
    mkdir -p restart_mid_output restart_resume_output restart_cont_output restart_legacy_output
  )
}

clean_bw_case() {
  (
    cd "$BW_CASE" || exit 1
    ./clean.sh
    mkdir -p bw_output bw_event_output bw_tpht_output
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
    cp -p fw_tracking/*.ids "$dest"/ 2>/dev/null || true
    cp -p tpht_analysis/* "$dest"/ 2>/dev/null || true
  )
}

run_fw_case() {
  local case_name=$1
  local ranks=$2
  local init_conf=$3
  local run_conf=$4
  clean_fw_case
  (
    set -e
    cd "$FW_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n "$ranks" ./scale-rm_init "$init_conf"
    "${MPI_LAUNCHER[@]}" -n "$ranks" ./scale-rm "$run_conf"
  )
  local status=$?
  archive_case_outputs "$case_name" "$FW_CASE"
  return "$status"
}

run_fw_model_only() {
  local ranks=$1
  local run_conf=$2
  (
    set -e
    cd "$FW_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n "$ranks" ./scale-rm "$run_conf"
  )
}

save_current_ordinary_outputs() {
  local dest=$1
  mkdir -p "$dest"
  cp -p "$FW_CASE"/SD_all_NetCDF_*.pe* "$dest"/ 2>/dev/null || true
  local count
  count=$(find "$dest" -maxdepth 1 -name 'SD_all_NetCDF_*.pe*' -type f | wc -l | tr -d ' ')
  echo "ordinary_output_files_saved=$count"
  if [ "$count" -eq 0 ]; then
    echo "ERROR: no ordinary SD_all_NetCDF files were saved to $dest" >&2
    return 1
  fi
  return 0
}

run_tpht_chain() {
  clean_fw_case
  clean_bw_case
  (
    set -e
    cd "$FW_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n 1 ./scale-rm_init init.conf
    "${MPI_LAUNCHER[@]}" -n 1 ./scale-rm run.conf
  ) || return 1

  "$PYTHON_BIN" "$SMOKE_DIR/validate_event_types.py" \
    --event-glob "$FW_CASE/SD_event_collision_NetCDF_*.pe*" \
    --require-trigger-code 2 \
    --require-trigger-code 3 || return 1

  "$PYTHON_BIN" "$SMOKE_DIR/validate_cold_tracking_schema.py" \
    --glob "$FW_CASE/SD_selected_NetCDF_*.pe*" \
    --label fw_tpht_selected || return 1

  "$PYTHON_BIN" "$SMOKE_DIR/merge_tracking_interest_ids.py" \
    --input-glob "$FW_CASE/fw_tracking/tracking_interest_ids.pe*.ids" \
    --output "$FW_CASE/fw_tracking/tracking_interest_ids_merged.ids" \
    --rank-bucket-basename "$FW_CASE/fw_tracking/tracking_interest_ids_dedup" || return 1

  (
    set -e
    cd "$BW_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n 1 ../fw_event_type_smoke/scale-rm_init init.conf
    "${MPI_LAUNCHER[@]}" -n 1 ../fw_event_type_smoke/scale-rm run_tpht.conf
  ) || return 1

  "$PYTHON_BIN" "$SMOKE_DIR/check_tpht_consistency.py" \
    --fw-ids "$FW_CASE/fw_tracking/tracking_interest_ids_merged.ids" \
    --bw-glob "$BW_CASE/SD_selected_NetCDF_*.pe*" || return 1

  "$PYTHON_BIN" "$SMOKE_DIR/validate_cold_tracking_schema.py" \
    --glob "$BW_CASE/SD_selected_NetCDF_*.pe*" \
    --label bw_tpht_selected || return 1

  "$PYTHON_BIN" "$SMOKE_DIR/analyze_tpht_tracks.py" \
    --bw-glob "$BW_CASE/SD_selected_NetCDF_*.pe*" \
    --event-glob "$BW_CASE/SD_event_collision_NetCDF_*.pe*" \
    --output-dir "$BW_CASE/tpht_analysis" || return 1

  archive_case_outputs "abc_fw_tpht" "$FW_CASE"
  archive_case_outputs "abc_bw_tpht" "$BW_CASE"
  return 0
}

record "# Cold TPHT v1.3 ABC Gap Validation"
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
run_cmd "build FW cold smoke" bash -lc "cd '$FW_CASE' && env SCALE_ENABLE_SDM=T make -j4"

run_cmd "A: FW/BW/TPHT consistency chain" run_tpht_chain
run_cmd "B: FW collision schema sweep" "$PYTHON_BIN" "$SMOKE_DIR/validate_v12_schema_sweep.py" \
  --glob "$FW_CASE/SD_event_collision_NetCDF_*.pe*" \
  --require-trigger-code 2 \
  --require-trigger-code 3
run_cmd "B: BW collision schema sweep" "$PYTHON_BIN" "$SMOKE_DIR/validate_v12_schema_sweep.py" \
  --glob "$BW_CASE/SD_event_collision_NetCDF_*.pe*" \
  --require-trigger-code 2 \
  --require-trigger-code 3

run_cmd "B: significant collision smoke" run_fw_case abc_collision_significant 1 init.conf run_collision_significant.conf
run_cmd "B: significant collision schema sweep" "$PYTHON_BIN" "$SMOKE_DIR/validate_v12_schema_sweep.py" \
  --glob "$FW_CASE/SD_event_collision_NetCDF_*.pe*" \
  --require-trigger-pair 2:2 \
  --require-trigger-pair 3:2
run_cmd "B: significant singleproc smoke" run_fw_case abc_singleproc_significant 1 init.conf run_singleproc_significant.conf
run_cmd "B: significant singleproc schema sweep" "$PYTHON_BIN" "$SMOKE_DIR/validate_v12_schema_sweep.py" \
  --glob "$FW_CASE/SD_event_singleproc_NetCDF_*.pe*" \
  --require-trigger-pair 8:2 \
  --require-trigger-pair 9:2 \
  --require-trigger-pair 4:2 \
  --require-trigger-pair 5:2
run_cmd "B: diagnostic smoke" run_fw_case abc_diag 1 init.conf run_diag.conf
run_cmd "B: diagnostic schema sweep" "$PYTHON_BIN" "$SMOKE_DIR/validate_v12_schema_sweep.py" \
  --glob "$FW_CASE/SD_event_diag_NetCDF_*.pe*" \
  --require-trigger-code 101 \
  --require-trigger-code 102 \
  --require-trigger-code 104 \
  --require-trigger-code 105 \
  --require-trigger-code 106 \
  --require-trigger-code 107

run_cmd "C: rank-2 spatial MPI smoke" run_fw_case abc_spatial_rank2 2 init_rank2.conf run_v12_spatial_rank2.conf
run_cmd "C: rank-2 spatial flag validator" "$PYTHON_BIN" "$SMOKE_DIR/validate_spatial_visit_flag.py" \
  --glob "$FW_CASE/SD_selected_NetCDF_*.pe*" \
  --expect-any
run_cmd "C: rank-2 spatial selected schema" "$PYTHON_BIN" "$SMOKE_DIR/validate_cold_tracking_schema.py" \
  --glob "$FW_CASE/SD_selected_NetCDF_*.pe*" \
  --label spatial_rank2_selected

run_cmd "C: rank-2 restart continuous smoke" run_fw_case abc_restart_rank2_continuous 2 init_rank2.conf run_v12_restart_rank2_continuous.conf
run_cmd "C: archive rank-2 restart continuous ordinary output" save_current_ordinary_outputs "$RESULT_DIR/restart_rank2_continuous_exact"
run_cmd "C: rank-2 restart part1 smoke" run_fw_case abc_restart_rank2_part1 2 init_rank2.conf run_v12_restart_rank2_part1.conf
run_cmd "C: rank-2 restart part2 smoke" run_fw_model_only 2 run_v12_restart_rank2_part2.conf
archive_case_outputs "abc_restart_rank2_part2" "$FW_CASE"
run_cmd "C: rank-2 restart exact validator" "$PYTHON_BIN" "$SMOKE_DIR/compare_cold_interval_fields.py" \
  --left-glob "$RESULT_DIR/restart_rank2_continuous_exact/SD_all_NetCDF_*.pe*" \
  --right-glob "$FW_CASE/SD_all_NetCDF_*.pe*"

record ""
record "## Documented Coverage Gaps"
record ""
record "- Coalescence trigger_level coverage is handled by the dedicated v1.3 coalescence seed in the full validation runner; this ABC schema/TPHT run does not require trigger_code=1 level=1/2."
record "- Phase-state runtime coverage for 0, 11, and 99 is not required unless the physical smoke produces dry-aerosol, mixed, or missing/invalid event contexts."
record "- Full production-scale I/O and lifecycle invalid-slot tests are tracked outside this ABC schema/TPHT/MPI gap run."
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
