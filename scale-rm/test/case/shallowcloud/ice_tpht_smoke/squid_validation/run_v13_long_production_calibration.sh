#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SMOKE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO_ROOT=$(cd "$SMOKE_DIR/../../../../.." && pwd)
COLD_CASE="$SMOKE_DIR/fw_event_type_smoke"
STAMP=${COLD_TPHT_VALIDATION_STAMP:-v13_long_prod_$(date +%Y%m%dT%H%M%S)}
RESULT_DIR="$SCRIPT_DIR/results/$STAMP"
LOG="$RESULT_DIR/run.log"
REPORT="$RESULT_DIR/long_production_calibration_report.md"
SUMMARY_JSON="$RESULT_DIR/long_production_summary.json"
PYTHON_BIN=${COLD_TPHT_PYTHON:-python3.11}
RANKS=${COLD_TPHT_LONG_RANKS:-${COLD_TPHT_PROD_RANKS:-4}}
FAILURES=0
MPI_LAUNCHER_TEXT=${COLD_TPHT_MPI_LAUNCHER:-"mpiexec --oversubscribe"}
read -r -a MPI_LAUNCHER <<< "$MPI_LAUNCHER_TEXT"

case "$RANKS" in
  1) INIT_CONF=init.conf ;;
  2) INIT_CONF=init_rank2.conf ;;
  4) INIT_CONF=init_rank4.conf ;;
  *)
    echo "unsupported COLD_TPHT_LONG_RANKS=$RANKS; supported: 1, 2, 4" >&2
    exit 2
    ;;
esac

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
  )
}

run_cold_case() {
  local case_name=$1
  local init_conf=$2
  local run_conf=$3
  clean_cold_case
  (
    set -e
    cd "$COLD_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n "$RANKS" ./scale-rm_init "$init_conf"
    local start
    local end
    start=$(date +%s)
    "${MPI_LAUNCHER[@]}" -n "$RANKS" ./scale-rm "$run_conf"
    end=$(date +%s)
    echo "$((end - start))" > "$RESULT_DIR/${case_name}_wallclock_sec.txt"
  )
  local status=$?
  archive_case_outputs "$case_name" "$COLD_CASE"
  return "$status"
}

record "# Cold TPHT v1.3 Longer Production Calibration"
record ""
record "- timestamp: $STAMP"
record "- workdir: $(pwd)"
record "- git_commit: $(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
record "- python: $($PYTHON_BIN --version 2>&1)"
record "- ranks: $RANKS"
record "- init_conf: $INIT_CONF"
record "- long_duration_sec: ${COLD_TPHT_LONG_DURATION_SEC:-180.0}"
record "- low_diag_duration_sec: ${COLD_TPHT_LONG_LOW_DIAG_DURATION_SEC:-30.0}"
record ""
record "## Command Status"
record ""
record "| step | exit_code |"
record "|---|---:|"

run_cmd "generate v1.3 longer production configs" "$PYTHON_BIN" "$SCRIPT_DIR/make_v13_long_production_configs.py"
if ! run_cmd "build cold smoke" bash -lc "cd '$COLD_CASE' && env SCALE_ENABLE_SDM=T make -j4"; then
  record ""
  record "Build failed; longer production calibration cases were not run."
  exit 1
fi

run_cmd "long conservative case" run_cold_case long_conservative "$INIT_CONF" run_v13_long_conservative.conf
run_cmd "long low significant case" run_cold_case long_low_significant "$INIT_CONF" run_v13_long_low_significant.conf
run_cmd "long moderate diagnostic case" run_cold_case long_moderate_diag "$INIT_CONF" run_v13_long_moderate_diag.conf
run_cmd "capped low diagnostic case" run_cold_case long_low_diag_cap "$INIT_CONF" run_v13_long_low_diag_cap.conf
run_cmd "long optional context case" run_cold_case long_context_on "$INIT_CONF" run_v13_long_context_on.conf

CASE_DIRS=(
  "$RESULT_DIR/long_conservative"
  "$RESULT_DIR/long_low_significant"
  "$RESULT_DIR/long_moderate_diag"
  "$RESULT_DIR/long_low_diag_cap"
  "$RESULT_DIR/long_context_on"
)

run_cmd "v1.3 longer I/O and threshold summary" "$PYTHON_BIN" "$SMOKE_DIR/summarize_cold_tpht_v13_long_production.py" \
  "${CASE_DIRS[@]}" \
  --json-out "$SUMMARY_JSON"

record ""
record "## Summary Artifacts"
record ""
record "- JSON summary: $SUMMARY_JSON"
for case_name in long_conservative long_low_significant long_moderate_diag long_low_diag_cap long_context_on; do
  if [ -f "$RESULT_DIR/${case_name}_wallclock_sec.txt" ]; then
    record "- ${case_name} wallclock_sec: $(cat "$RESULT_DIR/${case_name}_wallclock_sec.txt")"
  fi
done
record ""
record "## Calibration Notes"
record ""
record "- v1.3 process event files use trigger_code as process identity and trigger_level as occurrence/significant."
record "- The low-diagnostic case is duration-capped by default to avoid repeating the known near-zero diagnostic I/O explosion at full longer duration."
record "- The context case enables Kohler, aerosol, and thermodynamic context groups to estimate optional context overhead."
record "- Use trigger_pair_counts, mask bit counts, denominator-filtered mass-change metrics, absolute mass-change metrics, and Kohler margin statistics when available."
record ""
record "## Final Status"
record ""
record "- required_failures: $FAILURES"

echo
echo "LONG_PRODUCTION_REPORT=$REPORT"
echo "LONG_PRODUCTION_SUMMARY_JSON=$SUMMARY_JSON"
echo "REQUIRED_FAILURES=$FAILURES"

if [ "$FAILURES" -ne 0 ]; then
  exit 1
fi
exit 0
