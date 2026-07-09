#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SMOKE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO_ROOT=$(cd "$SMOKE_DIR/../../../../.." && pwd)
COLD_CASE="$SMOKE_DIR/fw_event_type_smoke"
WARM_CASE="$SMOKE_DIR/../tpht_test/ft_interest_id_baseline"
STAMP=${COLD_TPHT_VALIDATION_STAMP:-$(date +%Y%m%dT%H%M%S)}
RESULT_DIR="$SCRIPT_DIR/results/$STAMP"
LOG="$RESULT_DIR/run.log"
REPORT="$RESULT_DIR/validation_report.md"
PYTHON_BIN=${COLD_TPHT_PYTHON:-python3.11}
FAILURES=0
OPTIONAL_FAILURES=0
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

run_optional_cmd() {
  local label=$1
  shift
  echo
  echo "### optional: $label"
  echo "CMD: $*"
  "$@"
  local status=$?
  echo "OPTIONAL_EXIT: $status"
  record "| optional: $label | $status |"
  if [ "$status" -ne 0 ]; then
    OPTIONAL_FAILURES=$((OPTIONAL_FAILURES + 1))
  fi
  return 0
}

clean_cold_case() {
  (
    cd "$COLD_CASE" || exit 1
    ./clean.sh
    rm -f SD_all_NetCDF_* SD_all_history.pe* SD_selected_history.pe* SD_lifecycle_NetCDF_*
    mkdir -p fw_tracking fw_output fw_event_output
    mkdir -p restart_mid_output restart_resume_output restart_cont_output restart_legacy_output
  )
}

clean_warm_case() {
  (
    cd "$WARM_CASE" || exit 1
    ./clean.sh
    mkdir -p fw_tracking fw_output fw_event_output
  )
}

clean_stale_sdm_build_artifacts() {
  (
    set -e
    cd "$REPO_ROOT" || exit 1
    find contrib/SDM -maxdepth 1 \( -name "*.o" -o -name "*.mod" -o -name "*.a" \) -exec rm -f {} +
    rm -f lib/libsdm.a
    rm -rf "$COLD_CASE/.libs" "$WARM_CASE/.libs"
    rm -f "$COLD_CASE/scale-rm" "$COLD_CASE/scale-rm_init" "$COLD_CASE/scale-rm_pp"
    rm -f "$WARM_CASE/scale-rm" "$WARM_CASE/scale-rm_init" "$WARM_CASE/scale-rm_pp"
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
    "${MPI_LAUNCHER[@]}" -n "$ranks" ./scale-rm_init "$init_conf"
    "${MPI_LAUNCHER[@]}" -n "$ranks" ./scale-rm "$run_conf"
  )
  local status=$?
  archive_case_outputs "$case_name" "$COLD_CASE"
  return "$status"
}

run_warm_case() {
  clean_warm_case
  (
    set -e
    cd "$WARM_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n 1 ./scale-rm_init init_rank1_warm.conf
    "${MPI_LAUNCHER[@]}" -n 1 ./scale-rm run_rank1_warm.conf
  )
  local status=$?
  archive_case_outputs "warm_baseline" "$WARM_CASE"
  return "$status"
}

run_cold_model_only() {
  local ranks=$1
  local run_conf=$2
  (
    set -e
    cd "$COLD_CASE" || exit 1
    "${MPI_LAUNCHER[@]}" -n "$ranks" ./scale-rm "$run_conf"
  )
}

count_trigger_codes() {
  local label=$1
  local pattern=$2
  echo
  echo "### trigger summary: $label"
  "$PYTHON_BIN" - "$pattern" <<'PY'
import glob
import re
import subprocess
import sys
from collections import Counter
paths = sorted(glob.glob(sys.argv[1]))
counts = Counter()
for path in paths:
    out = subprocess.run(["ncdump", "-v", "trigger_code", path], text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if out.returncode != 0:
        continue
    data = out.stdout.split("data:", 1)[-1]
    counts.update(int(v) for v in re.findall(r"-?\d+", data))
print(f"files={len(paths)}")
print("trigger_code_counts=" + ",".join(f"{k}:{v}" for k, v in sorted(counts.items())))
PY
}

count_phase_states() {
  local label=$1
  local pattern=$2
  echo
  echo "### phase-state summary: $label"
  "$PYTHON_BIN" - "$pattern" <<'PY'
import glob
import re
import subprocess
import sys
from collections import Counter
fields = [
    "phase_state_pre", "phase_state_post",
    "phase_state1_pre", "phase_state1_post",
    "phase_state2_pre", "phase_state2_post",
]
paths = sorted(glob.glob(sys.argv[1]))
counts = Counter()
for path in paths:
    for field in fields:
        out = subprocess.run(["ncdump", "-v", field, path], text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if out.returncode != 0:
            continue
        data = out.stdout.split("data:", 1)[-1]
        counts.update(int(v) for v in re.findall(r"-?\d+", data))
print(f"files={len(paths)}")
print("phase_state_counts=" + ",".join(f"{k}:{v}" for k, v in sorted(counts.items())))
PY
}

record "# Cold TPHT v1.3 Validation Run"
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
run_cmd "clean stale SDM build artifacts" clean_stale_sdm_build_artifacts
if ! run_cmd "build cold smoke" bash -lc "cd '$COLD_CASE' && env SCALE_ENABLE_SDM=T make -j4"; then
  record ""
  record "Build failed; smoke validation was not run."
  exit 1
fi
if ! run_cmd "build warm baseline" bash -lc "cd '$WARM_CASE' && env SCALE_ENABLE_SDM=T make -j4"; then
  record ""
  record "Warm build failed; smoke validation was not run."
  exit 1
fi

run_cmd "warm baseline regression" run_warm_case
run_cmd "warm schema validator" bash -lc "cd '$WARM_CASE' && '$PYTHON_BIN' '$SMOKE_DIR/validate_warm_legacy_schema.py' --coal-glob 'SD_coal_output_NetCDF_*.pe*' --ordinary-glob 'SD_all_NetCDF_*.pe*'"

run_cmd "cold collision occurrence smoke" run_cold_case collision_occurrence 1 init.conf run.conf
run_cmd "cold collision occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_types.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --require-trigger-pair 2:1 --require-trigger-pair 3:1 --forbid-optional-groups"
count_trigger_codes collision_occurrence "$COLD_CASE/SD_event_collision_NetCDF_*.pe*"
count_phase_states collision_occurrence "$COLD_CASE/SD_event_collision_NetCDF_*.pe*"

run_cmd "cold collision significant smoke" run_cold_case collision_significant 1 init.conf run_collision_significant.conf
run_cmd "cold collision significant validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_types.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --require-trigger-pair 2:2 --require-trigger-pair 3:2"
count_trigger_codes collision_significant "$COLD_CASE/SD_event_collision_NetCDF_*.pe*"

run_cmd "cold coalescence occurrence seed" run_cold_case coalescence_occurrence 1 init.conf run_v13_coalescence_seed_occurrence.conf
run_cmd "cold coalescence occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_cold_liq_liq_seed.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-level 1"

run_cmd "cold coalescence significant seed" run_cold_case coalescence_significant 1 init.conf run_v13_coalescence_seed_significant.conf
run_cmd "cold coalescence significant validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_cold_liq_liq_seed.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-level 2"
run_cmd "cold coalescence significant mask validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_mask_bits.py --glob 'SD_selected_NetCDF_*.pe*' --mask-variable sd_event_sig_mask --require-bit 1"

run_cmd "cold singleproc occurrence smoke" run_cold_case singleproc_occurrence 1 init.conf run_singleproc.conf
run_cmd "cold singleproc occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_singleproc_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --require-trigger-pair 4:1 --require-trigger-pair 5:1 --forbid-extended-geometry"

run_cmd "cold singleproc significant smoke" run_cold_case singleproc_significant 1 init.conf run_singleproc_significant.conf
run_cmd "cold singleproc significant validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_singleproc_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --require-trigger-pair 8:2 --require-trigger-pair 9:2 --require-trigger-pair 4:2 --require-trigger-pair 5:2 --forbid-extended-geometry"
run_cmd "cold significant event mask validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_mask_bits.py --glob 'SD_selected_NetCDF_*.pe*' --mask-variable sd_event_sig_mask --require-bit 128 --require-bit 256 --require-bit 8 --require-bit 16"
count_trigger_codes singleproc_significant "$COLD_CASE/SD_event_singleproc_NetCDF_*.pe*"
count_phase_states singleproc_significant "$COLD_CASE/SD_event_singleproc_NetCDF_*.pe*"

run_cmd "vapor-growth occurrence smoke" run_cold_case vapor_occurrence 1 init.conf run_v12_vapor_occurrence.conf
run_cmd "vapor-growth occurrence validator: deposition covered" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_singleproc_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --require-trigger-pair 6:1"
run_cmd "vapor sublimation occurrence seed" run_cold_case vapor_sublimation_occurrence 1 init.conf run_vapor_occurrence_sublimation_seed.conf
run_cmd "vapor sublimation occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_vapor_occurrence_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-sublimation-occurrence --check-event-mask --forbid-lifecycle-records"

run_cmd "Kohler activation occurrence smoke" run_cold_case kohler_activation_occurrence 1 init.conf run_kohler_activation_occurrence.conf
run_cmd "Kohler activation occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_kohler_activation_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-activation --require-kohler-context --forbid-significant-activation --forbid-lifecycle-activation --check-event-mask --check-no-phase-change-forcing --require-condensation-overlap"
run_cmd "vapor condensation occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_vapor_occurrence_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-condensation-occurrence --check-event-mask --forbid-lifecycle-records"

run_cmd "Kohler activation significant smoke" run_cold_case kohler_activation_significant 1 init.conf run_kohler_activation_significant.conf
run_cmd "Kohler activation significant validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_kohler_activation_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-significant-activation --require-kohler-context --forbid-lifecycle-activation --check-event-mask --check-no-phase-change-forcing --require-condensation-overlap"

run_cmd "Kohler deactivation occurrence smoke" run_cold_case kohler_deactivation_occurrence 1 init.conf run_kohler_deactivation_occurrence.conf
run_cmd "Kohler deactivation occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_kohler_activation_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-deactivation --require-kohler-context --forbid-significant-activation --forbid-lifecycle-activation --check-event-mask --check-no-phase-change-forcing --require-evaporation-overlap"
run_cmd "vapor evaporation occurrence validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_vapor_occurrence_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-evaporation-occurrence --check-event-mask --forbid-lifecycle-records"

run_cmd "Kohler deactivation significant smoke" run_cold_case kohler_deactivation_significant 1 init.conf run_kohler_deactivation_significant.conf
run_cmd "Kohler deactivation significant validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_kohler_activation_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --ordinary-glob 'SD_selected_NetCDF_*.pe*' --require-significant-deactivation --require-kohler-context --forbid-lifecycle-activation --check-event-mask --check-no-phase-change-forcing --require-evaporation-overlap"

run_cmd "diagnostic threshold smoke" run_cold_case diagnostic_threshold 1 init.conf run_diag.conf
run_cmd "diagnostic threshold validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_diag_events.py --event-glob 'SD_event_diag_NetCDF_*.pe*' --require-trigger-code 101 --require-trigger-code 102 --require-trigger-code 104 --require-trigger-code 105 --require-trigger-code 106 --require-trigger-code 107"
count_trigger_codes diagnostic_threshold "$COLD_CASE/SD_event_diag_NetCDF_*.pe*"

run_cmd "ordinary schema smoke" run_cold_case ordinary_schema 1 init.conf run_ordinary.conf
run_cmd "ordinary schema validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_cold_tracking_schema.py --glob 'SD_all_NetCDF_*.pe*' --label ordinary"
run_cmd "spatial disabled validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_spatial_visit_flag.py --glob 'SD_all_NetCDF_*.pe*' --expect-zero"

run_cmd "history schema smoke" run_cold_case history_schema 1 init.conf run_hist.conf
run_cmd "history schema validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_cold_tracking_hist_schema.py --glob 'SD_all_history.pe*' --label history"

run_cmd "selected-history schema smoke" run_cold_case selected_history_schema 1 init.conf run_hist_selected.conf
run_cmd "selected-history schema validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_cold_tracking_hist_schema.py --glob 'SD_selected_history.pe*' --label selected_history"

run_cmd "ice geometry output group smoke" run_cold_case ice_geometry_group 1 init.conf run_v12_ice_geometry.conf
run_cmd "ice geometry group validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_singleproc_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --require-trigger-code 8 --require-extended-geometry"

run_cmd "rime morphology output group smoke" run_cold_case rime_morphology_group 1 init.conf run_v12_rime_morphology.conf
run_cmd "rime morphology group validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_types.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --require-trigger-code 2 --require-rime-morphology"

run_cmd "legacy extended geometry alias smoke" run_cold_case extended_geometry_alias 1 init.conf run_singleproc_extended.conf
run_cmd "legacy extended geometry alias validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_singleproc_events.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --require-trigger-code 6 --require-extended-geometry"

run_cmd "aerosol context output group smoke" run_cold_case aerosol_context_group 1 init.conf run_v13_aerosol_context.conf
run_cmd "aerosol context output group validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_optional_context_groups.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --require-aerosol --forbid-thermo --require-kohler"

run_cmd "thermo context output group smoke" run_cold_case thermo_context_group 1 init.conf run_v13_thermo_context.conf
run_cmd "thermo context output group validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_optional_context_groups.py --event-glob 'SD_event_singleproc_NetCDF_*.pe*' --forbid-aerosol --require-thermo --require-kohler"

run_cmd "spatial enabled smoke" run_cold_case spatial_enabled 1 init.conf run_spatial.conf
run_cmd "spatial enabled validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_spatial_visit_flag.py --glob 'SD_all_NetCDF_*.pe*' --expect-any"

run_cmd "spatial rank-pruning baseline rank2" run_cold_case spatial_rank_pruning_off_rank2 2 init_rank2.conf run_v13_spatial_rank_pruning_off_rank2.conf
run_cmd "spatial rank-pruning enabled rank2" run_cold_case spatial_rank_pruning_on_rank2 2 init_rank2.conf run_v13_spatial_rank_pruning_on_rank2.conf
run_cmd "spatial rank-pruning validator" "$PYTHON_BIN" "$SMOKE_DIR/validate_spatial_rank_pruning.py" --baseline-glob "$RESULT_DIR/spatial_rank_pruning_off_rank2/SD_all_NetCDF_*.pe*" --candidate-glob "$RESULT_DIR/spatial_rank_pruning_on_rank2/SD_all_NetCDF_*.pe*" --expect-any --log-glob "$LOG" --expect-active-ranks 1 --expect-inactive-ranks 1

run_cmd "restart continuous smoke" run_cold_case restart_continuous 1 init.conf run_restart_continuous.conf
mkdir -p "$RESULT_DIR/restart_continuous_exact"
cp -p "$COLD_CASE"/SD_selected_NetCDF_*.pe* "$RESULT_DIR/restart_continuous_exact"/ 2>/dev/null || true
cp -p "$COLD_CASE"/SD_all_NetCDF_*.pe* "$RESULT_DIR/restart_continuous_exact"/ 2>/dev/null || true
run_cmd "restart part1 smoke" run_cold_case restart_part1 1 init.conf run_restart_part1.conf
run_cmd "restart part2 smoke" run_cold_model_only 1 run_restart_part2.conf
archive_case_outputs "restart_part2" "$COLD_CASE"
run_cmd "restart exact validator" "$PYTHON_BIN" "$SMOKE_DIR/compare_cold_interval_fields.py" --left-glob "$RESULT_DIR/restart_continuous_exact/SD_all_NetCDF_*.pe*" --right-glob "$COLD_CASE/SD_all_NetCDF_*.pe*"

run_cmd "spatial restart continuous smoke" run_cold_case spatial_restart_continuous 1 init.conf run_v12_spatial_restart_continuous.conf
mkdir -p "$RESULT_DIR/spatial_restart_continuous_exact"
cp -p "$COLD_CASE"/SD_selected_NetCDF_*.pe* "$RESULT_DIR/spatial_restart_continuous_exact"/ 2>/dev/null || true
cp -p "$COLD_CASE"/SD_all_NetCDF_*.pe* "$RESULT_DIR/spatial_restart_continuous_exact"/ 2>/dev/null || true
run_cmd "spatial restart part1 smoke" run_cold_case spatial_restart_part1 1 init.conf run_v12_spatial_restart_part1.conf
run_cmd "spatial restart part2 smoke" run_cold_model_only 1 run_v12_spatial_restart_part2.conf
archive_case_outputs "spatial_restart_part2" "$COLD_CASE"
run_cmd "spatial restart exact validator" "$PYTHON_BIN" "$SMOKE_DIR/compare_cold_interval_fields.py" --left-glob "$RESULT_DIR/spatial_restart_continuous_exact/SD_all_NetCDF_*.pe*" --right-glob "$COLD_CASE/SD_all_NetCDF_*.pe*"
run_cmd "spatial restart nonzero validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_spatial_visit_flag.py --glob 'SD_all_NetCDF_*.pe*' --expect-any"

run_cmd "2-rank cold MPI smoke" run_cold_case mpi_rank2 2 init_rank2.conf run_rank2.conf
run_cmd "2-rank cold MPI validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_types.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --require-trigger-code 2 --require-trigger-code 3"

run_optional_cmd "4-rank cold MPI smoke" run_cold_case mpi_rank4 4 init_rank4.conf run_rank4.conf
run_optional_cmd "4-rank cold MPI validator" bash -lc "cd '$COLD_CASE' && '$PYTHON_BIN' ../validate_event_types.py --event-glob 'SD_event_collision_NetCDF_*.pe*' --require-trigger-code 2 --require-trigger-code 3"

record ""
record "## Coverage Notes"
record ""
record "- v1.3 coalescence coverage uses trigger_code=1 with trigger_level=1/2. The old v1.2 trigger_code=14 significant-coalescence code is intentionally inactive."
record "- Vapor-growth occurrence coverage uses controlled seeds for deposition/sublimation/condensation/evaporation (trigger_code=6/7/8/9, trigger_level=1). The older physical-case-dependent sub/cond/evap probe is superseded and no longer contributes optional failures."
record "- Follow-up coverage includes significant activation/deactivation (trigger_code=10/11, trigger_level=2), spatial rank-pruning on/off comparison, and aerosol/thermo/Kohler optional context schema checks."
record "- Diagnostic trigger 103 depends on mixed SDs and may remain uncovered if the smoke does not produce STAT_MIX."
record "- Phase-state codes 0, 11, and 99 are reported only if physically produced by the run; the suite does not synthesize dry-aerosol, mixed, or invalid event records."
record "- TPHT FW/BW consistency is tracked separately from this single-case schema suite."

record ""
record "## Final Status"
record ""
record "- required_failures: $FAILURES"
record "- optional_failures: $OPTIONAL_FAILURES"

echo
echo "VALIDATION_REPORT=$REPORT"
echo "REQUIRED_FAILURES=$FAILURES"
echo "OPTIONAL_FAILURES=$OPTIONAL_FAILURES"

if [ "$FAILURES" -ne 0 ]; then
  exit 1
fi
exit 0
