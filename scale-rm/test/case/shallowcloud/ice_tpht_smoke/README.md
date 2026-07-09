# Ice-Phase FW/BW/TPHT Smoke Cases

This directory contains short smoke cases for validating the FW/BW/TPHT tracking
path with `sdm_cold = .true.`. The default cases are rank-1; `fw_event_type_smoke`
also includes a rank-2 variant for local MPI boundary/copy coverage.

## Cases

- `fw_ice_smoke`: forward tracking with selected-SD output, TPHT interest-ID
  output, cold-SDM fields, and coalescence-event output enabled. It uses
  `tracking_selection_mode = "random"` with `tracking_fraction = 1.0d0` so the
  full valid SD set is tagged while `SD_selected_NetCDF_*` is still emitted for
  smoke inspection.
- `bw_ice_smoke`: backward tracking that consumes the deduplicated FW
  `tracking_interest_ids_dedup.peXXXXXX.ids` handoff.
- `fw_event_type_smoke`: forward tracking with a default-off synthetic cold SDM
  event seed enabled in the case namelist. It keeps only a small `1 liquid + 3
  ice` SD set active so cold coalescence writes both `trigger_code=2` riming and
  `trigger_code=3` aggregation records.
- `bw_event_type_smoke`: backward tracking variants for the same synthetic cold
  event seed. `run.conf` is a standalone BW event-type smoke. `run_tpht.conf`
  consumes the deduplicated FW event IDs and is used for TPHT consistency.
- `fw_event_type_smoke/run_singleproc.conf`: targeted freezing, melting, and
  significant deposition/sublimation process-hit smoke.
- `fw_event_type_smoke/run_singleproc_extended.conf`: same targeted
  single-process smoke with `TRACK_COLD_EVENT_EXTENDED_GEOMETRY=.true.` to
  validate optional single-process `ice_re/rp/rho` pre/post geometry.
- `fw_event_type_smoke/run_singleproc_significant.conf`: targeted significant
  single-process smoke with low thresholds for significant condensation,
  evaporation, freezing, and melting.
- `fw_event_type_smoke/run_collision_significant.conf`: targeted collision
  smoke with low thresholds for significant riming, aggregation, and
  liquid-liquid coalescence. The default ice seed reliably produces significant
  riming and aggregation records.
- `fw_event_type_smoke/run_cold_liq_liq_seed.conf`: dedicated cold coalescence
  seed. It keeps the cold collision path active while forcing a liquid/liquid
  participant pair. v1.3 validation uses generated occurrence/significant
  variants and checks `trigger_code=1` with `trigger_level=1/2`,
  `EVENT_COALESCENCE` in `sd_event_mask`, the matching `sd_event_sig_mask` bit
  for significant records, and collision-family `event_multiplicity`.
- `fw_event_type_smoke/run_spatial.conf`: v1.2+ spatial-visit smoke with a
  broad box and ordinary `SD_all_NetCDF_*` output. The current implementation
  uses Level-3 segment-box intersection plus the retained Level-2
  current-position check; this case validates `sd_spatial_visit_flag` schema and
  requires at least one visit hit.
- `fw_event_type_smoke/run_v12_spatial_level3.conf`: v1.2+ controlled
  Level-3 spatial crossing smoke. A smoke-only displacement crosses a thin
  spatial box by segment intersection, while endpoint-only point sampling would
  be insufficient for this test geometry.
- `fw_event_type_smoke/run_v13_spatial_rank_pruning_off_rank2.conf` and
  `run_v13_spatial_rank_pruning_on_rank2.conf`: v1.3+ rank-2 spatial pruning
  baseline and default-off pruning-on subset checks. The pruning-on case should
  preserve `sd_spatial_visit_flag` counts while one non-overlap rank is marked
  inactive.
- `fw_event_type_smoke/run_v13_spatial_rank_pruning_all_overlap_rank2.conf` and
  `run_v13_spatial_rank_pruning_margin_rank2.conf`: v1.3+ active-rank coverage
  checks for all-overlap and safety-margin behavior.
- `fw_event_type_smoke/run_vapor_occurrence_sublimation_seed.conf`: controlled
  v1.3+ sublimation occurrence seed. It lowers water vapor before the real
  `sdm_subldep` call so the process computes sublimation from its normal
  ice-mass update; it does not write event records directly. The expected
  process pair is `trigger_code=7`, `trigger_level=1`.
- `fw_event_type_smoke/run_diag.conf`: low-threshold diagnostic event smoke with
  extended collision geometry enabled.
- `fw_event_type_smoke/run_ordinary.conf`: ordinary `SD_all_NetCDF_*` cold
  tracking schema smoke.
- `fw_event_type_smoke/run_hist.conf`: history-style `SD_all_history.peXXXXXX`
  cold tracking schema smoke for `sdm_outnetcdf_hist`.
- `fw_event_type_smoke/run_hist_selected.conf`: selected history-style
  `SD_selected_history.peXXXXXX` cold tracking schema smoke. This is the path to
  use when the third `sdm_dmpvar` digit selects the sampled/selected SD subset.
- `fw_event_type_smoke/run_restart_*.conf`: continuous, split-restart, and legacy
  fallback restart checks for cold interval masks, flags, and maxima.
- `fw_event_type_smoke/run_terminal_precip_cleanup3.conf`: controlled
  driver-level terminal precipitation cleanup smoke. A test-only flag forces
  one live liquid/ice SD into the precipitation sentinel path before
  `sdm_sd2prec`, and `validate_invalid_cleanup.py` verifies that newly invalid
  slots are reset.
- `fw_event_type_smoke/run_vertical_outflow_cleanup3.conf`: controlled
  driver-level vertical outflow cleanup smoke. A test-only flag forces one live
  SD beyond the upper boundary before `sdm_jdginvdv`, and
  `validate_invalid_cleanup.py` verifies newly invalid slot reset.
- `fw_event_type_smoke/run_kohler_activation_occurrence.conf`: controlled
  v1.3+ derived Kohler activation crossing smoke. It enables
  `TRACK_COLD_OUTPUT_KOHLER_CONTEXT`, `tracking_evt_activation_enable`, and a
  test-only seed that starts one liquid SD below the selected-output critical
  radius before `sdm_condevp`; the expected event is `trigger_code=10`,
  `trigger_level=1`.
- `fw_event_type_smoke/run_kohler_deactivation_occurrence.conf`: controlled
  v1.3+ derived Kohler deactivation crossing smoke. It enables
  `TRACK_COLD_OUTPUT_KOHLER_CONTEXT`, `tracking_evt_deactivation_enable`, and a
  test-only seed that starts one liquid SD above the selected-output critical
  radius before `sdm_condevp`; the expected event is `trigger_code=11`,
  `trigger_level=1`.
- `fw_event_type_smoke/run_kohler_activation_significant.conf`: controlled
  v1.3+ significant Kohler activation crossing smoke. It enables both
  occurrence and significant activation controls with a zero absolute
  post-crossing margin threshold; the expected event is `trigger_code=10`,
  `trigger_level=2`, with the `EVENT_ACTIVATION` bit set in both
  `sd_event_mask` and `sd_event_sig_mask`.
- `fw_event_type_smoke/run_kohler_deactivation_significant.conf`: controlled
  v1.3+ significant Kohler deactivation crossing smoke. It enables both
  occurrence and significant deactivation controls with a zero absolute
  post-crossing negative-margin threshold; the expected event is
  `trigger_code=11`, `trigger_level=2`, with the `EVENT_DEACTIVATION` bit set
  in both `sd_event_mask` and `sd_event_sig_mask`.
- `fw_event_type_smoke/run_v13_aerosol_context.conf`: v1.3+ optional aerosol
  context smoke. It writes single-process `aerosol_total_mass_pre/post` and
  `aerosol_kohler_solute_pre/post` when
  `TRACK_COLD_OUTPUT_AEROSOL_CONTEXT=.true.`.
- `fw_event_type_smoke/run_v13_thermo_context.conf`: v1.3+ optional
  thermodynamic context smoke. It writes single-process `air_temperature`,
  `air_pressure`, and `water_vapor_mixing_ratio` when
  `TRACK_COLD_OUTPUT_THERMO_CONTEXT=.true.`.

The cases are intentionally short so they can be used after source edits as a
quick regression check.

## Local Smoke Workflow

From this directory:

```bash
cd fw_ice_smoke
../../../../../../bin/scale-rm_init init.conf
mkdir -p fw_output fw_tracking
../../../../../../bin/scale-rm run.conf
cd ..

python3 merge_tracking_interest_ids.py \
  --input-glob "./fw_ice_smoke/fw_tracking/tracking_interest_ids.pe*.ids" \
  --output "./fw_ice_smoke/fw_tracking/tracking_interest_ids_merged.ids" \
  --rank-bucket-basename "./fw_ice_smoke/fw_tracking/tracking_interest_ids_dedup"

cd bw_ice_smoke
../../../../../../bin/scale-rm_init init.conf
mkdir -p bw_output
../../../../../../bin/scale-rm run.conf
cd ..

python3 check_tpht_consistency.py \
  --fw-ids "./fw_ice_smoke/fw_tracking/tracking_interest_ids_merged.ids" \
  --bw-glob "./bw_ice_smoke/SD_selected_NetCDF_*.pe*"

python3 analyze_tpht_tracks.py \
  --bw-glob "./bw_ice_smoke/SD_selected_NetCDF_*.pe*" \
  --event-glob "./bw_ice_smoke/SD_event_collision_NetCDF_*.pe*" \
  --output-dir "./tpht_analysis"
```

## Targeted Ice Event-Type Workflow

Use these cases when a change must prove that cold coalescence writes riming
and aggregation trigger codes and that TPHT selects exactly the FW-discovered
particles.

```bash
cd fw_event_type_smoke
SCALE_ENABLE_SDM=T make -j4
./scale-rm_init init.conf
mkdir -p fw_event_output fw_tracking
./scale-rm run.conf
cd ..

python3 validate_event_types.py \
  --event-glob "./fw_event_type_smoke/SD_event_collision_NetCDF_*.pe*" \
  --require-trigger-pair 2:1 \
  --require-trigger-pair 3:1

python3 merge_tracking_interest_ids.py \
  --input-glob "./fw_event_type_smoke/fw_tracking/tracking_interest_ids.pe*.ids" \
  --output "./fw_event_type_smoke/fw_tracking/tracking_interest_ids_merged.ids" \
  --rank-bucket-basename "./fw_event_type_smoke/fw_tracking/tracking_interest_ids_dedup"

cd bw_event_type_smoke
../fw_event_type_smoke/scale-rm_init init.conf
mkdir -p bw_event_output
../fw_event_type_smoke/scale-rm run.conf
mkdir -p bw_tpht_output
../fw_event_type_smoke/scale-rm run_tpht.conf
cd ..

python3 validate_event_types.py \
  --event-glob "./bw_event_type_smoke/SD_event_collision_NetCDF_*.pe*" \
  --require-trigger-pair 2:1 \
  --require-trigger-pair 3:1

python3 check_tpht_consistency.py \
  --fw-ids "./fw_event_type_smoke/fw_tracking/tracking_interest_ids_merged.ids" \
  --bw-glob "./bw_event_type_smoke/SD_selected_NetCDF_*.pe*"
```

For targeted single-process events:

```bash
cd fw_event_type_smoke
rm -f SD_event_singleproc_NetCDF_*.pe*
./scale-rm run_singleproc.conf
cd ..

python3 validate_singleproc_events.py \
  --event-glob "./fw_event_type_smoke/SD_event_singleproc_NetCDF_*.pe*" \
  --require-trigger-pair 4:1 \
  --require-trigger-pair 5:1 \
  --forbid-extended-geometry

cd fw_event_type_smoke
rm -f SD_event_singleproc_NetCDF_*.pe*
./scale-rm run_singleproc_extended.conf
cd ..

python3 validate_singleproc_events.py \
  --event-glob "./fw_event_type_smoke/SD_event_singleproc_NetCDF_*.pe*" \
  --require-trigger-pair 4:1 \
  --require-trigger-pair 5:1 \
  --require-extended-geometry
```

For targeted v1.3 significant process-hit events:

```bash
cd fw_event_type_smoke
rm -f SD_event_singleproc_NetCDF_*.pe* SD_selected_NetCDF_*.pe*
./scale-rm run_singleproc_significant.conf
cd ..

python3 validate_singleproc_events.py \
  --event-glob "./fw_event_type_smoke/SD_event_singleproc_NetCDF_*.pe*" \
  --require-trigger-pair 8:2 \
  --require-trigger-pair 9:2 \
  --require-trigger-pair 4:2 \
  --require-trigger-pair 5:2 \
  --forbid-extended-geometry

python3 validate_event_mask_bits.py \
  --glob "./fw_event_type_smoke/SD_selected_NetCDF_*.pe*" \
  --mask-variable sd_event_sig_mask \
  --require-bit 128 \
  --require-bit 256 \
  --require-bit 8 \
  --require-bit 16

cd fw_event_type_smoke
rm -f SD_event_collision_NetCDF_*.pe* SD_selected_NetCDF_*.pe*
./scale-rm run_collision_significant.conf
cd ..

python3 validate_event_types.py \
  --event-glob "./fw_event_type_smoke/SD_event_collision_NetCDF_*.pe*" \
  --require-trigger-pair 2:2 \
  --require-trigger-pair 3:2

python3 validate_event_mask_bits.py \
  --glob "./fw_event_type_smoke/SD_selected_NetCDF_*.pe*" \
  --mask-variable sd_event_sig_mask \
  --require-bit 2 \
  --require-bit 4
```

For v1.3+ derived Kohler activation/deactivation occurrence events:

```bash
cd fw_event_type_smoke
./clean.sh
mkdir -p fw_tracking fw_output fw_event_output
./scale-rm_init init.conf
./scale-rm run_kohler_activation_occurrence.conf
cd ..

python3 validate_kohler_activation_events.py \
  --event-glob "./fw_event_type_smoke/SD_event_singleproc_NetCDF_*.pe*" \
  --ordinary-glob "./fw_event_type_smoke/SD_selected_NetCDF_*.pe*" \
  --require-activation \
  --require-kohler-context \
  --forbid-significant-activation \
  --forbid-lifecycle-activation \
  --check-event-mask \
  --check-no-phase-change-forcing \
  --require-condensation-overlap

cd fw_event_type_smoke
./clean.sh
mkdir -p fw_tracking fw_output fw_event_output
./scale-rm_init init.conf
./scale-rm run_kohler_deactivation_occurrence.conf
cd ..

python3 validate_kohler_activation_events.py \
  --event-glob "./fw_event_type_smoke/SD_event_singleproc_NetCDF_*.pe*" \
  --ordinary-glob "./fw_event_type_smoke/SD_selected_NetCDF_*.pe*" \
  --require-deactivation \
  --require-kohler-context \
  --forbid-significant-activation \
  --forbid-lifecycle-activation \
  --check-event-mask \
  --check-no-phase-change-forcing \
  --require-evaporation-overlap
```

For v1.2 spatial-visit output:

```bash
cd fw_event_type_smoke
./scale-rm run_spatial.conf
cd ..

python3 validate_cold_tracking_schema.py \
  --glob "./fw_event_type_smoke/SD_all_NetCDF_*" \
  --label spatial_cold_all

python3 validate_spatial_visit_flag.py \
  --glob "./fw_event_type_smoke/SD_all_NetCDF_*" \
  --expect-any
```

For the controlled Level-3 thin-box crossing smoke:

```bash
cd fw_event_type_smoke
./scale-rm run_v12_spatial_level3.conf
cd ..

python3 validate_spatial_visit_flag.py \
  --glob "./fw_event_type_smoke/SD_all_NetCDF_*" \
  --expect-any
```

For the dedicated cold coalescence coverage seed:

```bash
python3 squid_validation/make_v12_validation_configs.py

cd fw_event_type_smoke
rm -f SD_event_collision_NetCDF_*.pe* SD_selected_NetCDF_*.pe*
./scale-rm run_v13_coalescence_seed_occurrence.conf
cd ..

python3 validate_cold_liq_liq_seed.py \
  --event-glob "./fw_event_type_smoke/SD_event_collision_NetCDF_*" \
  --ordinary-glob "./fw_event_type_smoke/SD_selected_NetCDF_*" \
  --require-level 1

cd fw_event_type_smoke
rm -f SD_event_collision_NetCDF_*.pe* SD_selected_NetCDF_*.pe*
./scale-rm run_v13_coalescence_seed_significant.conf
cd ..

python3 validate_cold_liq_liq_seed.py \
  --event-glob "./fw_event_type_smoke/SD_event_collision_NetCDF_*" \
  --ordinary-glob "./fw_event_type_smoke/SD_selected_NetCDF_*" \
  --require-level 2
```

For driver-level invalid-slot cleanup smokes:

```bash
cd fw_event_type_smoke
./clean.sh
mkdir -p fw_tracking fw_output fw_event_output restart_mid_output restart_resume_output restart_cont_output restart_legacy_output
mpiexec --oversubscribe --map-by slot -n 1 ./scale-rm_init init.conf
mpiexec --oversubscribe --map-by slot -n 1 ./scale-rm run_terminal_precip_cleanup3.conf \
  > /private/tmp/cold_tpht_cleanup_smoke/terminal_cleanup3.log 2>&1
python3 ../validate_invalid_cleanup.py \
  --log /private/tmp/cold_tpht_cleanup_smoke/terminal_cleanup3.log \
  --require-label terminal_precip

./clean.sh
mkdir -p fw_tracking fw_output fw_event_output restart_mid_output restart_resume_output restart_cont_output restart_legacy_output
mpiexec --oversubscribe --map-by slot -n 1 ./scale-rm_init init.conf
mpiexec --oversubscribe --map-by slot -n 1 ./scale-rm run_vertical_outflow_cleanup3.conf \
  > /private/tmp/cold_tpht_cleanup_smoke/vertical_cleanup3.log 2>&1
python3 ../validate_invalid_cleanup.py \
  --log /private/tmp/cold_tpht_cleanup_smoke/vertical_cleanup3.log \
  --require-label vertical_outflow
```

The cleanup smoke flags are disabled by default and are not production
controls:

- `tracking_cleanup_debug_enable` prints compact `COLD_TPHT_CLEANUP` counters.
- `tracking_cleanup_terminal_smoke_enable` forces one SD through the terminal
  precipitation conversion path for validation.
- `tracking_cleanup_vertical_smoke_enable` forces one SD through the upper
  vertical outflow invalidation path for validation.

`tracking_event_type_smoke_enable` is disabled by default in the model. It is
only enabled in the targeted event smoke configs and now validates cold
`trigger_code` records; it should not be used for production or science runs.

Cold TPHT v1.2 changes occurrence-level defaults: all `tracking_evt_*_enable`
flags default to `.false.`. Smoke cases that need the v1.1 occurrence records
explicitly enable the relevant flags in `&PARAM_ATMOS_PHY_MP_SDM_TRACKING`.
Continuous vapor-growth occurrence records for deposition, sublimation,
condensation, and evaporation remain available but default disabled.

On macOS rank-1 local runs, OpenMPI may print a `bind() failed: Operation not
permitted (1)` warning from the TCP transport. Treat it as an environment
warning when `scale-rm_init` and `scale-rm` both exit with status 0.

## Expected Checks

- FW writes `SD_selected_NetCDF_*` with `sd_id`, `dm_id`, `sd_liqice`, and
  cold-SDM `sdi_*` fields.
- FW writes raw TPHT interest IDs under `fw_ice_smoke/fw_tracking/`.
- BW writes `SD_selected_NetCDF_*` with `pre_sdid`, `pre_dmid`,
  `sd_event_mask`, `sd_event_sig_mask`, `sd_liqice`, and cold-SDM `sdi_*`
  fields.
- Cold collision event files use `SD_event_collision_NetCDF_*` and include
  `trigger_code`, `trigger_level`, `target_reason_mask`, participant IDs,
  `event_multiplicity`, `phase_state1_pre/post`, `phase_state2_pre/post`,
  pair-midpoint `x/y/z`, participant coordinates
  `x1/y1/z1` and `x2/y2/z2`, `sd_n` pre/post, and phase-aware
  `hydro_radius`/`hydro_mass` pre/post. They do not write legacy `event_type`,
  `num_col`, `sd_liqice1`, or `sd_liqice2`.
- Cold event-file phase states use one `sd_liqice`-compatible code set:
  `0` dry aerosol / aerosol-only, `1` liquid, `10` ice, `11` mixed, and `99`
  none / invalid / missing / not applicable. Validators reject obsolete
  canonical ice/mixed codes `2` and `3`, reject `transition_code*`, and do not
  expect `raw_model_phase_state_*` or `canonical_phase_state_*` variables.
- Dynamic tracking IDs are supported by validators and sidecars. `sd_id=-999`
  and `dm_id=-999` are invalid / missing; non-negative `sd_id` values and
  `sd_id<=-1000` are valid. `validate_dynamic_id_predicates.py` checks that
  active TPHT paths no longer use generic `sd_id < 0`,
  `sd_id <= INVALID_i4`, or `sd_id > INVALID_i4` as ID predicates.
- `SD_lifecycle_NetCDF_*` is available as a lightweight lifecycle sidecar
  skeleton for `LIFE_SDADD_SPLIT=1` and
  `LIFE_SDREMOVE_INVALIDATION=2`. The current implementation is compile-safe
  tracking-aware support only: `sdm_sdremove` / `sdm_sdadd` carry tracking
  arrays and lifecycle sidecar calls, but the active smoke path does not
  exercise them because the production driver does not call `sdm_adjsdnum`.
  Runtime lifecycle tests should use a dedicated test-only harness rather than
  activating the production SD-number adjustment path.
- The driver-level terminal precipitation and vertical outflow cleanup smokes
  exercise `LIFE_SDREMOVE_INVALIDATION=2` records for newly invalid slots. They
  validate the approved driver cleanup paths, not inactive `sdm_sdremove` /
  `sdm_sdadd` runtime behavior.
- `validate_lifecycle_harness.py` is a test-only lifecycle harness. It validates
  reset, remove-record, split-record, parentless-entry, dynamic child-ID, and
  fixed-width lifecycle schema semantics without activating production
  `sdm_adjsdnum` or `sdm_aslform`.
- `validate_lifecycle_schema_static.py` checks the fixed-width
  `SD_lifecycle_NetCDF_*` writer schema, lifecycle code constants, dynamic-ID
  restart counter handling, and that `sdm_adjsdnum` remains inactive.
- `validate_kohler_activation_helpers.py` checks the shared Kohler
  critical-radius helper, selected-output `sdtype='activated'` refactor,
  condensation/evaporation wrapper pre/post activation-state capture, optional
  Kohler context schema strings, and the rule that activation/deactivation are
  not lifecycle records.
- `SD_event_singleproc_NetCDF_*` writes freezing/melting occurrence records and
  enabled significant deposition, sublimation, condensation, evaporation,
  freezing, and melting process-hit records with `phase_state_pre`,
  `phase_state_post`, `hydro_radius_pre/post`, `liq_mass_pre/post`, and
  `ice_mass_pre/post`, but without `event_multiplicity`. When
  `TRACK_COLD_OUTPUT_ICE_GEOMETRY=.true.` or the backward-compatible
  `TRACK_COLD_EVENT_EXTENDED_GEOMETRY=.true.`, it also writes optional
  `ice_re_pre/post`, `ice_rp_pre/post`, and `ice_rho_pre/post`, corresponding
  to `sdi%re` ice-crystal equatorial radius, `sdi%rp` ice-crystal polar radius,
  and `sdi%rho` ice-crystal density.
- v1.3 process event files use `trigger_code` for process identity and
  `trigger_level` for occurrence (`1`) versus significant (`2`) records.
  Significant deposition/sublimation use relative ice-mass change; significant
  condensation/evaporation use relative liquid-mass change; significant
  freezing/melting use relative component mass transfer; significant
  collision-family records use pre-event participant mass ratios and keep
  `event_multiplicity=num_col`.
- v1.3 occurrence triggers are disabled by default and require
  process-specific `tracking_evt_*_enable` namelist entries. Collision-family
  occurrence records keep `event_multiplicity=num_col`; single-process
  occurrence records do not write `event_multiplicity`. If a significant gate
  hits for the same process update, only the `trigger_level=2` record is
  written.
- v1.3+ activation/deactivation are derived process events based on event-time
  Kohler critical-radius crossing around `sdm_condevp`. They write
  `SD_event_singleproc_NetCDF_*` occurrence records with `trigger_code=10/11`
  and `trigger_level=1` when `tracking_evt_activation_enable` or
  `tracking_evt_deactivation_enable` is enabled. They write default-off
  significant records with `trigger_level=2` when the corresponding
  `tracking_sig_*_enable` flag is enabled and the absolute post-crossing
  Kohler margin threshold is reached.
- v1.3+ SQUID full validation passed as `832697.sqd` with
  `required_failures=0`; that completed job recorded vapor-growth occurrence
  `7:1`, `8:1`, and `9:1` as optional-not-covered at the time. A later local
  follow-up added controlled seeds and `validate_vapor_occurrence_events.py`
  coverage for all three pairs. v1.3+ SQUID production I/O / threshold
  calibration passed as `832700.sqd` with
  `required_failures=0`. The low-diagnostic stress case produced
  `75,828,480` diagnostic records and should be interpreted as an I/O stress
  result, not a runtime failure. Longer production calibration scripts are
  prepared as `run_v13_long_production_calibration.sh` and
  `squid_run_v13_long_production_calibration.sh`, but they are not validation
  evidence until a new run completes with `required_failures=0`.
- v1.2 optional output groups are controlled by
  `TRACK_COLD_OUTPUT_ICE_GEOMETRY` and
  `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY`. The legacy
  `TRACK_COLD_EVENT_EXTENDED_GEOMETRY` switch remains a broad alias for both
  groups. `TRACK_COLD_OUTPUT_KOHLER_CONTEXT` enables optional single-process
  `kohler_rcrit_pre/post`, `kohler_margin_pre/post`, and
  `activated_state_pre/post` fields. `TRACK_COLD_OUTPUT_AEROSOL_CONTEXT`
  writes optional single-process `aerosol_total_mass_pre/post` and
  `aerosol_kohler_solute_pre/post`; `TRACK_COLD_OUTPUT_THERMO_CONTEXT` writes
  optional single-process `air_temperature`, `air_pressure`, and
  `water_vapor_mixing_ratio`. Collision aerosol/thermo context remains future
  work.
- `SD_event_diag_NetCDF_*` writes one record per SD per diagnostic bit at the SD
  output time, including `x/y/z`, masks/flags, and the seven interval maxima.
  `sd_diag_mask` is evaluated from the maintained interval maxima; `x/y/z` are
  output-time coordinates, not the maximum-occurrence coordinates.
- Cold ordinary `SD_all_NetCDF_*` and selected output write `sd_event_mask`,
  `sd_diag_mask`, `sd_phase_change_flag`, `sd_spatial_visit_flag`, and seven
  interval maxima; cold ordinary output does not write legacy `if_coal` by
  default.
- Cold history `SD_all_history.peXXXXXX` output writes the same cold tracking
  fields with per-output-time suffixes and does not write legacy `if_coal_*` by
  default.
- Cold selected-history `SD_selected_history.peXXXXXX` output writes the same
  fields for the selected/sampled subset when selected SD output is requested.
- New-format cold SD restarts preserve interval masks, flags, spatial-visit
  state, and maxima across a mid-interval restart. Legacy cold SD restarts
  without the optional tracking block still read successfully and initialize the
  cold interval state.
- `analyze_tpht_tracks.py` writes cold trigger-code counts into
  `tpht_analysis/tpht_analysis_summary.json` and
  `tpht_analysis/tpht_collision_time_summary.csv`.
