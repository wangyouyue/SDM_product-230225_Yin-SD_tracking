# Introduction to the Super-Droplet Method (SDM) in SCALE-SDM

## What is SDM?
The **Super-Droplet Method (SDM)**, originally introduced by **Shima et al. (2009)**, represents a significant advancement in Lagrangian cloud microphysics simulations. Unlike traditional Eulerian microphysics schemes, SDM employs **super-droplets (SDs)**, which are computational particles representing a large number of real droplets, aerosols, or precipitation particles with similar properties. This framework enables explicit treatment of key cloud-microphysical processes, including condensation, evaporation, and coalescence, thereby improving the physical interpretability of cloud formation and precipitation evolution.

## Integration with SCALE-SDM
The **SCALE (Scalable Computing for Advanced Library and Environment)** framework, developed with co-design by computational and computer science researchers, provides a robust, scalable platform for high-resolution atmospheric simulations (Nishizawa et al., 2015; Sato et al., 2015). SCALE-SDM integrates the SDM within this framework, offering a powerful tool for simulating cloud microphysics and dynamics. For more information on the SCALE platform, please visit the [SCALE's official website](https://scale.riken.jp/) (last access: 6 April 2026).

This repository integrates SDM into SCALE version 5.2.6, leveraging both SDM’s microphysical precision and SCALE’s computational scalability to facilitate simulations of cloud systems that capture the complex interactions between microphysics and atmospheric dynamics. Further details about this version of SCALE can be explored in the [SCALE version 5.2.6 archives](https://scale.riken.jp/archives/5.2.6/) (last access: 6 April 2026).

# Super-Droplet Tracking (Merged Forward/Backward Version)

## Contents
- [1. Overview of New Features (vs. SCALE-SDM)](#1-overview-of-new-features-vs-scale-sdm)
- [2. Core Tracking State Variables (`sd_id/dm_id`, `pre_sdid/pre_dmid`, `if_coal`)](#2-core-tracking-state-variables-sd_iddm_id-pre_sdidpre_dmid-if_coal)
- [3. Sampling Algorithms](#3-sampling-algorithms)
- [4. Forward Tracking Algorithm](#4-forward-tracking-algorithm)
- [5. Backward Tracking Algorithm](#5-backward-tracking-algorithm)
- [6. Collision-Coalescence Event Recording](#6-collision-coalescence-event-recording)
- [7. SD Output, `sdm_dmpvar`, and File Naming Rules](#7-sd-output-sdm_dmpvar-and-file-naming-rules)
- [8. Namelist Parameters Added/Extended](#8-namelist-parameters-addedextended)
- [9. Random Perturbations in SD Motion](#9-random-perturbations-in-sd-motion)
- [10. Installation and Usage](#10-installation-and-usage)
- [11. Representativeness Tests](#11-representativeness-tests)
- [12. Cold TPHT Current Base and Future Porting Plan](#12-cold-tpht-current-base-and-future-porting-plan)

## 1. Overview of New Features (vs. SCALE-SDM)
Compared with the original SCALE-SDM branch, this merged branch introduces a unified, sampling-aware tracking framework with the following characteristics:
- **Forward tracking (FW)** and **backward tracking (BW)** are supported within the same executable and are controlled by a unified switch (`tracking_mode`) while preserving legacy compatibility.
- A common sampling subsystem is shared by FW and BW, with two selectable algorithms: `random` and `stratified`.
- A new **Two-Pass Hybrid Tracking (TPHT)** workflow is supported: FW can first scan the whole SD population, accumulate only the SD IDs that satisfy user-defined interest conditions, and BW can then reinitialize from this cumulative ID set instead of from a fresh random/stratified sample.
- Direction-dependent tracking identifiers are incorporated into SD NetCDF outputs: FW writes `sd_id/dm_id`, BW writes `pre_sdid/pre_dmid`, warm ordinary output can write `if_coal`, and cold ordinary output writes `sd_event_mask`/`sd_diag_mask`/`sd_phase_change_flag` when cold tracking arrays are present. Together these fields support trajectory-chain reconstruction in both directions during post-processing.
- Collision-event outputs are appended to time-bucketed NetCDF files: warm mode keeps `SD_coal_output_NetCDF_*`, while cold mode writes `SD_event_collision_NetCDF_*`.
- Random perturbation of SD motion is exposed through the namelist (`random_perturbation_enable`, `random_perturbation_amp`) and is physically inactive when the amplitude is zero.
- Representativeness test suites are provided for both FW and BW to quantify sampling-induced bias and chain-level fidelity.

### 1.1 Two-Pass Hybrid Tracking (TPHT)
This update adds a new workflow intended for “find first, trace later” studies:
- **Pass 1 (FW discovery pass):** run forward tracking on the full valid SD population and append rank-local raw ID records for SDs that satisfy one or more user-defined interest conditions.
- **Pass 2 (BW reconstruction pass):** offline-deduplicate the FW ID records, repartition the deduplicated target set by `dm_id`, and then run backward tracking from the per-rank FW ID buckets so that BW follows exactly the SD set discovered in the FW scan.

In this README, this workflow is referred to as **Two-Pass Hybrid Tracking**, abbreviated as **TPHT**. The word “hybrid” emphasizes that the method combines:
- a forward discovery pass,
- an interest-based ID reduction step,
- and a backward reconstruction pass.

The currently implemented interest conditions are:
- liquid-radius threshold: `sd_r >= tracking_interest_radius_threshold`,
- recorded collision participation: warm `if_coal > 0`, or cold
  `sd_event_mask` with a collision-family bit set,
- ice radius threshold, using the phase-aware volume-equivalent `hydro_radius`
  helper when cold-SDM fields are available,
- ice or mixed-phase status via `tracking_interest_ice_phase_enable`,
- rimed-mass threshold, using `sdi%mrime >= tracking_interest_rime_mass_threshold` when cold-SDM fields are available.

These FW `.ids` interest filters are separate from the cold event-file
`target_reason_mask`. Cold TPHT v1 does not define `TARGET_BY_STATE` or
`TRIG_STATE_*`; phase state is used as an internal gate, a radius-definition
selector, and diagnostic context, not as a separate event-file target-reason
category.

When `tracking_id_output_basename` is enabled together with at least one interest condition, FW no longer relies on `SD_selected_NetCDF_*` as the carrier of the target set. Instead, it appends raw `(dm_id, sd_id)` pairs to rank-local `.ids` files for later offline merge/dedup. Each append block also writes a `# TPHT_META PRC_NUM_X PRC_NUM_Y PRC_nprocs` header so that BW can validate MPI-decomposition compatibility before accepting the handoff file.

### 1.2 Cold SD event tracking
When `sdm_cold = .false.`, ordinary warm SD NetCDF files can still write
`if_coal`, and warm collision-event files keep the `SD_coal_output_NetCDF_*`
name. Collision multiplicity uses the output variable name
`event_multiplicity`; the internal coalescence code still uses the local
variable name `num_col`. This is an intentional warm event-output schema rename:
`event_multiplicity` in `SD_coal_output_NetCDF_*` is equivalent to the internal
`num_col` value, and the event file no longer writes a `num_col` variable.

When `sdm_cold = .true.`, cold TPHT uses cold-specific event and diagnostic
variables:
- Ordinary SD NetCDF output writes `sd_event_mask`, `sd_event_sig_mask`,
  `sd_diag_mask`, `sd_phase_change_flag`, `sd_spatial_visit_flag`, and interval
  maxima for liquid radius, ice volume-equivalent radius, mixed-phase
  volume-equivalent radius, rimed mass, rime fraction, monomer count, and aspect
  ratio. It does not write legacy `if_coal` by default.
- v1.3 process-event files use `trigger_code` for process identity and
  `trigger_level` for occurrence/significance. `trigger_level=1` means an
  occurrence-level record; `trigger_level=2` means a significant-level record.
  For one physical process update, significant overrides occurrence: if the
  significant gate is enabled and hit, only one significant-level record is
  written; otherwise an occurrence-level record is written only when the
  matching occurrence enable flag is true.
- `sd_event_mask` is an interval bitmask for processes that occurred at any
  tracked level during the last SD output interval. `sd_event_sig_mask` uses the
  same bit table, but records only processes that reached significant level.
  Bits are `1` coalescence, `2` riming, `4` aggregation, `8` freezing, `16`
  melting, `32` deposition, `64` sublimation, `128` condensation, `256`
  evaporation, `512` activation, and `1024` deactivation. Activation and
  deactivation are derived event-time process events based on Kohler
  critical-radius crossing around the condensation/evaporation update. Multiple
  processes are represented by the sum of their bits; this is not an event
  count.
- `SD_event_collision_NetCDF_*` replaces cold `SD_coal_output_NetCDF_*`.
  It writes process `trigger_code` values `1` coalescence, `2` riming, and `3`
  aggregation. The canonical user-visible name is `coalescence`; in the current
  cold collision path it refers to liquid-liquid coalescence.
  `event_multiplicity` is the output name for the internal `num_col` value and
  applies only to collision-family records. Core variables include
  `phase_state1_pre`, `phase_state1_post`, `phase_state2_pre`,
  `phase_state2_post`, pair-midpoint `x/y/z`, participant coordinates
  `x1/y1/z1` and `x2/y2/z2`, participant `sd_n*_pre/post`, phase-aware
  `hydro_radius*_pre/post`, and `hydro_mass*_pre/post`. Participant identifiers
  are `sd_id*`/`dm_id*` in FW mode and `pre_sdid*`/`pre_dmid*` in BW mode.
  The collision phase-state fields are context variables for interpreting
  `hydro_radius*_pre/post` and debugging collision outcomes. They are not, by
  themselves, physical phase-change target records; apparent riming-related
  phase-state changes can reflect multiplicity/slot bookkeeping and do not
  automatically set `TARGET_BY_TRANSITION` or create a transition-only record.
- `SD_event_singleproc_NetCDF_*` stores single-SD process records. Process
  `trigger_code` values are `4` freezing, `5` melting, `6` deposition, `7`
  sublimation, `8` condensation, `9` evaporation, `10` activation, and `11`
  deactivation. Activation/deactivation are derived from pre/post Kohler
  activated-state crossing and are written only as occurrence-level records in
  the current implementation. Significant deposition/sublimation
  are relative ice-mass-change process hits; significant
  condensation/evaporation are relative liquid-mass-change process hits;
  significant freezing/melting are relative component-mass-transfer hits.
  Significant activation/deactivation are default-off Kohler margin hits and
  write `trigger_level=2` records when enabled and thresholded. All
  occurrence and significant modes are disabled by default and require
  process-specific enable flags; significant modes also require thresholds.
  Core variables include `phase_state_pre`, `phase_state_post`, output-time
  `x/y/z`, `hydro_radius_pre`, `hydro_radius_post`, `liq_mass_pre/post`, and
  `ice_mass_pre/post`; single-process files do not write `event_multiplicity`.
  When `TRACK_COLD_OUTPUT_ICE_GEOMETRY = .true.`, single-process files also
  write optional ice geometry fields `ice_re_pre/post`, `ice_rp_pre/post`, and
  `ice_rho_pre/post`. These expose `sdi%re` (ice-crystal equatorial radius),
  `sdi%rp` (ice-crystal polar radius), and `sdi%rho` (ice-crystal density)
  before and after the process update. They are optional morphology context;
  `hydro_radius_pre/post` remains the mandatory compact phase-aware size
  diagnostic.
  When `TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.`, single-process files also
  write optional derived Kohler context fields `kohler_rcrit_pre/post`,
  `kohler_margin_pre/post`, and `activated_state_pre/post`. These fields are
  optional context for interpreting activation/deactivation crossing and are not
  lifecycle variables.
  Local v1.3+ controlled smokes cover activation occurrence
  (`trigger_code=10`, `trigger_level=1`) and deactivation occurrence
  (`trigger_code=11`, `trigger_level=1`) with `TRACK_COLD_OUTPUT_KOHLER_CONTEXT`
  enabled. Follow-up controlled smokes cover significant activation
  (`trigger_code=10`, `trigger_level=2`) and significant deactivation
  (`trigger_code=11`, `trigger_level=2`). Validators also confirm these records
  do not write lifecycle records and do not force `sd_phase_change_flag`.
  Freezing/melting process triggers do not necessarily imply
  `phase_state_pre != phase_state_post`; for example, a future mixed SD can
  freeze internally and remain in the `PHASE_MIXED` category.
- `SD_event_diag_NetCDF_*` is an output-interval-end diagnostic target file: each
  record corresponds to one SD and one diagnostic bit evaluated at the SD output
  time from maintained interval maximum diagnostics. `sd_diag_mask` records
  whether each diagnostic threshold was reached by the maintained interval
  maximum during the previous SD output interval. It does not identify the exact
  microphysical timestep when the maximum occurred. The default file includes
  `x/y/z`, the three interval masks/flags, and the interval maxima; `x/y/z` are
  the SD position at output time, not the position where an interval maximum
  occurred.
- `target_reason_mask` is a bitmask describing why a record was selected; in v1
  it uses event, diagnostic, and phase-state-category-change reason bits.
  `trigger_code` remains categorical and should not be decoded as a bitmask.
- `sd_spatial_visit_flag` is a v1.2+ ordinary/selected/history interval flag.
  When `tracking_spatial_visit_enable = .true.`, the model checks one
  configured spatial box. The current implementation includes Level-3
  segment-box intersection after motion/advection: the previous and current SD
  positions form a straight segment, and an intersection sets the flag to `1`
  for the current SD output interval. The Level-2 current-position-inside-box
  test remains as a fallback and output-time check. Level 3 is stricter than
  point sampling, but it still follows the model-time discretized straight
  segment and does not reconstruct a continuous curved trajectory within one
  timestep. v1.3+ adds default-off rank-subdomain pruning through
  `tracking_spatial_rank_pruning_enable`; when disabled, spatial tracking keeps
  the conservative all-rank-active behavior.
  timestep.
- Cold TPHT tracking IDs use a lifecycle-safe namespace. `sd_id = -999` and
  `dm_id = -999` are the invalid / missing sentinels. Non-negative `sd_id`
  values remain valid initial/static IDs, and `sd_id <= -1000` is reserved for
  valid dynamic IDs assigned by lifecycle-enabled paths such as future
  reserve/global-halo/seeding entry or `sdm_sdadd` split support. `dm_id`
  remains a non-negative domain ID or `-999`.
- FW initializes IDs as before and preserves them across MPI migration. BW
  treats `sd_id/dm_id` as previous-output-time identity and must allow targets
  that first appear after model initialization. TPHT histories therefore may
  begin at the target's first valid output time unless a lifecycle parent link
  allows continued backtracking.
- `SD_lifecycle_NetCDF_*` is a lightweight lifecycle sidecar skeleton for
  parent/child identity records. The current schema writes `time`,
  `lifecycle_code`, `n_parent`, `n_child`, two parent ID slots, four child ID
  slots, and `x/y/z`. `LIFE_SDADD_SPLIT=1` and
  `LIFE_SDREMOVE_INVALIDATION=2` are active if those currently inactive
  routines are called; domain/global-halo/seeding/aslform entry codes are
  reserved inactive.
- Cold TPHT now uses one phase-state code in event files. The only output
  variables are `phase_state_pre/post` and `phase_state1/2_pre/post`; there is
  no separate `raw_model_phase_state_*` / `canonical_phase_state_*` schema.
  The code follows the model `sd_liqice` convention where possible: `0`
  `PHASE_DRY_AEROSOL` for dry aerosol / aerosol-only, `1` `PHASE_LIQUID`,
  `10` `PHASE_ICE`, `11` `PHASE_MIXED`, and `99` `PHASE_NONE` for
  none / invalid / missing / not-applicable context. Do not read `0` as none.
- `hydro_radius` means a phase-aware volume-equivalent radius. Liquid SDs use
  the liquid radius, ice SDs use an ice volume-equivalent radius, and mixed SDs
  use a combined liquid-plus-ice volume-equivalent radius.
- `coalescence_output_enable` is retained as the legacy master switch. In cold
  mode it controls the new `SD_event_*` files; in warm mode it still controls
  `SD_coal_output_NetCDF_*`.
- Cold tracking interval state is restart-exact for new-format cold SD restarts:
  when `sdm_cold = .true.`, restart output appends a
  cold tracking block containing `sd_event_mask`, `sd_event_sig_mask`,
  `sd_diag_mask`, `sd_phase_change_flag`, `sd_spatial_visit_flag`, and the
  seven interval maxima. The current v4 block also saves and restores
  `tracking_next_dynamic_sd_id`.
- Older cold SD restarts without this optional block are still accepted as legacy
  restarts; their cold tracking interval state is initialized to zero and
  `tracking_next_dynamic_sd_id` is initialized to `TRACK_SD_ID_DYNAMIC_START`
  (`-1000`). This counter fallback is safe for native older restarts because
  they cannot contain dynamic SD IDs; it would need a compatibility shim only if
  dynamic-ID production had been externally backported into an older restart
  format. A present but incompatible or incomplete cold tracking restart block
  is treated as an error instead of being partially recovered.
- `sd_phase_change_flag` is an ordinary/selected/history output interval flag:
  `0` means no phase-state category change during the previous SD output
  interval, and `1` means at least one phase-state category change occurred.
  The flag should not be read as a freezing/melting count or as proof that every
  freezing/melting event changed the SD phase category.
- `sdm_aslform` remains future lifecycle-tracking work. `sdm_sdremove` /
  `sdm_sdadd` have compile-safe tracking-aware support: their signatures carry
  tracking arrays, dynamic child ID assignment, interval-state copy/reset, and
  lifecycle sidecar calls. Runtime lifecycle validation is still blocked because
  the active driver does not call `sdm_adjsdnum`; runtime tests should use a
  test-only lifecycle harness rather than activating the production SD-number
  adjustment path.
- The test-only lifecycle harness validates helper-level reset, remove, split,
  parentless-entry, dynamic-ID, and fixed-width sidecar schema semantics. It
  does not imply production `sdm_sdremove` / `sdm_sdadd` runtime validation.
- Driver-level terminal precipitation and vertical outflow cleanup paths have
  controlled local runtime-smoke coverage. These smokes use test-only
  `tracking_cleanup_*` namelist flags that default to `.false.` and are not
  production controls.
- Occurrence-level and significant-level controls are independent, but the
  event stream is mutually exclusive for a single process update. When a
  significant mode is enabled and its threshold is exceeded, the significant
  record uses the same process `trigger_code` with `trigger_level=2` and
  overrides the occurrence-level record for that update. If the significant
  gate is not hit and the occurrence control is enabled, one occurrence record
  is written with `trigger_level=1`.
- New cold tracking comments follow the local Fortran comment style of
  neighboring SDM files. They document mask semantics, interval
  update/reset/restart behavior, categorical-vs-bitmask code usage,
  collision-only `event_multiplicity`, phase-state context interpretation,
  phase-aware `hydro_radius`, thresholded significant dep/sub process hits, and
  diagnostic output-time coordinates without adding long prose blocks.

### 1.3 Cold tracking variable and code types

| Variable or setting | Current type / NetCDF convention | Meaning |
|---|---|---|
| `sd_event_mask` | Fortran default `integer`; NetCDF `NF90_INT` | Ordinary-output interval event bitmask. |
| `sd_diag_mask` | Fortran default `integer`; NetCDF `NF90_INT` | Ordinary-output interval diagnostic bitmask evaluated from interval maxima. |
| `sd_phase_change_flag` | Fortran default `integer`; NetCDF `NF90_INT` | 0/1 interval flag for phase-state category change. |
| `sd_spatial_visit_flag` | Fortran default `integer`; NetCDF `NF90_INT` | v1.2+ 0/1 interval flag for one configured spatial box; current code uses Level-3 segment-box intersection plus Level-2 point checks. |
| `trigger_code` | Fortran default `integer`; NetCDF `NF90_INT` | Categorical record trigger code, not a bitmask. |
| `target_reason_mask` | Fortran default `integer`; NetCDF `NF90_INT` | Broad reason bitmask: event, diagnostic, phase-state category change. No v1 state-target bit. |
| `phase_state_pre/post` / `phase_state1/2_pre/post` | Fortran default `integer`; NetCDF `NF90_INT` | Single event-context phase-state code used to interpret phase-aware radius fields: `0` `PHASE_DRY_AEROSOL`, `1` `PHASE_LIQUID`, `10` `PHASE_ICE`, `11` `PHASE_MIXED`, `99` `PHASE_NONE` / invalid / missing. |
| `event_multiplicity` | Fortran default `integer`; NetCDF `NF90_INT` | Collision-family multiplicity; warm/cold event-output name for internal `num_col`. |
| `sd_n`, `sd_n*_pre/post` | Existing `integer(DP)` source type; event files use `NF90_INT64` | Super-droplet multiplicity. |
| radii, masses, interval maxima, `x/y/z` | `real(RP)`; NetCDF follows `RP` (`NF90_DOUBLE` when `RP == DP`) | Physical coordinates and diagnostics. |
| `TRACK_COLD_EVENT_EXTENDED_GEOMETRY` | namelist `logical` | Backward-compatible broad alias; when true, enables `TRACK_COLD_OUTPUT_ICE_GEOMETRY` and `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY`. |
| `TRACK_COLD_OUTPUT_ICE_GEOMETRY` | namelist `logical` | Enables optional event-file ice geometry fields. |
| `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY` | namelist `logical` | Enables optional collision rime and morphology fields. |
| `TRACK_COLD_OUTPUT_AEROSOL_CONTEXT` | namelist `logical` | Enables optional single-process aerosol context fields: `aerosol_total_mass_pre/post` and `aerosol_kohler_solute_pre/post`. |
| `TRACK_COLD_OUTPUT_THERMO_CONTEXT` | namelist `logical` | Enables optional single-process thermodynamic context fields: `air_temperature`, `air_pressure`, and `water_vapor_mixing_ratio`. |
| `TRACK_COLD_OUTPUT_KOHLER_CONTEXT` | namelist `logical` | Enables optional single-process derived Kohler context fields for activation/deactivation interpretation. |
| `tracking_cleanup_debug_enable` | namelist `logical` | Test-only cleanup smoke diagnostic; prints compact `COLD_TPHT_CLEANUP` counters when enabled. |
| `tracking_cleanup_terminal_smoke_enable` | namelist `logical` | Test-only terminal precipitation cleanup injector; default `.false.`, not for production runs. |
| `tracking_cleanup_vertical_smoke_enable` | namelist `logical` | Test-only vertical outflow cleanup injector; default `.false.`, not for production runs. |

Masks and categorical codes do not need 64-bit integers in v1. Only `sd_n` and,
depending on future source ranges, collision multiplicity may need wider integer
types; the current event writers keep `event_multiplicity` as `NF90_INT`.

### 1.4 Cold TPHT design boundaries and future plan

The detailed near-term and long-term design notes are maintained in
[`docs/cold_tpht_v1_future_development_plan.md`](docs/cold_tpht_v1_future_development_plan.md).
That document separates current code facts from v1.1 candidates and long-term
future work. The main boundaries are:

- v1.3 distinguishes occurrence-level records from thresholded significant
  process-hit records with `trigger_level`. `trigger_code` now identifies the
  process, `sd_event_mask` records any tracked level, and
  `sd_event_sig_mask` records the significant subset. All occurrence-level
  records are controlled by process-specific `tracking_evt_*_enable` flags and
  default to disabled. All significant modes remain controlled by
  process-specific enable flags and thresholds and default to disabled.
- `tracking_interest_ice_phase_enable` is the current coarse FW interest
  selector for `STAT_ICE` or `STAT_MIX` SDs. It is not a generic phase gate and
  does not introduce `TARGET_BY_STATE` or `TRIG_STATE_*`.
- Current valid tracking IDs are local positive slot indices with `dm_id=mype`;
  `INVALID_i4=-999` is the invalid sentinel. MPI boundary migration preserves
  IDs and cold interval state. No active per-timestep global-halo ID assignment
  path is confirmed in the current driver.
- v1.2+ implements driver-level cleanup for terminal precipitation conversion
  after `sdm_sd2prec` and for vertical outflow after vertical boundary
  processing. These paths snapshot live slots before the operation, detect
  newly invalid terminal slots afterward, optionally write
  `LIFE_SDREMOVE_INVALIDATION`, and call the centralized tracking cleanup
  helper. Horizontal MPI source invalidation and collision invalidation remain
  untouched because their event/send-buffer ordering needs a separate
  send-buffer-safe design. Controlled local runtime smokes validate the
  approved terminal and vertical cleanup paths; `sdm_sdremove` / `sdm_sdadd`
  production runtime validation remains blocked because production
  `sdm_adjsdnum` is inactive. The test-only lifecycle harness validates helper
  semantics only.
- v1.2 activates optional output-group switches:
  `TRACK_COLD_OUTPUT_ICE_GEOMETRY`,
  `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY`,
  `TRACK_COLD_OUTPUT_AEROSOL_CONTEXT`,
  `TRACK_COLD_OUTPUT_THERMO_CONTEXT`, and
  `TRACK_COLD_OUTPUT_KOHLER_CONTEXT`. Ice geometry, rime morphology, and Kohler
  context control real optional output groups. The aerosol and thermodynamic
  switches currently write safe single-process context variables; collision
  aerosol and thermodynamic context remain future work.
  `TRACK_COLD_EVENT_EXTENDED_GEOMETRY` remains as a backward-compatible broad
  alias for the first two groups.
- v1.2+ adds spatial-visit tracking with `sd_spatial_visit_flag`; the current
  version uses Level-3 segment-box intersection plus Level-2 point checks.
  v1.3+ adds default-off rank-subdomain pruning; enabled spatial tracking uses
  conservative all-active ranks unless
  `tracking_spatial_rank_pruning_enable = .true.`.
- Long-term work includes deliquescence/efflorescence, dry/wet aerosol
  phase-state expansion, aerosol-core and water-film morphology, lifecycle
  tracking, primary ice nucleation source tracking, secondary ice /
  fragmentation, multi-region spatial-visit metadata, hail wet/dry growth
  attribution, and OpenMP reproducibility characterization.

## 2. Core Tracking State Variables (`sd_id/dm_id`, `pre_sdid/pre_dmid`, `if_coal`)
The merged implementation retains both identifier systems because FW and BW encode different trajectory semantics:
- `sd_id`, `dm_id`: forward-identity pair.
- `pre_sdid`, `pre_dmid`: backward predecessor-link pair.
- `if_coal`: coalescence marker in one SD-output interval (`0`: no recorded coalescence since last SD dump, `1`: at least one coalescence was recorded).

### 2.1 Assignment timing and rules
1. **Runtime storage**
   - The model keeps one in-memory identifier pair, `sd_id/dm_id`, for both FW and BW.
   - BW does not maintain a different in-memory initialization for `pre_sdid/pre_dmid`; instead, the BW writers output the same runtime arrays under the variable names `pre_sdid/pre_dmid`.
2. **Initial selection and subset refresh**
   - Tracking IDs are assigned in `sdm_select_stratified_random_particles`.
   - For each selected SD, `sd_id = local_index` and `dm_id = mype`.
   - For each unselected or invalid SD, `sd_id = INVALID_i4` and `dm_id = INVALID_i4`.
   - `if_coal` is set to `0` for all SDs whenever the tracked subset is initialized or refreshed.
   - If `tracking_fraction <= 0`, no SD is tracked.
   - If `tracking_fraction >= 1` and `max_tracked_sds = 0`, all valid SDs are tracked.
   - Therefore, under otherwise identical settings, the first BW output writes `pre_sdid/pre_dmid` with the same numerical values that FW writes as `sd_id/dm_id` for the same SD.
3. **Microphysical coalescence stage**
   - When coalescence-event recording is enabled, every recorded coalescence event sets `if_coal = 1` for both involved SDs during the current SD-output interval.
   - When coalescence-event recording is disabled, `if_coal` is neither defined in SD snapshot files nor updated by the coalescence-event writer.
4. **SD snapshot output stage**
   - FW writes the runtime arrays as `sd_id/dm_id`; BW writes the same runtime arrays as `pre_sdid/pre_dmid`.
   - In `sdm_outnetcdf`, `if_coal` is written and then reset to `0`, so each SD dump stores whether that SD participated in at least one recorded coalescence since the previous SD dump.
   - In history-style output (`sdm_outnetcdf_hist`), BW refreshes the tracked subset after writing through `sdm_assign_tracking_subset`; this does not draw a new random subset, but remaps the currently retained tracked SDs to their current local indices for the next output window and simultaneously zeros `if_coal`.

### 2.2 Forward identifier semantics (aligned with `SDM_selected_SD`)
- **Forward tracking (`tracking_mode = 1`)**:
  - `dm_id`: domain identifier (`mype`) of the tracked SD.
  - `sd_id`: SD identifier within that domain for the tracked SD.
  - In forward chain analysis, `(dm_id, sd_id)` is the identity key carried along the trajectory.
  - Normal NetCDF writes `sd_id`/`dm_id`; history NetCDF writes `sd_id_####`/`dm_id_####`.
  - Coalescence files store `sd_id1/2`, `dm_id1/2`.

### 2.3 Backward identifier semantics (aligned with `SDM_SD_tracking`)
- **Backward tracking (`tracking_mode = 2`)**:
  - `pre_dmid`: domain identifier of the predecessor SD at the previous output interval.
  - `pre_sdid`: SD identifier of the predecessor SD at the previous output interval.
  - In backward chain analysis, `(pre_dmid, pre_sdid)` links the current record to earlier output levels.
  - Normal NetCDF writes `pre_sdid`/`pre_dmid`; history NetCDF writes `pre_sdid_####`/`pre_dmid_####`.
  - Coalescence files store `pre_sdid1/2`, `pre_dmid1/2`.

### 2.4 Why both name sets exist
- `sd_id/dm_id` is the forward identity key.
- `pre_sdid/pre_dmid` is the backward predecessor key.
- They have different meanings and must be interpreted with the corresponding tracking direction.

## 3. Sampling Algorithms
Sampling is executed by `sdm_select_stratified_random_particles` and is shared by forward/backward tracking.

### 3.1 Random sampling
- Candidate set: all valid SDs (`sd_rk > VALID2INVALID`).
- Randomization rule: for each candidate, the code draws from the tracking-specific RNG initialized by `tracking_sampling_seed + rank` and keeps the SD when `rand_tracking <= tracking_fraction`.
- Cap control: if `max_tracked_sds > 0`, accepted count is hard-capped; once the cap is reached, later candidates in the loop are rejected.
- Corner behavior:
  - `tracking_fraction <= 0`: no SD is tracked.
  - `tracking_fraction >= 1` with `max_tracked_sds = 0`: all valid SDs are tracked.

### 3.2 Stratified sampling
Stratified sampling is activated when:
- `tracking_selection_mode = "stratified"`,
- `tracking_nz_bin > 0`, `tracking_nr_bin > 0`,
- `tracking_height_max > tracking_height_min`.

Algorithm steps:
1. Build candidate set constrained by height/radius:
   - `tracking_height_min <= z <= tracking_height_max`,
   - `r >= tracking_radius_min`,
   - and if `tracking_radius_max > tracking_radius_min`, additionally `r <= tracking_radius_max`; otherwise upper radius bound is taken from runtime candidate maximum.
2. Compute target size:
   - `target_cnt = round(candidate_cnt * tracking_fraction)`,
   - lower bounded to 1 when candidates exist and fraction is positive,
   - upper bounded by `candidate_cnt` and `max_tracked_sds` (if set).
3. Build 2D bins (`z`, `log(r)`) with `tracking_nz_bin × tracking_nr_bin`.
4. Assign per-bin quota proportional to bin population, so each bin is sampled at nearly the same fraction of its local population.
5. Enforce minimum quota in non-empty bins using `tracking_min_per_bin` under global budget consistency.
6. After integer truncation, distribute any remaining quota across bins that still have spare candidates; if quota must be reduced, remove excess while respecting enforced minimum counts.
7. Perform sequential stochastic selection within each bin using dynamic acceptance probability:
   - `needed_in_bin / remaining_in_bin`.

This design makes the stratified sampler a random equal-proportion selector: the target is not an equal count from each bin, but an approximately equal retained fraction across the stratified population, subject to integer rounding and optional minimum-per-bin constraints.

### 3.3 Stratified constraints and fallback
- `tracking_min_per_bin` is not unconditional; it is effectively clipped by global target size and non-empty bin count.
- If stratified selection cannot be applied:
  - with `tracking_fallback_to_random = .true.`: fallback to random sampling;
  - with `tracking_fallback_to_random = .false.`: no SD is selected.

## 4. Forward Tracking Algorithm
Forward mode (`tracking_mode = 1`) records trajectories from earlier to later simulation time:
1. Initialize tracked subset by the selected sampling algorithm.
2. Propagate `(dm_id, sd_id)` together with SD motion and MPI boundary exchanges.
3. Record selected SD states to `SD_selected_NetCDF_*` when selected-output digit of `sdm_dmpvar` is enabled.
4. Record collision events into `SD_coal_output_NetCDF_*` when coalescence event output is enabled.
5. When `tracking_id_output_basename` is set and at least one interest-condition switch is enabled, FW additionally appends a raw per-rank interest-ID stream to `tracking_id_output_basename.peXXXXXX.ids`; duplicate removal is deferred to offline post-processing.

In TPHT mode, the FW initialization keeps the original sampling interface intact, but the effective target subset for BW is no longer defined by the one-time FW sampling output. Instead, the final BW target set is the union of all FW SDs that ever satisfied the interest condition(s) during the run.

Forward outputs are directly interpreted as current-time SD identities.

## 5. Backward Tracking Algorithm
Backward mode (`tracking_mode = 2`) uses the same runtime fields and output infrastructure, but the trajectory interpretation is reversed during analysis:
1. Initialize a target subset with the same sampling subsystem as forward mode; numerically, the initial BW identifiers are the same runtime `sd_id/dm_id` values that FW uses, but BW writes them to file as `pre_sdid/pre_dmid`.
2. At each saved time, use predecessor IDs (`pre_dmid`, `pre_sdid`) to connect an SD to earlier outputs.
3. In history-style BW output, after one dump is written, `sdm_assign_tracking_subset` updates the in-memory identifier pair of the currently retained tracked SDs to their current local `(dm_id, sd_id)` values for the next output window; therefore BW and FW are identical at initialization, but BW later rewrites predecessor links between successive outputs.
4. Use `if_coal` and coalescence event files to detect branch/merge points in ancestry chains.

The BW initialization now supports two priority levels:
- **Priority 1:** when `tracking_id_input_basename` is non-empty, BW reads the target SD set from that input file.
- **Priority 2:** only when no input basename is supplied does BW fall back to the configured random/stratified sampling settings.

This priority rule is what makes TPHT operational: BW can consume the cumulative FW-discovered target set without changing the legacy sampling interface used by ordinary BW experiments.

Runtime implementation for SD dumps and coalescence-event dumps is shared with forward mode; the key difference is chain reconstruction direction during post-processing.

When tracking is enabled and at least one SD dump is written in a time step, `LOG.pe000000` records `tracking_chain_count`, memory estimates for the tracking-ID arrays and `if_coal`, and average timings for ID assignment and lateral-boundary tracking operations. Here `tracking_chain_count` is the number of SDs whose identifier pair is currently valid in the dump step, i.e., the number of SD records that remain on an active tracking chain at that output.

## 6. Collision-Coalescence Event Recording
Collision events are evaluated at each microphysical coalescence sub-step. The effective writing logic is as follows:
- When coalescence-event recording is enabled, every recorded coalescence event is written to `SD_coal_output_NetCDF_*`, regardless of whether FW or BW is enabled.
- When FW or BW is enabled, the event file additionally stores the corresponding tracking identifiers.
- When both FW and BW are disabled, coalescence events are still recorded, but identifier variables are omitted because the event partners have no tracking IDs.
- When the microphysical coalescence process itself is disabled, coalescence-event recording is forcibly disabled, regardless of the namelist setting.

### 6.1 Forward vs backward: what is different and what is not
- **Runtime event detection and physical event storage:** identical for FW and BW.
- **Difference in post-processing interpretation:**
  - FW treats events as downstream trajectory evolution.
  - BW treats events as ancestry constraints for predecessor reconstruction.

### 6.2 Tracking-state dependence of event output
- **FW or BW enabled:** event files contain both physical variables and the corresponding tracking identifiers.
- **FW and BW both disabled:** event files still contain physical event variables, but do not define `sd_id*`, `dm_id*`, `pre_sdid*`, or `pre_dmid*`.
- **Coalescence-event recording disabled:** no event file is produced, and `if_coal` is not written to SD snapshot output.

### 6.3 Event variables and identifier semantics
- `event_time`: model time (seconds).
- **Forward mode (`tracking_mode = 1`)**: `sd_id1`, `dm_id1`, `sd_r1`, `sd_n1` and `sd_id2`, `dm_id2`, `sd_r2`, `sd_n2`.
- **Backward mode (`tracking_mode = 2`)**: `pre_sdid1`, `pre_dmid1`, `sd_r1`, `sd_n1` and `pre_sdid2`, `pre_dmid2`, `sd_r2`, `sd_n2`.
- **No tracking mode (`tracking_mode = 0`, or both FW/BW disabled)**: collision files contain only physical event fields (`event_time`, `sd_r1`, `sd_n1`, `sd_r2`, `sd_n2`, `event_multiplicity`, and the common pre/post diagnostics) and omit tracking identifiers.
- The difference between FW and BW is restricted to identifier semantics (`current identity` versus `predecessor identity`); physical fields are shared.
- `event_multiplicity`: output variable for collision multiplicity/count. The
  warm and cold coalescence routines keep the internal local variable name
  `num_col`, but this name is not written to the current event files. The warm
  `SD_coal_output_NetCDF_*` rename is intentional: `event_multiplicity` is the
  event-output schema name for internal `num_col`.
- **Warm SDM (`sdm_cold = .false.`)**: legacy collision records are still
  written to `SD_coal_output_NetCDF_*`. They write `event_multiplicity`,
  pre-event midpoint `x/y/z`, pre/post `sd_n` for both participants, and
  pre/post warm `hydro_radius`/`hydro_mass`. In warm mode `hydro_radius` is the
  liquid droplet radius `sd_r`, and `hydro_mass` is `(4/3) * pi * CONST_DWATR *
  sd_r**3`.
- **Cold SDM (`sdm_cold = .true.`)**: cold collision records are written to
  `SD_event_collision_NetCDF_*`. They use categorical process `trigger_code`
  (`1`: coalescence, `2`: riming, `3`: aggregation) and `trigger_level`
  (`1`: occurrence, `2`: significant). In the current cold collision path,
  coalescence means liquid-liquid coalescence.
  `event_multiplicity` for the internal `num_col` value,
  `phase_state1_pre/post`, `phase_state2_pre/post`, participant IDs,
  pair-midpoint `x/y/z`, participant
  coordinates `x1/y1/z1` and `x2/y2/z2`, and pre/post `sd_n`,
  `hydro_radius`, and `hydro_mass` fields.

### 6.4 Time-step behavior and file append rule
- Coalescence condition is evaluated every microphysics coalescence step.
- `SD_coal_output_NetCDF_*` is append-based along unlimited `event` dimension.
- Files are split by output-time bucket and MPI rank.

## 7. SD Output, `sdm_dmpvar`, and File Naming Rules
`sdm_dmpvar` is a digit-coded output switch:
- 1st digit (`10^0`): ASCII SD output (`sdm_dmpitva`, unsupported in this workflow).
- 2nd digit (`10^1`): all-SD NetCDF output (`sdm_dmpitvb`).
- 3rd digit (`10^2`): selected-SD NetCDF output (`sdm_dmpitvl`).
- Digit value meaning: `0` off, `1` normal NetCDF (`sdm_outnetcdf`), `2` history-style NetCDF (`sdm_outnetcdf_hist`).

### 7.1 Typical settings and outputs
- `sdm_dmpvar = 010`: output all SDs only (`SD_all_NetCDF_*`).
- `sdm_dmpvar = 100`: output selected SDs only (`SD_selected_NetCDF_*`).
- `sdm_dmpvar = 020`: use history-style NetCDF writer for all SDs
  (`SD_all_history.peXXXXXX`, with suffixed variables such as
  `sd_event_mask_0001`).
- `sdm_dmpvar = 200`: use history-style NetCDF writer for selected SDs
  (`SD_selected_history.peXXXXXX`, with suffixed variables).

### 7.2 Effect of sampling on selected output
- **Sampling enabled and partial (`0 < tracking_fraction < 1`):**
  - `SD_selected_NetCDF_*` contains only selected tracked subset.
- **Sampling effectively full (`tracking_fraction = 1`, `max_tracked_sds = 0`):**
  - all valid SDs are tagged as tracked, so `SD_selected_NetCDF_*` contains the full valid-SD set whenever the selected-output digit of `sdm_dmpvar` is enabled.
  - The 2nd digit of `sdm_dmpvar` may remain `0`; selected output depends on the 3rd digit and `sdm_dmpitvl`, not on all-SD output being enabled.
- **TPHT interest-ID discovery (`tracking_id_output_basename` enabled with interest conditions):**
  - FW may track the full valid-SD set for ID propagation while writing the BW target set to a separate rank-local `.ids` handoff stream.
  - In the recommended TPHT configuration, `sdm_dmpvar = 000`, so `SD_selected_NetCDF_*` is not used as the carrier of the BW target definition.
- **Sampling disabled or empty selection:**
  - selected-output pipeline still exists when enabled by `sdm_dmpvar`, but selected count can be zero.
- **Coalescence-event recording enabled:**
  - `if_coal` is written to SD snapshot output irrespective of whether FW or BW is enabled.
- **Coalescence-event recording disabled:**
  - `if_coal` is not defined in SD snapshot files.

### 7.3 File naming patterns
- Default SD output without file tag: `SD_output_NetCDF_YYYYMMDD-hhmmss.mmm.peXXXXXX`
- All-SD output: `SD_all_NetCDF_YYYYMMDD-hhmmss.mmm.peXXXXXX`
- Selected-SD output: `SD_selected_NetCDF_YYYYMMDD-hhmmss.mmm.peXXXXXX`
- All-SD history output: `SD_all_history.peXXXXXX`
- Selected-SD history output: `SD_selected_history.peXXXXXX`
- Coalescence events: `SD_coal_output_NetCDF_YYYYMMDD-hhmmss.mmm.peXXXXXX`
- TPHT rank-local interest-ID output: `tracking_id_output_basename.peXXXXXX.ids`
- In this document and scripts, `YYYYMMDD` is the real simulation date token and is not fixed to `00000101`.

### 7.4 Supported/unsupported `sdm_dmpvar` modes in this tracking workflow
For the current SD-tracking workflow and provided test/tutorial cases, the following modes are supported and validated:
- `000`: no SD snapshot output; used in TPHT FW discovery cases where the BW target set is written through `tracking_id_output_basename`
- `010`: all-SD NetCDF output (`SD_all_NetCDF_*`)
- `100`: selected-SD NetCDF output (`SD_selected_NetCDF_*`)
- `020`: all-SD history NetCDF output (`SD_all_history.peXXXXXX` with suffixed history variables)
- `200`: selected-SD history NetCDF output (`SD_selected_history.peXXXXXX` with suffixed history variables)

The file basename is determined by both the writer tag (`all` or `selected`)
and whether the normal or history-style NetCDF writer is selected. Normal output
uses the time-stamped `SD_all_NetCDF_*` / `SD_selected_NetCDF_*` names, whereas
history-style output uses persistent `SD_all_history.peXXXXXX` /
`SD_selected_history.peXXXXXX` files with per-output-time variable suffixes.

The following modes are not supported in this workflow:
- ASCII-only mode (`001`)
- large-droplet-only mode (third digit as large-droplet binary option in legacy SDM usage)
- simultaneous all+selected output mode (`110`)

## 8. Namelist Parameters Added/Extended
Main parameters under `&PARAM_ATMOS_PHY_MP_SDM` are summarized below with code defaults.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tracking_mode` | integer | `0` | Unified tracking switch: `0` no tracking, `1` forward, `2` backward. |
| `tracking_selection_mode` | string | `"random"` | Sampling algorithm selector: `"random"` applies Bernoulli sampling, `"stratified"` allocates per-bin quotas, and `"none"` means no additional random/stratified sampling mode is requested. |
| `tracking_fraction` | real | `1.0` | Fraction of the candidate set retained after initialization. `1.0` means “do not downsample”. |
| `max_tracked_sds` | integer | `0` | Hard cap for tracked SD count (`0` means unlimited). |
| `tracking_sampling_seed` | integer | `0` | Dedicated seed for random/stratified tracking sampling. It initializes a per-rank tracking RNG as `tracking_sampling_seed + rank` and does not change SDM initialization, coalescence, or the main SDM RNG stream. |
| `tracking_height_min` | real [m] | `400.0` | Lower height bound for stratified candidate filter. |
| `tracking_height_max` | real [m] | `800.0` | Upper height bound for stratified candidate filter. |
| `tracking_radius_min` | real [m] | `1.0e-6` | Lower radius bound for stratified candidate filter. |
| `tracking_radius_max` | real [m] | `0.0` | Upper radius bound; if `<= tracking_radius_min`, runtime candidate maximum radius is used. |
| `tracking_nz_bin` | integer | `8` | Number of height bins in stratified sampling. |
| `tracking_nr_bin` | integer | `10` | Number of radius bins in stratified sampling. |
| `tracking_min_per_bin` | integer | `1` | Requested minimum tracked SDs in each non-empty bin (subject to global budget). |
| `tracking_fallback_to_random` | logical | `.true.` | Fallback from stratified to random when stratified is not feasible. |
| `tracking_id_input_basename` | string | `''` | Optional BW target-set input. If non-empty, BW reads tracked IDs from this basename and overrides random/stratified initialization. The runtime appends `.peXXXXXX.ids` when needed, and an exact `.nc` path is automatically mapped to its sibling `.ids` file. The recommended TPHT handoff is a deduplicated per-rank basename produced by the offline post-processing utility. BW also validates the embedded `TPHT_META` decomposition header and aborts if `PRC_NUM_X`, `PRC_NUM_Y`, or total MPI rank count do not match the FW handoff file. |
| `tracking_id_output_basename` | string | `''` | Optional FW interest-ID output basename. When enabled together with at least one interest condition, FW appends selected `(dm_id, sd_id)` records to `tracking_id_output_basename.peXXXXXX.ids`. Duplicate removal is performed offline, and each append block writes `# TPHT_META PRC_NUM_X PRC_NUM_Y PRC_nprocs` so that BW and merge tools can verify decomposition compatibility. |
| `tracking_interest_radius_enable` | logical | `.false.` | Enable radius-threshold interest detection during FW discovery. |
| `tracking_interest_radius_threshold` | real [m] | `0.0` | Radius threshold used when `tracking_interest_radius_enable = .true.`; SDs with `sd_r >= tracking_interest_radius_threshold` are treated as interest targets. |
| `tracking_interest_coalescence_enable` | logical | `.false.` | Enable interest detection for SDs that participated in at least one recorded coalescence event in the current SD-output interval. |
| `tracking_interest_ice_radius_enable` | logical | `.false.` | Enable cold-SDM ice-size interest detection during FW discovery. This uses the phase-aware volume-equivalent `hydro_radius` helper when `sdm_cold = .true.`. |
| `tracking_interest_ice_radius_threshold` | real [m] | `0.0` | Phase-aware volume-equivalent radius threshold used when `tracking_interest_ice_radius_enable = .true.`. |
| `tracking_interest_ice_phase_enable` | logical | `.false.` | Enable cold-SDM phase interest detection for SDs whose `sd_liqice` status is ice or mixed phase. |
| `tracking_interest_rime_mass_enable` | logical | `.false.` | Enable cold-SDM rimed-mass interest detection during FW discovery. |
| `tracking_interest_rime_mass_threshold` | real [kg] | `0.0` | Rimed-mass threshold used when `tracking_interest_rime_mass_enable = .true.`. |
| `coalescence_output_enable` | logical | `.true.` | Master switch for writing `SD_coal_output_NetCDF_*`. |
| `random_perturbation_enable` | logical | `.false.` | Master switch of SD motion perturbation. |
| `random_perturbation_amp` | real [m^1.5 s^-0.5] | `0.0` | User-facing perturbation amplitude; the configuration routine copies it to the internal runtime variable `sdm_noise_amp`. |

`tracking_sample_initialized` is an internal runtime state variable rather than a user-facing namelist control, so it is intentionally omitted from the table above.

### 8.1 Cold tracking namelist

Cold TPHT v1 adds `&PARAM_ATMOS_PHY_MP_SDM_TRACKING`. If the namelist is absent,
the code uses the defaults below.

| Parameter | Type | Default | Controls |
|---|---|---|---|
| `TRACK_COLD_EVENT_EXTENDED_GEOMETRY` | logical | `.false.` | Backward-compatible broad alias. If true, setup enables `TRACK_COLD_OUTPUT_ICE_GEOMETRY` and `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY`. Participant coordinates `x1/y1/z1` and `x2/y2/z2` are core fields and do not require this switch. |
| `TRACK_COLD_OUTPUT_ICE_GEOMETRY` | logical | `.false.` | Optional `ice_re*`, `ice_rp*`, and `ice_rho*` pre/post fields for collision and single-process event files. |
| `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY` | logical | `.false.` | Optional collision-event `rime_mass*`, `rime_frac*`, and `aspect_ratio*` pre/post fields. Ordinary-output interval maxima are core fields and are not controlled by this switch. |
| `TRACK_COLD_OUTPUT_AEROSOL_CONTEXT` | logical | `.false.` | Optional `SD_event_singleproc_NetCDF_*` aerosol context fields: `aerosol_total_mass_pre/post` and `aerosol_kohler_solute_pre/post`. Collision aerosol context remains future work. |
| `TRACK_COLD_OUTPUT_THERMO_CONTEXT` | logical | `.false.` | Optional `SD_event_singleproc_NetCDF_*` thermodynamic context fields: `air_temperature`, `air_pressure`, and `water_vapor_mixing_ratio`. Collision thermodynamic context remains future work. |
| `TRACK_COLD_OUTPUT_KOHLER_CONTEXT` | logical | `.false.` | Optional single-process fields `kohler_rcrit_pre/post`, `kohler_margin_pre/post`, and `activated_state_pre/post` for derived Kohler activation/deactivation records. |
| `tracking_evt_coalescence_enable` | logical | `.false.` | Enables coalescence occurrence records, `trigger_code=1`, `trigger_level=1`, `EVENT_COALESCENCE`. In the current cold collision path, coalescence means liquid-liquid coalescence. |
| `tracking_evt_liq_liq_coal_enable` | logical | `.false.` | Deprecated alias for `tracking_evt_coalescence_enable`; compatibility is conservative OR-style, so do not set old and new names inconsistently. |
| `tracking_evt_riming_enable` | logical | `.false.` | Enables riming occurrence records, `trigger_code=2`, `EVENT_RIMING`. |
| `tracking_evt_aggregation_enable` | logical | `.false.` | Enables aggregation occurrence records, `trigger_code=3`, `EVENT_AGGREGATION`. |
| `tracking_evt_freezing_enable` | logical | `.false.` | Enables freezing occurrence records, `trigger_code=4`, `EVENT_FREEZING`. |
| `tracking_evt_melting_enable` | logical | `.false.` | Enables melting occurrence records, `trigger_code=5`, `EVENT_MELTING`. |
| `tracking_evt_deposition_enable` | logical | `.false.` | Enables deposition occurrence records, `trigger_code=6`, `trigger_level=1`, `EVENT_DEPOSITION`; continuous vapor-growth occurrence output should remain opt-in. |
| `tracking_evt_sublimation_enable` | logical | `.false.` | Enables sublimation occurrence records, `trigger_code=7`, `trigger_level=1`, `EVENT_SUBLIMATION`; continuous vapor-growth occurrence output should remain opt-in. |
| `tracking_evt_condensation_enable` | logical | `.false.` | Enables condensation occurrence records, `trigger_code=8`, `trigger_level=1`, `EVENT_CONDENSATION`; continuous vapor-growth occurrence output should remain opt-in. |
| `tracking_evt_evaporation_enable` | logical | `.false.` | Enables evaporation occurrence records, `trigger_code=9`, `trigger_level=1`, `EVENT_EVAPORATION`; continuous vapor-growth occurrence output should remain opt-in. |
| `tracking_evt_activation_enable` | logical | `.false.` | Enables derived occurrence-level activation crossing records, `trigger_code=10`, `trigger_level=1`, using the shared Kohler selected-output critical-radius criterion around `sdm_condevp`. |
| `tracking_evt_deactivation_enable` | logical | `.false.` | Enables derived occurrence-level deactivation crossing records, `trigger_code=11`, `trigger_level=1`, using the shared Kohler selected-output critical-radius criterion around `sdm_condevp`. |
| `tracking_sig_deposition_enable` | logical | `.false.` | Enables significant-deposition single-process hit records, `trigger_code=6`, `trigger_level=2`, `EVENT_DEPOSITION` in `sd_event_mask` and `sd_event_sig_mask`. |
| `tracking_sig_deposition_relmass_threshold` | real | `huge(1.0_RP)` | Relative ice-mass-change threshold for significant deposition. |
| `tracking_sig_sublimation_enable` | logical | `.false.` | Enables significant-sublimation single-process hit records, `trigger_code=7`, `trigger_level=2`. |
| `tracking_sig_sublimation_relmass_threshold` | real | `huge(1.0_RP)` | Relative ice-mass-change threshold for significant sublimation. |
| `tracking_sig_condensation_enable` | logical | `.false.` | Enables significant-condensation single-process hit records, `trigger_code=8`, `trigger_level=2`. |
| `tracking_sig_condensation_relmass_threshold` | real | `huge(1.0_RP)` | Relative liquid-mass-increase threshold for significant condensation. |
| `tracking_sig_evaporation_enable` | logical | `.false.` | Enables significant-evaporation single-process hit records, `trigger_code=9`, `trigger_level=2`. |
| `tracking_sig_evaporation_relmass_threshold` | real | `huge(1.0_RP)` | Relative liquid-mass-decrease threshold for significant evaporation. |
| `tracking_sig_freezing_enable` | logical | `.false.` | Enables significant-freezing single-process hit records, `trigger_code=4`, `trigger_level=2`. |
| `tracking_sig_freezing_relmass_threshold` | real | `huge(1.0_RP)` | Relative component-mass-transfer threshold for significant freezing. |
| `tracking_sig_melting_enable` | logical | `.false.` | Enables significant-melting single-process hit records, `trigger_code=5`, `trigger_level=2`. |
| `tracking_sig_melting_relmass_threshold` | real | `huge(1.0_RP)` | Relative component-mass-transfer threshold for significant melting. |
| `tracking_sig_riming_enable` | logical | `.false.` | Enables significant-riming collision-family records, `trigger_code=2`, `trigger_level=2`. |
| `tracking_sig_riming_relmass_threshold` | real | `huge(1.0_RP)` | Relative collected-liquid/collector-mass threshold for significant riming. |
| `tracking_sig_aggregation_enable` | logical | `.false.` | Enables significant-aggregation collision-family records, `trigger_code=3`, `trigger_level=2`. |
| `tracking_sig_aggregation_relmass_threshold` | real | `huge(1.0_RP)` | Relative smaller-ice-mass/larger-ice-mass threshold for significant aggregation. |
| `tracking_sig_coalescence_enable` | logical | `.false.` | Enables significant coalescence records, `trigger_code=1`, `trigger_level=2`. |
| `tracking_sig_coalescence_relmass_threshold` | real | `huge(1.0_RP)` | Relative smaller-liquid-mass/larger-liquid-mass threshold for significant coalescence. |
| `tracking_sig_liq_liq_coal_enable` | logical | `.false.` | Deprecated alias for `tracking_sig_coalescence_enable`; compatibility is conservative OR-style. |
| `tracking_sig_liq_liq_coal_relmass_threshold` | real | `huge(1.0_RP)` | Deprecated alias threshold for `tracking_sig_coalescence_relmass_threshold`. |
| `tracking_sig_activation_enable` | logical | `.false.` | Enables significant derived Kohler activation records, `trigger_code=10`, `trigger_level=2`. |
| `tracking_sig_activation_radius_threshold` | real [m] | `huge(1.0_RP)` | Absolute post-crossing Kohler margin threshold: `sd_r_post - r_crit_post`. |
| `tracking_sig_deactivation_enable` | logical | `.false.` | Enables significant derived Kohler deactivation records, `trigger_code=11`, `trigger_level=2`. |
| `tracking_sig_deactivation_radius_threshold` | real [m] | `huge(1.0_RP)` | Absolute post-crossing Kohler negative-margin threshold: `r_crit_post - sd_r_post`. |
| `tracking_diag_liq_radius_enable` | logical | `.false.` | Enables liquid-radius diagnostic bit `DIAG_LIQ_RADIUS_LARGE`, `trigger_code=101`. |
| `tracking_diag_liq_radius_threshold` | real [m] | `huge(1.0_RP)` | Interval maximum liquid-radius threshold. |
| `tracking_diag_ice_rvol_enable` | logical | `.false.` | Enables ice volume-equivalent radius diagnostic bit `DIAG_ICE_RVOL_LARGE`, `trigger_code=102`. |
| `tracking_diag_ice_rvol_threshold` | real [m] | `huge(1.0_RP)` | Interval maximum ice volume-equivalent radius threshold. |
| `tracking_diag_mixed_rvol_enable` | logical | `.false.` | Enables mixed-phase volume-equivalent radius diagnostic bit `DIAG_MIXED_RVOL_LARGE`, `trigger_code=103`. |
| `tracking_diag_mixed_rvol_threshold` | real [m] | `huge(1.0_RP)` | Interval maximum mixed-phase volume-equivalent radius threshold. |
| `tracking_diag_rime_mass_enable` | logical | `.false.` | Enables rime-mass diagnostic bit `DIAG_RIME_MASS_LARGE`, `trigger_code=104`. |
| `tracking_diag_rime_mass_threshold` | real [kg] | `huge(1.0_RP)` | Interval maximum rime-mass threshold. |
| `tracking_diag_rime_frac_enable` | logical | `.false.` | Enables rime-fraction diagnostic bit `DIAG_RIME_FRAC_LARGE`, `trigger_code=105`. |
| `tracking_diag_rime_frac_threshold` | real | `huge(1.0_RP)` | Interval maximum rime-mass-fraction threshold. |
| `tracking_diag_nmono_enable` | logical | `.false.` | Enables monomer-count diagnostic bit `DIAG_NMONO_LARGE`, `trigger_code=106`. |
| `tracking_diag_nmono_threshold` | real | `huge(1.0_RP)` | Interval maximum monomer-count threshold. |
| `tracking_diag_aspect_ratio_enable` | logical | `.false.` | Enables aspect-ratio diagnostic bit `DIAG_ASPECT_RATIO_MATCH`, `trigger_code=107`. |
| `tracking_diag_aspect_ratio_threshold` | real | `huge(1.0_RP)` | Interval maximum ice-aspect-ratio threshold. |
| `tracking_spatial_visit_enable` | logical | `.false.` | Enables v1.2+ spatial-visit interval tracking for one box. Current code uses Level-3 segment-box intersection plus Level-2 point checks. |
| `tracking_spatial_x_min` / `tracking_spatial_x_max` | real [m] | `huge` / `-huge` | Spatial-visit x bounds. Bounds must be valid when spatial tracking is enabled. |
| `tracking_spatial_y_min` / `tracking_spatial_y_max` | real [m] | `huge` / `-huge` | Spatial-visit y bounds. |
| `tracking_spatial_z_min` / `tracking_spatial_z_max` | real [m] | `huge` / `-huge` | Spatial-visit z bounds. |
| `tracking_spatial_rank_pruning_enable` | logical | `.false.` | Enables v1.3+ rank-subdomain pruning. Default off preserves conservative all-rank-active behavior. |
| `tracking_spatial_rank_margin` | real [m] | `0.0` | Safety margin used when rank-subdomain pruning compares the configured spatial box with each rank-local bbox. |

Default strategy: all occurrence-level records are disabled; all significant
process-hit triggers are disabled; diagnostic triggers are disabled; spatial
visit is disabled; all thresholds default to a huge/sentinel value. Smoke
configs explicitly enable the records they validate.

| Process | Occurrence enable | Process trigger_code | Occurrence level | Significant enable | Significant level | Significant threshold | Output file | Default behavior |
|---|---|---:|---|---|---:|---|---|---|
| coalescence | `tracking_evt_coalescence_enable` | 1 | `trigger_level=1` | `tracking_sig_coalescence_enable` | `trigger_level=2` | relative smaller/larger liquid mass | `SD_event_collision_NetCDF_*` | both off |
| riming | `tracking_evt_riming_enable` | 2 | `trigger_level=1` | `tracking_sig_riming_enable` | `trigger_level=2` | relative collected-liquid / collector mass | `SD_event_collision_NetCDF_*` | both off |
| aggregation | `tracking_evt_aggregation_enable` | 3 | `trigger_level=1` | `tracking_sig_aggregation_enable` | `trigger_level=2` | relative smaller/larger ice mass | `SD_event_collision_NetCDF_*` | both off |
| freezing | `tracking_evt_freezing_enable` | 4 | `trigger_level=1` | `tracking_sig_freezing_enable` | `trigger_level=2` | relative component mass transfer | `SD_event_singleproc_NetCDF_*` | both off |
| melting | `tracking_evt_melting_enable` | 5 | `trigger_level=1` | `tracking_sig_melting_enable` | `trigger_level=2` | relative component mass transfer | `SD_event_singleproc_NetCDF_*` | both off |
| deposition | `tracking_evt_deposition_enable` | 6 | `trigger_level=1` | `tracking_sig_deposition_enable` | `trigger_level=2` | relative ice-mass increase | `SD_event_singleproc_NetCDF_*` | both off |
| sublimation | `tracking_evt_sublimation_enable` | 7 | `trigger_level=1` | `tracking_sig_sublimation_enable` | `trigger_level=2` | relative ice-mass decrease | `SD_event_singleproc_NetCDF_*` | both off |
| condensation | `tracking_evt_condensation_enable` | 8 | `trigger_level=1` | `tracking_sig_condensation_enable` | `trigger_level=2` | relative liquid-mass increase | `SD_event_singleproc_NetCDF_*` | both off |
| evaporation | `tracking_evt_evaporation_enable` | 9 | `trigger_level=1` | `tracking_sig_evaporation_enable` | `trigger_level=2` | relative liquid-mass decrease | `SD_event_singleproc_NetCDF_*` | both off |
| activation | `tracking_evt_activation_enable` | 10 | `trigger_level=1`, derived Kohler crossing | `tracking_sig_activation_enable` | `trigger_level=2` | absolute post-crossing Kohler margin | `SD_event_singleproc_NetCDF_*` | both levels off by default |
| deactivation | `tracking_evt_deactivation_enable` | 11 | `trigger_level=1`, derived Kohler crossing | `tracking_sig_deactivation_enable` | `trigger_level=2` | absolute post-crossing Kohler negative margin | `SD_event_singleproc_NetCDF_*` | both levels off by default |

Reserved v1.2 trigger-code constants are documented in the helper module but
inactive: `20-23` for primary ice nucleation source candidates, `30-32` for
secondary ice / breakup / fragmentation candidates, and `40-41` for
deliquescence/efflorescence candidates. They are not written by current
process hooks and are not active `sd_event_mask` bits.

### 8.2 Namelist notes specific to TPHT
- **FW discovery pass:** set `tracking_mode = 1`, set `tracking_id_output_basename`, enable one or both interest-condition switches, and typically use `sdm_dmpvar = 000`.
- **BW reconstruction pass:** set `tracking_mode = 2` and point `tracking_id_input_basename` to the FW deduplicated per-rank ID basename.
- **Priority rule:** `tracking_id_input_basename` takes precedence over `tracking_fraction`, `tracking_selection_mode`, and the stratified bounds during BW initialization.
- **Common recommendation:** in TPHT FW cases, use `tracking_selection_mode = "none"`, `tracking_fraction = 1.0`, and `max_tracked_sds = 0` so that no extra random/stratified downsampling is introduced before the interest-condition reduction step.
- **MPI decomposition must match between FW and BW:** `PRC_NUM_X`, `PRC_NUM_Y`, and total MPI rank count must match the FW run that produced the handoff IDs. BW now checks the `TPHT_META` header in the input `.ids` file and aborts on mismatch.
- **Interest-condition combination:** warm and cold interest switches can each be used independently. If multiple switches are enabled, the current implementation uses logical OR; logical AND is not implemented in the present code.

## 9. Random Perturbations in SD Motion
The motion update uses a single effective perturbation amplitude:
- `random_perturbation_amp` is the user-set namelist parameter in `run.conf`;
- `sdm_noise_amp` is the internal runtime variable used in the motion routine;
- when `random_perturbation_enable = .true.`, the configuration routine sets `sdm_noise_amp = random_perturbation_amp`;
- when `random_perturbation_enable = .false.`, it sets `sdm_noise_amp = 0`.

Therefore, users should control the perturbation amplitude through `random_perturbation_amp`; `sdm_noise_amp` is only the internal variable that actually appears in the update formula.

In the motion routine, independent random numbers are generated for the `x`, `y`, and vertical directions, and the perturbations are added as:

$$\Delta x = (\xi_x - 0.5) \times \sqrt{\frac{\Delta t}{r}} \times \text{sdm\\_noise\\_amp}$$

$$\Delta y = (\xi_y - 0.5) \times \sqrt{\frac{\Delta t}{r}} \times \text{sdm\\_noise\\_amp}$$

$$\Delta k = (\xi_k - 0.5) \times \sqrt{\frac{\Delta t}{r}} \times \frac{1}{DZ} \times \text{sdm\\_noise\\_amp}$$

where:
- $\xi_x$ and $\xi_y$ are independent random numbers for the horizontal perturbations in the physical coordinates `x` and `y`;
- $\xi_k$ is an independent random number for the perturbation of the vertical SD coordinate `sd_rk`, which is a grid-index-like vertical coordinate rather than the physical height itself;
- the vertical formula therefore includes `1/DZ`, because the model first evaluates a vertical displacement amplitude in physical space and then converts it into the increment of the vertical grid coordinate;
- $\Delta t$ is the SD advection time step;
- $r$ is SD radius;
- `sdm_noise_amp` is the internal runtime copy of the user parameter `random_perturbation_amp`.

## 10. Installation and Usage
For more details, please see [the official user guide of SCALE](https://scale.riken.jp/archives/scale_users_guide_En.v5.2.6.pdf) (last access: 6 April 2026).

### Prerequisites
- **Compilers:** Fortran and C compilers are required.
- **Libraries:** MPI, NetCDF4, and HDF5 libraries must be installed.
- **Post-processing:** Python 3. `netCDF4` is preferred for NetCDF reads; the TPHT smoke checkers and summaries can fall back to the system `ncdump` executable. Plotting workflows still use packages such as `numpy` and `matplotlib`.

### Environment preparation
```bash
export SCALE_SYS=Linux64-intel-impi
```

### Lightweight FW/BW/TPHT test examples
This branch keeps the routine FW/BW/TPHT validation cases under `scale-rm/test/case/shallowcloud/` and avoids depending on manuscript-specific benchmark suites or generated analysis-output trees.

Use the following cases as the baseline tutorial/test matrix:
- **(a) 2D forward, with coalescence output:** `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward`
- **(b) 2D backward, with coalescence output:** `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward`
- **(c) 2D backward, no coalescence output and no large-scale horizontal wind:** `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward_no_coal_no_uv`
- **(d) 3D forward, with stratified sampling and coalescence output:** `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline`
- **(e) 3D backward, with stratified sampling and coalescence output:** `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline`
- **(f) Forward and backward representativeness ensembles:** `scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests` and `scale-rm/test/case/shallowcloud/backward_tracking_representativeness_tests`

For a no-coalescence-output variant of the sampling cases, copy case (d) or (e), then set `coalescence_output_enable = .false.` in the copied `run.conf`.

The retained TPHT tutorial pair is:
- **FW TPHT discovery case**  
  `scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline`
- **BW TPHT reconstruction case**  
  `scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline`

The retained ice-phase smoke pair is:
- **FW ice TPHT smoke case**:
  `scale-rm/test/case/shallowcloud/ice_tpht_smoke/fw_ice_smoke`
- **BW ice TPHT smoke case**:
  `scale-rm/test/case/shallowcloud/ice_tpht_smoke/bw_ice_smoke`

Use `scale-rm/test/case/shallowcloud/ice_tpht_smoke/README.md` for the
rank-1 local workflow that checks cold selected-output fields, cold TPHT ID
handoff, and cold event metadata schema.

### 10.1 Practical tutorial: Two-Pass Hybrid Tracking (TPHT)
This subsection gives a concrete, end-to-end workflow for running TPHT with the retained shallowcloud tutorial cases.

#### Step 1. Use the dedicated tutorial cases
The tutorial pair is:
- FW discovery case: `scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline`
- BW reconstruction case: `scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline`

These two tutorial cases are placed under the same `tpht_test/` parent directory so that the FW-to-BW handoff path stays short and FW/BW restart outputs can be isolated with `fw_output/` and `bw_output/`.

#### Step 2. Check the FW discovery namelist
Enter the FW tutorial directory first:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline
vi ./run.conf
```

Key TPHT FW settings are:

| Parameter | Recommended tutorial value | Purpose |
|---|---|---|
| `tracking_mode` | `1` | Enable FW discovery pass. |
| `tracking_id_output_basename` | `"./fw_tracking/tracking_interest_ids"` | Basename of the raw per-rank interest-ID output stream that is later merged and deduplicated offline. |
| `SD_OUT_BASENAME` | `"./fw_output/superdroplet_restart"` | Keep FW restart output separate from BW when both cases live under `tpht_test/`. |
| `tracking_interest_radius_enable` | `.true.` | Enable radius-based discovery. |
| `tracking_interest_radius_threshold` | `1.0d-6` | Mark SDs with radius greater than or equal to this value as “interest” targets. |
| `tracking_interest_coalescence_enable` | `.true.` | Also include SDs flagged by coalescence participation. |
| `tracking_selection_mode` | `"none"` | Do not request any additional random/stratified sampling mode in the FW discovery pass. |
| `tracking_fraction` | `1.0d0` | Disable downsampling so the effective FW initialization reduces to the configured height/radius window. |
| `max_tracked_sds` | `0` | No artificial cap on discovery candidates. |
| `sdm_dmpvar` | `000` | Disable `SD_selected_NetCDF_*`; use the separate raw interest-ID stream instead. |
| `HISTORY_DEFAULT_TINTERVAL` | `0.0D0` | Disable Eulerian history output in the FW discovery pass; the TPHT handoff artifact is `tracking_interest_ids.pe*.ids`. |

Notes:
- `tracking_mode` already drives FW/BW enablement, so `forward_tracking_enable` / `backward_tracking_enable` are not user-facing namelist keys in the current TPHT tutorial cases; they remain only as internal runtime flags derived from `tracking_mode`.
- `tracking_selection_mode="none"` means no additional random/stratified sampling mode is requested.
- `tracking_fraction=1.0d0` means the post-initialization downsampling step is disabled. In the TPHT FW tutorial, that leaves the configured height/radius window as the effective initialization filter.
- For the shallowcloud tracking cases in this repository, `ATMOS_DYN_TYPE` is set to `HEVI`, `TIME_DT = 0.1D0`, `TIME_DT_ATMOS_DYN = 0.05D0`, other atmospheric physics time steps are `0.1D0`, the sponge layer starts at `ATMOS_DYN_wdamp_height = 1400.D0`, and `ATMOS_DYN_wdamp_tau = 0.5D0`.
- `ATMOS_PHY_TB_SMG_horizontal` remains `.false.` here because these tutorial cases use `ATMOS_PHY_TB_TYPE = 'SMAGORINSKY'`, not the HYBRID gray-zone/RANS setting that requires horizontal LES viscosity.

#### Step 3. Build the executable
From the same FW TPHT case directory, rebuild SCALE-RM with SDM enabled. A generic example is:

```bash
module purge
module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1
make allclean
time make -j SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
ln -fsv `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
```

#### Step 4. Run the FW discovery pass
Still in the FW TPHT case directory, make sure the FW output directories exist, and submit the job:

```bash
mkdir -p fw_output fw_tracking
qsub UoH_run.pbs
```

or on SQUID:

```bash
mkdir -p fw_output fw_tracking
qsub squid_run.sh
```

For the Cold TPHT v1.3+ smoke suite, the prepared SQUID wrapper is:

```bash
cd scale-rm/test/case/shallowcloud/ice_tpht_smoke/squid_validation
qsub squid_run_v13_kohler_validation.sh
```

The v1.3+ SQUID full validation passed as `832697.sqd` with
`required_failures=0`. The v1.3+ production I/O calibration wrapper is
`squid_run_v13_production_io_threshold.sh`; it passed as `832700.sqd` with
`required_failures=0`. A later local follow-up added controlled vapor occurrence
coverage for `7:1`, `8:1`, and `9:1`; the completed SQUID job remains historical
for that optional gap. Older v1.2 SQUID production I/O results are historical
evidence only and should not be reported as v1.3 validation.

Longer v1.3 production calibration scripts are prepared but not submitted:
`squid_run_v13_long_production_calibration.sh` and
`run_v13_long_production_calibration.sh`. Do not report longer production
validation until a new job finishes with `required_failures=0`.

If you rerun the FW case, remove old `fw_tracking/` outputs first so that previously appended `tracking_interest_ids.pe*.ids` files do not contaminate the new TPHT handoff set.

#### Step 5. Confirm that rank-local interest-ID files were generated
After FW finishes, check that files like the following exist:

```bash
ls ./fw_tracking/tracking_interest_ids.pe*.ids
```

Each rank-local `.ids` file stores raw discovered `(dm_id, sd_id)` records in append order. Each append block begins with:
- `# TPHT_META PRC_NUM_X PRC_NUM_Y PRC_nprocs`

This metadata is later preserved by the merge utility and checked again by BW before it accepts the handoff file.

FW interest-ID output is evaluated every microphysics step so that transient coalescence events between regular SD dump times are also captured. Multiple records for the same pair can appear across output times by design; the later offline post-processing step removes duplicates and repartitions the unique pairs into BW-ready per-rank handoff files.

#### Step 6. Deduplicate the FW rank-local ID files and repartition them by `dm_id`
Run the dedicated post-processing utility after FW finishes:

```bash
cd ..
python merge_tracking_interest_ids.py \
  --input-glob "./ft_interest_id_baseline/fw_tracking/tracking_interest_ids.pe*.ids" \
  --output "./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids" \
  --rank-bucket-basename "./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_dedup"
```

This creates:
- `tracking_interest_ids_merged.ids` as a global deduplicated summary file
- `tracking_interest_ids_dedup.peXXXXXX.ids` as the BW-ready per-rank handoff files
- a preserved `# TPHT_META PRC_NUM_X PRC_NUM_Y PRC_nprocs` header at the top of every output `.ids` file

The post-processing step removes duplicate `(dm_id, sd_id)` pairs across both rank-local files and repeated output times while preserving the FW decomposition metadata required by BW runtime validation. It then repartitions the unique pairs by `dm_id` so that each BW rank reads only its own target subset.

#### Step 7. Check the BW reconstruction namelist
Still in the `tpht_test/` parent directory, view:

```bash
vi ./bt_interest_id_baseline/run.conf
```

Key BW TPHT settings are:

| Parameter | Recommended tutorial value | Purpose |
|---|---|---|
| `tracking_mode` | `2` | Enable BW reconstruction pass. |
| `SD_OUT_BASENAME` | `"./bw_output/superdroplet_restart"` | Keep BW restart output separate from FW under `tpht_test/`. |
| `tracking_id_input_basename` | `"../ft_interest_id_baseline/fw_tracking/tracking_interest_ids_dedup"` | Read BW targets from the deduplicated per-rank FW handoff basename. The runtime appends `.peXXXXXX.ids` per rank. |
| `sdm_dmpitvl` | `5.0D0` | Output `SD_selected_NetCDF_*` every 5 s for BW trajectory analysis. |
| `HISTORY_DEFAULT_TINTERVAL` | `5.0D0` | Keep Eulerian history enabled in BW and align its interval with `sdm_dmpitvl` for synchronized analysis. |

In other words, BW no longer starts from a newly sampled subset; it starts from the FW-discovered target set. In the updated tutorial, the BW job reads the deduplicated per-rank handoff basename directly, so each rank consumes only its own target bucket.

`tracking_interest_ids*.ids` is sufficient only for defining **which SDs belong to the BW target set**. It is not a full SD restart. The actual SD state evolution still comes from the same initial physical state used by FW, while the ID file only filters the tracked subset. BW must not read any FW discovery `superdroplet_restart` written at a later time, because those evolved SD IDs no longer represent the initial-time ID mapping that BW needs for reconstruction. The FW/BW `superdroplet_restart` files kept under `fw_output/` and `bw_output/` are output archives for restart/debugging, not the TPHT handoff definition itself.

#### Step 8. Prepare the BW executable and run the reconstruction pass
No rebuild is required if the BW case uses the same source tree and the same build options as the FW case, because the TPHT FW→BW handoff changes only the namelist input and the deduplicated tracking-ID files. Before submission, place the executable in the BW case directory by linking or copying the already built FW executable. A minimal example is:

```bash
cd ./bt_interest_id_baseline
mkdir -p bw_output
ln -fsv ../ft_interest_id_baseline/scale-rm* .
qsub UoH_run.pbs
```

or on SQUID:

```bash
cd ./bt_interest_id_baseline
mkdir -p bw_output
ln -fsv ../ft_interest_id_baseline/scale-rm* .
qsub squid_run.sh
```

If you changed the source code or the build flags after finishing the FW pass, rebuild once and refresh the executable link before launching BW.

#### Step 9. Verify the handoff and analyze results
The minimum practical validation checks are:
1. FW produced `fw_tracking/tracking_interest_ids.peXXXXXX.ids`.
2. The offline post-processing step produced both `fw_tracking/tracking_interest_ids_merged.ids` and `fw_tracking/tracking_interest_ids_dedup.peXXXXXX.ids`.
3. BW started successfully without falling back to random/stratified sampling.
4. The BW outputs contain only SD chains whose `(dm_id, sd_id)` appeared in the deduplicated per-rank FW handoff files.

A simple manual check is:

```bash
cd ..
wc -l ./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids

python - <<'PY'
from pathlib import Path
path = Path('./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids')
with path.open() as f:
    for _ in range(5):
        line = f.readline()
        if not line:
            break
        print(line.rstrip())
PY
```

An automated strict FW/BW consistency check is also available. It compares the deduplicated per-rank FW handoff set with the predecessor-ID pairs stored in the earliest BW `SD_selected_NetCDF_*` output group:

```bash
python verify_tpht_restart_consistency.py \
  --ids-basename "./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_dedup" \
  --bw-glob "./bt_interest_id_baseline/SD_selected_NetCDF_*"
```

The script reports:
- `BW_MATCH rank=...`: per-rank exact-set comparison against the deduplicated FW handoff bucket
- `fw_unique_pairs` / `bw_unique_pairs`: global unique `(dm_id, sd_id)` counts
- `RESULT: PASS/FAIL`: whether the earliest BW output reproduces the TPHT handoff exactly

An optional TPHT trajectory/collision-history summary can then be generated from the BW outputs:

```bash
python analyze_tpht_tracks.py \
  --bw-glob "./bt_interest_id_baseline/SD_selected_NetCDF_*" \
  --coal-glob "./bt_interest_id_baseline/SD_coal_output_NetCDF_*" \
  --output-dir "./tpht_analysis"
```

This writes:
- `tpht_time_summary.csv`: per-output-time tracked-count and spatial/radius summary
- `tpht_pair_timeseries.csv`: full TPHT pair-wise trajectory table keyed by `(pre_dmid, pre_sdid)`
- `tpht_pair_overview.csv`: per-pair lifetime, radius envelope, and collision-summary table
- `tpht_collision_summary.csv`: non-zero collision history summary, using `SD_coal_output_NetCDF_*` when available and `if_coal` flags otherwise
- `tpht_trajectory_samples.csv`: a compact sample of representative TPHT trajectories
- `tpht_analysis_summary.json`: machine-readable summary of the whole analysis run

For the strict consistency check above, the expected result is `RESULT: PASS`.

For trajectory reconstruction and visualization, you can then reuse the existing Python post-processing scripts in each case's `results/` directory.

Concrete files to edit/check for these settings:
- Main model switch file (all cases): `run.conf`
- Forward 3D baseline config: `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline/run.conf`
- Backward 3D baseline config: `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline/run.conf`
- Forward TPHT discovery config: `scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline/run.conf`
- Backward TPHT reconstruction config: `scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline/run.conf`
- Ice-phase FW TPHT smoke config: `scale-rm/test/case/shallowcloud/ice_tpht_smoke/fw_ice_smoke/run.conf`
- Ice-phase BW TPHT smoke config: `scale-rm/test/case/shallowcloud/ice_tpht_smoke/bw_ice_smoke/run.conf`
- TPHT merge utility: `scale-rm/test/case/shallowcloud/tpht_test/merge_tracking_interest_ids.py`
- TPHT consistency checker: `scale-rm/test/case/shallowcloud/tpht_test/check_tpht_consistency.py`
- Ice-phase TPHT smoke README: `scale-rm/test/case/shallowcloud/ice_tpht_smoke/README.md`
- 2D forward case config: `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward/run.conf`
- 2D backward case config: `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward/run.conf`
- 2D backward no-coalescence/no-UV config: `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward_no_coal_no_uv/run.conf`
- Trajectory post-processing script (forward): `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline/results/sd_output.py`
- Trajectory post-processing script (backward): `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline/results/sd_output.py`
- Generic smoke post-processing helper: `scale-rm/test/case/shallowcloud/tracking_postprocess/summarize_tracking_smoke.py`
- Representativeness evaluation (forward): `scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests/evaluate_representativeness.py`
- Representativeness evaluation (backward): `scale-rm/test/case/shallowcloud/backward_tracking_representativeness_tests/evaluate_representativeness.py`

The generic smoke helper is intended for quick local or CI validation after a
short FW/BW/TPHT run has produced outputs under one parent directory. With the
default temporary layout used by the examples, run:

```bash
cd scale-rm/test/case/shallowcloud/tracking_postprocess
python3 summarize_tracking_smoke.py \
  --smoke-root ../codex_smoke \
  --output-base ../codex_smoke/postprocess/tracking_smoke_summary
```

It writes `tracking_smoke_summary.{csv,json,md,tex}` and checks:
- FW selected-output files contain `sd_id/dm_id`.
- BW selected-output files contain `pre_sdid/pre_dmid`.
- TPHT BW selected-output IDs exactly reproduce the deduplicated FW handoff set.

For the rank-1 ice-phase smoke pair, pass the case paths explicitly:

```bash
cd scale-rm/test/case/shallowcloud
python3 tracking_postprocess/summarize_tracking_smoke.py \
  --fw-case ./ice_tpht_smoke/fw_ice_smoke \
  --bw-case ./ice_tpht_smoke/bw_ice_smoke \
  --tpht-fw-case ./ice_tpht_smoke/fw_ice_smoke \
  --tpht-bw-case ./ice_tpht_smoke/bw_ice_smoke \
  --tpht-fw-ids ./ice_tpht_smoke/fw_ice_smoke/fw_tracking/tracking_interest_ids_merged.ids \
  --output-base ./ice_tpht_smoke/tpht_analysis/tracking_smoke_summary
```

The helper prefers Python `netCDF4` when available and falls back to the system
`ncdump` executable, which is useful on lightweight local environments where
Python NetCDF bindings are not installed.

### How to configure “no coalescence output” in (d)/(f)
In the copied case directory:
1. Open `run.conf`.
2. Set:
   ```fortran
   coalescence_output_enable = .false.
   ```
3. Keep tracking/sampling settings unchanged unless your experiment requires otherwise.

### Build and run workflow (what each step does)
1. **Enter target case directory** (isolates all case-specific `run.conf`, job scripts, and outputs).
2. **Load compiler/library modules** (ensures ABI-compatible MPI/NetCDF/HDF5 toolchain).
3. **Run `make allclean` twice with SDM flags** (removes stale objects and stale dependency artifacts).
4. **Compile SCALE-RM with SDM flags** (builds executable for SDM-enabled physics).
5. **Link executable to case directory** (ensures local run script/job script can find binaries).
6. **Submit simulation job** (launches model integration on scheduler).

Rokko (University of Hyogo) example (`dycoms2_rf02_sdm_2D_forward`):
```bash
cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward
module purge
module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1
make allclean
make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
make SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
ln -fsv `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
qsub UoH_run.pbs
```

SQUID (Osaka University) example (`dycoms2_rf02_sdm_2D_forward`):
```bash
cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward
source /etc/profile.d/modules.sh
module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1
make allclean
make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
make SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
ln -fsv `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
qsub squid_run.sh
```

### Job scripts and environment/module settings (Rokko & SQUID)
- For each concrete case directory in (a)-(e), and for generated or retained ensemble case directories under (f):
  - Rokko job script: `UoH_run.pbs`
  - SQUID job script: `squid_run.sh`
- Rokko script module section (`UoH_run.pbs`) uses:
  - `module purge`
  - `module use --append /home/s274s011/modules`
  - `module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1`
- SQUID script module section (`squid_run.sh`) uses:
  - `source /etc/profile.d/modules.sh`
  - `module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1`

### How to disable Large-scale horizontal winds
When a case is marked `no_uv`, the large-scale horizontal wind must be set to zero consistently in both the initialization procedure and the forcing specification.

Required modifications are:
1. In `scale-rm/src/preprocess/mod_mkinit.f90` for the DYCOMS RF02 initialization, set:
   ```fortran
   velx(k,i,j) = 0.0_RP
   vely(k,i,j) = 0.0_RP
   ```
2. In `code/mod_user.f90` of the target `no_uv` case, within the `USER_LS_FLG == 1` branch, set:
   ```fortran
   U_GEOS(k) = 0.0_RP
   V_GEOS(k) = 0.0_RP
   ```
3. Rebuild the case executable after the source modification.

Notes:
- This modification should be applied only to cases explicitly marked `no_uv`.
- In this repository, the corresponding case is `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward_no_coal_no_uv`.

### Post-processing workflow
- **Python trajectory extraction (`results/sd_output.py`)**: reads SD NetCDF snapshots and reconstructs trajectory chains.
- **Python plotting (`results/traj_plot.py`)**: visualizes sampled trajectory paths for sanity checks.
- The old per-case `random_traj.py` and `particle_tracer_opt.f90` helper copies have been removed from the shallowcloud test tree. Use the maintained scripts under `scale-rm/test/case/shallowcloud/tpht_test/` for TPHT handoff checks and the per-case `results/` scripts for trajectory extraction and plotting.

2D backward no-coalescence post-processing example:
```bash
cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward_no_coal_no_uv/results
qsub run_py.pbs
python traj_plot.py
```

2D forward post-processing example:
```bash
cd scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward/results
qsub run_py.pbs
```

### Current-version namelist settings used by these cases
- `tracking_mode = 1` for forward cases and `tracking_mode = 2` for backward cases.
- `tracking_selection_mode`, `tracking_fraction`, `max_tracked_sds`, and stratified bounds (`tracking_height_*`, `tracking_radius_*`) are used instead of legacy `num_selected/height_min/radius_min`.
- `tracking_selection_mode="none"` means no additional random/stratified sampling mode is requested, while `tracking_fraction=1.0` disables any further downsampling.
- `tracking_id_output_basename` enables FW raw per-rank interest-ID output, and `tracking_id_input_basename` lets BW read a basename-expanded per-rank `.ids` handoff set directly. The recommended TPHT workflow is to deduplicate the FW records offline and then let each BW rank read its own deduplicated `.peXXXXXX.ids` file.
- `tracking_interest_radius_enable`, `tracking_interest_radius_threshold`, and `tracking_interest_coalescence_enable` define the warm/liquid TPHT interest filter, with the radius condition evaluated as `sd_r >= tracking_interest_radius_threshold`.
- `tracking_interest_ice_radius_enable`, `tracking_interest_ice_radius_threshold`, `tracking_interest_ice_phase_enable`, `tracking_interest_rime_mass_enable`, and `tracking_interest_rime_mass_threshold` extend the TPHT interest filter to cold-SDM fields when `sdm_cold = .true.`.
- If multiple warm or cold interest switches are enabled, the current implementation uses logical OR.
- `coalescence_output_enable` controls writing `SD_coal_output_NetCDF_*`, but it is forced to `.false.` when the microphysical coalescence process is disabled.
- Cold-SDM collision files use `SD_event_collision_NetCDF_*` with
  `trigger_code`, `target_reason_mask`, participant IDs, `event_multiplicity`,
  `phase_state*_pre/post`, and pre/post hydrometeor diagnostics. Warm mode
  keeps `SD_coal_output_NetCDF_*` and warm ordinary `if_coal`, while using
  `event_multiplicity` in the warm collision-event file.
- `random_perturbation_amp` is the user-facing amplitude, and `sdm_noise_amp` is its internal runtime copy.
- `sdm_dmpvar = 100` writes `SD_selected_NetCDF_*`, while `sdm_dmpvar = 200`
  writes `SD_selected_history.peXXXXXX`.
- `sdm_dmpvar = 010` writes `SD_all_NetCDF_*`, while `sdm_dmpvar = 020` writes
  `SD_all_history.peXXXXXX`.
- `sdm_dmpvar = 000` is the recommended FW setting for TPHT discovery cases that output the BW target set only through `tracking_id_output_basename`.
- Current shallowcloud examples now use `ATMOS_DYN_TYPE = "HEVI"`, `TIME_DT = 0.1D0`, `TIME_DT_ATMOS_DYN = 0.05D0`, other atmospheric physics time steps `= 0.1D0`, and `ATMOS_DYN_wdamp_tau = 0.5D0`.

## 11. Representativeness Tests
This section is placed after installation intentionally, because representativeness evaluation depends on successful build-and-run workflows.

### 11.1 Scientific objective
Representativeness tests quantify how well a sampled tracked subset reproduces the full SD population statistics and trajectory-chain properties. The target is to measure sampling-induced bias and uncertainty in both Eulerian-like moments and Lagrangian chain diagnostics.

### 11.2 Test design
Directories:
- `scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests`
- `scale-rm/test/case/shallowcloud/backward_tracking_representativeness_tests`

Case generation script:
- Each directory has its own `generate_cases.py`, which creates:
  - two reference full-tracking cases (`fraction = 1.0`) for `stratified` and `random`,
  - ensemble sampled cases for fractions `0.05`, `0.10`, `0.20`,
  - 10 seeds per (mode, fraction) group,
  - deterministic tracking-sampling seeding via `tracking_sampling_seed`.

### 11.3 Execution steps
1. Build binaries for target architecture.
2. Generate/refresh test cases (if needed) with `generate_cases.py`.
3. Submit all simulation jobs (`submit_all_uoh.sh` or `submit_all_squid.sh`).
4. Wait for all outputs (`SD_selected_NetCDF_*` or fallback `SD_all_NetCDF_*`) to complete.
5. Run post-processing in the direction-specific directory:
```bash
cd scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests
python evaluate_representativeness.py

cd ../backward_tracking_representativeness_tests
python evaluate_representativeness.py
```

### 11.4 Post-processing workflow
`evaluate_representativeness.py` performs:
1. Parse per-case metadata from `run.conf` (mode, fraction, binning, bounds, seed).
2. Read SD files grouped by output time and MPI rank.
3. Compute weighted microphysical metrics from `sd_r`, `sd_n`, `sd_z`, and `sd_asl` (moments, quantiles, concentration proxies).
4. Build 2D (`z`, `r`) weighted distributions and derive distribution-distance diagnostics.
5. Reconstruct backward chain-like statistics from identifier evolution and `if_coal`, including chain length and collision-chain fraction.
6. Export numerical summary to `representativeness_metrics.csv` and figures to `figures/`.

### 11.5 Interpretation guideline
- Compare sampled cases against full-tracking references within the same selection mode.
- Evaluate both scalar metrics and distribution metrics; agreement in one does not guarantee agreement in the other.
- Use multi-seed spread as an uncertainty estimate for sampling robustness.

## 12. Cold TPHT Current Base and Future Porting Plan

Cold TPHT extends the FW/BW/TPHT tracking workflow from warm liquid
coalescence histories to cold and mixed-phase SDM process histories. The
current v1.3+ implementation is intended as a research-development branch with
validated smoke, schema, restart, MPI, FW/BW/TPHT, and SQUID calibration
coverage. It should not be read as full real-science production validation.
When `sdm_cold = .false.`, the warm tracking boundary is preserved: legacy
ordinary `if_coal`, warm `SD_coal_output_NetCDF_*`, and collision
`event_multiplicity` behavior remain the compatibility path.

For cold mode, the main user-facing additions are:
- process event files with `trigger_code` as process identity and
  `trigger_level=1/2` for occurrence/significant records;
- interval masks and flags in ordinary, selected, and history output:
  `sd_event_mask`, `sd_event_sig_mask`, `sd_diag_mask`,
  `sd_phase_change_flag`, `sd_spatial_visit_flag`, and seven maintained
  interval maxima;
- cold event streams `SD_event_collision_NetCDF_*`,
  `SD_event_singleproc_NetCDF_*`, `SD_event_diag_NetCDF_*`, and the lightweight
  `SD_lifecycle_NetCDF_*` sidecar skeleton;
- optional output groups for ice geometry, rime morphology, single-process
  aerosol context, single-process thermodynamic context, and derived Kohler
  activation/deactivation context;
- Level-3 spatial-visit tracking with a Level-2 current-position fallback and a
  default-off rank-subdomain pruning option.

The implementation intentionally leaves several extension interfaces in place
without claiming that all future runtime paths are already active:
- process-specific occurrence controls (`tracking_evt_*_enable`) and
  significant controls (`tracking_sig_*_enable` plus thresholds);
- shared process bit mapping for `sd_event_mask` and `sd_event_sig_mask`;
- optional context switches:
  `TRACK_COLD_OUTPUT_ICE_GEOMETRY`,
  `TRACK_COLD_OUTPUT_RIME_MORPHOLOGY`,
  `TRACK_COLD_OUTPUT_AEROSOL_CONTEXT`,
  `TRACK_COLD_OUTPUT_THERMO_CONTEXT`,
  `TRACK_COLD_OUTPUT_KOHLER_CONTEXT`, and the legacy alias
  `TRACK_COLD_EVENT_EXTENDED_GEOMETRY`;
- lifecycle-safe ID predicates, a dynamic-ID namespace where
  `sd_id <= -1000` is valid and `-999` is the only invalid sentinel, an exact
  `(sd_id, dm_id)` lookup path, and restart storage for the dynamic-ID counter;
- a fixed-width `SD_lifecycle_NetCDF_*` schema for future SD creation, removal,
  split, domain-entry, global-halo, seeding, and `sdm_aslform`-like lifecycle
  records;
- test-only lifecycle harnesses for helper semantics, without activating the
  production `sdm_adjsdnum` path.

Before using this branch as a base for new science cases, the recommended
verification ladder is:
1. build the model and run the warm regression;
2. run cold schema smokes for ordinary, selected, history, collision,
   single-process, diagnostic, and lifecycle files;
3. cover process triggers and `trigger_level` combinations required by the
   study;
4. run restart-exact checks, including `sd_event_sig_mask` and spatial flags;
5. run MPI rank-2/rank-4 smokes and selected-ID copy checks;
6. run FW/BW/TPHT consistency checks;
7. run SQUID full validation if the target machine is SQUID;
8. run production I/O and threshold calibration for the actual science case.

The next porting step is to move the Cold TPHT v1.3+ patch set to the latest
SCALE-SDM code base instead of assuming this development branch already tracks
upstream. High-risk files and modules for that port include
`contrib/SDM/scale_atmos_phy_mp_sdm.F90`, the mirrored
`scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_sdm.F90`,
`sdm_common.f90`, `sdm_tracking_cold.f90`, `sdm_io.f90`, `sdm_boundary.f90`,
`sdm_idutil.f90`, `sdm_coalescence_cold.f90`, `sdm_condensation_water.f90`,
`sdm_meltfreeze.f90`, `sdm_subldep.f90`, `sdm_memmgr.f90`, build-system
dependencies, restart read/write blocks, validators, and smoke-case configs.
The port should be treated as a staged integration: first isolate the patch
series and namelist/build changes, then reconnect interval state arrays,
restart/MPI/selected-copy paths, event writers, physical process hooks,
spatial/lifecycle helpers, validators, and finally SQUID/production I/O
calibration.

Longer-term Cold TPHT development should focus on:
- collision aerosol and thermodynamic context after a separate collision-path
  variable audit;
- active global-halo, non-periodic boundary-entry, seeding, and `sdm_aslform`
  lifecycle hooks;
- production `sdm_adjsdnum`, `sdm_sdadd`, and `sdm_sdremove` runtime lifecycle
  validation;
- horizontal MPI source cleanup and collision invalidation cleanup, only after
  send-buffer and event-capture ordering are explicitly validated;
- primary ice nucleation source split, secondary ice production, breakup,
  fragmentation, and SIP process tracking;
- deliquescence, efflorescence, dry-core/water-film aerosol morphology, and
  expanded dry/wet aerosol phase-state handling;
- possible `sd_id` / `dm_id` migration to 64-bit storage if future dynamic-ID
  production approaches 32-bit limits.

Final reminder: the current Cold TPHT v1.3+ implementation in this repository
is developed on the SCALE 5.2.6 / SCALE-SDM 5.2.6-2.3.1 code base. A future
port should update the implementation to the latest SCALE-SDM version and rerun
the full validation ladder before treating it as an upstream-ready or
production-ready feature.

## Acknowledgements
I would like to sincerely thank my advisor, Prof. Shin-ichiro Shima, for his invaluable suggestions on the code and algorithms, his scientific and technical guidance, and his generous support in providing computational resources. I would also like to express my special gratitude to my Ph.D. supervisor, Prof. Chunsong Lu, for his mentorship and cultivation throughout my doctoral studies. I would also like to thank Mikito Toda for his generous support and informative discussions.

## Support and Community
Questions, issues, and discussions about SCALE-SDM can be directed here. Contributions and feedback are highly encouraged to enhance the model's capabilities and user experience. Please feel free to contact me: yinchongzhi@gmail.com. :grin:

## Reference
*Shima, S.-i., Kusano, K., Kawano, A., Sugiyama, T., and Kawahara, S.: The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model, Quarterly Journal of the Royal Meteorological Society, 135, 1307-1320, [https://doi.org/10.1002/qj.441](https://doi.org/10.1002/qj.441), 2009.*

*Nishizawa, S., Yashiro, H., Sato, Y., Miyamoto, Y., and Tomita, H.: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations, Geoscientific Model Development, 8, 3393-3419, [https://doi.org/10.5194/gmd-8-33932015](https://doi.org/10.5194/gmd-8-33932015), 2015.*

*Sato, Y., Nishizawa, S., Yashiro, H., Miyamoto, Y., Kajikawa, Y., and Tomita, H.: Impacts of cloud microphysics on trade wind cumulus: which cloud microphysics processes contribute to the diversity in a large eddy simulation?, Progress in Earth and Planetary Science, 2, 1-16, [https://doi.org/10.1186/s40645-015-0053-6](https://doi.org/10.1186/s40645-015-0053-6), 2015.*

*Shima, S., Sato, Y., Hashimoto, A., and Misumi, R.: Predicting the morphology of ice particles in deep convection using the super-droplet method: development and evaluation of SCALE-SDM 0.2.5-2.2.0, -2.2.1, and -2.2.2, Geosci. Model Dev., 13, 4107-4157, [https://doi.org/10.5194/gmd-13-4107-2020](https://doi.org/10.5194/gmd-13-4107-2020), 2020.*

*Yin, C., Shima, S., Xue, L., Lu, C., et al.: Simulation of marine stratocumulus using the super-droplet method: numerical convergence and comparison to a double-moment bulk scheme using SCALE-SDM 5.2.6-2.3.1, Geoscientific Model Development, 17, 5167-5189, [https://doi.org/10.5194/gmd-17-5167-2024](https://doi.org/10.5194/gmd-17-5167-2024), 2024.*

*Yin, C., Shima, Si., Lu, C. et al.: Resolving Entrainment–Mixing in Marine Stratocumulus: The Role of LES Grid Resolution and Super-Droplet Number, Advances in Atmospheric Sciences, 43, 845-860, [https://doi.org/10.1007/s00376-025-5043-z](https://doi.org/10.1007/s00376-025-5043-z), 2026.*
