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

## 1. Overview of New Features (vs. SCALE-SDM)
Compared with the original SCALE-SDM branch, this merged branch introduces a unified, sampling-aware tracking framework with the following characteristics:
- **Forward tracking (FW)** and **backward tracking (BW)** are supported within the same executable and are controlled by a unified switch (`tracking_mode`) while preserving legacy compatibility.
- A common sampling subsystem is shared by FW and BW, with two selectable algorithms: `random` and `stratified`.
- A new **Two-Pass Hybrid Tracking (TPHT)** workflow is supported: FW can first scan the whole SD population, accumulate only the SD IDs that satisfy user-defined interest conditions, and BW can then reinitialize from this cumulative ID set instead of from a fresh random/stratified sample.
- Direction-dependent tracking identifiers are incorporated into SD NetCDF outputs: FW writes `sd_id/dm_id`, BW writes `pre_sdid/pre_dmid`, and `if_coal` is written when coalescence-event output is enabled. Together these fields support trajectory-chain reconstruction in both directions during post-processing.
- Collision–coalescence event outputs are written through a unified interface and appended to time-bucketed NetCDF files (`SD_coal_output_NetCDF_*`).
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
- radius threshold: `sd_r >= tracking_interest_radius_threshold`,
- recorded coalescence participation: `if_coal > 0`.

When `tracking_id_output_basename` is enabled together with at least one interest condition, FW no longer relies on `SD_selected_NetCDF_*` as the carrier of the target set. Instead, it appends raw `(dm_id, sd_id)` pairs to rank-local `.ids` files for later offline merge/dedup. Each append block also writes a `# TPHT_META PRC_NUM_X PRC_NUM_Y PRC_nprocs` header so that BW can validate MPI-decomposition compatibility before accepting the handoff file.

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
- Randomization rule: for each candidate, the code calls `random_number(rand_tracking)` and keeps the SD when `rand_tracking <= tracking_fraction`.
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
- **No tracking mode (`tracking_mode = 0`, or both FW/BW disabled)**: collision files contain only physical event fields (`event_time`, `sd_r1`, `sd_n1`, `sd_r2`, `sd_n2`, `num_col`) and omit tracking identifiers.
- The difference between FW and BW is restricted to identifier semantics (`current identity` versus `predecessor identity`); physical fields are shared.
- `num_col`: collision multiplicity/count used in SDM coalescence update.

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
- `sdm_dmpvar = 020/200`: use history-style NetCDF writer for all-SD or selected-SD group.

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
- Coalescence events: `SD_coal_output_NetCDF_YYYYMMDD-hhmmss.mmm.peXXXXXX`
- TPHT rank-local interest-ID output: `tracking_id_output_basename.peXXXXXX.ids`
- In this document and scripts, `YYYYMMDD` is the real simulation date token and is not fixed to `00000101`.

### 7.4 Supported/unsupported `sdm_dmpvar` modes in this tracking workflow
For the current SD-tracking workflow and provided test/tutorial cases, the following modes are supported and validated:
- `000`: no SD snapshot output; used in TPHT FW discovery cases where the BW target set is written through `tracking_id_output_basename`
- `010`: all-SD NetCDF output (`SD_all_NetCDF_*`)
- `100`: selected-SD NetCDF output (`SD_selected_NetCDF_*`)
- `020`: all-SD history NetCDF output (`SD_all_NetCDF_*` with suffixed history variables)
- `200`: selected-SD history NetCDF output (`SD_selected_NetCDF_*` with suffixed history variables)

The file basename is determined by the writer tag (`all` or `selected`), not by whether the writer is normal or history-style. Therefore, `010` and `020` both produce `SD_all_NetCDF_*`, whereas `100` and `200` both produce `SD_selected_NetCDF_*`; the difference is the internal variable layout.

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
| `coalescence_output_enable` | logical | `.true.` | Master switch for writing `SD_coal_output_NetCDF_*`. |
| `random_perturbation_enable` | logical | `.false.` | Master switch of SD motion perturbation. |
| `random_perturbation_amp` | real [m^1.5 s^-0.5] | `0.0` | User-facing perturbation amplitude; the configuration routine copies it to the internal runtime variable `sdm_noise_amp`. |

`tracking_sample_initialized` is an internal runtime state variable rather than a user-facing namelist control, so it is intentionally omitted from the table above.

### 8.1 Namelist notes specific to TPHT
- **FW discovery pass:** set `tracking_mode = 1`, set `tracking_id_output_basename`, enable one or both interest-condition switches, and typically use `sdm_dmpvar = 000`.
- **BW reconstruction pass:** set `tracking_mode = 2` and point `tracking_id_input_basename` to the FW deduplicated per-rank ID basename.
- **Priority rule:** `tracking_id_input_basename` takes precedence over `tracking_fraction`, `tracking_selection_mode`, and the stratified bounds during BW initialization.
- **Common recommendation:** in TPHT FW cases, use `tracking_selection_mode = "none"`, `tracking_fraction = 1.0`, and `max_tracked_sds = 0` so that no extra random/stratified downsampling is introduced before the interest-condition reduction step.
- **MPI decomposition must match between FW and BW:** `PRC_NUM_X`, `PRC_NUM_Y`, and total MPI rank count must match the FW run that produced the handoff IDs. BW now checks the `TPHT_META` header in the input `.ids` file and aborts on mismatch.
- **Interest-condition combination:** `tracking_interest_radius_enable` and `tracking_interest_coalescence_enable` can each be used independently. If both are enabled, the current implementation uses logical OR; logical AND is not implemented in the present code.

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
- **Post-processing:** Python 3 with `numpy`, `netCDF4`, and `matplotlib`.

### Environment preparation
```bash
export SCALE_SYS=Linux64-intel-impi
```

### Test examples used in this repository
The recommended baseline test matrix is:
- **(a) 2D backward, with coalescence output**  
  `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward`
- **(b) 2D backward, no coalescence output**  
  `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward_no_coal_no_uv`
- **(c) 3D backward, with sampling and coalescence output**  
  `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline`
- **(d) 3D backward, with sampling and no coalescence output**  
  Copy case (c), then set `coalescence_output_enable = .false.` in `run.conf`.
- **(e) 3D forward, with sampling and coalescence output**  
  `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline`
- **(f) 3D forward, with sampling and no coalescence output**  
  Copy case (e), then set `coalescence_output_enable = .false.` in `run.conf`.

Additional TPHT tutorial pair introduced in this update:
- **FW TPHT discovery case**  
  `scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline`
- **BW TPHT reconstruction case**  
  `scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline`

### 10.1 Practical tutorial: Two-Pass Hybrid Tracking (TPHT)
This subsection gives a concrete, end-to-end workflow for running TPHT with the shallowcloud test cases added in this update.

#### Step 1. Use the dedicated tutorial cases
The tutorial pair is:
- FW discovery case: `scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline`
- BW reconstruction case: `scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline`

These two tutorial cases are placed under the same `tpht_test/` parent directory so that the FW-to-BW handoff path stays short and FW/BW restart outputs can be isolated with `fw_output/` and `bw_output/`.

#### Step 2. Check the FW discovery namelist
View:
- `vi scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline/run.conf`

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
| `HISTORY_DEFAULT_TINTERVAL` | `0.0D0` | FW discovery can disable history output because the TPHT handoff artifact is `tracking_interest_ids.pe*.ids`, not the Eulerian history stream. |

Notes:
- `tracking_mode` already drives FW/BW enablement, so `forward_tracking_enable` / `backward_tracking_enable` are not user-facing namelist keys in the current TPHT tutorial cases; they remain only as internal runtime flags derived from `tracking_mode`.
- `tracking_selection_mode="none"` means no additional random/stratified sampling mode is requested.
- `tracking_fraction=1.0d0` means the post-initialization downsampling step is disabled. In the TPHT FW tutorial, that leaves the configured height/radius window as the effective initialization filter.
- For the shallowcloud tracking cases in this repository, `ATMOS_DYN_TYPE` is set to `HEVE`, the sponge layer starts at `ATMOS_DYN_wdamp_height = 1400.D0`, and `ATMOS_DYN_wdamp_tau = 0.2D0` matches `10 * TIME_DT_ATMOS_DYN` because `TIME_DT_ATMOS_DYN = 0.02D0`.
- `ATMOS_PHY_TB_SMG_horizontal` remains `.false.` here because these tutorial cases use `ATMOS_PHY_TB_TYPE = 'SMAGORINSKY'`, not the HYBRID gray-zone/RANS setting that requires horizontal LES viscosity.

#### Step 3. Build the executable
From the FW TPHT case directory, rebuild SCALE-RM with SDM enabled. A generic example is:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline
module purge
module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1
make allclean
time make -j SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
ln -fsv `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
```

#### Step 4. Run the FW discovery pass
Enter the FW TPHT case directory, make sure the FW output directories exist, and submit the job:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline
mkdir -p fw_output fw_tracking
qsub UoH_run.pbs
```

or on SQUID:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline
mkdir -p fw_output fw_tracking
qsub squid_run.sh
```

If you rerun the FW case, remove old `fw_tracking/` outputs first so that previously appended `tracking_interest_ids.pe*.ids` files do not contaminate the new TPHT handoff set.

#### Step 5. Confirm that rank-local interest-ID files were generated
After FW finishes, check that files like the following exist:

```bash
ls scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline/fw_tracking/tracking_interest_ids.pe*.ids
```

Each rank-local `.ids` file stores raw discovered `(dm_id, sd_id)` records in append order. Each append block begins with:
- `# TPHT_META PRC_NUM_X PRC_NUM_Y PRC_nprocs`

This metadata is later preserved by the merge utility and checked again by BW before it accepts the handoff file.

FW interest-ID output is evaluated every microphysics step so that transient coalescence events between regular SD dump times are also captured. Multiple records for the same pair can appear across output times by design; the later offline post-processing step removes duplicates and repartitions the unique pairs into BW-ready per-rank handoff files.

#### Step 6. Deduplicate the FW rank-local ID files and repartition them by `dm_id`
Run the dedicated post-processing utility after FW finishes:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test
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
View:
- `vi scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline/run.conf`

Key BW TPHT settings are:

| Parameter | Recommended tutorial value | Purpose |
|---|---|---|
| `tracking_mode` | `2` | Enable BW reconstruction pass. |
| `SD_OUT_BASENAME` | `"./bw_output/superdroplet_restart"` | Keep BW restart output separate from FW under `tpht_test/`. |
| `tracking_id_input_basename` | `"../ft_interest_id_baseline/fw_tracking/tracking_interest_ids_dedup"` | Read BW targets from the deduplicated per-rank FW handoff basename. The runtime appends `.peXXXXXX.ids` per rank. |
| `sdm_dmpitvl` | `5.0D0` | Output `SD_selected_NetCDF_*` every 5 s for BW trajectory analysis. |
| `HISTORY_DEFAULT_TINTERVAL` | `5.0D0` | Keep Eulerian history enabled in BW and align its interval with `sdm_dmpitvl` for synchronized analysis. |

In other words, BW no longer starts from a newly sampled subset; it starts from the FW-discovered target set. In the updated tutorial, the BW job reads the deduplicated per-rank handoff basename directly, so each rank consumes only its own target bucket.

`tracking_interest_ids*.ids` is sufficient only for defining **which SDs belong to the BW target set**. It is not a full SD restart. The actual SD state evolution still comes from the normal model initialization/restart inputs, while the ID file only filters the tracked subset. In the current TPHT tutorial, BW reruns from the same atmospheric/SD initial state as FW and reuses the FW-discovered ID list to initialize backward tracking; it does not read FW `superdroplet_restart` as an input requirement. The FW/BW `superdroplet_restart` files kept under `fw_output/` and `bw_output/` are output archives for restart/debugging, not the TPHT handoff definition itself.

#### Step 8. Prepare the BW executable and run the reconstruction pass
No rebuild is required if the BW case uses the same source tree and the same build options as the FW case, because the TPHT FW→BW handoff changes only the namelist input and the deduplicated tracking-ID files. Before submission, place the executable in the BW case directory by linking or copying the already built FW executable. A minimal example is:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline
mkdir -p bw_output
ln -fsv ../ft_interest_id_baseline/scale-rm* .
qsub UoH_run.pbs
```

or on SQUID:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline
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
wc -l scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids

python - <<'PY'
from pathlib import Path
path = Path('scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids')
with path.open() as f:
    for _ in range(5):
        line = f.readline()
        if not line:
            break
        print(line.rstrip())
PY
```

An automated FW/BW consistency check is also available. It compares the merged FW handoff set with the union of predecessor IDs stored in the earliest BW `SD_selected_NetCDF_*` output group:

```bash
cd scale-rm/test/case/shallowcloud/tpht_test
python check_tpht_consistency.py \
  --fw-ids "./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids" \
  --bw-glob "./bt_interest_id_baseline/bw_output/SD_selected_NetCDF_*.pe*.nc"
```

The script reports:
- `fw_unique_pairs`: unique `(dm_id, sd_id)` pairs in the merged FW handoff file
- `bw_valid_records`: total valid BW predecessor records found in the earliest BW output-time group
- `bw_unique_pairs`: unique BW predecessor pairs after rank-wise union
- `missing_in_bw` / `extra_in_bw`: set differences between FW and BW

The expected result is `missing_in_bw=0` and `extra_in_bw=0`.

For trajectory reconstruction and visualization, you can then reuse the existing Python post-processing scripts in each case's `results/` directory.

Concrete files to edit/check for these settings:
- Main model switch file (all cases): `run.conf`
- Forward 3D baseline config: `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline/run.conf`
- Backward 3D baseline config: `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline/run.conf`
- Forward TPHT discovery config: `scale-rm/test/case/shallowcloud/tpht_test/ft_interest_id_baseline/run.conf`
- Backward TPHT reconstruction config: `scale-rm/test/case/shallowcloud/tpht_test/bt_interest_id_baseline/run.conf`
- TPHT merge utility: `scale-rm/test/case/shallowcloud/tpht_test/merge_tracking_interest_ids.py`
- TPHT consistency checker: `scale-rm/test/case/shallowcloud/tpht_test/check_tpht_consistency.py`
- 2D forward case config: `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_forward/run.conf`
- 2D backward case config: `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward/run.conf`
- 2D backward no-coalescence/no-UV config: `scale-rm/test/case/shallowcloud/dycoms2_rf02_sdm_2D_backward_no_coal_no_uv/run.conf`
- Trajectory post-processing script (forward): `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline/results/sd_output.py`
- Trajectory post-processing script (backward): `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline/results/sd_output.py`
- Representativeness evaluation (forward): `scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests/evaluate_representativeness.py`
- Representativeness evaluation (backward): `scale-rm/test/case/shallowcloud/backward_tracking_representativeness_tests/evaluate_representativeness.py`

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
- For each case directory in (a)–(f):
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
- **Python plotting (`results/traj_plot.py` or `random_traj.py`)**: visualizes sampled trajectory paths for sanity checks.
- **Fortran tracer (`particle_tracer_opt.f90`, available in 2D cases)**: fast binary NetCDF traversal for large trajectory datasets.

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
- `tracking_interest_radius_enable`, `tracking_interest_radius_threshold`, and `tracking_interest_coalescence_enable` define the TPHT interest filter, with the radius condition evaluated as `sd_r >= tracking_interest_radius_threshold`; if both switches are enabled, the current implementation uses logical OR.
- `coalescence_output_enable` controls writing `SD_coal_output_NetCDF_*`, but it is forced to `.false.` when the microphysical coalescence process is disabled.
- `random_perturbation_amp` is the user-facing amplitude, and `sdm_noise_amp` is its internal runtime copy.
- `sdm_dmpvar = 100/200` writes `SD_selected_NetCDF_*`; `sdm_dmpvar = 010/020` writes `SD_all_NetCDF_*`.
- `sdm_dmpvar = 000` is the recommended FW setting for TPHT discovery cases that output the BW target set only through `tracking_id_output_basename`.
- Current shallowcloud examples now use `ATMOS_DYN_TYPE = "HEVI"`, `TIME_DT = 0.1D0`, `TIME_DT_ATMOS_DYN = 0.05D0`, other atmospheric physics time steps `= 0.1D0`, and `ATMOS_DYN_wdamp_tau = 0.5D0`.

## 11. Representativeness Tests
This section is placed after installation intentionally, because representativeness evaluation depends on successful build-and-run workflows.

### 11.1 Scientific objective
Representativeness tests quantify how well a sampled tracked subset reproduces the full SD population statistics and trajectory-chain properties. The target is to measure sampling-induced bias and uncertainty in both Eulerian-like moments and Lagrangian chain diagnostics.

### 11.2 Test design
Directory:
- `scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests`

Case generation script:
- `generate_cases.py` creates:
  - two reference full-tracking cases (`fraction = 1.0`) for `stratified` and `random`,
  - ensemble sampled cases for fractions `0.05`, `0.10`, `0.20`,
  - 10 seeds per (mode, fraction) group,
  - deterministic namelist seeding via `RANDOM_SEED_SCALE`.

### 11.3 Execution steps
1. Build binaries for target architecture.
2. Generate/refresh test cases (if needed) with `generate_cases.py`.
3. Submit all simulation jobs (`submit_all_uoh.sh` or `submit_all_squid.sh`).
4. Wait for all outputs (`SD_selected_NetCDF_*` or fallback `SD_all_NetCDF_*`) to complete.
5. Run post-processing:
```bash
cd scale-rm/test/case/shallowcloud/forward_tracking_representativeness_tests
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

## Acknowledgements
I would like to sincerely thank my advisor, Prof. Shin-ichiro Shima, for his invaluable suggestions on the code and algorithms, his scientific and technical guidance, and his generous support in providing computational resources. I would also like to express my special gratitude to my Ph.D. supervisor, Prof. Chunsong Lu, for his mentorship and cultivation throughout my doctoral studies. I would also like to thank Mikito Toda for his generous support and informative discussions.

## Support and Community
Questions, issues, and discussions about SCALE-SDM can be directed here. Contributions and feedback are highly encouraged to enhance the model's capabilities and user experience. Please feel free to contact me: yinchongzhi@gmail.com. :grin:

## Reference
*Nishizawa, S., Yashiro, H., Sato, Y., Miyamoto, Y., and Tomita, H.: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations, Geoscientific Model Development, 8, 3393-3419, [https://doi.org/10.5194/gmd-8-33932015](https://doi.org/10.5194/gmd-8-33932015), 2015.*

*Sato, Y., Nishizawa, S., Yashiro, H., Miyamoto, Y., Kajikawa, Y., and Tomita, H.: Impacts of cloud microphysics on trade wind cumulus: which cloud microphysics processes contribute to the diversity in a large eddy simulation?, Progress in Earth and Planetary Science, 2, 1-16, [https://doi.org/10.1186/s40645-015-0053-6](https://doi.org/10.1186/s40645-015-0053-6), 2015.*

*Shima, S.-i., Kusano, K., Kawano, A., Sugiyama, T., and Kawahara, S.: The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model, Quarterly Journal of the Royal Meteorological Society, 135, 1307-1320, [https://doi.org/10.1002/qj.441](https://doi.org/10.1002/qj.441), 2009.*

*Yin, C., Shima, Si., Lu, C. et al.: Resolving Entrainment–Mixing in Marine Stratocumulus: The Role of LES Grid Resolution and Super-Droplet Number, Advances in Atmospheric Sciences, 43, 845-860, [https://doi.org/10.1007/s00376-025-5043-z](https://doi.org/10.1007/s00376-025-5043-z), 2026.*

*Yin, C.\*, Shima, S., Xue, L., Lu, C., et al.: Simulation of marine stratocumulus using the super-droplet method: numerical convergence and comparison to a double-moment bulk scheme using SCALE-SDM 5.2.6-2.3.1, Geoscientific Model Development, 17, 5167-5189, [https://doi.org/10.5194/gmd-17-5167-2024](https://doi.org/10.5194/gmd-17-5167-2024), 2024.*
