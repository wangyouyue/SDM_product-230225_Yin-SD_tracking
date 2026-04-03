# Introduction to the Super-Droplet Method (SDM) in SCALE-SDM

## What is SDM?
The **Super-Droplet Method (SDM)**, originally introduced by **Shima et al. (2009)**, represents a significant advancement in Lagrangian cloud microphysics simulations. Unlike traditional Eulerian microphysics schemes, SDM employs **super-droplets (SDs)**, which are computational particles representing a large number of real droplets, aerosols, or precipitation particles with similar properties. This framework enables explicit treatment of key cloud-microphysical processes, including condensation, evaporation, and coalescence, thereby improving the physical interpretability of cloud formation and precipitation evolution.

## Integration with SCALE-SDM
The **SCALE (Scalable Computing for Advanced Library and Environment)** framework, developed with co-design by computational and computer science researchers, provides a robust, scalable platform for high-resolution atmospheric simulations (Nishizawa et al., 2015; Sato et al., 2015). SCALE-SDM integrates the SDM within this framework, offering a powerful tool for simulating cloud microphysics and dynamics. For more information on the SCALE platform, please visit the [SCALE's official website](http://scale.aics.riken.jp/).

This repository integrates SDM into SCALE version 5.2.6, leveraging both SDM’s microphysical precision and SCALE’s computational scalability to facilitate simulations of cloud systems that capture the complex interactions between microphysics and atmospheric dynamics. Further details about this version of SCALE can be explored in the [SCALE version 5.2.6 archives](https://scale.riken.jp/archives/5.2.6/).

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
- Direction-dependent tracking identifiers are incorporated into SD NetCDF outputs: FW writes `sd_id/dm_id`, BW writes `pre_sdid/pre_dmid`, and `if_coal` is written when coalescence-event output is enabled. Together these fields support trajectory-chain reconstruction in both directions during post-processing.
- Collision–coalescence event outputs are written through a unified interface and appended to time-bucketed NetCDF files (`SD_coal_output_NetCDF_*`).
- Random perturbation of SD motion is exposed through the namelist (`random_perturbation_enable`, `random_perturbation_amp`) and is physically inactive when the amplitude is zero.
- Representativeness test suites are provided for both FW and BW to quantify sampling-induced bias and chain-level fidelity.

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

Forward outputs are directly interpreted as current-time SD identities.

## 5. Backward Tracking Algorithm
Backward mode (`tracking_mode = 2`) uses the same runtime fields and output infrastructure, but the trajectory interpretation is reversed during analysis:
1. Initialize a target subset with the same sampling subsystem as forward mode; numerically, the initial BW identifiers are the same runtime `sd_id/dm_id` values that FW uses, but BW writes them to file as `pre_sdid/pre_dmid`.
2. At each saved time, use predecessor IDs (`pre_dmid`, `pre_sdid`) to connect an SD to earlier outputs.
3. In history-style BW output, after one dump is written, `sdm_assign_tracking_subset` updates the in-memory identifier pair of the currently retained tracked SDs to their current local `(dm_id, sd_id)` values for the next output window; therefore BW and FW are identical at initialization, but BW later rewrites predecessor links between successive outputs.
4. Use `if_coal` and coalescence event files to detect branch/merge points in ancestry chains.

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
- In this document and scripts, `YYYYMMDD` is the real simulation date token and is not fixed to `00000101`.

### 7.4 Supported/unsupported `sdm_dmpvar` modes in this tracking workflow
For the current SD-tracking workflow and provided test/tutorial cases, the following modes are supported and validated:
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
| `tracking_mode` | integer | `-1` | Unified tracking switch: `0` no tracking, `1` forward, `2` backward, `-1` compatibility mode using legacy logical switches. |
| `forward_tracking_enable` | logical | `.false.` | Legacy forward input. When `tracking_mode >= 0`, it is overwritten during configuration and then used internally by I/O, coalescence, and boundary routines. |
| `backward_tracking_enable` | logical | `.true.` | Legacy backward input. When `tracking_mode >= 0`, it is overwritten during configuration and then used internally by I/O, coalescence, and boundary routines. |
| `tracking_selection_mode` | string | `"random"` | Sampling algorithm: `"random"` or `"stratified"`. |
| `tracking_fraction` | real | `1.0` | Fraction of candidate SDs requested for tracking. |
| `max_tracked_sds` | integer | `0` | Hard cap for tracked SD count (`0` means unlimited). |
| `tracking_height_min` | real [m] | `400.0` | Lower height bound for stratified candidate filter. |
| `tracking_height_max` | real [m] | `800.0` | Upper height bound for stratified candidate filter. |
| `tracking_radius_min` | real [m] | `1.0e-6` | Lower radius bound for stratified candidate filter. |
| `tracking_radius_max` | real [m] | `0.0` | Upper radius bound; if `<= tracking_radius_min`, runtime candidate maximum radius is used. |
| `tracking_nz_bin` | integer | `8` | Number of height bins in stratified sampling. |
| `tracking_nr_bin` | integer | `10` | Number of radius bins in stratified sampling. |
| `tracking_min_per_bin` | integer | `1` | Requested minimum tracked SDs in each non-empty bin (subject to global budget). |
| `tracking_fallback_to_random` | logical | `.true.` | Fallback from stratified to random when stratified is not feasible. |
| `tracking_sample_initialized` | logical | `.false.` | Internal state indicating whether tracked subset is already initialized. |
| `coalescence_output_enable` | logical | `.true.` | Master switch for writing `SD_coal_output_NetCDF_*`. |
| `random_perturbation_enable` | logical | `.false.` | Master switch of SD motion perturbation. |
| `random_perturbation_amp` | real [m^1.5 s^-0.5] | `0.0` | User-facing perturbation amplitude; the configuration routine copies it to the internal runtime variable `sdm_noise_amp`. |

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
For more details, please see [the official user guide of SCALE](https://scale.riken.jp/archives/scale_users_guide_En.v5.2.6.pdf).

### Prerequisites
- **Compilers:** Fortran and C compilers are required.
- **Libraries:** MPI, NetCDF4, and HDF5 libraries must be installed.
- **Post-processing:** Python 3 with `numpy`, `netCDF4`, and `matplotlib`.

### Environment preparation
```bash
export SCALE_SYS=Linux64-intel-impi
```

### Test examples used in this repository (exact six categories)
The recommended test matrix is:
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

Concrete files to edit/check for these settings:
- Main model switch file (all cases): `run.conf`
- Forward 3D baseline config: `scale-rm/test/case/shallowcloud/forward_tracking_sampling_tests/ft_stratified_baseline/run.conf`
- Backward 3D baseline config: `scale-rm/test/case/shallowcloud/backward_tracking_sampling_tests/bt_stratified_baseline/run.conf`
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
- `forward_tracking_enable` and `backward_tracking_enable` remain in the namelist for legacy compatibility and as internal runtime logicals after `tracking_mode` is decoded; they are not an additional independent mode selector when `tracking_mode >= 0`.
- `tracking_selection_mode`, `tracking_fraction`, `max_tracked_sds`, and stratified bounds (`tracking_height_*`, `tracking_radius_*`) are used instead of legacy `num_selected/height_min/radius_min`.
- `coalescence_output_enable` controls writing `SD_coal_output_NetCDF_*`, but it is forced to `.false.` when the microphysical coalescence process is disabled.
- `random_perturbation_amp` is the user-facing amplitude, and `sdm_noise_amp` is its internal runtime copy.
- `sdm_dmpvar = 100/200` writes `SD_selected_NetCDF_*`; `sdm_dmpvar = 010/020` writes `SD_all_NetCDF_*`.

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

## Support and Community
Questions, issues, and discussions about SCALE-SDM can be directed here. Contributions and feedback are highly encouraged to enhance the model's capabilities and user experience. Please feel free to contact me: yinchongzhi@gmail.com. :grin:

## Reference
*Nishizawa, S., Yashiro, H., Sato, Y., Miyamoto, Y., and Tomita, H.: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations, Geoscientific Model Development, 8, 3393-3419, [https://doi.org/10.5194/gmd-8-33932015](https://doi.org/10.5194/gmd-8-33932015), 2015.*

*Sato, Y., Nishizawa, S., Yashiro, H., Miyamoto, Y., Kajikawa, Y., and Tomita, H.: Impacts of cloud microphysics on trade wind cumulus: which cloud microphysics processes contribute to the diversity in a large eddy simulation?, Progress in Earth and Planetary Science, 2, 1-16, [https://doi.org/10.1186/s40645-015-0053-6](https://doi.org/10.1186/s40645-015-0053-6), 2015.*

*Shima, S.-i., Kusano, K., Kawano, A., Sugiyama, T., and Kawahara, S.: The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model, Quarterly Journal of the Royal Meteorological Society, 135, 1307-1320, [https://doi.org/10.1002/qj.441](https://doi.org/10.1002/qj.441), 2009.*

*Yin, C., Shima, Si., Lu, C. et al.: Resolving Entrainment–Mixing in Marine Stratocumulus: The Role of LES Grid Resolution and Super-Droplet Number, Advances in Atmospheric Sciences, 43, 845-860, [https://doi.org/10.1007/s00376-025-5043-z](https://doi.org/10.1007/s00376-025-5043-z), 2026.*

*Yin, C.\*, Shima, S., Xue, L., Lu, C., et al.: Simulation of marine stratocumulus using the super-droplet method: numerical convergence and comparison to a double-moment bulk scheme using SCALE-SDM 5.2.6-2.3.1, Geoscientific Model Development, 17, 5167-5189, [https://doi.org/10.5194/gmd-17-5167-2024](https://doi.org/10.5194/gmd-17-5167-2024), 2024.*
