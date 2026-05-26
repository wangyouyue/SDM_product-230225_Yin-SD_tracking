# GMD2026 Shallow-Cloud Benchmark and TPHT Experiment Suite

This suite supports the revised GMD manuscript response with controlled DYCOMS-II RF02 SCALE-SDM computational evidence and a practical 3D TPHT workflow for interest-restricted backward reconstruction.

## Groups

Required groups are `00_base_restart_3d`, `01_bench_3d_samp_30min`, `02_tpht_3d_interest_70min`, and `03_fw_rep_2d_600s`. Optional groups are `04_sdnc_scaling_lite_30min` and `05_outint_io_lite_10min`.

`03_fw_rep_2d_600s` is a 600 s sampling-procedure verification test. It verifies the sampling implementation only and is not a proof of long-time or 3D spatial representativeness.

## Build and Makefile Check

From `scale-rm/test/case/shallowcloud/GMD2026`:

```bash
bash common/check_makefiles.sh
bash run_all_static_audits.sh
```

From any model run directory:

```bash
../../common/build_squid.sh
```

The local `Makefile` keeps `TOPDIR` as the SCALE repository root required by `scale-rm/test/Makefile.common`, and also defines `SCALE_RM_DIR := $(abspath ../../../../../..)` for the GMD2026 two-level case layout.

`check_makefiles.sh` and `build_squid.sh` load the SQUID build modules internally and set `SCALE_SYS=Linux64-intel-impi` when it is not already set. `run_all_static_audits.sh` loads the SQUID Python module through `common/load_basepy_2026_quiet.sh` before running Python audit scripts. Before the first build, `check_makefiles.sh` may report that `${TOPDIR}/bin` is not present; this is informational and is expected until `build_squid.sh` creates the SCALE executables.

`run_all_static_audits.sh` additionally checks the expected GMD2026 namelist settings, physical restart chains, and benchmark-diagnostic availability. It writes `config_audit.*`, `restart_audit.*`, and `benchmark_diagnostics_availability.*` in this directory.

`build_squid.sh` copies the exact executables `scale-rm` and `scale-rm_init` into the current case directory. It intentionally avoids symbolic links and `scale-rm*` wildcards.

## Submit

```bash
cd 00_base_restart_3d && bash submit_base_restart_squid.sh
cd ../01_bench_3d_samp_30min && bash submit_bench_3d_squid.sh
cd ../02_tpht_3d_interest_70min && bash submit_tpht_squid.sh
cd ../03_fw_rep_2d_600s && bash submit_fw_rep_2d_squid.sh
cd ../04_sdnc_scaling_lite_30min && bash submit_sdnc_scaling_squid.sh
cd ../05_outint_io_lite_10min && bash submit_outint_io_squid.sh
```

Batch scripts use `qsub --after` so that post-processing is submitted after the model jobs.
On the current SQUID queue, the suite uses a single-dependency job chain rather than colon-separated dependency lists, because this environment rejects multi-ID `--after` strings.
`common/submit_utils.sh` tries both SQUID dependency forms, `qsub --after <jobid> <script>` and `qsub <script> --after <jobid>`, to handle site-specific NQSV argument parsing.
All post-processing scripts submitted through `qsub --after` include the SQUID `#PBS --group=hp250136` header, so model jobs and post-processing jobs follow the same group-account requirement.
Every model job writes a per-case status JSON into `GMD2026/job_status/` and appends one line to `GMD2026/job_status/job_events.tsv`. The batch scripts also submit `common/serial_collect_failures.sh`, which collects those files into `job_failure_summary.csv`, `job_failure_summary.md`, and `job_failure_summary.json`.

`common/squid_model_job.sh` creates parent directories declared by output basenames in `init.conf` and `run.conf`, including `restart_output/`, `fw_output/`, `bw_output/`, History basenames, and TPHT ID output directories. It also checks every required `tracking_id_input_basename.peXXXXXX.ids` file before BW reconstruction starts, so a failed TPHT merge is reported as a clear preflight failure instead of an MPI abort inside SDM initialization.

`04_sdnc_scaling_lite_30min/submit_sdnc_scaling_squid.sh` submits the SDNC = 10, 20, 40, and 80 NT/FW/BW cases as a dependency chain, then runs `postprocess/serial_postprocess.sh`. That post-processing job runs `collect_sdnc_scaling.py` and `plot_sdnc_scaling.py`, then the common failure-summary collector.

`05_outint_io_lite_10min/submit_outint_io_squid.sh` submits `nt_nolog_10min` plus the 30, 60, and 120 s FW/BW output-interval cases as a dependency chain, then runs `postprocess/serial_postprocess.sh`. That post-processing job runs `collect_outint_io.py` and `plot_outint_io.py`, then the common failure-summary collector.

## Restart Policy

`00_base_restart_3d/base_restart_3d` is retained as a common 3D spin-up/restart example, but the current benchmark groups run from the initial DYCOMS-II RF02 state instead of reading that 60 min restart. This avoids inconsistent tracking state when the base restart was generated with `tracking_mode = 0`.

`01_bench_3d_samp_30min` runs 0-30 min from cold start. `02_tpht_3d_interest_70min` runs FW discovery and BW reconstruction from the same initial physical state for 0-70 min; BW additionally reads deduplicated `.ids` from `fw_discovery/fw_tracking/` and must not read any FW discovery `fw_output/superdroplet_restart` written after the TPHT initial time. `03_fw_rep_2d_600s`, `04_sdnc_scaling_lite_30min`, and `05_outint_io_lite_10min` also run from cold start.

Do not rerun TPHT `bw_reconstruction` until `02_tpht_3d_interest_70min/fw_discovery/fw_tracking/tracking_interest_ids_dedup.pe000000.ids` and the remaining per-rank deduplicated ID files exist. If they are missing, rerun `postprocess/serial_merge_ids.sh` first and inspect `postprocess/logs/merge.log` plus `gmd02_merge.o*` / `gmd02_merge.e*`.

Within each controlled group, SQUID scripts use shared initialization where it is scientifically valid:

- `01_bench_3d_samp_30min`: `nt_nolog` runs `scale-rm_init`; `nt_coallog`, `fw005_coallog`, and `bw005_coallog` copy its `init_dycom_rf02_*` and `random_number_init.pe*`.
- `02_tpht_3d_interest_70min`: `fw_discovery` provides the initial physical/random state; `bw_reconstruction` copies that initial state and reads only the deduplicated TPHT `.ids` as the target handoff.
- `03_fw_rep_2d_600s`: `ref_full_stratified` runs initialization once, and all other cases copy that initial state. Seed-to-seed variation is controlled only by `&PARAM_ATMOS_PHY_MP_SDM tracking_sampling_seed` in `run.conf`; this seed drives tracking-sampling RNG state and does not change SDM initialization or coalescence RNG streams.
- `04_sdnc_scaling_lite_30min`: only cases with the same SDNC share init; for example, `sdnc20_fw005` and `sdnc20_bw005` copy from `sdnc20_nt`.
- `05_outint_io_lite_10min`: `nt_nolog_10min` provides the common initial state for all output-interval FW/BW cases.

The default SQUID layout uses 64 MPI ranks for 3D cases (`PRC_NUM_X = 8`, `PRC_NUM_Y = 8`, local `IMAX = JMAX = 5`) and 20 MPI ranks for quasi-2D representativeness cases (`PRC_NUM_X = 20`, `PRC_NUM_Y = 1`). The 3D choice preserves the global 40 x 40 x 300 grid while using most of one 76-core SQUID node without forcing inefficient 16-rank runs.

## Post-processing

Benchmark:

```bash
cd 01_bench_3d_samp_30min
python postprocess/collect_bench_3d_table.py
python postprocess/plot_bench_3d_overhead.py
```

TPHT:

```bash
cd 02_tpht_3d_interest_70min
python postprocess/check_tpht_consistency.py
python postprocess/analyze_tpht_summary.py
python postprocess/plot_tpht_targets.py
```

Representativeness:

```bash
cd 03_fw_rep_2d_600s
python postprocess/evaluate_representativeness.py
python postprocess/plot_representativeness.py
```

SDNC scaling:

```bash
cd 04_sdnc_scaling_lite_30min
python postprocess/collect_sdnc_scaling.py
python postprocess/plot_sdnc_scaling.py
```

This writes `sdnc_scaling_30min_summary.csv`, `sdnc_scaling_30min_summary.md`, `sdnc_scaling_30min_summary.tex`, `sdnc_scaling_30min_summary.json`, and `sdnc_scaling_30min.{pdf,svg,png}`.

Output-interval I/O scaling:

```bash
cd 05_outint_io_lite_10min
python postprocess/collect_outint_io.py
python postprocess/plot_outint_io.py
```

This writes `outint_io_10min_summary.csv`, `outint_io_10min_summary.md`, `outint_io_10min_summary.tex`, `outint_io_10min_summary.json`, and `outint_io_10min.{pdf,svg,png}`.

Plot scripts save `.pdf`, `.svg`, and 600 dpi `.png` using `common/nature_style.py`.

## TPHT Consistency

`missing_in_bw` counts FW target `(dm_id, sd_id)` pairs that are not found in the earliest BW selected output. `extra_in_bw` counts BW selected pairs that were not in the FW handoff. A strict TPHT handoff should give `missing_in_bw = 0` and `extra_in_bw = 0`. This check is meaningful only when BW starts from the same TPHT initial physical state as FW; using a later FW restart changes the ID/state mapping and invalidates the comparison.

## Manuscript Interpretation

The revised manuscript does not present full-population 3D backward tracking as the recommended practical configuration. Instead, it provides a controlled 3D benchmark for practical sampled FW/BW tracking and a TPHT experiment for interest-restricted 3D backward reconstruction. This addresses the storage limitation of brute-force full-lineage backward tracking while preserving the O(1) direct-lookup advantage for reconstructed target chains.

TPHT FW discovery uses `HISTORY_DEFAULT_TINTERVAL = 0.0D0` and `sdm_dmpvar = 000` so that Eulerian history and `SD_selected_NetCDF_*` output are suppressed. The TPHT handoff is the `.ids` stream written through `tracking_id_output_basename`.

TPHT BW reconstruction writes Eulerian History output every 10 s, matching `sdm_dmpitvl = 10.0D0`, so reconstructed target trajectories can be interpreted together with the meteorological fields. `HISTORY_OUTPUT_STEP0 = .true.` is used for History-output cases; History-disabled cases keep `HISTORY_OUTPUT_STEP0 = .false.`.

## Benchmark Diagnostics

The SQUID job scripts provide `job_metrics.json`, and the post-processing scripts can compute file sizes and coalescence-event counts after model output exists. Existing SCALE-SDM logs also expose tracking-chain counts and tracking-ID memory estimates when tracking output diagnostics are written.

Rank-level peak memory uses binary MiB fields (`peak_memory_rank_max_mib`, `peak_memory_rank_sum_mib`) from `/proc/self/status` `VmHWM`. I/O timing is split into total and mean fields, for example `sd_output_write_time_total_s` and `sd_output_write_time_mean_s`; the hook smoke-test checker verifies that last-write timing is not used as total timing.

The Fortran hooks are compiled into the source but remain disabled unless a case sets:

```fortran
gmd_benchmark_diag_enable = .true.
```

This keeps ordinary production runs clean while allowing GMD2026 cases to opt in explicitly. Groups `00`-`05` set this switch in `run.conf`, and their table scripts parse the same `GMD_BENCH_DIAG` and `GMD_IO_DIAG` fields that are smoke-tested in group `06`.

For `04_sdnc_scaling_lite_30min`, NT cases intentionally set `tracking_selection_mode = "none"` and do not write SD output. They are compute baselines for SDNC scaling, not I/O-equivalent output baselines. FW/BW cases use `tracking_selection_mode = "stratified"` with `tracking_fraction = 0.05D0`.
