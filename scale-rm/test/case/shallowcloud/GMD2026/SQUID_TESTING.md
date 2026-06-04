# SQUID Testing Guide for the GMD2026 Suite

This guide assumes the archive was unpacked on SQUID and the current directory is the repository root.

## 1. Static checks

```bash
cd scale-rm/test/case/shallowcloud/GMD2026
bash common/check_makefiles.sh
bash run_all_static_audits.sh
```

`check_makefiles.sh` sources `/etc/profile.d/modules.sh`, loads the SQUID build modules, and sets `SCALE_SYS=Linux64-intel-impi` automatically. `run_all_static_audits.sh` sources `common/load_basepy_2026_quiet.sh`, so it uses `BasePy/2026` or `${HOME}/venvs/gmd2026-netcdf4/bin/python` when that optional venv exists. If `${TOPDIR}/bin` is absent before the first build, the script prints an informational message only; it is not a failed check.

Expected result:

- `check_makefiles.sh` must not fail.
- `run_all_static_audits.sh` must end with `GMD2026 static audit completed successfully.`
- `FAIL: 0 cases` is required for both config and restart audits.
- `restart_audit.md` includes the SQUID shared-init policy for each case.

The current benchmark and TPHT groups intentionally run from the initial DYCOMS-II RF02 state. Therefore first-hour supersaturation cap evidence is not a restart-time warning for those groups; it only documents that the cold-start runs include the first-hour regime.

## 1.1 Python environment

Post-processing job scripts load `BasePy/2026` through `common/load_basepy_2026_quiet.sh`. If you need NetCDF variable checks, create the optional venv once:

```bash
module purge
module load BasePy/2026
python3 -m venv "$HOME/venvs/gmd2026-netcdf4"
source "$HOME/venvs/gmd2026-netcdf4/bin/activate"
python -m pip install --upgrade pip setuptools wheel
python -m pip install netCDF4 numpy matplotlib
python -c "import netCDF4; print(netCDF4.__version__)"
```

The helper uses `$HOME/venvs/gmd2026-netcdf4/bin/python` automatically when it exists. It does not write persistent environment dump logs during normal runs.

## 2. Build

Build from one model run directory:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/00_base_restart_3d/base_restart_3d
../../common/build_squid.sh
```

If the build succeeds, the case directory should contain copied executables `scale-rm` and `scale-rm_init`. The build helper deliberately avoids `scale-rm*` wildcards and symbolic links.

## 3. Optional base 3D restart

Submit the 60 min spin-up restart example if you need a restart sanity check:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/00_base_restart_3d
bash submit_base_restart_squid.sh
```

After the job finishes, verify:

```bash
cd base_restart_3d
ls restart_output/base_restart_00000101-010000.000.pe*.nc
grep -n "TIME:" LOG.pe000000 | tail
cat job_metrics.json
```

The restart basename must encode the 60 min state: `00000101-010000.000`. This restart is not used by the current 30 min benchmark, 70 min TPHT, or 10 min output-interval groups because those groups need tracking state initialized from model start.
The common model runner creates `restart_output/` from `ATMOS_RESTART_OUT_BASENAME` before the model starts. If an older run reached 60 min but failed with `failed to open file :./restart_output/base_restart_...`, rerun `base_restart_3d` after updating the scripts or manually create `restart_output/` before rerunning.

## 4. Controlled 3D benchmark

Submit:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/01_bench_3d_samp_30min
bash submit_bench_3d_squid.sh
```

The submission script enters each case directory before `qsub` and chains model jobs with one `--after` dependency at a time:

```text
nt_nolog -> nt_coallog -> fw005_coallog -> bw005_coallog
```

This is slower than submitting all four model jobs concurrently, but it is robust on the SQUID configuration that rejects colon-separated multi-job dependencies.
Use the suite-level `analysis/` jobs for tables and figures after the model outputs exist.

After jobs finish, generate benchmark tables and figures through the unified analysis workflow from the suite root. For a group-specific run:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026
source common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" analysis/01_benchmark/analyze_01_benchmark.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs"
"${GMD2026_PYTHON}" analysis/01_benchmark/plot_01_benchmark.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs"
ls analysis_outputs/tables/01_benchmark_summary.*
```

Check that benchmark runs start from the initial physical state and use the 64-rank decomposition:

```bash
for case in nt_nolog nt_coallog fw005_coallog bw005_coallog; do
  echo "== $case =="
  grep -n "Open restart file" "$case/LOG.pe000000" | head
  grep -n "TIME:" "$case/LOG.pe000000" | head
  grep -n "PRC_NUM_X\\|PRC_NUM_Y" "$case/run.conf"
done
```

Expected: `ATMOS_RESTART_IN_BASENAME` points to `./init_dycom_rf02_00000101-000000.000`, model time starts at 0 s, and `PRC_NUM_X = PRC_NUM_Y = 8`.

The benchmark cases share the same initial files inside this group: `nt_nolog` runs `scale-rm_init`; the other three cases copy `init_dycom_rf02_*` and `random_number_init.pe*` from `../nt_nolog`.

## 5. TPHT 3D experiment

Submit the complete TPHT chain:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/02_tpht_3d_interest_70min
bash submit_tpht_squid.sh
```

The intended chain is:

1. `fw_discovery`: starts from the initial physical state and writes interest `.ids`.
2. `postprocess/serial_merge_ids.sh`: merges and deduplicates FW `.ids`.
3. `bw_reconstruction`: starts from the same initial physical state as FW plus the deduplicated `.ids` handoff.
4. `postprocess/serial_check_tpht.sh`: checks strict FW/BW ID consistency.

Important TPHT state rule:

- BW must not read `fw_discovery/fw_output/superdroplet_restart*`.
- BW must read `./init_dycom_rf02_00000101-000000.000`, matching FW discovery.
- BW must read target IDs from `../fw_discovery/fw_tracking/tracking_interest_ids_dedup`.

Before rerunning BW manually, verify that the merge really produced all per-rank handoff files:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/02_tpht_3d_interest_70min
ls fw_discovery/fw_tracking/tracking_interest_ids_dedup.pe000000.ids
ls fw_discovery/fw_tracking/tracking_interest_ids_dedup.pe*.ids | wc -l
awk '$1 == "TPROC" && $2 == "=" {print $3}' bw_reconstruction/Makefile
```

The file count should match `TPROC`. If `tracking_interest_ids_dedup.pe000000.ids` is missing, do not rerun `bw_reconstruction` directly. Rerun the merge step first:

```bash
cd postprocess
qsub serial_merge_ids.sh
# or, for an interactive check after loading the Python environment:
bash serial_merge_ids.sh
tail -n 40 logs/merge_tracking_interest_ids.log
```

`common/squid_model_job.sh` also performs this same `.ids` preflight before starting BW. This is intentional: SQUID `qsub --after` starts dependent jobs after the previous job finishes, not necessarily only after it succeeds.

The merge job uses `serial_merge_ids.sh`, but it is no longer single-core: the SQUID script requests 8 cores, defaults `TPHT_MERGE_WORKERS=8`, and writes its Python log to `postprocess/logs/merge_tracking_interest_ids.log`.

Verify the rule before submission:

```bash
grep -n "ATMOS_RESTART_IN_BASENAME\\|SD_IN_BASENAME\\|SD_OUT_BASENAME\\|tracking_id_input_basename" \
  bw_reconstruction/run.conf
```

Expected:

- `ATMOS_RESTART_IN_BASENAME = "./init_dycom_rf02_00000101-000000.000"`.
- `SD_IN_BASENAME` is absent or empty.
- `SD_OUT_BASENAME = "./bw_output/superdroplet_restart"`.
- `tracking_id_input_basename = "../fw_discovery/fw_tracking/tracking_interest_ids_dedup"`.
- `fw_discovery` has `HISTORY_DEFAULT_TINTERVAL = 0.0D0`.
- `bw_reconstruction` has `HISTORY_DEFAULT_TINTERVAL = 10.0D0`, matching `sdm_dmpitvl = 10.0d0`.
- `sdm_dmpitvl = 10.0d0` in both TPHT run directories.
- `HISTORY_OUTPUT_STEP0 = .true.` in BW because History output is enabled.

After the TPHT chain finishes:

```bash
ls fw_discovery/fw_tracking/tracking_interest_ids.pe*.ids
ls fw_discovery/fw_tracking/tracking_interest_ids_dedup.pe*.ids
python postprocess/check_tpht_consistency.py
```

Expected strict consistency:

```text
missing_in_bw = 0
extra_in_bw = 0
```

If `missing_in_bw` or `extra_in_bw` is non-zero, first check whether BW accidentally used a later FW physical restart, a different initial state, or a different MPI decomposition.

BW shares only FW's initial files. It must not use FW output at a non-initial time as a physical restart.

Full TPHT summary tables and plots are generated from the suite-level `analysis/` workflow, not from `02_tpht_3d_interest_70min/postprocess/`. For large outputs, prefer:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026
export GMD2026_ROOT=$(pwd)
export GMD2026_OUTDIR=$(pwd)/analysis_outputs
qsub analysis/job_scripts/submit_02_tpht_summary_parallel.sh
qsub analysis/job_scripts/submit_02_tpht_heavy_parallel.sh
```

## 6. 2D representativeness test

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/03_fw_rep_2d_600s
bash submit_fw_rep_2d_squid.sh
```

The script chains all generated 2D cases. Use the suite-level `analysis/03_sampling/` scripts for tables and figures after the model outputs exist.

After completion:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026
source common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" analysis/03_sampling/analyze_03_sampling.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs"
"${GMD2026_PYTHON}" analysis/03_sampling/plot_03_sampling.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs"
ls analysis_outputs/tables/03_sampling_metrics.*
```

Interpret this group only as a short sampling-procedure verification test, not as proof of long-time 3D representativeness.
The 60 s SD output interval over a 600 s quasi-2D run gives 11 output times, which is more useful for checking whether sampling errors are stable in time. The analysis uses `sd_n` multiplicity-weighted radius distributions and compares each sampling design against its own full-reference output. It is still not sufficient to claim 3D spatial representativeness. Use the 3D benchmark/TPHT groups, or add a dedicated optional 3D sampling test, for manuscript-level spatial-representativeness claims.

## 7. Optional groups

SDNC scaling:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/04_sdnc_scaling_lite_30min
bash submit_sdnc_scaling_squid.sh
```

This group runs SDNC = 10, 20, 40, and 80 for NT/FW005/BW005, all as 30 min cold-start performance cases without Eulerian history output.
The submission order is the shell-sorted case list. Generate tables and figures through `analysis/04_sdnc_scaling/` after the model jobs finish.
Within each SDNC value, the NT case provides the shared init for FW/BW. NT uses `tracking_selection_mode = "none"` and writes no SD output; it is a compute baseline, not an I/O-equivalent output baseline.

Output-interval I/O scaling:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/05_outint_io_lite_10min
bash submit_outint_io_squid.sh
```

This group runs 10 min cold-start FW005/BW005 cases for 30, 60, and 120 s SD output intervals plus `nt_nolog_10min`.
The submission order is `nt_nolog_10min`, then the shell-sorted FW/BW output-interval cases. Generate tables and figures through `analysis/05_outint_io/` after the model jobs finish.
`nt_nolog_10min` provides the shared init for the FW/BW output-interval cases.

## 8. Central failure summary

Each model job writes a status file under:

```text
scale-rm/test/case/shallowcloud/GMD2026/job_status/
```

The immediately updated single-file log is:

```text
scale-rm/test/case/shallowcloud/GMD2026/job_status/job_events.tsv
```

The batch scripts submit `common/serial_collect_failures.sh` at the end of each chain. If a chain stops before the final collector can run, inspect `job_status/job_events.tsv` first, then run the collector manually:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026
source common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" common/collect_job_failures.py
ls job_failure_summary.*
```

Inspect `job_failure_summary.md` first when a chain fails.

## 9. Final files to collect

Collect these files for the manuscript-response benchmark archive after `analysis_outputs/` and, if needed, `analysis_outputs_for_GMD/` are current:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026
ls config_audit.* restart_audit.* benchmark_diagnostics_availability.*
ls job_failure_summary.*
ls analysis_outputs/tables/01_benchmark_summary.*
ls analysis_outputs/tables/02_tpht_summary.*
ls analysis_outputs/tables/03_sampling_metrics.*
ls analysis_outputs/tables/04_sdnc_scaling_summary.*
ls analysis_outputs/tables/05_outint_io_summary.*
ls analysis_outputs_for_GMD/README_GMD_CANDIDATES.md
find . -name "*.pdf" -o -name "*.svg" -o -name "*.png"
```
