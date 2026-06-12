# Fortran Diagnostics Hook Smoke Test on SQUID

This directory provides SQUID-first smoke tests for the optional GMD benchmark
diagnostics hook.  These cases are engineering checks only; they are not
manuscript science experiments.

The diagnostics are disabled by default in source code.  A case emits
`GMD_BENCH_DIAG` and `GMD_IO_DIAG` only when:

```fortran
gmd_benchmark_diag_enable = .true.
```

## Diagnostics Semantics

`GMD_BENCH_DIAG` includes:

- rank-level memory from Linux `/proc/self/status`;
- tracking ID array memory estimate;
- `if_coal` memory estimate;
- ID-assignment timing;
- lateral-boundary tracking timing.

Memory fields use binary MiB.  `peak_memory_rank_max_mib` and
`peak_memory_rank_sum_mib` are converted from `/proc/self/status` `VmHWM` kB by
dividing by 1024.  They must not be mixed with scheduler-level job memory.

`GMD_IO_DIAG` uses explicit last/total/count fields:

- `*_write_time_last_s` is the maximum rank-local time for the latest write call.
- `*_write_time_total_s` is the MPI-summed cumulative rank-local write time.
- `*_write_count` is the MPI-summed number of completed rank-local write calls.
- `tpht_id_records_written_total` is the MPI-summed number of TPHT ID records written.

Do not interpret a `*_last_s` field as a total or mean cost.

## Smoke-Test Cases

`prepare_hook_smoke_cases.py` creates:

| Case | Purpose |
| --- | --- |
| `shared_init_hook_2min` | Runs `scale-rm_init` once and provides common `init_dycom_rf02_*` and `random_number_init.pe*` files. |
| `nt_nolog_hook_2min` | No tracking, no coalescence logging; checks zero-output diagnostics. |
| `fw_all_hook_2min` | FW `tracking_selection_mode="none"` with `sdm_dmpvar=010`; writes `SD_all_NetCDF_*`. |
| `fw_all_unhooked_2min` | Same FW all-output case with `gmd_benchmark_diag_enable=.false.`. |
| `bw_selected_hook_2min` | BW stratified sampled tracking; writes `SD_selected_NetCDF_*`. |
| `tpht_ids_hook_2min` | TPHT FW interest-ID stream with `sdm_dmpvar=000`; writes `.ids` only. |
| `tpht_bw_hook_2min` | TPHT BW reconstruction after merge/dedup; writes selected target output. |
| `tpht_bw_missing_ids_failfast` | Negative control; missing `.ids` input must fail instead of falling back to sampling. |

All model runs use the shared initialization.  This is required for exact
hooked-vs-unhooked comparison: independent `scale-rm_init` runs can produce
different initial SD populations even if the namelist is the same.

`tracking_selection_mode="none"` means all valid SDs are tracked.  For ordinary
all-SD output, the corresponding smoke-test configuration uses
`sdm_dmpvar=010`, so the output file family is `SD_all_NetCDF_*`.  The source
also contains a compatibility guard that routes `none + sdm_dmpvar=100` to
all-SD output when no TPHT `.ids` input is active.  This guard does not affect
random, stratified, or TPHT BW `.ids` target reconstruction.

The `sdm_sdnmlvol`, `sdm_inisdnc`, and `sdm_zupper` namelist values are kept in
`init.conf` intentionally.  They are used by SDM initialization and particle
number calculations; the source defaults do not match this DYCOMS smoke-test
configuration.

## Python Environment

On SQUID, create the recommended netCDF4 venv once:

```bash
module purge
module load BasePy/2026
python3 -m venv "$HOME/venvs/gmd2026-netcdf4"
source "$HOME/venvs/gmd2026-netcdf4/bin/activate"
python -m pip install --upgrade pip setuptools wheel
python -m pip install netCDF4
python -c "import netCDF4; print(netCDF4.__version__)"
```

The helper `../common/load_basepy_2026_quiet.sh` uses this interpreter if it
exists.  It also suppresses module pager output that can otherwise require
manual `q` input.

## SQUID Procedure

From the unpacked repository:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/06_fortran_diagnostics_hook_smoke_squid
source ../common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" prepare_hook_smoke_cases.py
bash build_hook_smoke_squid.sh
bash submit_hook_smoke_squid.sh
```

`build_hook_smoke_squid.sh` regenerates the cases, builds once, and copies the
exact executables `scale-rm` and `scale-rm_init` into every generated case
directory.  It does not use symbolic links or `scale-rm*` wildcards.

`submit_hook_smoke_squid.sh` enters each case directory before calling `qsub`
and uses a single-dependency chain because this SQUID environment does not
accept colon-separated `--after` dependency lists.

The submission chain is:

```text
shared_init_hook_2min
  -> nt_nolog_hook_2min
  -> fw_all_hook_2min
  -> fw_all_unhooked_2min
  -> bw_selected_hook_2min
  -> tpht_ids_hook_2min
  -> tpht_bw_missing_ids_failfast
  -> merge/dedup
  -> tpht_bw_hook_2min
  -> postprocess/serial_check_hook.sh
  -> ../common/serial_collect_failures.sh
```

To run the checker manually after all jobs finish:

```bash
source ../common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" postprocess/check_hook_diagnostics.py
```

The checker writes:

```text
postprocess/hook_diagnostics_check.csv
postprocess/hook_diagnostics_check.md
postprocess/hook_diagnostics_check.json
postprocess/logs/check_hook_diagnostics.log
```

The final failure-summary job updates `../job_failure_summary.*`. These files
are shared by groups `00`-`06` and are the first place to check for failed or
not-yet-run model jobs.

For Mac or other non-Linux inspection only, missing `/proc/self/status` can be
allowed explicitly:

```bash
HOOK_DIAG_ALLOW_NONLINUX=1 "${GMD2026_PYTHON}" postprocess/check_hook_diagnostics.py
```
