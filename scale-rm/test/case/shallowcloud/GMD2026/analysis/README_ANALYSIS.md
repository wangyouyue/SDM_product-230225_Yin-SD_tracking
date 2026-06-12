# GMD2026 Analysis Workflow

This directory contains command-line analysis and plotting scripts for the GMD2026 SCALE-SDM manuscript revision. The scripts target completed outputs on SQUID and do not assume that simulation outputs exist on the local workstation.

The analysis separates physical sanity checks, tracking overhead, coalescence-log overhead, SD output I/O cost, memory cost, SDNC scaling, output-frequency sensitivity, and TPHT target-restricted backward reconstruction efficiency.

## Run on SQUID

The completed GMD2026 simulation outputs are expected to already exist under the suite root. From that root, run post-processing in this order:

1. Load the same Python environment helper used by the SQUID job scripts.
2. Run the cost estimator and optional foreground metadata checks.
3. Submit full analysis jobs with qsub for table generation.
4. Submit plotting only after the analysis tables are complete.
5. Inspect warnings, `02_tpht_analysis_feature_status.*`, and `GMD2026_reviewer_coverage.md` before moving numbers into the manuscript.

```bash
cd /path/to/scale-rm/test/case/shallowcloud/GMD2026
source common/load_basepy_2026_quiet.sh
echo "$GMD2026_PYTHON_ENV"
"${GMD2026_PYTHON}" -c "import netCDF4; print(netCDF4.__version__)"

"${GMD2026_PYTHON}" analysis/common/estimate_analysis_cost.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs"
```

Use `--strict` to stop on missing cases, missing LOG.pe* diagnostics, failed TPHT consistency, or failed QC checks. Without `--strict`, missing diagnostics are written as `NA`, warnings are emitted, and the remaining tables continue to be generated. Use `--dry-run` to verify command routing without writing outputs.

## Foreground vs qsub execution

Start with the cost estimator after outputs exist on SQUID:

```bash
source common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" analysis/common/estimate_analysis_cost.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs"
```

This writes `tables/analysis_cost_estimate.csv`, `.md`, and `.json` with `foreground_ok`, `qsub_recommended`, or `qsub_required` classifications.

Foreground-safe commands:

- `--dry-run`
- `--quick`
- `--metadata-only`
- `00_restart` summary
- small log-only summaries
- plotting from already generated tables

Use qsub for:

- full `run_all_analysis.sh`
- `01_benchmark` full NetCDF scan
- `02_tpht` heavy target histories and chain validity
- `03_sampling` ensemble analysis
- `04_sdnc_scaling` full scan
- `05_outint_io` full scan
- `run_all_plots.sh` if generating all figures on the login node is discouraged

Foreground mode flags:

- `--quick`: minimal check and small summary.
- `--metadata-only`: parse logs, `job_metrics.json`, `time_scale_rm_main.log`, file sizes, and NetCDF headers without reading large variables.
- `--skip-heavy-netcdf`: skip DSD, z-r histograms, target histories, and chain reconstruction.
- `--max-files N`: process at most N NetCDF files for debugging.
- `--max-records N`: process at most N SD records per file when chunked reading is available.
- `--workers N`: pass a worker count to scripts that implement parallel processing.
- `--chunk-size N`: NetCDF records per chunk for memory-safe variable reading.

Quick local check:

```bash
bash analysis/run_all_analysis.sh \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs_quick" \
  --dry-run
```

Metadata-only foreground check:

```bash
source common/load_basepy_2026_quiet.sh
"${GMD2026_PYTHON}" analysis/01_benchmark/analyze_01_benchmark.py \
  --root "$(pwd)" \
  --outdir "$(pwd)/analysis_outputs_quick" \
  --metadata-only
```

Full analysis and plotting via qsub:

```bash
export GMD2026_ROOT=$(pwd)
export GMD2026_OUTDIR=$(pwd)/analysis_outputs
qsub analysis/job_scripts/submit_analysis_then_plots.sh
```

The qsub scripts default to `GMD2026_ROOT=$(pwd)` and `GMD2026_OUTDIR=$(pwd)/analysis_outputs` when these variables are not set. They source the suite Python helper (`common/load_basepy_2026_quiet.sh`) before calling `${GMD2026_PYTHON}`, set `PYTHONUNBUFFERED=1`, write logs under `analysis_outputs/logs/`, and emit `*.analysis_job_metrics.json`. The helper prefers the conda environment `sdm_env`; set `GMD2026_CONDA_ENV` before submission if the environment name differs. The shared qsub wrapper records `python_env`, `python_executable`, and the `netCDF4` import status in each job log.

TPHT is split into light and heavy execution. The light part reads logs, `.ids` metadata, handoff integrity, and inexpensive rank-load summaries. The heavy part reads BW selected-output NetCDF variables for target categories, histories, chain-validity diagnostics, target trajectories, if_coal timelines, interval diagnostics, and target-linked event summaries.

The serial TPHT jobs remain available:

```bash
qsub analysis/job_scripts/submit_02_tpht_summary.sh
qsub analysis/job_scripts/submit_02_tpht_heavy.sh
```

For large TPHT outputs, prefer the parallel jobs. They request 8 cores by default and use `TPHT_LIGHT_WORKERS` or `TPHT_HEAVY_WORKERS` to control the parser worker count:

```bash
qsub analysis/job_scripts/submit_02_tpht_summary_parallel.sh
qsub analysis/job_scripts/submit_02_tpht_heavy_parallel.sh
```

`submit_analysis_then_plots.sh` first tries `qsub --after "$job00:$job01:$job02l:$job02h:$job03:$job04:$job05"`. If SQUID rejects colon-separated dependencies, it creates a conservative polling aggregator job that waits for all analysis job IDs before running the plotting job.

## Manual qsub Workflow

Use this sequence when you do not want to submit the whole dependency workflow at once. The 00 restart sanity check is small, but submitting it keeps all logs and metrics in the same `analysis_outputs/logs/` structure.

```bash
cd /path/to/scale-rm/test/case/shallowcloud/GMD2026
export GMD2026_ROOT=$(pwd)
export GMD2026_OUTDIR=$(pwd)/analysis_outputs

qsub analysis/job_scripts/submit_00_restart.sh
qsub analysis/job_scripts/submit_01_benchmark.sh
qsub analysis/job_scripts/submit_02_tpht_summary_parallel.sh
```

After `submit_02_tpht_summary_parallel.sh` finishes, inspect `analysis_outputs/tables/02_tpht_summary.*`, `analysis_outputs/tables/02_tpht_consistency.*`, and `analysis_outputs/tables/02_tpht_analysis_feature_status.*`. If the TPHT handoff is valid and the BW selected output exists, submit the heavy TPHT job:

```bash
qsub analysis/job_scripts/submit_02_tpht_heavy_parallel.sh
```

Use the non-parallel scripts only when you intentionally want a single-worker diagnostic run or are debugging worker-specific parser behavior.

Then submit the remaining full analyses:

```bash
qsub analysis/job_scripts/submit_03_sampling.sh
qsub analysis/job_scripts/submit_04_sdnc_scaling.sh
qsub analysis/job_scripts/submit_05_outint_io.sh
```

Submit plotting only after all relevant analysis jobs have finished. `submit_all_plots.sh` also refreshes the combined manuscript summary before plotting:

```bash
qsub analysis/job_scripts/submit_all_plots.sh
```

Regenerate the manuscript-candidate GMD figure set as a separate background job after `analysis_outputs/` is current:

```bash
export GMD2026_GMD_OUTDIR=$(pwd)/analysis_outputs_for_GMD
qsub analysis/job_scripts/submit_regenerate_gmd_figures.sh
```

A monolithic table-generation job is also available, but it is `qsub_required` and should not be run in the login foreground:

```bash
qsub analysis/job_scripts/submit_all_analysis.sh
```

Use `qstat` and the files in `analysis_outputs/logs/` to check progress and failures. Re-run a failed group-specific submit script after fixing the underlying input or environment issue; missing diagnostics remain `NA` unless `--strict` is used.

## Outputs

Tables are written as `.csv`, `.md`, `.tex`, and `.json` under `analysis_outputs/tables/`. Figures are written as `.pdf`, `.svg`, and 600 dpi `.png` under `analysis_outputs/figures/`.

The cost-estimator table is written as `.csv`, `.md`, and `.json`. The manuscript key-results and reviewer-coverage notes are Markdown-only summary documents.

## Analysis Script Outputs

| Script or group | Output tables / data files | Main data contained | Output figures |
|---|---|---|---|
| `common/estimate_analysis_cost.py` | `tables/analysis_cost_estimate.csv`, `.md`, `.json` | per-group file counts, NetCDF/IDS/log counts, byte totals, estimated runtime class, memory class, and recommended foreground/qsub mode | none |
| `00_restart/analyze_00_restart.py` | `tables/00_restart_summary.csv`, `.md`, `.tex`, `.json` | restart-file presence, restart bytes, wallclock/core-hours, model-time coverage, exit status, first-hour cap flag, notes | none |
| `00_restart/plot_00_restart_timeline.py` | none | reads `00_restart_summary.csv` | `figures/00_restart_timeline.pdf`, `.svg`, `.png` |
| `01_benchmark/analyze_01_benchmark.py` | `tables/01_benchmark_summary.csv`, `.md`, `.tex`, `.json` | no-tracking baseline, coalescence-log overhead, FW/BW tracking overhead, wallclock/core-hours, memory, tracking chain count, I/O timing, output sizes, coalescence-event count, QC warnings | none |
| `01_benchmark/plot_01_benchmark.py` | none | reads `01_benchmark_summary.csv` | `figures/01_wallclock_relative.*`, `01_core_hours.*`, `01_output_size.*`, `01_peak_memory.*`, `01_runtime_components.*`, `01_tracking_chain_count.*` |
| `02_tpht/analyze_02_tpht.py` | `tables/02_tpht_summary.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_rank_load_balance.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_analysis_feature_status.csv`, `.md`, `.tex`, `.json` | FW/BW wallclock and core-hours, memory, ID handoff counts, missing/extra BW IDs, target reduction, storage estimates, rank load balance, implementation status of advanced TPHT diagnostics | none |
| `02_tpht/analyze_02_tpht_light_parallel.py` | `tables/02_tpht_summary.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_rank_load_balance.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_analysis_feature_status.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_consistency.csv`, `.md`, `.tex`, `.json` | parallel light TPHT driver for summary, rank-load, feature-status, and consistency tables; uses `--workers` or `TPHT_LIGHT_WORKERS` for `.ids`/rank-file parsing | none |
| `02_tpht/check_02_tpht_consistency.py` | `tables/02_tpht_consistency.csv`, `.md`, `.tex`, `.json` | TPHT consistency pass/fail fields: missing/extra BW IDs, epoch match, fallback-to-sampling status, decomposition match | none |
| `02_tpht/analyze_02_tpht_targets.py` | `tables/02_tpht_target_categories.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_discovery_time.csv`, `.md`, `.tex`, `.json` | interest-condition decomposition (`radius_only`, `coal_only`, `both`, `unknown`) from stepwise BW predecessor-link chains; target first/last discovery time, new targets per time, target records per time | none |
| `02_tpht/analyze_02_tpht_chains.py` | `tables/02_tpht_chain_validity.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_target_histories.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_chain_links_by_time.csv`, `.md`, `.tex`, `.json` | adjacent-output predecessor-link validity, reconstructed time coverage, memory-safe mean radius/height histories, compact chain-length distribution, and per-output link fractions | none |
| `02_tpht/analyze_02_tpht_science.py` | `tables/02_tpht_science_summary.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_science_pathways.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_science_formation_height_bins.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_science_coalescence_counts.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_science_time_series.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_science_target_summary.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_target_occurrence_zt.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_chain_links_by_time.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_target_trajectory_records.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_target_ifcoal_timeline.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_target_interval_diagnostics.csv`, `.md`, `.tex`, `.json`; `tables/02_tpht_target_event_links.csv`, `.md`, `.tex`, `.json` | TPHT diagnostic characterization from stepwise BW predecessor-link chains; first threshold crossing, if_coal occurrence proxy, target height-time occurrence, target trajectory records, interval diagnostics, and linked coalescence-event counters when available | none |
| `02_tpht/analyze_02_tpht_heavy_parallel.py` | heavy TPHT target, science, and chain tables listed above | parallel heavy TPHT driver for science, target, and chain diagnostics; uses `--workers` or `TPHT_HEAVY_WORKERS` for NetCDF/coalescence readers | none |
| `02_tpht/plot_02_tpht_science.py` | none | reads TPHT science diagnostic tables | `figures/02_tpht_science_formation_pathways.*`, `02_tpht_science_formation_height.*`, `02_tpht_science_coalescence_counts.*`, `02_tpht_science_history.*`, `02_tpht_science_condition_composition.*` |
| `02_tpht/plot_02_tpht_predecessor_tree.py` | `tables/02_tpht_predecessor_tree_examples.csv`, `.md`, `.tex`, `.json` when the cache is absent | reads or creates compact representative TPHT predecessor-tree example rows; branches are linked if_coal occurrence proxies using coalescence-log IDs and target height at the next selected-output level | `figures/02_tpht_predecessor_tree_examples.*` |
| `02_tpht/plot_02_tpht_predecessor_tree_final_window.py` | `analysis_outputs_for_GMD/tables/supplement_candidates/supp_candidate_TPHT_predecessor_tree_examples_final10_selected_targets.*` | final-window supplementary predecessor-tree diagnostic examples using target trajectory and target-linked event tables; writes into the manuscript-candidate output tree | `analysis_outputs_for_GMD/figures/supplement_candidates/supp_candidate_TPHT_predecessor_tree_examples_*.{pdf,svg,png}` |
| `02_tpht/plot_02_tpht.py` | none | reads the TPHT summary, rank-load, category, discovery-time, chain-validity, and target-history tables | `figures/02_tpht_workflow.*`, `02_tpht_cost.*`, `02_tpht_storage.*`, `02_tpht_id_counts.*`, `02_tpht_target_categories.*`, `02_tpht_discovery_time.*`, `02_tpht_rank_load_balance.*`, `02_tpht_chain_validity.*`, `02_tpht_target_histories.*` |
| `03_sampling/analyze_03_sampling.py` | `tables/03_sampling_metrics.csv`, `.md`, `.tex`, `.json`; `tables/03_sampling_seed_statistics.csv`, `.md`, `.tex`, `.json` | short 2D sampling-verification metrics by time, sample mode, fraction, and seed; `sd_n` multiplicity-weighted radius statistics; weighted radius-distribution L1/L2 errors vs full reference; weighted threshold fractions; seed counts and aggregate errors | none |
| `03_sampling/plot_03_sampling.py` | none | reads `03_sampling_metrics.csv` and `03_sampling_seed_statistics.csv` | `figures/03_sampling_error_boxplot.*`, `03_sampling_error_vs_fraction.*`, `03_sampling_radius_threshold_error.*`, `03_sampling_dsd_example.*` |
| `04_sdnc_scaling/analyze_04_sdnc_scaling.py` | `tables/04_sdnc_scaling_summary.csv`, `.md`, `.tex`, `.json`; `tables/04_scaling_slopes.csv`, `.md`, `.tex`, `.json` | SDNC-dependent wallclock/core-hours, memory, output size, chain count, FW/BW overhead ratios, and log-log scaling slopes for wallclock, core-hours, memory, chain count, and output size | none |
| `04_sdnc_scaling/plot_04_sdnc_scaling.py` | none | reads `04_sdnc_scaling_summary.csv` | `figures/04_wallclock_vs_sdnc.*`, `04_core_hours_vs_sdnc.*`, `04_memory_vs_sdnc.*`, `04_tracking_overhead_ratio.*`, `04_chain_count_vs_sdnc.*` |
| `05_outint_io/analyze_05_outint_io.py` | `tables/05_outint_io_summary.csv`, `.md`, `.tex`, `.json` | output-interval sensitivity, wallclock overhead vs `nt_nolog_10min`, selected-output bytes, file counts, write-time totals/means, bytes/write-time per event or tracked chain | none |
| `05_outint_io/plot_05_outint_io.py` | none | reads `05_outint_io_summary.csv` | `figures/05_output_size_vs_interval.*`, `05_write_time_vs_interval.*`, `05_file_count_vs_interval.*`, `05_wallclock_overhead_vs_interval.*`, `05_output_efficiency.*` |
| `create_combined_summary.py` | `tables/GMD2026_combined_performance_summary.csv`, `.md`, `.tex`, `.json`; `tables/GMD2026_key_results_for_manuscript.md`; `tables/GMD2026_reviewer_coverage.md` | cross-group manuscript summary, concise key-result bullets, and reviewer-concern coverage map with remaining manuscript-only issues | none |
| `run_all_analysis.sh` | all analysis tables above, including `analysis_cost_estimate.*` and combined manuscript summaries | orchestrates all analysis scripts; use qsub for full runs | none |
| `run_all_plots.sh` | none | reads already generated tables | all figures listed above |

## Interpretation Guardrails

`00_base_restart_3d` is a restart sanity check only and is not included in tracking-performance comparisons.

`01_bench_3d_samp_30min` is a controlled cold-start 3D computational benchmark. It should not be described as a mature-cloud physical benchmark.

`02_tpht_3d_interest_70min` is the key 3D TPHT demonstration. It analyzes FW discovery, raw `.ids`, merge/dedup, BW reconstruction, target-set handoff, consistency, storage, rank_load_balance, chain_validity, and target-history diagnostics.

`03_fw_rep_2d_600s` is a short quasi-2D sampling-procedure verification. It is not proof of long-time 3D statistical representativeness.

## Diagnostics and Keywords

The parsers recognize `GMD_BENCH_DIAG`, `GMD_IO_DIAG`, `GMD_CAP_DIAG`, `TPHT_ID_DIAG`, `job_metrics.json`, `time_scale_rm_main.log`, `LOG.pe*`, `SD_selected_NetCDF_*`, `SD_coal_output_NetCDF_*`, `tracking_interest_ids`, `tracking_interest_ids_merged.ids`, and `tracking_interest_ids_dedup`.

Reported manuscript metrics include `wallclock_s`, `core_hours`, `peak_memory_rank_max_mib`, `peak_memory_rank_sum_mib`, `tracking_chain_count`, `tracking_id_memory_bytes`, `if_coal_memory_bytes`, `id_assignment_time_s`, `boundary_tracking_time_s`, `sd_output_write_time_total_s`, `coalescence_output_write_time_total_s`, `tpht_id_write_time_total_s`, `target_reduction_ratio`, `dedup_reduction_ratio`, `rank_load_balance`, `chain_validity`, sampling error, random sampling, stratified sampling, SDNC scaling, and output interval I/O scaling.

Plotting uses publication-quality Matplotlib only, with the Okabe-Ito palette.
