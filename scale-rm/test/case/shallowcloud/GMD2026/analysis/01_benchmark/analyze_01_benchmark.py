#!/usr/bin/env python3
"""Analyze the controlled cold-start 3D GMD2026 sampled benchmark."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.parse_netcdf import validate_selected_output  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import append_warning, safe_float, safe_ratio, safe_subtract, write_table_bundle  # noqa: E402


COLUMNS = [
    "case_name",
    "tracking_label",
    "tracking_mode",
    "wallclock_s",
    "wallclock_relative_to_nt_nolog",
    "core_hours",
    "core_hours_relative_to_nt_nolog",
    "mpi_ranks",
    "omp_threads",
    "node_count",
    "peak_memory_rank_max_mib",
    "peak_memory_rank_sum_mib",
    "tracking_chain_count",
    "tracking_id_memory_bytes",
    "if_coal_memory_bytes",
    "id_assignment_time_s",
    "boundary_tracking_time_s",
    "sd_output_write_time_total_s",
    "sd_output_write_time_mean_s",
    "sd_output_write_count",
    "coalescence_output_write_time_total_s",
    "coalescence_output_write_time_mean_s",
    "coalescence_output_write_count",
    "total_output_bytes",
    "scientific_output_bytes",
    "scientific_output_file_count",
    "auxiliary_file_bytes",
    "auxiliary_file_count",
    "sd_selected_output_bytes",
    "sd_selected_file_count",
    "sd_all_output_bytes",
    "sd_all_file_count",
    "coalescence_log_bytes",
    "coalescence_log_file_count",
    "history_output_bytes",
    "history_output_file_count",
    "number_of_output_files",
    "mean_file_size_bytes",
    "coalescence_event_count",
    "bytes_per_coalescence_event",
    "first_hour_cap_active_fraction",
    "coalescence_logging_overhead_s",
    "fw_tracking_overhead_s",
    "bw_tracking_overhead_s",
    "output_bytes_gib",
    "selected_output_gib",
    "coalescence_log_gib",
    "warnings",
]


def _check_outputs(row: dict, files_by_metric: dict, warnings: list[str]) -> None:
    """Run group-01 output-presence and NetCDF variable checks."""
    case_name = row.get("case_name")
    selected_bytes = safe_float(row.get("sd_selected_output_bytes")) or 0.0
    coal_bytes = safe_float(row.get("coalescence_log_bytes")) or 0.0
    if case_name == "nt_nolog":
        if selected_bytes > 0.0:
            append_warning(warnings, "nt_nolog should not have SD_selected output")
        if coal_bytes > 0.0:
            append_warning(warnings, "nt_nolog should not have coalescence-log output")
    if case_name == "nt_coallog":
        if selected_bytes > 0.0:
            append_warning(warnings, "nt_coallog should not have SD_selected output")
        if coal_bytes <= 0.0:
            append_warning(warnings, "nt_coallog should have coalescence-log output")
            row["coalescence_log_bytes"] = None
            row["coalescence_log_gib"] = None
    if case_name in {"fw005_coallog", "bw005_coallog"}:
        if selected_bytes <= 0.0:
            append_warning(warnings, f"{case_name} should have selected SD output")
            row["sd_selected_output_bytes"] = None
            row["selected_output_gib"] = None
        if coal_bytes <= 0.0:
            append_warning(warnings, f"{case_name} should have coalescence-log output")
            row["coalescence_log_bytes"] = None
            row["coalescence_log_gib"] = None
        mode = "FW" if case_name.startswith("fw") else "BW"
        validation, validation_warnings = validate_selected_output(files_by_metric.get("sd_selected_output_bytes", []), mode)
        row.update(validation)
        warnings.extend(validation_warnings)
    event_count = safe_float(row.get("coalescence_event_count"))
    if coal_bytes > 0.0 and (event_count is None or event_count < 10):
        append_warning(warnings, "Coalescence-log overhead is weakly constrained because event count is small.")


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write the group-01 benchmark summary."""
    config = load_config()
    rows: list[dict] = []
    warnings_by_case: dict[str, list[str]] = {}
    for case in iter_cases(config, "01_bench_3d_samp_30min"):
        row, files_by_metric, _records, warnings = collect_case_metrics(root, case, options)
        _check_outputs(row, files_by_metric, warnings)
        row["warnings"] = warning_text(warnings)
        rows.append(row)
        warnings_by_case[row["case_name"]] = warnings

    by_case = {row["case_name"]: row for row in rows}
    baseline_wall = by_case.get("nt_nolog", {}).get("wallclock_s")
    baseline_core = by_case.get("nt_nolog", {}).get("core_hours")
    coallog_overhead = safe_subtract(by_case.get("nt_coallog", {}).get("wallclock_s"), baseline_wall)
    fw_overhead = safe_subtract(by_case.get("fw005_coallog", {}).get("wallclock_s"), by_case.get("nt_coallog", {}).get("wallclock_s"))
    bw_overhead = safe_subtract(by_case.get("bw005_coallog", {}).get("wallclock_s"), by_case.get("nt_coallog", {}).get("wallclock_s"))

    for row in rows:
        row["wallclock_relative_to_nt_nolog"] = safe_ratio(row.get("wallclock_s"), baseline_wall)
        row["core_hours_relative_to_nt_nolog"] = safe_ratio(row.get("core_hours"), baseline_core)
        row["coalescence_logging_overhead_s"] = coallog_overhead
        row["fw_tracking_overhead_s"] = fw_overhead
        row["bw_tracking_overhead_s"] = bw_overhead

    all_warnings = [message for warnings in warnings_by_case.values() for message in warnings]
    if strict and all_warnings:
        raise RuntimeError("; ".join(all_warnings))
    write_table_bundle(rows, outdir / "tables" / "01_benchmark_summary", COLUMNS)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["tables/01_benchmark_summary.{csv,md,tex,json}"])
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
