#!/usr/bin/env python3
"""Analyze GMD2026 selected-output interval I/O sensitivity."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import safe_float, safe_ratio, safe_subtract, write_table_bundle  # noqa: E402


COLUMNS = [
    "case_name",
    "tracking_label",
    "output_interval_s",
    "wallclock_s",
    "wallclock_overhead_vs_nt_s",
    "wallclock_ratio_vs_nt",
    "core_hours",
    "peak_memory_rank_max_mib",
    "peak_memory_rank_sum_mib",
    "tracking_chain_count",
    "sd_output_write_time_total_s",
    "sd_output_write_time_mean_s",
    "sd_output_write_count",
    "total_output_bytes",
    "scientific_output_bytes",
    "scientific_output_file_count",
    "sd_selected_output_bytes",
    "sd_selected_file_count",
    "sd_all_output_bytes",
    "sd_all_file_count",
    "sd_particle_output_bytes",
    "sd_particle_file_count",
    "history_output_bytes",
    "history_output_file_count",
    "auxiliary_file_bytes",
    "auxiliary_file_count",
    "number_of_output_files",
    "mean_file_size_bytes",
    "output_bytes_per_event",
    "output_bytes_per_tracked_chain",
    "write_time_per_event_s",
    "write_time_per_tracked_chain_s",
    "warnings",
]


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write output-interval I/O summary tables."""
    config = load_config()
    rows: list[dict] = []
    warnings: list[str] = []
    options = options or {}
    cases = iter_cases(config, "05_outint_io_lite_10min")
    if options.get("quick"):
        cases = [case for case in cases if case["case_name"] in {"nt_nolog_10min", "out60_fw005", "out60_bw005"}]
    for case in cases:
        row, _files, _records, current_warnings = collect_case_metrics(root, case, options)
        if case.get("sd_output_enabled") and (safe_float(row.get("sd_selected_output_bytes")) or 0.0) <= 0.0:
            current_warnings.append(f"{case['case_name']} should have selected SD output")
            row["sd_selected_output_bytes"] = None
        row["warnings"] = warning_text(current_warnings)
        warnings.extend(current_warnings)
        rows.append(row)

    by_case = {row["case_name"]: row for row in rows}
    baseline = by_case.get("nt_nolog_10min", {})
    baseline_wall = baseline.get("wallclock_s")
    observed = {(row.get("tracking_label"), safe_float(row.get("output_interval_s"))) for row in rows if row.get("tracking_label") != "NT"}
    for row in rows:
        chain_count = safe_float(row.get("tracking_chain_count"))
        write_count = safe_float(row.get("sd_output_write_count"))
        selected_bytes = safe_float(row.get("sd_selected_output_bytes"))
        write_total = safe_float(row.get("sd_output_write_time_total_s"))
        row["wallclock_overhead_vs_nt_s"] = safe_subtract(row.get("wallclock_s"), baseline_wall)
        row["wallclock_ratio_vs_nt"] = safe_ratio(row.get("wallclock_s"), baseline_wall)
        row["output_bytes_per_event"] = safe_ratio(selected_bytes, write_count)
        row["output_bytes_per_tracked_chain"] = safe_ratio(selected_bytes, chain_count)
        row["write_time_per_event_s"] = safe_ratio(write_total, write_count)
        row["write_time_per_tracked_chain_s"] = safe_ratio(write_total, chain_count)

    expected = {("FW005", 30.0), ("FW005", 60.0), ("FW005", 120.0), ("BW005", 30.0), ("BW005", 60.0), ("BW005", 120.0)}
    missing = sorted(expected - observed)
    if missing:
        warnings.append(f"Output intervals must include 30, 60, 120 s for FW and BW; missing {missing}")
    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle(rows, outdir / "tables" / "05_outint_io_summary", COLUMNS)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["tables/05_outint_io_summary.{csv,md,tex,json}"])
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
