#!/usr/bin/env python3
"""Analyze GMD2026 scaling with super-droplet number (SDNC)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import safe_float, safe_ratio, safe_subtract, write_table_bundle  # noqa: E402


COLUMNS = [
    "sdnc",
    "tracking_label",
    "case_name",
    "wallclock_s",
    "core_hours",
    "wallclock_relative_to_sdnc10_same_mode",
    "wallclock_relative_to_nt_same_sdnc",
    "core_hours_relative_to_nt_same_sdnc",
    "peak_memory_rank_max_mib",
    "peak_memory_rank_sum_mib",
    "tracking_chain_count",
    "tracking_id_memory_bytes",
    "if_coal_memory_bytes",
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
    "sd_output_write_time_total_s",
    "id_assignment_time_s",
    "boundary_tracking_time_s",
    "fw_overhead_vs_nt_s",
    "bw_overhead_vs_nt_s",
    "fw_overhead_ratio",
    "bw_overhead_ratio",
    "warnings",
]


def _log_slope(rows: list[dict], metric: str) -> float | None:
    """Fit log(metric) = a + b log(SDNC) and return b."""
    points: list[tuple[float, float]] = []
    for row in rows:
        sdnc = safe_float(row.get("sdnc"))
        value = safe_float(row.get(metric))
        if sdnc is not None and value is not None and sdnc > 0.0 and value > 0.0:
            points.append((math.log(sdnc), math.log(value)))
    if len(points) < 2:
        return None
    mean_x = sum(point[0] for point in points) / len(points)
    mean_y = sum(point[1] for point in points) / len(points)
    variance = sum((point[0] - mean_x) ** 2 for point in points)
    if variance == 0.0:
        return None
    return sum((point[0] - mean_x) * (point[1] - mean_y) for point in points) / variance


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write SDNC scaling summary and slope tables."""
    config = load_config()
    rows: list[dict] = []
    warnings: list[str] = []
    options = options or {}
    cases = iter_cases(config, "04_sdnc_scaling_lite_30min")
    if options.get("quick"):
        cases = [case for case in cases if case["case_name"] in {"sdnc10_nt", "sdnc10_fw005", "sdnc10_bw005"}]
    for case in cases:
        row, _files, _records, current_warnings = collect_case_metrics(root, case, options)
        if case.get("sd_output_enabled") and (safe_float(row.get("sd_selected_output_bytes")) or 0.0) <= 0.0:
            current_warnings.append(f"{case['case_name']} should have selected SD output")
            row["sd_selected_output_bytes"] = None
        row["warnings"] = warning_text(current_warnings)
        warnings.extend(current_warnings)
        rows.append(row)

    by_key = {(safe_float(row.get("sdnc")), row.get("tracking_label")): row for row in rows}
    for row in rows:
        sdnc = safe_float(row.get("sdnc"))
        label = row.get("tracking_label")
        baseline_same_mode = by_key.get((10.0, label), {})
        nt_same_sdnc = by_key.get((sdnc, "NT"), {})
        row["wallclock_relative_to_sdnc10_same_mode"] = safe_ratio(row.get("wallclock_s"), baseline_same_mode.get("wallclock_s"))
        row["wallclock_relative_to_nt_same_sdnc"] = safe_ratio(row.get("wallclock_s"), nt_same_sdnc.get("wallclock_s"))
        row["core_hours_relative_to_nt_same_sdnc"] = safe_ratio(row.get("core_hours"), nt_same_sdnc.get("core_hours"))
        fw_row = by_key.get((sdnc, "FW005"), {})
        bw_row = by_key.get((sdnc, "BW005"), {})
        nt_wall = nt_same_sdnc.get("wallclock_s")
        row["fw_overhead_vs_nt_s"] = safe_subtract(fw_row.get("wallclock_s"), nt_wall)
        row["bw_overhead_vs_nt_s"] = safe_subtract(bw_row.get("wallclock_s"), nt_wall)
        row["fw_overhead_ratio"] = safe_ratio(fw_row.get("wallclock_s"), nt_wall)
        row["bw_overhead_ratio"] = safe_ratio(bw_row.get("wallclock_s"), nt_wall)

    slope_rows = []
    metrics = {
        "wallclock_s": "wallclock",
        "core_hours": "core-hours",
        "peak_memory_rank_max_mib": "memory",
        "tracking_chain_count": "tracking_chain_count",
        "scientific_output_bytes": "scientific_output_size",
        "sd_selected_output_bytes": "selected_sd_output_size",
    }
    for label in ("NT", "FW005", "BW005"):
        label_rows = [row for row in rows if row.get("tracking_label") == label]
        for metric, metric_label in metrics.items():
            slope_rows.append(
                {
                    "tracking_label": label,
                    "metric": metric_label,
                    "log_log_slope_b": _log_slope(label_rows, metric),
                    "fit_model": "log(metric) = a + b log(SDNC)",
                    "positive_finite_points": sum(
                        1
                        for row in label_rows
                        if safe_float(row.get("sdnc")) is not None
                        and safe_float(row.get(metric)) is not None
                        and (safe_float(row.get(metric)) or 0.0) > 0.0
                    ),
                }
            )

    observed_sdnc = {int(row["sdnc"]) for row in rows if safe_float(row.get("sdnc")) is not None}
    if observed_sdnc != {10, 20, 40, 80}:
        warnings.append(f"SDNC values should include 10, 20, 40, 80; found {sorted(observed_sdnc)}")
    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle(rows, outdir / "tables" / "04_sdnc_scaling_summary", COLUMNS)
    write_table_bundle(slope_rows, outdir / "tables" / "04_scaling_slopes")


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(
            Path(__file__).name,
            root,
            outdir,
            ["tables/04_sdnc_scaling_summary.{csv,md,tex,json}", "tables/04_scaling_slopes.{csv,md,tex,json}"],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
