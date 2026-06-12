#!/usr/bin/env python3
"""Analyze TPHT target diagnostic histories from BW selected output."""

from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import safe_float, safe_ratio, write_table_bundle  # noqa: E402
from tpht_chain_utils import (  # noqa: E402
    EVENT_LINK_COLUMNS,
    IFCOAL_TIMELINE_COLUMNS,
    INTERVAL_COLUMNS,
    LINK_COLUMNS,
    TRAJECTORY_COLUMNS,
    build_optional_source_tables,
    build_stepwise_bw_diagnostics,
)


HEIGHT_BINS_M = (0.0, 250.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, math.inf)

SUMMARY_COLUMNS = [
    "target_count",
    "record_count",
    "large_target_count",
    "coal_target_count",
    "large_and_coal_target_count",
    "large_without_coal_target_count",
    "coal_before_large_count",
    "large_before_coal_count",
    "simultaneous_large_and_coal_count",
    "timing_unknown_count",
    "median_first_large_time_s",
    "median_first_coal_time_s",
    "median_first_large_height_m",
    "median_first_large_radius_um",
    "median_final_radius_um",
    "median_max_radius_um",
    "median_coal_record_count",
    "median_coal_episode_count_proxy",
    "coal_event_target_count",
    "coal_event_count_total",
    "coal_num_col_sum_total",
    "warnings",
]

PATHWAY_COLUMNS = ["pathway", "target_count", "target_fraction", "warnings"]
HEIGHT_BIN_COLUMNS = ["height_bin_m", "target_count", "target_fraction", "warnings"]
COAL_COLUMNS = ["coal_episode_count_proxy", "target_count", "target_fraction", "warnings"]
TIME_COLUMNS = [
    "time_s",
    "record_count",
    "mean_radius_um",
    "mean_height_m",
    "large_record_fraction",
    "coal_record_fraction",
    "warnings",
]
TARGET_SUMMARY_COLUMNS = [
    "target_id",
    "dm_id",
    "sd_id",
    "record_count",
    "first_time_s",
    "last_time_s",
    "first_selected_category",
    "category",
    "first_selected_time_s",
    "first_large_time_s",
    "first_ifcoal_time_s",
    "first_coal_event_time_s",
    "first_large_height_m",
    "first_large_radius_um",
    "first_large_x_m",
    "first_large_y_m",
    "final_radius_um",
    "max_radius_um",
    "duration_ge_15_s",
    "duration_ge_20_s",
    "coal_record_count",
    "coal_episode_count_proxy",
    "coal_event_count",
    "coal_num_col_sum",
    "diagnostic_pathway",
    "warnings",
]
OCCURRENCE_ZT_COLUMNS = [
    "time_s",
    "height_bin_m",
    "target_count",
    "large_target_count",
    "ifcoal_target_count",
    "warnings",
]
EXAMPLE_COLUMNS = [
    "rank",
    "dm_id",
    "sd_id",
    "record_count",
    "first_time_s",
    "last_time_s",
    "first_large_time_s",
    "first_coal_time_s",
    "first_coal_event_time_s",
    "first_large_height_m",
    "first_large_x_m",
    "first_large_y_m",
    "final_radius_um",
    "max_radius_um",
    "coal_record_count",
    "coal_episode_count_proxy",
    "coal_event_count",
    "coal_num_col_sum",
    "formation_pathway",
    "warnings",
]


def _log(message: str) -> None:
    """Write a flushed progress message for long SQUID TPHT analysis jobs."""
    print(f"[02_tpht_science] {message}", flush=True)


def _emit_warnings(warnings: list[str], limit: int = 20) -> None:
    """Write a bounded warning list to stderr without flooding batch logs."""
    if not warnings:
        return
    print(f"[02_tpht_science] warnings={len(warnings)}", file=sys.stderr, flush=True)
    for warning in warnings[:limit]:
        print(f"[02_tpht_science][warning] {warning}", file=sys.stderr, flush=True)
    if len(warnings) > limit:
        print(f"[02_tpht_science][warning] ... {len(warnings) - limit} more warnings omitted", file=sys.stderr, flush=True)


def _missing_rows(columns: list[str], warnings: list[str]) -> list[dict[str, Any]]:
    """Return one explicit NA row with warnings for an unavailable heavy table."""
    row = {column: None for column in columns}
    row["warnings"] = warning_text(warnings)
    return [row]


def _quantile(values: list[float], fraction: float) -> float | None:
    """Return a linear quantile for finite values."""
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return None
    position = (len(finite) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(finite) - 1)
    weight = position - lower
    return finite[lower] * (1.0 - weight) + finite[upper] * weight


def _pathway(state: dict[str, Any]) -> str:
    """Classify the order of large-radius and coalescence diagnostics."""
    has_large = bool(state.get("has_large"))
    has_coal = bool(state.get("has_coal"))
    if has_large and not has_coal:
        return "large_without_coal"
    if has_coal and not has_large:
        return "coal_without_large"
    if not has_large and not has_coal:
        return "no_large_or_coal"
    large_time = state.get("first_large_time_s")
    coal_time = state.get("first_coal_time_s")
    if large_time is not None and coal_time is not None:
        if coal_time < large_time:
            return "coal_before_large"
        if large_time < coal_time:
            return "large_before_coal"
        return "simultaneous_large_and_coal"
    large_order = state.get("first_large_order")
    coal_order = state.get("first_coal_order")
    if large_order is not None and coal_order is not None:
        if coal_order < large_order:
            return "coal_before_large"
        if large_order < coal_order:
            return "large_before_coal"
        return "simultaneous_large_and_coal"
    return "timing_unknown"


def _height_bin(height_m: float | None) -> str:
    """Return a compact height-bin label."""
    if height_m is None:
        return "unknown"
    for lower, upper in zip(HEIGHT_BINS_M[:-1], HEIGHT_BINS_M[1:]):
        if lower <= height_m < upper:
            return f"{lower:.0f}-{upper:.0f}" if math.isfinite(upper) else f">={lower:.0f}"
    return "unknown"


def _build_rows(targets: dict[str, dict[str, Any]], by_time: dict[float, dict[str, float]], by_time_height: dict[tuple[float, float], dict[str, float]], warnings: list[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Build aggregate science-analysis output rows."""
    pathway_counts: Counter[str] = Counter()
    height_counts: Counter[str] = Counter()
    coal_counts: Counter[str] = Counter()
    first_large_times: list[float] = []
    first_coal_times: list[float] = []
    first_large_heights: list[float] = []
    first_large_radii_um: list[float] = []
    final_radii_um: list[float] = []
    max_radii_um: list[float] = []
    coal_record_counts: list[float] = []
    coal_episode_counts: list[float] = []
    target_summary_rows: list[dict[str, Any]] = []
    example_rows: list[dict[str, Any]] = []

    for rank, (chain_id, state) in enumerate(sorted(targets.items(), key=lambda item: safe_float(item[1].get("max_radius_m")) or -1.0, reverse=True)[:30], start=1):
        pathway = _pathway(state)
        example_rows.append(
            {
                "rank": rank,
                "dm_id": state.get("dm_id"),
                "sd_id": state.get("sd_id"),
                "record_count": state.get("record_count"),
                "first_time_s": state.get("first_time_s"),
                "last_time_s": state.get("last_time_s"),
                "first_large_time_s": state.get("first_large_time_s"),
                "first_coal_time_s": state.get("first_coal_time_s"),
                "first_coal_event_time_s": state.get("first_coal_event_time_s"),
                "first_large_height_m": state.get("first_large_height_m"),
                "first_large_x_m": state.get("first_large_x_m"),
                "first_large_y_m": state.get("first_large_y_m"),
                "final_radius_um": (state.get("final_radius_m") * 1.0e6) if state.get("final_radius_m") is not None else None,
                "max_radius_um": (state.get("max_radius_m") * 1.0e6) if state.get("max_radius_m") is not None else None,
                "coal_record_count": state.get("coal_record_count"),
                "coal_episode_count_proxy": state.get("coal_episode_count_proxy"),
                "coal_event_count": state.get("coal_event_count"),
                "coal_num_col_sum": state.get("coal_num_col_sum"),
                "formation_pathway": pathway,
                "warnings": None,
            }
        )

    for state in targets.values():
        pathway = _pathway(state)
        pathway_counts[pathway] += 1
        height_counts[_height_bin(state.get("first_large_height_m"))] += 1
        coal_counts[str(min(int(state.get("coal_episode_count_proxy") or 0), 5)) if (state.get("coal_episode_count_proxy") or 0) < 5 else ">=5"] += 1
        if state.get("first_large_time_s") is not None:
            first_large_times.append(float(state["first_large_time_s"]))
        if state.get("first_coal_time_s") is not None:
            first_coal_times.append(float(state["first_coal_time_s"]))
        if state.get("first_large_height_m") is not None:
            first_large_heights.append(float(state["first_large_height_m"]))
        if state.get("first_large_radius_m") is not None:
            first_large_radii_um.append(float(state["first_large_radius_m"]) * 1.0e6)
        if state.get("final_radius_m") is not None:
            final_radii_um.append(float(state["final_radius_m"]) * 1.0e6)
        if state.get("max_radius_m") is not None:
            max_radii_um.append(float(state["max_radius_m"]) * 1.0e6)
        coal_record_counts.append(float(state.get("coal_record_count") or 0.0))
        coal_episode_counts.append(float(state.get("coal_episode_count_proxy") or 0.0))

    for chain_id, state in sorted(targets.items(), key=lambda item: (safe_float(item[1].get("dm_id")) or -1.0, safe_float(item[1].get("sd_id")) or -1.0, item[0])):
        target_summary_rows.append(
            {
                "target_id": chain_id,
                "dm_id": state.get("dm_id"),
                "sd_id": state.get("sd_id"),
                "record_count": state.get("record_count"),
                "first_time_s": state.get("first_time_s"),
                "last_time_s": state.get("last_time_s"),
                "first_selected_category": state.get("first_selected_category") or "unknown",
                "category": state.get("first_selected_category") or "unknown",
                "first_selected_time_s": state.get("first_selected_time_s"),
                "first_large_time_s": state.get("first_large_time_s"),
                "first_ifcoal_time_s": state.get("first_coal_time_s"),
                "first_coal_event_time_s": state.get("first_coal_event_time_s"),
                "first_large_height_m": state.get("first_large_height_m"),
                "first_large_radius_um": (state.get("first_large_radius_m") * 1.0e6) if state.get("first_large_radius_m") is not None else None,
                "first_large_x_m": state.get("first_large_x_m"),
                "first_large_y_m": state.get("first_large_y_m"),
                "final_radius_um": (state.get("final_radius_m") * 1.0e6) if state.get("final_radius_m") is not None else None,
                "max_radius_um": (state.get("max_radius_m") * 1.0e6) if state.get("max_radius_m") is not None else None,
                "duration_ge_15_s": state.get("duration_ge_15_s"),
                "duration_ge_20_s": state.get("duration_ge_20_s"),
                "coal_record_count": state.get("coal_record_count"),
                "coal_episode_count_proxy": state.get("coal_episode_count_proxy"),
                "coal_event_count": state.get("coal_event_count"),
                "coal_num_col_sum": state.get("coal_num_col_sum"),
                "diagnostic_pathway": _pathway(state),
                "warnings": None,
            }
        )

    target_count = len(targets)
    summary_rows = [
        {
            "target_count": target_count or None,
            "record_count": sum(int(state.get("record_count") or 0) for state in targets.values()) if targets else None,
            "large_target_count": sum(1 for state in targets.values() if state.get("has_large")) if targets else None,
            "coal_target_count": sum(1 for state in targets.values() if state.get("has_coal")) if targets else None,
            "large_and_coal_target_count": sum(1 for state in targets.values() if state.get("has_large") and state.get("has_coal")) if targets else None,
            "large_without_coal_target_count": pathway_counts.get("large_without_coal") if targets else None,
            "coal_before_large_count": pathway_counts.get("coal_before_large") if targets else None,
            "large_before_coal_count": pathway_counts.get("large_before_coal") if targets else None,
            "simultaneous_large_and_coal_count": pathway_counts.get("simultaneous_large_and_coal") if targets else None,
            "timing_unknown_count": pathway_counts.get("timing_unknown") if targets else None,
            "median_first_large_time_s": _quantile(first_large_times, 0.50),
            "median_first_coal_time_s": _quantile(first_coal_times, 0.50),
            "median_first_large_height_m": _quantile(first_large_heights, 0.50),
            "median_first_large_radius_um": _quantile(first_large_radii_um, 0.50),
            "median_final_radius_um": _quantile(final_radii_um, 0.50),
            "median_max_radius_um": _quantile(max_radii_um, 0.50),
            "median_coal_record_count": _quantile(coal_record_counts, 0.50),
            "median_coal_episode_count_proxy": _quantile(coal_episode_counts, 0.50),
            "coal_event_target_count": sum(1 for state in targets.values() if state.get("has_coal_event")) if targets else None,
            "coal_event_count_total": sum(int(state.get("coal_event_count") or 0) for state in targets.values()) if targets else None,
            "coal_num_col_sum_total": sum(float(state.get("coal_num_col_sum") or 0.0) for state in targets.values()) if targets else None,
            "warnings": warning_text(warnings),
        }
    ]

    pathway_order = ["coal_before_large", "large_before_coal", "simultaneous_large_and_coal", "large_without_coal", "coal_without_large", "no_large_or_coal", "timing_unknown"]
    pathway_rows = [
        {"pathway": pathway, "target_count": pathway_counts.get(pathway), "target_fraction": safe_ratio(pathway_counts.get(pathway), target_count), "warnings": warning_text(warnings) if index == 0 and warnings else None}
        for index, pathway in enumerate(pathway_order)
        if pathway_counts.get(pathway)
    ] or [{"pathway": None, "target_count": None, "target_fraction": None, "warnings": warning_text(warnings)}]

    height_order = [f"{lower:.0f}-{upper:.0f}" if math.isfinite(upper) else f">={lower:.0f}" for lower, upper in zip(HEIGHT_BINS_M[:-1], HEIGHT_BINS_M[1:])] + ["unknown"]
    height_rows = [
        {"height_bin_m": label, "target_count": height_counts.get(label), "target_fraction": safe_ratio(height_counts.get(label), target_count), "warnings": warning_text(warnings) if index == 0 and warnings else None}
        for index, label in enumerate(height_order)
        if height_counts.get(label)
    ] or [{"height_bin_m": None, "target_count": None, "target_fraction": None, "warnings": warning_text(warnings)}]

    coal_rows = [
        {"coal_episode_count_proxy": label, "target_count": coal_counts.get(label), "target_fraction": safe_ratio(coal_counts.get(label), target_count), "warnings": warning_text(warnings) if index == 0 and warnings else None}
        for index, label in enumerate(["0", "1", "2", "3", "4", ">=5"])
        if coal_counts.get(label)
    ] or [{"coal_episode_count_proxy": None, "target_count": None, "target_fraction": None, "warnings": warning_text(warnings)}]

    time_rows = []
    for time_s, values in sorted(by_time.items()):
        count = values.get("count", 0.0)
        time_rows.append(
            {
                "time_s": time_s,
                "record_count": count,
                "mean_radius_um": safe_ratio(values.get("radius_sum", 0.0) * 1.0e6, count),
                "mean_height_m": safe_ratio(values.get("height_sum", 0.0), values.get("height_count", 0.0)),
                "large_record_fraction": safe_ratio(values.get("large_count", 0.0), count),
                "coal_record_fraction": safe_ratio(values.get("coal_count", 0.0), count),
                "warnings": None,
            }
        )
    if not time_rows:
        time_rows = [{"time_s": None, "record_count": None, "mean_radius_um": None, "mean_height_m": None, "large_record_fraction": None, "coal_record_fraction": None, "warnings": warning_text(warnings)}]

    occurrence_rows = []
    for (time_s, height_center), values in sorted(by_time_height.items()):
        occurrence_rows.append(
            {
                "time_s": time_s,
                "height_bin_m": height_center,
                "target_count": values.get("target_count"),
                "large_target_count": values.get("large_target_count"),
                "ifcoal_target_count": values.get("ifcoal_target_count"),
                "warnings": None,
            }
        )
    if not occurrence_rows:
        occurrence_rows = [{"time_s": None, "height_bin_m": None, "target_count": None, "large_target_count": None, "ifcoal_target_count": None, "warnings": warning_text(warnings)}]

    if not target_summary_rows:
        target_summary_rows = [{column: None for column in TARGET_SUMMARY_COLUMNS}]
        target_summary_rows[0]["warnings"] = warning_text(warnings)

    if not example_rows:
        example_rows = [{column: None for column in EXAMPLE_COLUMNS}]
        example_rows[0]["warnings"] = warning_text(warnings)

    return summary_rows, pathway_rows, height_rows, coal_rows, time_rows, target_summary_rows, occurrence_rows, example_rows


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write TPHT threshold-crossing diagnostic tables."""
    options = options or {}
    config = load_config()
    bw_case = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}["bw_reconstruction"]
    _bw_row, bw_files, _records, warnings = collect_case_metrics(root, bw_case, options)
    selected_files = bw_files.get("sd_selected_output_bytes", [])
    coal_files = bw_files.get("coalescence_log_bytes", [])
    _log(f"root={root}")
    _log(f"outdir={outdir}")
    _log(f"case_dir={_bw_row.get('case_dir')}")
    _log(f"selected_files={len(selected_files)} coal_files={len(coal_files)} options={options}")
    targets, by_time, by_time_height, link_rows, science_warnings = build_stepwise_bw_diagnostics(selected_files, coal_files, options)
    warnings.extend(science_warnings)
    _log(f"stepwise targets={len(targets)} by_time={len(by_time)} by_time_height={len(by_time_height)} link_rows={len(link_rows)}")
    optional_sources, optional_warnings = build_optional_source_tables(selected_files, coal_files, targets, options) if targets else ({}, [])
    warnings.extend(optional_warnings)
    optional_counts = {name: len(rows) for name, rows in optional_sources.items()}
    _log(f"optional_source_rows={optional_counts}")
    if targets and coal_files and not optional_sources.get("event_links"):
        warnings.append(
            "No target-linked coalescence event rows were produced although BW selected-output and coalescence-log files were found; "
            "check whether the job used --metadata-only/--skip-heavy-netcdf/--max-files, whether coalescence-log times overlap selected-output intervals, "
            "and whether event pre_dmid/pre_sdid match the following selected-output level."
        )
    if not targets:
        warnings.append("No BW target chains were reconstructed; heavy TPHT science tables contain NA placeholders.")
    summary_rows, pathway_rows, height_rows, coal_rows, time_rows, target_summary_rows, occurrence_rows, example_rows = _build_rows(targets, by_time, by_time_height, warnings)
    if not link_rows:
        link_rows = [{column: None for column in LINK_COLUMNS}]
        link_rows[0]["warnings"] = warning_text(warnings)
    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    _emit_warnings(warnings)
    write_table_bundle(summary_rows, outdir / "tables" / "02_tpht_science_summary", SUMMARY_COLUMNS)
    write_table_bundle(pathway_rows, outdir / "tables" / "02_tpht_science_pathways", PATHWAY_COLUMNS)
    write_table_bundle(height_rows, outdir / "tables" / "02_tpht_science_formation_height_bins", HEIGHT_BIN_COLUMNS)
    write_table_bundle(coal_rows, outdir / "tables" / "02_tpht_science_coalescence_counts", COAL_COLUMNS)
    write_table_bundle(time_rows, outdir / "tables" / "02_tpht_science_time_series", TIME_COLUMNS)
    write_table_bundle(target_summary_rows, outdir / "tables" / "02_tpht_science_target_summary", TARGET_SUMMARY_COLUMNS)
    write_table_bundle(occurrence_rows, outdir / "tables" / "02_tpht_target_occurrence_zt", OCCURRENCE_ZT_COLUMNS)
    write_table_bundle(example_rows, outdir / "tables" / "02_tpht_science_example_targets", EXAMPLE_COLUMNS)
    write_table_bundle(link_rows, outdir / "tables" / "02_tpht_chain_links_by_time", LINK_COLUMNS)
    write_table_bundle(optional_sources.get("trajectory") or _missing_rows(TRAJECTORY_COLUMNS, warnings), outdir / "tables" / "02_tpht_target_trajectory_records", TRAJECTORY_COLUMNS)
    write_table_bundle(optional_sources.get("ifcoal_timeline") or _missing_rows(IFCOAL_TIMELINE_COLUMNS, warnings), outdir / "tables" / "02_tpht_target_ifcoal_timeline", IFCOAL_TIMELINE_COLUMNS)
    write_table_bundle(optional_sources.get("interval") or _missing_rows(INTERVAL_COLUMNS, warnings), outdir / "tables" / "02_tpht_target_interval_diagnostics", INTERVAL_COLUMNS)
    write_table_bundle(optional_sources.get("event_links") or _missing_rows(EVENT_LINK_COLUMNS, warnings), outdir / "tables" / "02_tpht_target_event_links", EVENT_LINK_COLUMNS)
    _log("finished writing TPHT science tables")


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
            [
                "tables/02_tpht_science_summary.{csv,md,tex,json}",
                "tables/02_tpht_science_pathways.{csv,md,tex,json}",
                "tables/02_tpht_science_formation_height_bins.{csv,md,tex,json}",
                "tables/02_tpht_science_coalescence_counts.{csv,md,tex,json}",
                "tables/02_tpht_science_time_series.{csv,md,tex,json}",
                "tables/02_tpht_science_target_summary.{csv,md,tex,json}",
                "tables/02_tpht_target_occurrence_zt.{csv,md,tex,json}",
                "tables/02_tpht_science_example_targets.{csv,md,tex,json}",
                "tables/02_tpht_chain_links_by_time.{csv,md,tex,json}",
                "tables/02_tpht_target_trajectory_records.{csv,md,tex,json}",
                "tables/02_tpht_target_ifcoal_timeline.{csv,md,tex,json}",
                "tables/02_tpht_target_interval_diagnostics.{csv,md,tex,json}",
                "tables/02_tpht_target_event_links.{csv,md,tex,json}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
