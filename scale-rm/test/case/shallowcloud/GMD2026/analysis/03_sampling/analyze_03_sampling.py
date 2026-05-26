#!/usr/bin/env python3
"""Analyze the short 2D GMD2026 sampling-procedure verification."""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.parse_netcdf import distribution_summary, read_radius_records  # noqa: E402
from common.paths import analysis_options, build_arg_parser, case_path, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import safe_float, safe_ratio, write_table_bundle  # noqa: E402


COLUMNS = [
    "case_name",
    "sample_mode",
    "sampling_fraction",
    "seed",
    "time_s",
    "selected_record_count",
    "weighted_mean_radius",
    "weighted_median_radius",
    "radius_q25",
    "radius_q75",
    "radius_q95",
    "max_radius",
    "fraction_r_ge_10um",
    "fraction_r_ge_15um",
    "fraction_r_ge_20um",
    "dsd_l1_error_vs_full_reference",
    "dsd_l2_error_vs_full_reference",
    "zr_hist_l1_error_vs_full_reference",
    "if_coal_fraction",
    "postprocess_wallclock_s",
    "warnings",
]


def _mode_for_case(case: dict) -> str:
    """Infer sampling mode from metadata or case name."""
    if case.get("sample_mode"):
        return str(case["sample_mode"])
    if "stratified" in case["case_name"]:
        return "stratified"
    return "random"


def _histogram(values: list[float], bins: list[float]) -> list[float]:
    """Build a normalized histogram without depending on NumPy."""
    if not values:
        return [0.0] * (len(bins) - 1)
    counts = [0.0] * (len(bins) - 1)
    for value in values:
        for index in range(len(bins) - 1):
            if bins[index] <= value < bins[index + 1] or (index == len(bins) - 2 and value == bins[index + 1]):
                counts[index] += 1.0
                break
    total = sum(counts)
    return [count / total for count in counts] if total else counts


def _histogram_errors(values: list[float], reference: list[float]) -> tuple[float | None, float | None]:
    """Return L1 and L2 errors against a full-reference radius histogram."""
    if not values or not reference:
        return None, None
    maximum = max(max(values), max(reference), 25.0e-6)
    bins = [index * maximum / 30.0 for index in range(31)]
    hist = _histogram(values, bins)
    ref_hist = _histogram(reference, bins)
    diffs = [left - right for left, right in zip(hist, ref_hist)]
    l1_error = sum(abs(value) for value in diffs)
    l2_error = math.sqrt(sum(value * value for value in diffs))
    return l1_error, l2_error


def _rows_for_case(case: dict, root: Path, reference_values: list[float], elapsed_s: float, options: dict) -> tuple[list[dict], list[str]]:
    """Compute sampling rows for one case."""
    case_dir = case_path(root, case)
    records, warnings = read_radius_records(
        case_dir,
        max_files=options.get("max_files"),
        max_records=options.get("max_records"),
        chunk_size=options.get("chunk_size", 100000),
        metadata_only=options.get("metadata_only", False),
        skip_heavy_netcdf=options.get("skip_heavy_netcdf", False),
    )
    by_time: dict[float | None, list[dict]] = defaultdict(list)
    for record in records:
        by_time[safe_float(record.get("time_s"))].append(record)
    if not by_time:
        by_time[None] = []

    rows: list[dict] = []
    for time_s, records_at_time in sorted(by_time.items(), key=lambda item: -1.0 if item[0] is None else item[0]):
        radii = [safe_float(record.get("radius_m")) for record in records_at_time]
        radii = [value for value in radii if value is not None and value > 0.0]
        coal_flags = [safe_float(record.get("if_coal")) for record in records_at_time]
        coal_flags = [value for value in coal_flags if value is not None]
        summary = distribution_summary(radii)
        l1_error, l2_error = _histogram_errors(radii, reference_values)
        rows.append(
            {
                "case_name": case["case_name"],
                "sample_mode": _mode_for_case(case),
                "sampling_fraction": case.get("sampling_fraction"),
                "seed": case.get("seed"),
                "time_s": time_s,
                "selected_record_count": summary["count"],
                "weighted_mean_radius": summary["mean"],
                "weighted_median_radius": summary["median"],
                "radius_q25": summary["q25"],
                "radius_q75": summary["q75"],
                "radius_q95": summary["q95"],
                "max_radius": summary["maximum"],
                "fraction_r_ge_10um": safe_ratio(sum(1 for value in radii if value >= 10.0e-6), len(radii)),
                "fraction_r_ge_15um": safe_ratio(sum(1 for value in radii if value >= 15.0e-6), len(radii)),
                "fraction_r_ge_20um": safe_ratio(sum(1 for value in radii if value >= 20.0e-6), len(radii)),
                "dsd_l1_error_vs_full_reference": l1_error,
                "dsd_l2_error_vs_full_reference": l2_error,
                "zr_hist_l1_error_vs_full_reference": None,
                "if_coal_fraction": safe_ratio(sum(1 for value in coal_flags if value > 0.0), len(coal_flags)),
                "postprocess_wallclock_s": elapsed_s,
                "warnings": warning_text(warnings),
            }
        )
    return rows, warnings


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write sampling verification metrics and seed statistics."""
    import time

    start = time.time()
    options = options or {}
    config = load_config()
    cases = iter_cases(config, "03_fw_rep_2d_600s")
    if options.get("quick"):
        cases = [case for case in cases if case["case_name"] in {"random_full", "stratified_full", "random_f005_s00", "stratified_f005_s00"}]
    warnings: list[str] = []
    reference_by_mode: dict[str, list[float]] = {}
    for case in cases:
        if case["case_name"] not in {"random_full", "stratified_full"}:
            continue
        records, current_warnings = read_radius_records(
            case_path(root, case),
            max_files=options.get("max_files"),
            max_records=options.get("max_records"),
            chunk_size=options.get("chunk_size", 100000),
            metadata_only=options.get("metadata_only", False),
            skip_heavy_netcdf=options.get("skip_heavy_netcdf", False),
        )
        warnings.extend(current_warnings)
        reference_by_mode[_mode_for_case(case)] = [
            value for value in (safe_float(record.get("radius_m")) for record in records) if value is not None and value > 0.0
        ]
    if not reference_by_mode.get("random") and reference_by_mode.get("stratified"):
        reference_by_mode["random"] = reference_by_mode["stratified"]
        warnings.append("random_full missing; stratified_full used as the random-sampling reference")
    if not reference_by_mode.get("stratified") and reference_by_mode.get("random"):
        reference_by_mode["stratified"] = reference_by_mode["random"]
        warnings.append("stratified_full missing; random_full used as the stratified-sampling reference")

    rows: list[dict] = []
    for case in cases:
        row, _files, _records, metric_warnings = collect_case_metrics(root, case, options)
        case_warnings = list(metric_warnings)
        reference = reference_by_mode.get(_mode_for_case(case), [])
        case_rows, radius_warnings = _rows_for_case(case, root, reference, time.time() - start, options)
        case_warnings.extend(radius_warnings)
        for case_row in case_rows:
            if case_warnings and not case_row.get("warnings"):
                case_row["warnings"] = warning_text(case_warnings)
            rows.append(case_row)
        warnings.extend(case_warnings)

    seed_groups: dict[tuple[str, float], list[dict]] = defaultdict(list)
    for row in rows:
        seed = row.get("seed")
        fraction = safe_float(row.get("sampling_fraction"))
        if seed is None or fraction is None or fraction >= 1.0:
            continue
        seed_groups[(str(row.get("sample_mode")), fraction)].append(row)
    seed_rows = []
    for (mode, fraction), group_rows in sorted(seed_groups.items()):
        seeds = sorted({int(row["seed"]) for row in group_rows if row.get("seed") is not None})
        l1_values = [safe_float(row.get("dsd_l1_error_vs_full_reference")) for row in group_rows]
        l1_values = [value for value in l1_values if value is not None]
        note = None
        if len(seeds) != 10:
            note = f"Expected 10 seeds, found {len(seeds)}"
            warnings.append(f"{mode} fraction {fraction}: {note}")
        mean_l1 = sum(l1_values) / len(l1_values) if l1_values else None
        seed_rows.append(
            {
                "sample_mode": mode,
                "sampling_fraction": fraction,
                "seed_count": len(seeds),
                "expected_seed_count": 10,
                "mean_dsd_l1_error_vs_full_reference": mean_l1,
                "missing_seed_warning": note,
            }
        )

    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle(rows, outdir / "tables" / "03_sampling_metrics", COLUMNS)
    write_table_bundle(seed_rows, outdir / "tables" / "03_sampling_seed_statistics")


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
            ["tables/03_sampling_metrics.{csv,md,tex,json}", "tables/03_sampling_seed_statistics.{csv,md,tex,json}"],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
