#!/usr/bin/env python3
"""Check raw TPHT BW selected-output time order and chain traces.

This diagnostic script reads the raw ``SD_selected_NetCDF_*`` files from
``02_tpht_3d_interest_70min/bw_reconstruction`` one output time at a time.  It
does not rely on the already generated TPHT science tables for the direction
test.  The goal is to verify whether the raw ``pre_dmid/pre_sdid`` fields are
more consistent with ascending model-output time, descending model-output time,
or both.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from common.paths import analysis_options, build_arg_parser, case_path, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import read_csv_rows, safe_float, safe_ratio, write_table_bundle  # noqa: E402
from tpht_chain_utils import _group_paths_by_time, _read_selected_level  # noqa: E402


SUMMARY_COLUMNS = [
    "check_name",
    "case_dir",
    "time_group_count",
    "raw_selected_file_count",
    "ascending_mean_valid_link_fraction",
    "ascending_min_valid_link_fraction",
    "ascending_unmatched_total",
    "descending_mean_valid_link_fraction",
    "descending_min_valid_link_fraction",
    "descending_unmatched_total",
    "direction_inference",
    "selected_target_count",
    "target_ids",
    "warnings",
]

LEVEL_COLUMNS = [
    "direction",
    "level_index",
    "time_s",
    "record_count",
    "matched_records",
    "unmatched_records",
    "valid_link_fraction",
    "unique_output_ids",
    "duplicate_output_ids",
    "warnings",
]

SAMPLE_COLUMNS = [
    "target_id",
    "time_s",
    "time_min",
    "dm_id",
    "sd_id",
    "radius_um",
    "height_m",
    "if_coal",
    "source_file",
    "direction_assumption",
    "warnings",
]

CHAIN_SUMMARY_COLUMNS = [
    "target_id",
    "record_count",
    "first_time_s",
    "last_time_s",
    "first_radius_um",
    "last_radius_um",
    "peak_radius_um",
    "peak_time_s",
    "first_ge_5um_time_s",
    "first_ge_15um_time_s",
    "peak_in_first_10min",
    "last_radius_ge_5um",
    "last_radius_ge_15um",
    "warnings",
]


def _limited_time_groups(groups: list[tuple[float | None, list[Path]]], max_time_groups: int | None) -> list[tuple[float | None, list[Path]]]:
    """Limit output time groups while preserving model-time ordering."""
    if max_time_groups is None or max_time_groups < 0 or len(groups) <= max_time_groups:
        return groups
    return groups[:max_time_groups]


def _target_sort_key(row: dict[str, str]) -> tuple[float, float, float]:
    """Sort target-summary rows by diagnostic size and record count."""
    return (
        safe_float(row.get("max_radius_um")) or -math.inf,
        safe_float(row.get("record_count")) or -math.inf,
        -(safe_float(row.get("first_large_time_s")) or math.inf),
    )


def _choose_target_ids(outdir: Path, explicit_ids: list[str], top_n: int) -> tuple[list[str], list[str]]:
    """Choose chain IDs to trace from an existing science summary table."""
    if explicit_ids:
        return list(dict.fromkeys(explicit_ids)), []
    rows = read_csv_rows(outdir / "tables" / "02_tpht_science_target_summary.csv")
    candidates = [row for row in rows if row.get("target_id") not in (None, "", "NA")]
    if not candidates:
        return [], ["No target IDs were provided and 02_tpht_science_target_summary.csv was not available; first chains will be sampled"]
    ranked = sorted(candidates, key=_target_sort_key, reverse=True)
    return [str(row["target_id"]) for row in ranked[: max(1, top_n)]], []


def _direction_rows(
    groups: list[tuple[float | None, list[Path]]],
    direction: str,
    options: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Scan raw selected-output levels in one direction and summarize ID links."""
    previous_ids: set[tuple[int, int]] = set()
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    ordered_groups = groups if direction == "ascending" else list(reversed(groups))
    for level_index, (time_s, paths) in enumerate(ordered_groups):
        records, level_warnings = _read_selected_level(paths, time_s, options)
        warnings.extend(level_warnings)
        current_ids: set[tuple[int, int]] = set()
        duplicate_count = 0
        matched_count = 0
        for record in records:
            output_id = record["output_id"]
            if output_id in current_ids:
                duplicate_count += 1
            else:
                current_ids.add(output_id)
            if level_index > 0 and output_id in previous_ids:
                matched_count += 1
        unmatched_count = (len(records) - matched_count) if level_index > 0 else None
        rows.append(
            {
                "direction": direction,
                "level_index": level_index,
                "time_s": time_s,
                "record_count": len(records) if records else None,
                "matched_records": matched_count if level_index > 0 else None,
                "unmatched_records": unmatched_count,
                "valid_link_fraction": safe_ratio(matched_count, len(records)) if level_index > 0 else None,
                "unique_output_ids": len(current_ids) if current_ids else None,
                "duplicate_output_ids": duplicate_count if duplicate_count else 0,
                "warnings": warning_text(level_warnings),
            }
        )
        previous_ids = current_ids
    return rows, warnings


def _summarize_direction(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Return aggregate link metrics for one direction."""
    fractions = [safe_float(row.get("valid_link_fraction")) for row in rows]
    fractions = [value for value in fractions if value is not None]
    unmatched = [safe_float(row.get("unmatched_records")) for row in rows]
    unmatched = [value for value in unmatched if value is not None]
    return {
        "mean_valid": sum(fractions) / len(fractions) if fractions else None,
        "min_valid": min(fractions) if fractions else None,
        "unmatched_total": sum(unmatched) if unmatched else None,
    }


def _direction_inference(ascending: dict[str, Any], descending: dict[str, Any]) -> str:
    """Describe which raw-file ordering is better supported by the ID links."""
    asc = safe_float(ascending.get("mean_valid"))
    desc = safe_float(descending.get("mean_valid"))
    if asc is None and desc is None:
        return "no_direction_test_available"
    if asc is not None and desc is not None and abs(asc - desc) < 1.0e-9:
        return "both_directions_equally_match_ids_check_radius_traces_and_model_semantics"
    if asc is not None and (desc is None or asc > desc):
        return "ascending_model_time_has_higher_link_consistency"
    return "descending_model_time_has_higher_link_consistency"


def _trace_targets(
    groups: list[tuple[float | None, list[Path]]],
    target_ids: list[str],
    top_n: int,
    options: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    """Trace selected chain IDs in ascending model-output time."""
    warnings: list[str] = []
    selected_ids = list(dict.fromkeys(target_ids))
    auto_select = not selected_ids
    previous_id_to_chain: dict[tuple[int, int], str] = {}
    next_chain_index = 0
    sample_rows: list[dict[str, Any]] = []

    for _level_index, (time_s, paths) in enumerate(groups):
        records, level_warnings = _read_selected_level(paths, time_s, options)
        warnings.extend(level_warnings)
        current_id_to_chain: dict[tuple[int, int], str] = {}
        for record in records:
            output_id = record["output_id"]
            chain_id = previous_id_to_chain.get(output_id)
            if chain_id is None:
                chain_id = f"chain_{next_chain_index:08d}"
                next_chain_index += 1
            current_id_to_chain.setdefault(output_id, chain_id)
            if auto_select and len(selected_ids) < max(1, top_n):
                selected_ids.append(chain_id)
            if chain_id in selected_ids:
                radius_m = record.get("radius_m")
                sample_rows.append(
                    {
                        "target_id": chain_id,
                        "time_s": record.get("time_s"),
                        "time_min": (record.get("time_s") / 60.0) if record.get("time_s") is not None else None,
                        "dm_id": record.get("dm_id"),
                        "sd_id": record.get("sd_id"),
                        "radius_um": radius_m * 1.0e6 if radius_m is not None else None,
                        "height_m": record.get("height_m"),
                        "if_coal": record.get("if_coal"),
                        "source_file": record.get("source_file"),
                        "direction_assumption": "ascending_model_output_time",
                        "warnings": None,
                    }
                )
        previous_id_to_chain = current_id_to_chain
    if not sample_rows:
        warnings.append("No raw chain samples were collected for the selected target IDs")
    return sample_rows, selected_ids, warnings


def _chain_summaries(sample_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarize raw traced chain samples."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in sample_rows:
        target_id = row.get("target_id")
        if target_id not in (None, "", "NA"):
            grouped.setdefault(str(target_id), []).append(row)
    summaries: list[dict[str, Any]] = []
    for target_id, rows in sorted(grouped.items()):
        valid = [row for row in rows if safe_float(row.get("time_s")) is not None]
        valid.sort(key=lambda row: safe_float(row.get("time_s")) or -math.inf)
        radius_rows = [row for row in valid if safe_float(row.get("radius_um")) is not None]
        peak_row = max(radius_rows, key=lambda row: safe_float(row.get("radius_um")) or -math.inf) if radius_rows else {}
        first_ge_5 = next((row for row in valid if (safe_float(row.get("radius_um")) or -math.inf) >= 5.0), {})
        first_ge_15 = next((row for row in valid if (safe_float(row.get("radius_um")) or -math.inf) >= 15.0), {})
        first_row = valid[0] if valid else {}
        last_row = valid[-1] if valid else {}
        peak_time = safe_float(peak_row.get("time_s"))
        last_radius = safe_float(last_row.get("radius_um"))
        summaries.append(
            {
                "target_id": target_id,
                "record_count": len(valid) if valid else None,
                "first_time_s": safe_float(first_row.get("time_s")),
                "last_time_s": safe_float(last_row.get("time_s")),
                "first_radius_um": safe_float(first_row.get("radius_um")),
                "last_radius_um": last_radius,
                "peak_radius_um": safe_float(peak_row.get("radius_um")),
                "peak_time_s": peak_time,
                "first_ge_5um_time_s": safe_float(first_ge_5.get("time_s")),
                "first_ge_15um_time_s": safe_float(first_ge_15.get("time_s")),
                "peak_in_first_10min": (peak_time <= 600.0) if peak_time is not None else None,
                "last_radius_ge_5um": (last_radius >= 5.0) if last_radius is not None else None,
                "last_radius_ge_15um": (last_radius >= 15.0) if last_radius is not None else None,
                "warnings": None,
            }
        )
    return summaries or [{column: None for column in CHAIN_SUMMARY_COLUMNS}]


def _resolve_bw_case_dir(root: Path, case_dir: Path | None) -> Path:
    """Resolve the BW reconstruction directory."""
    if case_dir is not None:
        return case_dir.resolve()
    config = load_config()
    cases = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}
    return case_path(root, cases["bw_reconstruction"])


def analyze(root: Path, outdir: Path, args: Any, options: dict[str, Any]) -> None:
    """Run the raw time-order check and write diagnostic tables."""
    warnings: list[str] = []
    bw_dir = _resolve_bw_case_dir(root, args.case_dir)
    selected_files = sorted(path for path in bw_dir.rglob("SD_selected_NetCDF_*") if path.is_file())
    if not selected_files:
        warnings.append(f"No raw SD_selected_NetCDF_* files found below {bw_dir}")
    groups, group_warnings = _group_paths_by_time(selected_files, options.get("max_files"))
    warnings.extend(group_warnings)
    groups = _limited_time_groups(groups, args.max_time_groups)

    ascending_rows, ascending_warnings = _direction_rows(groups, "ascending", options) if groups else ([], [])
    descending_rows, descending_warnings = _direction_rows(groups, "descending", options) if groups else ([], [])
    warnings.extend(ascending_warnings)
    warnings.extend(descending_warnings)
    ascending_summary = _summarize_direction(ascending_rows)
    descending_summary = _summarize_direction(descending_rows)

    target_ids, target_warnings = _choose_target_ids(outdir, args.target_id, args.top_n)
    warnings.extend(target_warnings)
    sample_rows, selected_ids, sample_warnings = _trace_targets(groups, target_ids, args.top_n, options) if groups else ([], target_ids, [])
    warnings.extend(sample_warnings)
    chain_summary_rows = _chain_summaries(sample_rows)

    summary_row = {
        "check_name": "02_tpht_raw_time_order",
        "case_dir": bw_dir,
        "time_group_count": len(groups) if groups else None,
        "raw_selected_file_count": len(selected_files) if selected_files else None,
        "ascending_mean_valid_link_fraction": ascending_summary.get("mean_valid"),
        "ascending_min_valid_link_fraction": ascending_summary.get("min_valid"),
        "ascending_unmatched_total": ascending_summary.get("unmatched_total"),
        "descending_mean_valid_link_fraction": descending_summary.get("mean_valid"),
        "descending_min_valid_link_fraction": descending_summary.get("min_valid"),
        "descending_unmatched_total": descending_summary.get("unmatched_total"),
        "direction_inference": _direction_inference(ascending_summary, descending_summary),
        "selected_target_count": len(selected_ids) if selected_ids else None,
        "target_ids": ";".join(selected_ids) if selected_ids else None,
        "warnings": warning_text(warnings),
    }
    level_rows = ascending_rows + descending_rows
    if not level_rows:
        level_rows = [{column: None for column in LEVEL_COLUMNS}]
        level_rows[0]["warnings"] = warning_text(warnings)
    if not sample_rows:
        sample_rows = [{column: None for column in SAMPLE_COLUMNS}]
        sample_rows[0]["warnings"] = warning_text(warnings)
    write_table_bundle([summary_row], outdir / "tables" / "02_tpht_raw_time_order_check", SUMMARY_COLUMNS)
    write_table_bundle(level_rows, outdir / "tables" / "02_tpht_raw_time_order_by_level", LEVEL_COLUMNS)
    write_table_bundle(sample_rows, outdir / "tables" / "02_tpht_raw_chain_samples", SAMPLE_COLUMNS)
    write_table_bundle(chain_summary_rows, outdir / "tables" / "02_tpht_raw_chain_sample_summary", CHAIN_SUMMARY_COLUMNS)
    if args.strict and warnings:
        raise RuntimeError("; ".join(warnings))

    print("02_tpht_raw_time_order_check")
    print(f"case_dir={bw_dir}")
    print(f"time_group_count={summary_row['time_group_count']}")
    print(f"ascending_mean_valid_link_fraction={summary_row['ascending_mean_valid_link_fraction']}")
    print(f"descending_mean_valid_link_fraction={summary_row['descending_mean_valid_link_fraction']}")
    print(f"direction_inference={summary_row['direction_inference']}")
    print(f"target_ids={summary_row['target_ids']}")


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    parser.add_argument("--case-dir", type=Path, default=None, help="Optional explicit path to the bw_reconstruction directory.")
    parser.add_argument("--target-id", action="append", default=[], help="Target chain ID to trace, e.g. chain_00072917. May be repeated.")
    parser.add_argument("--top-n", type=int, default=12, help="Number of top-radius targets to trace when --target-id is not provided.")
    parser.add_argument("--max-time-groups", type=int, default=None, help="Debug limit on output time groups after grouping per-rank files.")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(
            Path(__file__).name,
            root,
            outdir,
            [
                "tables/02_tpht_raw_time_order_check.{csv,md,tex,json}",
                "tables/02_tpht_raw_time_order_by_level.{csv,md,tex,json}",
                "tables/02_tpht_raw_chain_samples.{csv,md,tex,json}",
                "tables/02_tpht_raw_chain_sample_summary.{csv,md,tex,json}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args, analysis_options(args))


if __name__ == "__main__":
    main()
