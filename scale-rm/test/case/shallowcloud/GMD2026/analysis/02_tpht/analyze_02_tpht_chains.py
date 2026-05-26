#!/usr/bin/env python3
"""Analyze TPHT chain-validity proxies and target-history diagnostics."""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.parse_netcdf import inspect_netcdf_files, read_selected_pairs  # noqa: E402
from common.paths import analysis_options, build_arg_parser, case_path, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import read_csv_rows, safe_float, safe_ratio, write_table_bundle  # noqa: E402
from tpht_chain_utils import LINK_COLUMNS, build_stepwise_bw_diagnostics, chain_lengths, parse_chain_distribution_text  # noqa: E402


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write chain-validity and target-history tables."""
    options = options or {}
    config = load_config()
    bw_case = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}["bw_reconstruction"]

    cached_link_rows = read_csv_rows(outdir / "tables" / "02_tpht_chain_links_by_time.csv")
    cached_time_rows = read_csv_rows(outdir / "tables" / "02_tpht_science_time_series.csv")
    if cached_link_rows and cached_time_rows and any(row.get("time_s") not in (None, "", "NA") for row in cached_link_rows):
        summary_row = (read_csv_rows(outdir / "tables" / "02_tpht_summary.csv") or [{}])[0]
        times = [safe_float(row.get("time_s")) for row in cached_link_rows]
        times = [time for time in times if time is not None]
        coverage = max(times) - min(times) if times else None
        matched_total = sum(int(safe_float(row.get("matched_parent_records")) or 0) for row in cached_link_rows)
        link_record_total = sum(int(safe_float(row.get("record_count")) or 0) for row in cached_link_rows if safe_float(row.get("matched_parent_records")) is not None)
        stepwise_valid_fraction = safe_ratio(matched_total, link_record_total)
        stepwise_invalid_fraction = 1.0 - stepwise_valid_fraction if stepwise_valid_fraction is not None else None
        chain_row = {
            "case_name": "bw_reconstruction",
            "valid_predecessor_fraction": stepwise_valid_fraction,
            "invalid_predecessor_fraction": stepwise_invalid_fraction,
            "chain_length_distribution": None,
            "reconstructed_time_coverage": coverage,
            "target_survival_fraction_by_time": None,
            "valid_link_fraction_by_time": stepwise_valid_fraction,
            "bw_valid_record_rows": summary_row.get("bw_valid_record_rows"),
            "bw_valid_records": summary_row.get("bw_valid_records"),
            "bw_unique_pairs": summary_row.get("bw_unique_pairs"),
            "warnings": "reused stepwise BW chain cache from 02_tpht_science outputs",
        }
        history_rows = []
        for row in cached_time_rows:
            mean_radius_um = safe_float(row.get("mean_radius_um"))
            history_rows.append(
                {
                    "time_s": safe_float(row.get("time_s")),
                    "record_count": safe_float(row.get("record_count")),
                    "mean_radius_history_m": mean_radius_um * 1.0e-6 if mean_radius_um is not None else None,
                    "median_radius_history_m": None,
                    "radius_q25_m": None,
                    "radius_q75_m": None,
                    "mean_height_history_m": safe_float(row.get("mean_height_m")),
                    "median_height_history_m": None,
                    "height_q25_m": None,
                    "height_q75_m": None,
                    "diagnostic_label": "TPHT diagnostic capability demonstration",
                    "warnings": "reused memory-safe means from stepwise BW science cache",
                }
            )
        write_table_bundle([chain_row], outdir / "tables" / "02_tpht_chain_validity")
        write_table_bundle(history_rows, outdir / "tables" / "02_tpht_target_histories")
        write_table_bundle(cached_link_rows, outdir / "tables" / "02_tpht_chain_links_by_time", LINK_COLUMNS)
        return

    bw_row, bw_files, _records, warnings = collect_case_metrics(root, bw_case, options)
    bw_dir = case_path(root, bw_case)
    selected_files = bw_files.get("sd_selected_output_bytes", [])
    pairs, valid_record_rows, pair_warnings = read_selected_pairs(
        selected_files,
        "BW",
        max_files=options.get("max_files"),
        max_records=options.get("max_records"),
        chunk_size=options.get("chunk_size", 100000),
        metadata_only=options.get("metadata_only", False) or options.get("skip_heavy_netcdf", False),
    )
    warnings.extend(pair_warnings)
    netcdf_stats, netcdf_warnings = inspect_netcdf_files(
        selected_files,
        max_files=options.get("max_files"),
        metadata_only=options.get("metadata_only", False),
    )
    warnings.extend(netcdf_warnings)
    selected_count = safe_float(netcdf_stats.get("selected_sd_record_count"))
    valid_fraction = safe_ratio(valid_record_rows, selected_count)
    invalid_fraction = 1.0 - valid_fraction if valid_fraction is not None else None
    unique_pair_count = len(pairs) if valid_record_rows is not None else None

    targets, by_time, _by_time_height, link_rows, chain_warnings = build_stepwise_bw_diagnostics(
        selected_files,
        bw_files.get("coalescence_log_bytes", []),
        options,
    )
    warnings.extend(chain_warnings)
    times = [safe_float(row.get("time_s")) for row in link_rows]
    times = [time for time in times if time is not None]
    coverage = max(times) - min(times) if times else None
    matched_total = sum(int(row.get("matched_parent_records") or 0) for row in link_rows)
    link_record_total = sum(int(row.get("record_count") or 0) for row in link_rows if row.get("matched_parent_records") is not None)
    stepwise_valid_fraction = safe_ratio(matched_total, link_record_total)
    stepwise_invalid_fraction = 1.0 - stepwise_valid_fraction if stepwise_valid_fraction is not None else None
    compact_lengths = parse_chain_distribution_text(chain_lengths(targets))

    chain_row = {
        "case_name": "bw_reconstruction",
        "valid_predecessor_fraction": stepwise_valid_fraction if stepwise_valid_fraction is not None else valid_fraction,
        "invalid_predecessor_fraction": stepwise_invalid_fraction if stepwise_invalid_fraction is not None else invalid_fraction,
        "chain_length_distribution": compact_lengths,
        "reconstructed_time_coverage": coverage,
        "target_survival_fraction_by_time": None,
        "valid_link_fraction_by_time": stepwise_valid_fraction,
        "bw_valid_record_rows": valid_record_rows,
        "bw_valid_records": unique_pair_count,
        "bw_unique_pairs": unique_pair_count,
        "warnings": warning_text(warnings),
    }

    history_rows = []
    for time_s, values in sorted(by_time.items(), key=lambda item: item[0]):
        count = safe_float(values.get("count"))
        height_count = safe_float(values.get("height_count"))
        history_rows.append(
            {
                "time_s": time_s,
                "record_count": count,
                "mean_radius_history_m": safe_ratio(values.get("radius_sum"), count),
                "median_radius_history_m": None,
                "radius_q25_m": None,
                "radius_q75_m": None,
                "mean_height_history_m": safe_ratio(values.get("height_sum"), height_count),
                "median_height_history_m": None,
                "height_q25_m": None,
                "height_q75_m": None,
                "diagnostic_label": "TPHT diagnostic capability demonstration",
                "warnings": "stepwise BW chain reconstruction; percentile histories are not computed in memory-safe mode",
            }
        )
    if not history_rows:
        history_rows = [
            {
                "time_s": None,
                "record_count": None,
                "median_radius_history_m": None,
                "radius_q25_m": None,
                "radius_q75_m": None,
                "median_height_history_m": None,
                "height_q25_m": None,
                "height_q75_m": None,
                "diagnostic_label": "TPHT diagnostic capability demonstration",
                "warnings": warning_text(warnings),
            }
        ]

    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle([chain_row], outdir / "tables" / "02_tpht_chain_validity")
    write_table_bundle(history_rows, outdir / "tables" / "02_tpht_target_histories")
    if not link_rows:
        link_rows = [{column: None for column in LINK_COLUMNS}]
        link_rows[0]["warnings"] = warning_text(warnings)
    write_table_bundle(link_rows, outdir / "tables" / "02_tpht_chain_links_by_time", LINK_COLUMNS)


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
                "tables/02_tpht_chain_validity.{csv,md,tex,json}",
                "tables/02_tpht_target_histories.{csv,md,tex,json}",
                "tables/02_tpht_chain_links_by_time.{csv,md,tex,json}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
