#!/usr/bin/env python3
"""Analyze TPHT target categories and target discovery timing."""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_02_tpht import build_summary  # noqa: E402
from common.parse_logs import collect_case_metrics  # noqa: E402
from common.parse_gmd_diag import parse_case_logs  # noqa: E402
from common.paths import analysis_options, build_arg_parser, case_path, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import read_csv_rows, safe_float, safe_ratio, write_table_bundle  # noqa: E402
from tpht_chain_utils import build_stepwise_bw_diagnostics  # noqa: E402


LINKAGE_NOTE = "Target categories are assigned from the first interest-condition record along stepwise BW predecessor-link chains."


def _is_earlier(time_s: float | None, order: float | None, ref_time: float | None, ref_order: float | None) -> bool:
    """Compare records by model time when available, otherwise by stream order."""
    if ref_order is None:
        return True
    if time_s is not None and ref_time is not None:
        return time_s < ref_time
    if order is None:
        return False
    return order < ref_order


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write target category and discovery-time tables."""
    config = load_config()
    cases = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}
    fw_dir = case_path(root, cases["fw_discovery"])
    options = options or {}
    summary, summary_warnings = build_summary(root, options)
    warnings = list(summary_warnings)

    cached_target_rows = read_csv_rows(outdir / "tables" / "02_tpht_science_target_summary.csv")
    counter: Counter[str] = Counter()
    if cached_target_rows and any(row.get("target_id") not in (None, "", "NA") for row in cached_target_rows):
        for row in cached_target_rows:
            category = row.get("first_selected_category") or row.get("category") or "unknown"
            counter[str(category)] += 1
    else:
        _bw_row, bw_files, _records, metric_warnings = collect_case_metrics(root, cases["bw_reconstruction"], options)
        targets, _by_time, _by_time_height, _link_rows, chain_warnings = build_stepwise_bw_diagnostics(
            bw_files.get("sd_selected_output_bytes", []),
            bw_files.get("coalescence_log_bytes", []),
            options,
        )
        warnings.extend(metric_warnings + chain_warnings)
        counter = Counter((state.get("first_selected_category") or "unknown") for state in targets.values())
    if not counter:
        counter["unknown"] = safe_float(summary.get("deduplicated_target_pairs")) or None

    total = sum(value for value in counter.values() if value is not None)
    category_rows = []
    for category in ("radius_only", "coal_only", "both", "unknown"):
        count = counter.get(category)
        category_rows.append(
            {
                "category": category,
                "target_count": count,
                "target_fraction": safe_ratio(count, total),
                "warnings": warning_text(warnings) if category == "unknown" and warnings else None,
            }
        )

    diag, diag_records, log_warnings = parse_case_logs(fw_dir)
    warnings.extend(log_warnings)
    discovery_records = [
        record for record in diag_records if record.get("record_type") == "TPHT_ID_DIAG" and safe_float(record.get("time_s")) is not None
    ]
    discovery_records.sort(key=lambda record: safe_float(record.get("time_s")) or 0.0)
    discovery_rows = []
    previous_total = 0.0
    first_detected = None
    last_detected = None
    for record in discovery_records:
        time_s = safe_float(record.get("time_s"))
        total_records = safe_float(record.get("tpht_id_records_written_total") or record.get("ids_loaded"))
        new_targets = total_records - previous_total if total_records is not None else None
        if total_records is not None:
            previous_total = total_records
        if new_targets is not None and new_targets > 0:
            first_detected = time_s if first_detected is None else first_detected
            last_detected = time_s
        discovery_rows.append(
            {
                "time_s": time_s,
                "new_targets_per_time": new_targets,
                "target_records_per_time": total_records,
                "target_first_detected_time_s": first_detected,
                "target_last_detected_time_s": last_detected,
                "first_hour_cap_region": time_s is not None and time_s <= 3600.0,
                "warnings": None,
            }
        )
    if not discovery_rows:
        discovery_rows = [
            {
                "time_s": None,
                "new_targets_per_time": None,
                "target_records_per_time": None,
                "target_first_detected_time_s": None,
                "target_last_detected_time_s": None,
                "first_hour_cap_region": None,
                "warnings": warning_text(warnings),
            }
        ]

    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle(category_rows, outdir / "tables" / "02_tpht_target_categories")
    write_table_bundle(discovery_rows, outdir / "tables" / "02_tpht_discovery_time")


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
            ["tables/02_tpht_target_categories.{csv,md,tex,json}", "tables/02_tpht_discovery_time.{csv,md,tex,json}"],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
