#!/usr/bin/env python3
"""Summarize FW, BW, and TPHT smoke-test outputs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any

from common.file_stats import collect_file_stats
from common.parse_ids import summarize_id_files
from common.parse_netcdf import inspect_netcdf_files, read_selected_pairs, validate_selected_output
from common.table_utils import write_table_bundle
from common.units import bytes_to_gib

DEFAULT_SMOKE_ROOT = Path(__file__).resolve().parents[1] / "codex_smoke"
CHAIN_COUNT_RE = re.compile(r"tracking_chain_count=\s*([0-9]+)")
DEFAULT_COLUMNS = [
    "smoke",
    "mode",
    "case_dir",
    "status",
    "selected_output_has_required_ids",
    "selected_output_id_variables",
    "selected_file_count",
    "latest_selected_valid_records",
    "latest_selected_unique_pairs",
    "tracking_chain_count_last",
    "tpht_fw_unique_pairs",
    "tpht_missing_in_bw",
    "tpht_extra_in_bw",
    "scientific_output_gib",
    "warnings",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke-root", type=Path, default=DEFAULT_SMOKE_ROOT)
    parser.add_argument("--fw-case", type=Path, default=None)
    parser.add_argument("--bw-case", type=Path, default=None)
    parser.add_argument("--tpht-fw-case", type=Path, default=None)
    parser.add_argument("--tpht-bw-case", type=Path, default=None)
    parser.add_argument("--tpht-fw-ids", type=Path, default=None)
    parser.add_argument("--output-base", type=Path, default=None)
    return parser.parse_args()


def resolve_cases(args: argparse.Namespace) -> dict[str, Path]:
    smoke_root = args.smoke_root.resolve()
    return {
        "fw": (args.fw_case or smoke_root / "fw").resolve(),
        "bw": (args.bw_case or smoke_root / "bw").resolve(),
        "tpht_fw": (args.tpht_fw_case or smoke_root / "tpht_fw").resolve(),
        "tpht_bw": (args.tpht_bw_case or smoke_root / "tpht_bw").resolve(),
    }


def parse_chain_counts(case_dir: Path) -> list[int]:
    """Read tracking-chain counts from SCALE log files."""
    counts: list[int] = []
    for log_path in sorted(case_dir.glob("LOG.pe*")):
        text = log_path.read_text(errors="ignore")
        counts.extend(int(match.group(1)) for match in CHAIN_COUNT_RE.finditer(text))
    return counts


def summarize_directional_case(name: str, mode: str, case_dir: Path) -> tuple[dict[str, Any], list[str]]:
    """Summarize one ordinary FW or BW selected-output smoke case."""
    stats, files_by_metric, warnings = collect_file_stats(case_dir)
    selected_files = files_by_metric.get("sd_selected_output_bytes", [])
    validation, validation_warnings = validate_selected_output(selected_files, mode)
    warnings.extend(validation_warnings)
    pairs, valid_records, pair_warnings = read_selected_pairs(selected_files, mode, latest_group_only=True)
    warnings.extend(pair_warnings)
    inspect_stats, inspect_warnings = inspect_netcdf_files(selected_files, max_files=1)
    warnings.extend(inspect_warnings)
    chain_counts = parse_chain_counts(case_dir)

    status = "pass"
    if validation.get("selected_output_has_required_ids") is not True or not pairs:
        status = "fail"

    row = {
        "smoke": name,
        "mode": mode,
        "case_dir": case_dir,
        "status": status,
        "selected_output_has_required_ids": validation.get("selected_output_has_required_ids"),
        "selected_output_id_variables": validation.get("selected_output_id_variables"),
        "selected_file_count": stats.get("sd_selected_file_count"),
        "latest_selected_valid_records": valid_records,
        "latest_selected_unique_pairs": len(pairs),
        "tracking_chain_count_last": chain_counts[-1] if chain_counts else None,
        "tpht_fw_unique_pairs": None,
        "tpht_missing_in_bw": None,
        "tpht_extra_in_bw": None,
        "scientific_output_gib": bytes_to_gib(stats.get("scientific_output_bytes")),
        "netcdf_dimensions": inspect_stats.get("netcdf_dimensions"),
        "warnings": "; ".join(warnings),
    }
    return row, warnings


def summarize_tpht_case(tpht_fw_case: Path, tpht_bw_case: Path, tpht_fw_ids: Path | None) -> tuple[dict[str, Any], list[str]]:
    """Summarize the TPHT FW-ID to BW-selected-output handoff."""
    warnings: list[str] = []
    stats_bw, files_by_metric_bw, stat_warnings_bw = collect_file_stats(tpht_bw_case)
    warnings.extend(stat_warnings_bw)
    selected_files_bw = files_by_metric_bw.get("sd_selected_output_bytes", [])

    id_path = tpht_fw_ids or tpht_fw_case / "fw_tracking" / "tracking_interest_ids_merged.ids"
    id_summary, fw_pairs_list, id_warnings = summarize_id_files([id_path])
    warnings.extend(id_warnings)
    fw_pairs = set(fw_pairs_list)

    validation, validation_warnings = validate_selected_output(selected_files_bw, "BW")
    warnings.extend(validation_warnings)
    bw_pairs, valid_records, pair_warnings = read_selected_pairs(selected_files_bw, "BW", latest_group_only=True)
    warnings.extend(pair_warnings)
    chain_counts = parse_chain_counts(tpht_bw_case)

    missing_in_bw = fw_pairs - bw_pairs
    extra_in_bw = bw_pairs - fw_pairs
    status = "pass"
    if not fw_pairs or missing_in_bw or extra_in_bw or validation.get("selected_output_has_required_ids") is not True:
        status = "fail"

    row = {
        "smoke": "tpht_handoff",
        "mode": "TPHT",
        "case_dir": tpht_bw_case,
        "status": status,
        "selected_output_has_required_ids": validation.get("selected_output_has_required_ids"),
        "selected_output_id_variables": validation.get("selected_output_id_variables"),
        "selected_file_count": stats_bw.get("sd_selected_file_count"),
        "latest_selected_valid_records": valid_records,
        "latest_selected_unique_pairs": len(bw_pairs),
        "tracking_chain_count_last": chain_counts[-1] if chain_counts else None,
        "tpht_fw_unique_pairs": id_summary.get("unique_pairs"),
        "tpht_missing_in_bw": len(missing_in_bw),
        "tpht_extra_in_bw": len(extra_in_bw),
        "scientific_output_gib": bytes_to_gib(stats_bw.get("scientific_output_bytes")),
        "warnings": "; ".join(warnings),
    }
    return row, warnings


def main() -> int:
    args = parse_args()
    cases = resolve_cases(args)
    output_base = args.output_base or args.smoke_root.resolve() / "postprocess" / "tracking_smoke_summary"

    rows: list[dict[str, Any]] = []
    rows.append(summarize_directional_case("fw_sampling", "FW", cases["fw"])[0])
    rows.append(summarize_directional_case("bw_sampling", "BW", cases["bw"])[0])
    rows.append(summarize_tpht_case(cases["tpht_fw"], cases["tpht_bw"], args.tpht_fw_ids)[0])

    write_table_bundle(rows, output_base, DEFAULT_COLUMNS)
    for row in rows:
        print(
            "{smoke}: status={status} selected_unique={latest_selected_unique_pairs} "
            "chain_last={tracking_chain_count_last} tpht_missing={tpht_missing_in_bw} "
            "tpht_extra={tpht_extra_in_bw}".format(**row)
        )

    failed = [row for row in rows if row.get("status") != "pass"]
    if failed:
        print("failing_smokes=" + ",".join(str(row["smoke"]) for row in failed), file=sys.stderr)
        return 1
    print(f"summary_base={output_base}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

