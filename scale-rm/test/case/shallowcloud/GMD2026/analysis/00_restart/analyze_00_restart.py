#!/usr/bin/env python3
"""Analyze the GMD2026 restart sanity check without using it for overhead claims."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import safe_float, write_table_bundle  # noqa: E402


COLUMNS = [
    "case_name",
    "wallclock_s",
    "core_hours",
    "restart_files_found",
    "restart_file_count",
    "restart_total_bytes",
    "first_model_time_s",
    "last_model_time_s",
    "expected_duration_s",
    "exit_status",
    "first_hour_cap_detected",
    "notes",
]


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write the group-00 restart sanity table."""
    config = load_config()
    case = iter_cases(config, "00_base_restart_3d")[0]
    row, _files, records, warnings = collect_case_metrics(root, case, options)
    case_dir = Path(row["case_dir"])

    restart_files = sorted(path for path in case_dir.rglob("*restart*") if path.is_file()) if case_dir.exists() else []
    times = [safe_float(record.get("time_s")) for record in records]
    times = [value for value in times if value is not None]
    cap_fraction = safe_float(row.get("first_hour_cap_active_fraction"))

    output_row = {
        "case_name": row.get("case_name"),
        "wallclock_s": row.get("wallclock_s"),
        "core_hours": row.get("core_hours"),
        "restart_files_found": bool(restart_files) if case_dir.exists() else None,
        "restart_file_count": len(restart_files) if case_dir.exists() else None,
        "restart_total_bytes": sum(path.stat().st_size for path in restart_files) if case_dir.exists() else None,
        "first_model_time_s": min(times) if times else None,
        "last_model_time_s": max(times) if times else None,
        "expected_duration_s": case.get("duration_s"),
        "exit_status": row.get("exit_status"),
        "first_hour_cap_detected": cap_fraction is not None and cap_fraction > 0.0,
        "notes": warning_text(warnings) or "Restart sanity only; excluded from tracking-performance comparisons.",
    }
    if strict and warnings:
        raise RuntimeError(output_row["notes"])
    write_table_bundle([output_row], outdir / "tables" / "00_restart_summary", COLUMNS)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["tables/00_restart_summary.{csv,md,tex,json}"])
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
