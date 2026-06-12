#!/usr/bin/env python3
"""Parallel lightweight TPHT summary and consistency driver for GMD2026."""

from __future__ import annotations

import os
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_02_tpht import COLUMNS, FEATURE_STATUS_COLUMNS, RANK_COLUMNS, build_feature_status_rows, build_rank_load_rows, build_summary  # noqa: E402
from check_02_tpht_consistency import COLUMNS as CONSISTENCY_COLUMNS  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.table_utils import write_table_bundle  # noqa: E402


def _parallel_options(args) -> dict:
    """Return analysis options with a useful TPHT light-worker default."""
    default_workers = max(1, int(os.environ.get("TPHT_LIGHT_WORKERS", "8")))
    if int(args.workers or 1) <= 1:
        args.workers = default_workers
    return analysis_options(args)


def main() -> None:
    """Run the light TPHT summary and consistency checks with parallel readers."""
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
                "tables/02_tpht_summary.{csv,md,tex,json}",
                "tables/02_tpht_rank_load_balance.{csv,md,tex,json}",
                "tables/02_tpht_analysis_feature_status.{csv,md,tex,json}",
                "tables/02_tpht_consistency.{csv,md,tex,json}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    options = _parallel_options(args)
    summary_cache: dict[str, dict] = {}
    row, warnings = build_summary(root, options, summary_cache)
    rank_rows, rank_warnings = build_rank_load_rows(
        root,
        options,
        summary_cache.get("raw_summary"),
        summary_cache.get("dedup_summary"),
    )
    warnings.extend(rank_warnings)
    consistency_row = {column: row.get(column) for column in CONSISTENCY_COLUMNS}
    if args.strict and (warnings or row.get("consistency_status") == "FAIL"):
        raise RuntimeError("; ".join(warnings) or "TPHT consistency failed")
    write_table_bundle([row], outdir / "tables" / "02_tpht_summary", COLUMNS)
    write_table_bundle(rank_rows, outdir / "tables" / "02_tpht_rank_load_balance", RANK_COLUMNS)
    write_table_bundle(build_feature_status_rows(), outdir / "tables" / "02_tpht_analysis_feature_status", FEATURE_STATUS_COLUMNS)
    write_table_bundle([consistency_row], outdir / "tables" / "02_tpht_consistency", CONSISTENCY_COLUMNS)


if __name__ == "__main__":
    main()
