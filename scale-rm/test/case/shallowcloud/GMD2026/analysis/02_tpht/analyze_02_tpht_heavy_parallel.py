#!/usr/bin/env python3
"""Parallel heavy TPHT target-history and chain-diagnostic driver for GMD2026."""

from __future__ import annotations

import os
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_02_tpht_chains import analyze as analyze_chains  # noqa: E402
from analyze_02_tpht_science import analyze as analyze_science  # noqa: E402
from analyze_02_tpht_targets import analyze as analyze_targets  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402


def _parallel_options(args) -> dict:
    """Return analysis options with a useful TPHT heavy-worker default."""
    default_workers = max(1, int(os.environ.get("TPHT_HEAVY_WORKERS", "8")))
    if int(args.workers or 1) <= 1:
        args.workers = default_workers
    return analysis_options(args)


def main() -> None:
    """Run heavy TPHT science, target, and chain diagnostics with parallel readers."""
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
                "tables/02_tpht_science_target_summary.{csv,md,tex,json}",
                "tables/02_tpht_target_trajectory_records.{csv,md,tex,json}",
                "tables/02_tpht_target_interval_diagnostics.{csv,md,tex,json}",
                "tables/02_tpht_target_event_links.{csv,md,tex,json}",
                "tables/02_tpht_target_categories.{csv,md,tex,json}",
                "tables/02_tpht_chain_validity.{csv,md,tex,json}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    options = _parallel_options(args)
    analyze_science(root, outdir, args.strict, options)
    analyze_targets(root, outdir, args.strict, options)
    analyze_chains(root, outdir, args.strict, options)


if __name__ == "__main__":
    main()
