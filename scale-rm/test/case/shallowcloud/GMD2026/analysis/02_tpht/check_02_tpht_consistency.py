#!/usr/bin/env python3
"""Check TPHT FW target-set handoff against BW reconstructed selected output."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_02_tpht import build_summary  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.table_utils import write_table_bundle  # noqa: E402


COLUMNS = [
    "missing_in_bw",
    "extra_in_bw",
    "fallback_to_sampling",
    "id_epoch_sec_fw",
    "id_epoch_sec_bw",
    "mpi_decomposition_match",
    "consistency_status",
    "warnings",
]


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write a compact consistency table."""
    row, warnings = build_summary(root, options)
    out = {column: row.get(column) for column in COLUMNS}
    if strict and (warnings or row.get("consistency_status") == "FAIL"):
        raise RuntimeError(out.get("warnings") or "TPHT consistency failed")
    write_table_bundle([out], outdir / "tables" / "02_tpht_consistency", COLUMNS)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["tables/02_tpht_consistency.{csv,md,tex,json}"])
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
