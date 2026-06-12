"""Path, configuration, and command-line helpers for GMD2026 analysis."""

from __future__ import annotations

import argparse
import json
import warnings as warning_module
from pathlib import Path
from typing import Any

from .table_utils import append_warning


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ANALYSIS_DIR / "analysis_config.json"


def build_arg_parser(description: str) -> argparse.ArgumentParser:
    """Create the common parser required by every analysis and plotting script."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--root", type=Path, default=ANALYSIS_DIR.parent, help="Path to the GMD2026 suite root.")
    parser.add_argument("--outdir", type=Path, default=None, help="Directory for generated tables and figures.")
    parser.add_argument("--strict", action="store_true", help="Treat missing diagnostics and failed QC checks as errors.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without writing outputs.")
    parser.add_argument("--quick", action="store_true", help="Run a minimal foreground-safe check and small summary.")
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Parse logs, job metrics, file sizes, and NetCDF headers without reading large variables.",
    )
    parser.add_argument(
        "--skip-heavy-netcdf",
        action="store_true",
        help="Skip DSD, z-r histograms, target histories, and chain reconstruction.",
    )
    parser.add_argument("--max-files", type=int, default=None, help="Process at most N NetCDF files for debugging.")
    parser.add_argument("--max-records", type=int, default=None, help="Process at most N SD records per NetCDF file when chunked reading is available.")
    parser.add_argument("--workers", type=int, default=1, help="Requested worker count for scripts that implement parallel processing.")
    parser.add_argument("--chunk-size", type=int, default=100000, help="Records per NetCDF variable chunk for memory-safe reading.")
    return parser


def analysis_options(args: argparse.Namespace) -> dict[str, Any]:
    """Return normalized execution-mode options from parsed arguments."""
    metadata_only = bool(args.metadata_only or args.quick)
    skip_heavy = bool(args.skip_heavy_netcdf or args.quick or metadata_only)
    max_files = args.max_files
    max_records = args.max_records
    if args.quick:
        max_files = 1 if max_files is None else min(max_files, 1)
        max_records = 1000 if max_records is None else min(max_records, 1000)
    return {
        "quick": bool(args.quick),
        "metadata_only": metadata_only,
        "skip_heavy_netcdf": skip_heavy,
        "max_files": max_files,
        "max_records": max_records,
        "workers": max(1, int(args.workers or 1)),
        "chunk_size": max(1, int(args.chunk_size or 100000)),
    }


def resolve_outdir(root: Path, outdir: Path | None) -> Path:
    """Return the output directory, defaulting below the GMD2026 root."""
    return (outdir or root / "analysis_outputs").resolve()


def load_config(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load the GMD2026 analysis configuration."""
    return json.loads(config_path.read_text())


def iter_cases(config: dict[str, Any], group_name: str | None = None) -> list[dict[str, Any]]:
    """Return flattened case metadata from the analysis configuration."""
    cases: list[dict[str, Any]] = []
    for current_group, group in config.get("groups", {}).items():
        if group_name is not None and current_group != group_name:
            continue
        for case in group.get("cases", []):
            row = dict(case)
            row.setdefault("group", current_group)
            cases.append(row)
    return cases


def get_group(config: dict[str, Any], group_name: str) -> dict[str, Any]:
    """Return one group configuration."""
    return config.get("groups", {}).get(group_name, {})


def case_path(root: Path, case: dict[str, Any]) -> Path:
    """Resolve a case path, using optional aliases when the canonical name was renamed."""
    root = root.resolve()
    group = case["group"]
    candidate_names = [case["case_name"]]
    candidate_names.extend(case.get("aliases", []))
    for name in candidate_names:
        candidate = root / group / name
        if candidate.exists():
            return candidate
    return root / group / case["case_name"]


def output_dirs(outdir: Path) -> tuple[Path, Path]:
    """Return table and figure directories."""
    return outdir / "tables", outdir / "figures"


def ensure_output_dirs(outdir: Path) -> None:
    """Create output subdirectories."""
    tables, figures = output_dirs(outdir)
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)


def warn_or_raise(warnings: list[str], message: str, strict: bool = False) -> None:
    """Record a warning and optionally raise in strict mode."""
    append_warning(warnings, message)
    warning_module.warn(message, RuntimeWarning, stacklevel=2)
    if strict:
        raise RuntimeError(message)


def warning_text(warnings: list[str]) -> str | None:
    """Serialize warnings for table output."""
    return "; ".join(warnings) if warnings else None


def dry_run_message(script_name: str, root: Path, outdir: Path, outputs: list[str]) -> None:
    """Print the actions a script would take during dry-run mode."""
    print(f"[dry-run] {script_name}")
    print(f"[dry-run] root={root}")
    print(f"[dry-run] outdir={outdir}")
    for output in outputs:
        print(f"[dry-run] would write {output}")
