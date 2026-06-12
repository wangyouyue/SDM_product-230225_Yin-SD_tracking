#!/usr/bin/env python3
"""Estimate foreground-vs-qsub execution mode for GMD2026 analysis scripts."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, get_group, iter_cases, load_config, resolve_outdir  # noqa: E402
from common.table_utils import clean_value  # noqa: E402


GIB = 1024**3
GROUPS = [
    "00_base_restart_3d",
    "01_bench_3d_samp_30min",
    "02_tpht_3d_interest_70min",
    "03_fw_rep_2d_600s",
    "04_sdnc_scaling_lite_30min",
    "05_outint_io_lite_10min",
    "run_all_analysis.sh",
    "run_all_plots.sh",
]


def _files(root: Path, group_name: str, pattern: str) -> list[Path]:
    """Return matching files for one group."""
    group_dir = root / group_name
    if not group_dir.exists():
        return []
    paths = sorted(path for path in group_dir.rglob(pattern) if path.is_file())
    if pattern.startswith(("SD_selected_NetCDF_", "SD_all_NetCDF_", "SD_coal_output_NetCDF_", "history.pe")):
        return [path for path in paths if path.suffix != ".ids"]
    return paths


def _sum(paths: list[Path]) -> int:
    """Return total size for existing paths."""
    return sum(path.stat().st_size for path in paths if path.exists())


def _classify(group_name: str, total_bytes: int, netcdf_count: int, sdnc80_bytes: int, bw_selected_bytes: int) -> tuple[str, str, str, str]:
    """Return runtime class, memory class, execution mode, and reason."""
    total_gib = total_bytes / GIB
    many_netcdf = netcdf_count >= 200
    runtime_class = "<5min"
    memory_class = "low"
    mode = "foreground_ok"
    reason = "small log/file-size scan"

    if total_gib > 10.0 or many_netcdf:
        runtime_class = ">30min"
        memory_class = "high"
        mode = "qsub_required"
        reason = "total data exceeds 10 GiB or many NetCDF files must be scanned"
    elif total_gib >= 1.0 or netcdf_count > 0:
        runtime_class = "5-30min"
        memory_class = "medium"
        mode = "qsub_recommended"
        reason = "NetCDF output exists or total data is 1-10 GiB"

    if group_name == "00_base_restart_3d" and total_gib < 1.0:
        return "<5min", "low", "foreground_ok", "restart sanity is log/file-size dominated"
    if group_name == "01_bench_3d_samp_30min" and netcdf_count > 0 and mode == "foreground_ok":
        return "5-30min", "medium", "qsub_recommended", "benchmark NetCDF output exists"
    if group_name == "02_tpht_3d_interest_70min":
        if bw_selected_bytes > 10 * GIB or bw_selected_bytes > 0 and netcdf_count >= 64:
            return ">30min", "high", "qsub_required", "TPHT heavy target histories and chain validity read BW selected output"
        return "5-30min", "medium", "qsub_recommended", "TPHT ID summary is moderate, heavy NetCDF analysis may be large"
    if group_name == "03_fw_rep_2d_600s":
        return ">30min", "high", "qsub_required", "many seeds and NetCDF distribution calculations"
    if group_name == "04_sdnc_scaling_lite_30min":
        if sdnc80_bytes > 5 * GIB or total_gib > 10.0:
            return ">30min", "high", "qsub_required", "SDNC80 output is large"
        return "5-30min", "medium", "qsub_recommended", "SDNC scaling scans multiple NetCDF-producing cases"
    if group_name == "05_outint_io_lite_10min":
        return "5-30min", "medium", "qsub_recommended", "output-interval I/O scan compares multiple output frequencies"
    if group_name == "run_all_analysis.sh":
        return ">30min", "high", "qsub_required", "full suite includes qsub-required sampling and TPHT heavy analysis"
    if group_name == "run_all_plots.sh":
        return "5-30min", "medium", "qsub_recommended", "all figures are generated; foreground is acceptable only when plotting from small tables is allowed"
    return runtime_class, memory_class, mode, reason


def estimate(root: Path) -> list[dict[str, Any]]:
    """Estimate execution cost for each analysis group and suite runner."""
    config = load_config()
    rows: list[dict[str, Any]] = []
    for group_name in GROUPS:
        if group_name.endswith(".sh"):
            case_count = None
            log_files: list[Path] = []
            netcdf_files: list[Path] = []
            ids_files: list[Path] = []
            selected_files: list[Path] = []
            coal_files: list[Path] = []
            sdnc80_bytes = 0
            bw_selected_bytes = 0
        else:
            group_dir = root / group_name
            group_config = get_group(config, group_name)
            case_count = len(group_config.get("cases", []))
            log_files = _files(root, group_name, "LOG.pe*")
            history_files = _files(root, group_name, "history.pe*")
            selected_files = _files(root, group_name, "SD_selected_NetCDF_*")
            all_sd_files = _files(root, group_name, "SD_all_NetCDF_*")
            coal_files = _files(root, group_name, "SD_coal_output_NetCDF_*")
            netcdf_files = history_files + selected_files + all_sd_files + coal_files
            netcdf_files = sorted(set(netcdf_files))
            ids_files = _files(root, group_name, "*.ids")
            sdnc80_bytes = _sum([path for path in selected_files + coal_files if "sdnc80" in str(path)])
            bw_selected_bytes = _sum([path for path in selected_files if "bw_reconstruction" in str(path)])
            if not group_dir.exists():
                log_files = []
                netcdf_files = []
                ids_files = []
                history_files = []
                selected_files = []
                all_sd_files = []
                coal_files = []
        all_read_files = sorted(set(log_files + netcdf_files + ids_files))
        total_bytes = _sum(all_read_files)
        runtime_class, memory_class, mode, reason = _classify(
            group_name,
            total_bytes,
            len(netcdf_files),
            sdnc80_bytes,
            bw_selected_bytes,
        )
        rows.append(
            {
                "group_name": group_name,
                "case_count": case_count,
                "log_file_count": len(log_files),
                "netcdf_file_count": len(netcdf_files),
                "ids_file_count": len(ids_files),
                "total_bytes": total_bytes,
                "history_output_bytes": _sum(history_files) if not group_name.endswith(".sh") else 0,
                "sd_selected_bytes": _sum(selected_files),
                "sd_all_bytes": _sum(all_sd_files) if not group_name.endswith(".sh") else 0,
                "coalescence_log_bytes": _sum(coal_files),
                "tpht_id_bytes": _sum(ids_files),
                "estimated_runtime_class": runtime_class,
                "estimated_memory_class": memory_class,
                "recommended_execution_mode": mode,
                "reason": reason,
            }
        )
    return rows


def write_outputs(rows: list[dict[str, Any]], outdir: Path) -> None:
    """Write CSV, Markdown, and JSON cost-estimate tables."""
    tables = outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0]) if rows else []
    cleaned = [{column: clean_value(row.get(column)) for column in columns} for row in rows]
    with (tables / "analysis_cost_estimate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(cleaned)
    (tables / "analysis_cost_estimate.json").write_text(json.dumps(cleaned, indent=2) + "\n")
    md_lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    md_lines.extend("| " + " | ".join(str(row.get(column, "NA")) for column in columns) + " |" for row in cleaned)
    (tables / "analysis_cost_estimate.md").write_text("\n".join(md_lines) + "\n")


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["tables/analysis_cost_estimate.{csv,md,json}"])
        return
    ensure_output_dirs(outdir)
    write_outputs(estimate(root), outdir)


if __name__ == "__main__":
    main()
