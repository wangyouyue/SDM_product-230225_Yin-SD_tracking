"""High-level GMD2026 case collectors built from common parsers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .file_stats import collect_file_stats
from .parse_gmd_diag import parse_case_logs
from .parse_job_metrics import parse_job_metrics
from .parse_netcdf import inspect_netcdf_files
from .parse_time_log import parse_time_log
from .paths import case_path
from .table_utils import safe_float


def collect_case_metrics(
    root: Path,
    case: dict[str, Any],
    options: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, list[Path]], list[dict[str, Any]], list[str]]:
    """Collect job, LOG, file-size, and NetCDF metrics for one case."""
    options = options or {}
    warnings: list[str] = []
    row = dict(case)
    current_case_path = case_path(root, case)
    row["case_dir"] = str(current_case_path)

    if not current_case_path.exists():
        warnings.append(f"case directory missing: {current_case_path}")
        return row, {}, [], warnings

    job = parse_job_metrics(current_case_path / "job_metrics.json")
    if job.get("wallclock_s") is None:
        time_metrics = parse_time_log(current_case_path / "time_scale_rm_main.log")
        if time_metrics.get("wallclock_s") is not None:
            job["wallclock_s"] = time_metrics.get("wallclock_s")
            job["exit_status"] = time_metrics.get("exit_status")
            job["max_resident_set_mib"] = time_metrics.get("max_resident_set_mib")
        else:
            warnings.append(f"job_metrics.json and usable time_scale_rm_main.log missing in {current_case_path}")
    row.update(job)

    diag, records, log_warnings = parse_case_logs(current_case_path)
    for key, value in diag.items():
        if key == "tracking_mode" and row.get("tracking_mode") is not None:
            row["diagnostic_tracking_mode"] = value
            continue
        if value is None and key in row:
            continue
        row[key] = value
    warnings.extend(log_warnings)

    stats, files_by_metric, stat_warnings = collect_file_stats(current_case_path)
    row.update(stats)
    warnings.extend(stat_warnings)

    coal_files = files_by_metric.get("coalescence_log_bytes", [])
    if coal_files:
        netcdf_stats, netcdf_warnings = inspect_netcdf_files(
            coal_files,
            max_files=options.get("max_files"),
            metadata_only=options.get("metadata_only", False) or options.get("quick", False),
        )
        row.update({key: value for key, value in netcdf_stats.items() if key == "coalescence_event_count"})
        warnings.extend(netcdf_warnings)
    else:
        row["coalescence_event_count"] = None

    chain_count = safe_float(row.get("tracking_chain_count"))
    selected_bytes = safe_float(row.get("sd_selected_output_bytes"))
    coal_bytes = safe_float(row.get("coalescence_log_bytes"))
    total_bytes = safe_float(row.get("total_output_bytes"))
    scientific_bytes = safe_float(row.get("scientific_output_bytes"))
    sd_particle_bytes = safe_float(row.get("sd_particle_output_bytes"))
    event_count = safe_float(row.get("coalescence_event_count"))
    row["output_bytes_gib"] = total_bytes / 1024.0**3 if total_bytes is not None else None
    row["scientific_output_gib"] = scientific_bytes / 1024.0**3 if scientific_bytes is not None else None
    row["sd_particle_output_gib"] = sd_particle_bytes / 1024.0**3 if sd_particle_bytes is not None else None
    row["selected_output_gib"] = selected_bytes / 1024.0**3 if selected_bytes is not None else None
    row["coalescence_log_gib"] = coal_bytes / 1024.0**3 if coal_bytes is not None else None
    row["bytes_per_coalescence_event"] = coal_bytes / event_count if coal_bytes is not None and event_count else None
    row["output_bytes_per_tracked_chain"] = (
        total_bytes / chain_count if total_bytes is not None and chain_count else None
    )
    return row, files_by_metric, records, warnings
