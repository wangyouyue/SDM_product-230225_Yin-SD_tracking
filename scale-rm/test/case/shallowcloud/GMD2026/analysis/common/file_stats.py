"""Collect GMD2026 output file-size metrics."""

from __future__ import annotations

from pathlib import Path
from typing import Any


PATTERNS = {
    "sd_selected_output_bytes": ["SD_selected_NetCDF_*"],
    "sd_all_output_bytes": ["SD_all_NetCDF_*"],
    "coalescence_log_bytes": ["SD_coal_output_NetCDF_*"],
    "tpht_id_bytes": ["*.ids"],
    "history_output_bytes": ["history.pe*"],
    "log_bytes": ["LOG.pe*"],
    "job_metric_bytes": ["job_metrics.json"],
    "time_log_bytes": ["time_scale_rm_main.log"],
}

SCIENTIFIC_OUTPUT_METRICS = {
    "sd_selected_output_bytes",
    "sd_all_output_bytes",
    "coalescence_log_bytes",
    "tpht_id_bytes",
    "history_output_bytes",
}
AUXILIARY_OUTPUT_METRICS = {"log_bytes", "job_metric_bytes", "time_log_bytes"}


def _matching_files(case_dir: Path, pattern: str) -> list[Path]:
    """Find files matching a pattern below a case directory."""
    paths = sorted(path for path in case_dir.rglob(pattern) if path.is_file())
    if pattern.startswith(("SD_selected_NetCDF_", "SD_all_NetCDF_", "SD_coal_output_NetCDF_", "history")):
        # SCALE-SDM sidecar ID files can share the NetCDF stem; keep them in
        # tpht_id_bytes rather than mixing them into NetCDF byte counts.
        return [path for path in paths if path.suffix != ".ids"]
    return paths


def collect_file_stats(case_dir: Path) -> tuple[dict[str, Any], dict[str, list[Path]], list[str]]:
    """Collect output sizes for a model case.

    Counts and byte sizes are true zeroes when a case directory exists and a
    pattern is intentionally absent. A missing case directory is reported with
    NA-like ``None`` values by callers.
    """
    warnings: list[str] = []
    files_by_metric: dict[str, list[Path]] = {}
    if not case_dir.exists():
        return {}, files_by_metric, [f"case directory missing: {case_dir}"]

    scientific_files: set[Path] = set()
    auxiliary_files: set[Path] = set()
    stats: dict[str, Any] = {}
    for metric, patterns in PATTERNS.items():
        paths: list[Path] = []
        for pattern in patterns:
            paths.extend(_matching_files(case_dir, pattern))
        unique_paths = sorted(set(paths))
        files_by_metric[metric] = unique_paths
        stats[metric] = sum(path.stat().st_size for path in unique_paths)
        if metric in SCIENTIFIC_OUTPUT_METRICS:
            scientific_files.update(unique_paths)
        elif metric in AUXILIARY_OUTPUT_METRICS:
            auxiliary_files.update(unique_paths)

    sd_particle_files = set(files_by_metric.get("sd_selected_output_bytes", [])) | set(files_by_metric.get("sd_all_output_bytes", []))
    stats["sd_selected_file_count"] = len(files_by_metric.get("sd_selected_output_bytes", []))
    stats["sd_all_file_count"] = len(files_by_metric.get("sd_all_output_bytes", []))
    stats["sd_particle_output_bytes"] = stats["sd_selected_output_bytes"] + stats["sd_all_output_bytes"]
    stats["sd_particle_file_count"] = len(sd_particle_files)
    stats["history_output_file_count"] = len(files_by_metric.get("history_output_bytes", []))
    stats["scientific_output_bytes"] = sum(path.stat().st_size for path in scientific_files)
    stats["scientific_output_file_count"] = len(scientific_files)
    stats["auxiliary_file_bytes"] = sum(path.stat().st_size for path in auxiliary_files)
    stats["auxiliary_file_count"] = len(auxiliary_files)
    # Backward-compatible names now use only scientific model outputs, not LOG.pe*
    # or job/time metadata. This keeps storage plots focused on data products.
    stats["total_output_bytes"] = stats["scientific_output_bytes"]
    stats["number_of_output_files"] = stats["scientific_output_file_count"]
    stats["mean_file_size_bytes"] = (
        stats["scientific_output_bytes"] / stats["scientific_output_file_count"]
        if stats["scientific_output_file_count"]
        else 0
    )
    stats["coalescence_log_file_count"] = len(files_by_metric.get("coalescence_log_bytes", []))
    stats["tpht_id_file_count"] = len(files_by_metric.get("tpht_id_bytes", []))
    return stats, files_by_metric, warnings
