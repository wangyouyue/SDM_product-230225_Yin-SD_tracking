#!/usr/bin/env python3
"""Classify availability of GMD2026 benchmark diagnostics."""

from __future__ import annotations

import csv
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
GMD_ROOT = SCRIPT_DIR.parent
REPO_ROOT = GMD_ROOT.parents[4]

CLASSIFICATIONS = {
    "wallclock_s": "available_from_job_metrics_json",
    "core_hours": "available_from_job_metrics_json",
    "mpi_ranks": "available_from_job_metrics_json",
    "omp_threads": "available_from_job_metrics_json",
    "node_count": "available_from_job_metrics_json",
    "peak_memory_job_mb": "available_after_model_run",
    "peak_memory_rank_max_mib": "requires_source_hook",
    "peak_memory_rank_sum_mib": "requires_source_hook",
    "tracking_chain_count": "available_after_model_run",
    "tracking_id_memory_bytes": "available_after_model_run",
    "if_coal_memory_bytes": "available_after_model_run",
    "id_assignment_time_s": "requires_source_hook",
    "boundary_tracking_time_s": "requires_source_hook",
    "sd_output_write_time_total_s": "requires_source_hook",
    "sd_output_write_time_mean_s": "requires_source_hook",
    "coalescence_output_write_time_total_s": "requires_source_hook",
    "coalescence_output_write_time_mean_s": "requires_source_hook",
    "tpht_id_write_time_total_s": "requires_source_hook",
    "tpht_id_write_time_mean_s": "requires_source_hook",
    "tpht_id_records_written_total": "requires_source_hook",
    "total_output_bytes": "available_after_model_run",
    "sd_selected_output_bytes": "available_after_model_run",
    "sd_all_output_bytes": "available_after_model_run",
    "coalescence_log_bytes": "available_after_model_run",
    "tpht_id_bytes": "available_after_model_run",
    "number_of_output_files": "available_after_model_run",
    "mean_file_size_bytes": "available_after_model_run",
    "coalescence_event_count": "available_after_model_run",
    "bytes_per_coalescence_event": "available_after_model_run",
    "postprocess_wallclock_s": "available_after_model_run",
}

NOTES = {
    "available_from_job_metrics_json": "Produced by squid_run.sh in job_metrics.json after a model job finishes.",
    "available_after_model_run": "Computed from model output, LOG.pe000000, time_scale_rm_main.log, or NetCDF files after a run.",
    "available_if_GMD_BENCH_DIAG_present": "Available only when grep-friendly GMD_BENCH_DIAG/GMD_IO_DIAG log hooks are present.",
    "available_from_scheduler_only": "Only available from scheduler accounting; do not mix with rank-level memory.",
    "requires_source_hook": "Requires additional Fortran instrumentation; table scripts must report NA until implemented.",
    "not_available": "No current source or post-processing path is known.",
}


def source_contains(pattern: str) -> bool:
    """Return True when a source pattern is present in the current tree."""
    roots = [REPO_ROOT / "contrib" / "SDM", REPO_ROOT / "scale-rm" / "src"]
    for root in roots:
        if not root.exists():
            continue
        for suffix in ("*.f90", "*.F90"):
            for path in root.rglob(suffix):
                try:
                    if pattern in path.read_text(errors="ignore"):
                        return True
                except OSError:
                    continue
    return False


def current_classifications() -> dict[str, str]:
    """Adjust classifications if optional source hooks are detected."""
    classes = dict(CLASSIFICATIONS)
    if source_contains("GMD_BENCH_DIAG"):
        classes["peak_memory_rank_max_mib"] = "available_if_GMD_BENCH_DIAG_present"
        classes["peak_memory_rank_sum_mib"] = "available_if_GMD_BENCH_DIAG_present"
        classes["id_assignment_time_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["boundary_tracking_time_s"] = "available_if_GMD_BENCH_DIAG_present"
    if source_contains("GMD_IO_DIAG"):
        classes["sd_output_write_time_total_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["sd_output_write_time_mean_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["coalescence_output_write_time_total_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["coalescence_output_write_time_mean_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["tpht_id_write_time_total_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["tpht_id_write_time_mean_s"] = "available_if_GMD_BENCH_DIAG_present"
        classes["tpht_id_records_written_total"] = "available_if_GMD_BENCH_DIAG_present"
    return classes


def write_reports(classes: dict[str, str]) -> None:
    """Write CSV and Markdown diagnostics-availability reports."""
    rows = [
        {
            "variable": variable,
            "availability": availability,
            "note": NOTES.get(availability, ""),
        }
        for variable, availability in classes.items()
    ]
    csv_path = GMD_ROOT / "benchmark_diagnostics_availability.csv"
    md_path = GMD_ROOT / "benchmark_diagnostics_availability.md"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variable", "availability", "note"])
        writer.writeheader()
        writer.writerows(rows)
    md_lines = [
        "| variable | availability | note |",
        "| --- | --- | --- |",
    ]
    for row in rows:
        md_lines.append(f"| {row['variable']} | {row['availability']} | {row['note']} |")
    md_path.write_text("\n".join(md_lines) + "\n")


def main() -> int:
    classes = current_classifications()
    write_reports(classes)
    for variable, availability in classes.items():
        print(f"{variable}: {availability}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
