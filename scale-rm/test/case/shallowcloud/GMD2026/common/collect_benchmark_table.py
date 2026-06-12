"""Collect GMD2026 benchmark diagnostics into readable tables."""

import csv
import json
import re
from pathlib import Path

from audit_benchmark_diagnostics import NOTES, current_classifications
from collect_file_sizes import collect_output_inventory
from parse_scale_logs import parse_job_metrics, parse_log


COLUMNS = [
    "case_name",
    "tracking_mode",
    "tracking_selection_mode",
    "tracking_fraction",
    "domain_km3",
    "resolution_m3",
    "SDNC",
    "start_time_min",
    "end_time_min",
    "output_interval_s",
    "coalescence_output_enable",
    "wallclock_s",
    "core_hours",
    "mpi_ranks",
    "omp_threads",
    "node_count",
    "peak_memory_job_mb",
    "peak_memory_rank_max_mib",
    "peak_memory_rank_sum_mib",
    "tracking_chain_count_mean",
    "tracking_chain_count_max",
    "tracking_id_memory_bytes",
    "if_coal_memory_bytes",
    "id_assignment_time_s",
    "boundary_tracking_time_s",
    "sd_output_write_time_total_s",
    "sd_output_write_time_mean_s",
    "coalescence_output_write_time_total_s",
    "coalescence_output_write_time_mean_s",
    "tpht_id_write_time_total_s",
    "tpht_id_write_time_mean_s",
    "tpht_id_records_written_total",
    "total_output_bytes",
    "sd_selected_output_bytes",
    "sd_all_output_bytes",
    "coalescence_log_bytes",
    "tpht_id_bytes",
    "number_of_output_files",
    "mean_file_size_bytes",
    "coalescence_event_count",
    "bytes_per_coalescence_event",
    "postprocess_wallclock_s",
]

ASSIGNMENT_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s*=\s*([^,!/]+)", re.MULTILINE)
UNAVAILABLE_CLASSES = {"requires_source_hook", "not_available"}


def parse_run_conf(path: Path) -> dict[str, str]:
    """Parse scalar namelist assignments from run.conf."""
    if not path.exists():
        return {}
    return {key: value.strip().strip('"\'') for key, value in ASSIGNMENT_RE.findall(path.read_text(errors="ignore"))}


def parse_float(conf: dict[str, str], key: str, default: float = 0.0) -> float:
    """Parse Fortran-style floating-point values."""
    try:
        return float(conf.get(key, "").replace("D", "E").replace("d", "e"))
    except Exception:
        return default


def count_coalescence_events(case_dir: Path) -> int:
    """Count coalescence records from NetCDF dimensions when netCDF4 is available."""
    try:
        from netCDF4 import Dataset
    except Exception:
        return 0

    total = 0
    for path in case_dir.glob("SD_coal_output_NetCDF_*.nc"):
        with Dataset(path) as nc:
            dimensions = [name for name, dim in nc.dimensions.items() if dim.isunlimited()]
            dimensions += [name for name in nc.dimensions if name.lower() in {"event", "events", "n_event", "coal_event"}]
            if dimensions:
                total += len(nc.dimensions[dimensions[0]])
            elif "time" in nc.variables:
                total += len(nc.variables["time"])
    return total


def parse_time_peak_memory_mb(case_dir: Path) -> float:
    """Read GNU time maximum resident set size when available."""
    path = case_dir / "time_scale_rm_main.log"
    if not path.exists():
        return -1
    for line in path.read_text(errors="ignore").splitlines():
        if "Maximum resident set size" not in line:
            continue
        try:
            return float(line.rsplit(":", 1)[1].strip()) / 1024.0
        except Exception:
            return -1
    return -1


def apply_diagnostic_availability(row: dict, availability: dict[str, str]) -> dict:
    """Replace unsupported benchmark diagnostics with NA."""
    aliases = {
        "tracking_chain_count": ["tracking_chain_count_mean", "tracking_chain_count_max"],
    }
    for metric, status in availability.items():
        for target in aliases.get(metric, [metric]):
            if target in row and status in UNAVAILABLE_CLASSES:
                row[target] = "NA"
    return row


def nonnegative_or_na(value) -> float | str:
    """Return NA for unavailable hook fields instead of reporting sentinel values."""
    try:
        number = float(value)
    except Exception:
        return "NA"
    return number if number >= 0.0 else "NA"


def mean_time(total, count) -> float | str:
    """Compute mean write time from total time and count."""
    total_value = nonnegative_or_na(total)
    count_value = nonnegative_or_na(count)
    if total_value == "NA" or count_value == "NA":
        return "NA"
    if count_value == 0:
        return 0.0
    return total_value / count_value


def case_row(case_dir: Path, start_time_min: float = 60.0) -> dict:
    """Collect one case row."""
    case_dir = Path(case_dir)
    conf = parse_run_conf(case_dir / "run.conf")
    job_path = case_dir / "job_metrics.json"
    log_path = case_dir / "LOG.pe000000"
    job = parse_job_metrics(job_path)
    log = parse_log(log_path)
    files = collect_output_inventory(case_dir)

    job_exists = job_path.exists()
    log_exists = log_path.exists()
    wallclock_s = float(job.get("wallclock_s", 0) or 0) if job_exists else 0.0
    mpi_ranks = int(job.get("mpi_ranks", 0) or 0) if job_exists else 0
    omp_threads = int(job.get("omp_threads", 1) or 1) if job_exists else 1
    coalescence_event_count = count_coalescence_events(case_dir)
    coalescence_log_bytes = files.get("coalescence_log_bytes", 0)

    row = {column: "" for column in COLUMNS}
    row.update(
        {
            "case_name": case_dir.name,
            "tracking_mode": conf.get("tracking_mode", ""),
            "tracking_selection_mode": conf.get("tracking_selection_mode", ""),
            "tracking_fraction": conf.get("tracking_fraction", ""),
            "domain_km3": "2 x 2 x 1.5",
            "resolution_m3": "50 x 50 x 5",
            "SDNC": conf.get("sdm_inisdnc", ""),
            "start_time_min": start_time_min,
            "end_time_min": start_time_min + parse_float(conf, "TIME_DURATION") / 60.0,
            "output_interval_s": conf.get("sdm_dmpitvl", conf.get("HISTORY_DEFAULT_TINTERVAL", "")),
            "coalescence_output_enable": conf.get("coalescence_output_enable", ""),
            "wallclock_s": wallclock_s if job_exists else "NA",
            "core_hours": wallclock_s * mpi_ranks * omp_threads / 3600.0 if job_exists and wallclock_s else "NA",
            "mpi_ranks": mpi_ranks if job_exists else "NA",
            "omp_threads": omp_threads if job_exists else "NA",
            "node_count": job.get("node_count", "NA") if job_exists else "NA",
            "coalescence_event_count": coalescence_event_count,
            "bytes_per_coalescence_event": coalescence_log_bytes / coalescence_event_count if coalescence_event_count else 0.0,
            "postprocess_wallclock_s": "NA",
        }
    )
    row.update(log)
    for key in (
        "peak_memory_rank_max_mib",
        "peak_memory_rank_sum_mib",
        "tracking_chain_count_mean",
        "tracking_chain_count_max",
        "tracking_id_memory_bytes",
        "if_coal_memory_bytes",
        "id_assignment_time_s",
        "boundary_tracking_time_s",
        "sd_output_write_time_total_s",
        "coalescence_output_write_time_total_s",
        "tpht_id_write_time_total_s",
        "tpht_id_records_written_total",
    ):
        row[key] = nonnegative_or_na(row.get(key, "NA"))
    row["sd_output_write_time_mean_s"] = mean_time(
        log.get("sd_output_write_time_total_s", -1),
        log.get("sd_output_write_count", -1),
    )
    row["coalescence_output_write_time_mean_s"] = mean_time(
        log.get("coalescence_output_write_time_total_s", -1),
        log.get("coalescence_output_write_count", -1),
    )
    row["tpht_id_write_time_mean_s"] = mean_time(
        log.get("tpht_id_write_time_total_s", -1),
        log.get("tpht_id_write_count", -1),
    )
    row.update(files)
    if not log_exists:
        for key in (
            "tracking_chain_count_mean",
            "tracking_chain_count_max",
            "tracking_id_memory_bytes",
            "if_coal_memory_bytes",
        ):
            row[key] = "NA"
    if float(row.get("peak_memory_job_mb") or -1) < 0:
        memory = parse_time_peak_memory_mb(case_dir)
        row["peak_memory_job_mb"] = memory if memory >= 0 else "NA"
    return row


def collect_cases(group_dir: Path, case_names: list[str], start_time_min: float = 60.0) -> list[dict]:
    """Collect benchmark rows for case_names."""
    availability = current_classifications()
    return [apply_diagnostic_availability(case_row(Path(group_dir) / name, start_time_min), availability) for name in case_names]


def availability_footnote(availability: dict[str, str]) -> str:
    """Build a compact source note for table fields."""
    grouped: dict[str, list[str]] = {}
    for metric, status in availability.items():
        grouped.setdefault(status, []).append(metric)
    parts = []
    for status, metrics in sorted(grouped.items()):
        parts.append(f"{status}: {', '.join(metrics)}. {NOTES.get(status, '')}")
    return " ".join(parts)


def write_outputs(rows: list[dict], output_base: Path) -> None:
    """Write CSV, Markdown, LaTeX, and JSON tables."""
    base = Path(output_base)
    base.parent.mkdir(parents=True, exist_ok=True)

    with base.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    base.with_suffix(".json").write_text(json.dumps(rows, indent=2))
    availability = current_classifications()
    md_table = (
        "| " + " | ".join(COLUMNS) + " |\n"
        + "| " + " | ".join(["---"] * len(COLUMNS)) + " |\n"
        + "".join("| " + " | ".join(str(row.get(column, "")) for column in COLUMNS) + " |\n" for row in rows)
    )
    md_table += "\nMetric availability note: " + availability_footnote(availability) + "\n"
    md_table += (
        "\nScheduler/job peak memory (`peak_memory_job_mb`) is not rank-summed MPI memory. "
        "`peak_memory_rank_max_mib` and `peak_memory_rank_sum_mib` are rank-level MiB values from "
        "`/proc/self/status` VmHWM via source-level hooks.\n"
    )
    base.with_suffix(".md").write_text(md_table)

    compact = [
        "case_name",
        "tracking_mode",
        "tracking_fraction",
        "coalescence_output_enable",
        "wallclock_s",
        "core_hours",
        "node_count",
        "total_output_bytes",
        "coalescence_event_count",
    ]
    lines = ["\\begin{tabular}{lllllllll}", "\\hline", " & ".join(compact).replace("_", "\\_") + " \\\\", "\\hline"]
    for row in rows:
        lines.append(" & ".join(str(row.get(column, "")) for column in compact).replace("_", "\\_") + " \\\\")
    lines.extend(
        [
            "\\hline",
            "\\end{tabular}",
            "",
            "% Metric availability note: " + availability_footnote(availability).replace("_", "\\_"),
            "% peak\\_memory\\_job\\_mb is scheduler/job-level memory, not rank-summed MPI memory.",
        ]
    )
    base.with_suffix(".tex").write_text("\n".join(lines) + "\n")
