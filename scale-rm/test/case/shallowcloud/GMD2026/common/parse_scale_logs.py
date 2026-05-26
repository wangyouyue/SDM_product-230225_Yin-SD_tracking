"""Parse SCALE-SDM benchmark diagnostics from LOG.pe000000."""

from __future__ import annotations

import json
import re
from pathlib import Path


DIAG_MARKERS = ("GMD_BENCH_DIAG", "GMD_IO_DIAG")
LOG_RECORD_STOP_PREFIXES = ("***", "######", "++++++")
KEY_VALUE_RE = re.compile(r"([A-Za-z0-9_]+)\s*=\s*([^\s,]+)")

DEFAULT_LOG_VALUES = {
    "tracking_chain_count_mean": -1,
    "tracking_chain_count_max": -1,
    "chain_count": -1,
    "full_scan_chain_count": -1,
    "tracking_id_memory_bytes": -1,
    "if_coal_memory_bytes": -1,
    "id_assignment_time_s": -1,
    "boundary_tracking_time_s": -1,
    "peak_memory_job_mb": -1,
    "peak_memory_rank_max_mib": -1,
    "peak_memory_rank_sum_mib": -1,
    "sd_output_write_time_last_s": -1,
    "sd_output_write_time_total_s": -1,
    "sd_output_write_count": -1,
    "coalescence_output_write_time_last_s": -1,
    "coalescence_output_write_time_total_s": -1,
    "coalescence_output_write_count": -1,
    "tpht_id_write_time_last_s": -1,
    "tpht_id_write_time_total_s": -1,
    "tpht_id_write_count": -1,
    "tpht_id_records_written_total": -1,
}


def parse_job_metrics(path: Path) -> dict:
    """Read job_metrics.json when present."""
    metric_path = Path(path)
    return json.loads(metric_path.read_text()) if metric_path.exists() else {}


def parse_number(value: str) -> float | None:
    """Parse a Fortran-style numeric value."""
    try:
        return float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def logical_diagnostic_records(text: str) -> list[str]:
    """Collect wrapped Fortran diagnostic writes into logical records."""
    records: list[str] = []
    current: list[str] = []

    def flush() -> None:
        if current:
            records.append(" ".join(part.strip() for part in current))
            current.clear()

    for line in text.splitlines():
        stripped = line.strip()
        if any(marker in line for marker in DIAG_MARKERS):
            flush()
            current.append(line)
            continue
        if not current:
            continue
        if not stripped or stripped.startswith(LOG_RECORD_STOP_PREFIXES):
            flush()
            continue
        current.append(line)
    flush()
    return records


def parse_legacy_tracking_lines(text: str, out: dict) -> None:
    """Parse pre-hook tracking memory and chain-count log lines if present."""
    chains: list[float] = []
    for line in text.splitlines():
        if "tracking_id_memory_estimate_MB" in line:
            match = re.search(r"=\s*([-+0-9.EeDd]+)", line)
            if match:
                value = parse_number(match.group(1))
                if value is not None:
                    out["tracking_id_memory_bytes"] = value * 1024.0 * 1024.0
        if "coal_flag_memory_estimate_MB" in line:
            match = re.search(r"=\s*([-+0-9.EeDd]+)", line)
            if match:
                value = parse_number(match.group(1))
                if value is not None:
                    out["if_coal_memory_bytes"] = value * 1024.0 * 1024.0
        if "tracking_chain_count" in line:
            match = re.search(r"tracking_chain_count[^0-9-]*(\d+)", line)
            if match:
                chains.append(float(match.group(1)))
    if chains and out.get("tracking_chain_count_mean", -1) < 0:
        out["tracking_chain_count_mean"] = sum(chains) / len(chains)
        out["tracking_chain_count_max"] = max(chains)


def parse_log(path: Path) -> dict:
    """Extract GMD_BENCH_DIAG, GMD_IO_DIAG, and legacy tracking diagnostics."""
    out = dict(DEFAULT_LOG_VALUES)
    log_path = Path(path)
    if not log_path.exists():
        return out

    text = log_path.read_text(errors="ignore")
    for record in logical_diagnostic_records(text):
        fields = dict(KEY_VALUE_RE.findall(" ".join(record.split())))
        for key, raw_value in fields.items():
            if key not in out:
                continue
            value = parse_number(raw_value)
            if value is not None:
                out[key] = value

    parse_legacy_tracking_lines(text, out)
    return out
