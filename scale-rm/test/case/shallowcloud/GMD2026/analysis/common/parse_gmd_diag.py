"""Parse GMD2026 SCALE-SDM diagnostics from LOG.pe* files."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

DIAG_MARKERS = ("GMD_BENCH_DIAG", "GMD_IO_DIAG", "GMD_CAP_DIAG", "TPHT_ID_DIAG")
STOP_PREFIXES = ("***", "######", "++++++")
KEY_VALUE_RE = re.compile(r"([A-Za-z0-9_]+)\s*=\s*(\"[^\"]*\"|'[^']*'|[^\s,]+)")
RANK_RE = re.compile(r"pe(\d{6})")

BENCH_FIELDS = [
    "time_s",
    "tracking_mode",
    "tracking_chain_count",
    "tracking_id_memory_bytes",
    "if_coal_memory_bytes",
    "peak_memory_rank_max_mib",
    "peak_memory_rank_sum_mib",
    "id_assignment_time_s",
    "boundary_tracking_time_s",
]
IO_FIELDS = [
    "sd_output_write_time_last_s",
    "sd_output_write_time_total_s",
    "sd_output_write_count",
    "sd_output_write_time_mean_s",
    "coalescence_output_write_time_last_s",
    "coalescence_output_write_time_total_s",
    "coalescence_output_write_count",
    "coalescence_output_write_time_mean_s",
    "tpht_id_write_time_last_s",
    "tpht_id_write_time_total_s",
    "tpht_id_write_count",
    "tpht_id_write_time_mean_s",
    "tpht_id_records_written_total",
]
TPHT_FIELDS = [
    "mode",
    "id_epoch_sec",
    "id_input_basename",
    "id_output_basename",
    "ids_loaded",
    "fallback_to_sampling",
]
OLD_IO_MAP = {
    "sd_output_write_time_s": "sd_output_write_time_last_s",
    "coalescence_output_write_time_s": "coalescence_output_write_time_last_s",
    "tpht_id_write_time_s": "tpht_id_write_time_last_s",
}


def parse_value(raw: str) -> Any:
    """Parse Fortran-style numbers, logicals, and quoted strings."""
    text = raw.strip().strip(",")
    if len(text) >= 2 and text[0] in {'"', "'"} and text[-1] == text[0]:
        return text[1:-1]
    lower = text.lower().strip(".")
    if lower in {"true", "t"}:
        return True
    if lower in {"false", "f"}:
        return False
    try:
        value = float(text.replace("D", "E").replace("d", "e"))
    except ValueError:
        return text
    return int(value) if value.is_integer() else value


def _logical_records(text: str) -> list[str]:
    """Collect wrapped diagnostic writes into logical one-line records."""
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
        if not stripped or stripped.startswith(STOP_PREFIXES):
            flush()
            continue
        current.append(line)
    flush()
    return records


def parse_log_file(path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    """Parse diagnostic records from one LOG.pe* file."""
    warnings: list[str] = []
    if not path.exists():
        return [], [f"LOG file missing: {path}"]

    rank_match = RANK_RE.search(path.name)
    rank = int(rank_match.group(1)) if rank_match else None
    records: list[dict[str, Any]] = []
    text = path.read_text(errors="ignore")
    for logical in _logical_records(text):
        marker = next((item for item in DIAG_MARKERS if item in logical), None)
        if marker is None:
            continue
        fields = {key: parse_value(value) for key, value in KEY_VALUE_RE.findall(logical)}
        fields["record_type"] = marker
        fields["source_log"] = path.name
        fields["rank"] = rank
        for old_key, new_key in OLD_IO_MAP.items():
            if old_key in fields and new_key not in fields:
                fields[new_key] = fields[old_key]
                warnings.append(
                    f"{path.name}: old {old_key} is treated as last-write timing only, not total write time"
                )
        records.append(fields)
    return records, warnings


def _numeric(value: Any) -> float | None:
    """Return a finite float or None."""
    if isinstance(value, bool):
        return None
    try:
        output = float(value)
    except (TypeError, ValueError):
        return None
    return output if math.isfinite(output) else None


def _latest_by_time(records: list[dict[str, Any]], key: str) -> Any:
    """Return the value from the latest record that contains key."""
    candidates = [record for record in records if key in record]
    if not candidates:
        return None
    candidates.sort(key=lambda record: (_numeric(record.get("time_s")) is not None, _numeric(record.get("time_s")) or -1.0))
    return candidates[-1].get(key)


def _max_numeric(records: list[dict[str, Any]], key: str) -> float | int | None:
    """Return the maximum finite numeric value for key."""
    values = [_numeric(record.get(key)) for record in records]
    finite = [value for value in values if value is not None]
    if not finite:
        return None
    output = max(finite)
    return int(output) if output.is_integer() else output


def summarize_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize LOG.pe* diagnostics into one table row."""
    summary: dict[str, Any] = {}
    for key in BENCH_FIELDS + IO_FIELDS + TPHT_FIELDS:
        summary[key] = None

    for key in ("peak_memory_rank_max_mib", "peak_memory_rank_sum_mib"):
        summary[key] = _max_numeric(records, key)
    for key in (
        "sd_output_write_time_total_s",
        "sd_output_write_count",
        "coalescence_output_write_time_total_s",
        "coalescence_output_write_count",
        "tpht_id_write_time_total_s",
        "tpht_id_write_count",
        "tpht_id_records_written_total",
        "tracking_chain_count",
    ):
        summary[key] = _max_numeric(records, key)
    for key in (
        "tracking_mode",
        "tracking_id_memory_bytes",
        "if_coal_memory_bytes",
        "id_assignment_time_s",
        "boundary_tracking_time_s",
        "sd_output_write_time_last_s",
        "sd_output_write_time_mean_s",
        "coalescence_output_write_time_last_s",
        "coalescence_output_write_time_mean_s",
        "tpht_id_write_time_last_s",
        "tpht_id_write_time_mean_s",
        "mode",
        "id_epoch_sec",
        "id_input_basename",
        "id_output_basename",
        "ids_loaded",
        "fallback_to_sampling",
    ):
        summary[key] = _latest_by_time(records, key)

    if summary["tracking_chain_count"] is None:
        for alias in ("chain_count", "tracking_chain_count_mean", "full_scan_chain_count"):
            value = _max_numeric(records, alias)
            if value is not None:
                summary["tracking_chain_count"] = value
                break

    cap_records = [record for record in records if record.get("record_type") == "GMD_CAP_DIAG"]
    active_values = [bool(record.get("supersaturation_cap_active")) for record in cap_records if "supersaturation_cap_active" in record]
    summary["supersaturation_cap_active_fraction"] = (
        sum(1.0 for active in active_values if active) / len(active_values) if active_values else None
    )
    first_hour_values = [
        bool(record.get("supersaturation_cap_active"))
        for record in cap_records
        if _numeric(record.get("time_s")) is not None
        and (_numeric(record.get("time_s")) or 0.0) <= 3600.0
        and "supersaturation_cap_active" in record
    ]
    summary["first_hour_cap_active_fraction"] = (
        sum(1.0 for active in first_hour_values if active) / len(first_hour_values) if first_hour_values else None
    )
    return summary


def parse_case_logs(case_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[str]]:
    """Parse all LOG.pe* files in a case directory."""
    warnings: list[str] = []
    log_files = sorted(case_dir.glob("LOG.pe*"))
    if not log_files:
        return summarize_records([]), [], [f"LOG.pe* missing in {case_dir}"]
    records: list[dict[str, Any]] = []
    for path in log_files:
        current, current_warnings = parse_log_file(path)
        records.extend(current)
        warnings.extend(current_warnings)
    return summarize_records(records), records, warnings

