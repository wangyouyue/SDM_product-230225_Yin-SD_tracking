"""Parse /usr/bin/time output from GMD2026 SQUID jobs."""

from __future__ import annotations

import re
from pathlib import Path

from .units import kb_to_mib


TIME_KEYS = {
    "User time (seconds)": "user_time_s",
    "System time (seconds)": "system_time_s",
    "Percent of CPU this job got": "cpu_percent",
    "Maximum resident set size (kbytes)": "max_resident_set_kb",
    "Major (requiring I/O) page faults": "major_page_faults",
    "Minor (reclaiming a frame) page faults": "minor_page_faults",
    "File system inputs": "file_system_inputs",
    "File system outputs": "file_system_outputs",
    "Exit status": "exit_status",
}


def parse_wallclock(value: str) -> float | None:
    """Parse wallclock strings in m:ss, h:mm:ss, or shell ``1m23s`` format."""
    text = value.strip()
    shell_match = re.match(r"^(?:(?P<minutes>\d+(?:\.\d+)?)m)?(?P<seconds>\d+(?:\.\d+)?)s$", text)
    if shell_match:
        minutes = float(shell_match.group("minutes") or 0.0)
        return minutes * 60.0 + float(shell_match.group("seconds"))

    parts = text.split(":")
    try:
        if len(parts) == 2:
            minutes, seconds = parts
            return float(minutes) * 60.0 + float(seconds)
        if len(parts) == 3:
            hours, minutes, seconds = parts
            return float(hours) * 3600.0 + float(minutes) * 60.0 + float(seconds)
    except ValueError:
        return None
    return None


def _parse_number(value: str) -> float | int | None:
    """Parse a numeric field from /usr/bin/time."""
    stripped = value.strip().rstrip("%")
    try:
        number = float(stripped)
    except ValueError:
        return None
    return int(number) if number.is_integer() else number


def parse_time_text(text: str) -> dict[str, float | int | None]:
    """Parse the text content of a time log."""
    out: dict[str, float | int | None] = {
        "wallclock_s": None,
        "user_time_s": None,
        "system_time_s": None,
        "cpu_percent": None,
        "max_resident_set_kb": None,
        "max_resident_set_mib": None,
        "major_page_faults": None,
        "minor_page_faults": None,
        "file_system_inputs": None,
        "file_system_outputs": None,
        "exit_status": None,
    }

    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line and ":" in line:
            match = re.search(r":\s*([0-9:.]+)\s*$", line)
            out["wallclock_s"] = parse_wallclock(match.group(1)) if match else None
            continue
        match = re.match(r"^\s*(real|user|sys)\s+(.+)$", line)
        if match:
            key = {"real": "wallclock_s", "user": "user_time_s", "sys": "system_time_s"}[match.group(1)]
            out[key] = parse_wallclock(match.group(2))
            continue
        if ":" not in line:
            continue
        key_text, raw_value = line.split(":", 1)
        key = TIME_KEYS.get(key_text.strip())
        if key:
            out[key] = _parse_number(raw_value)

    if out["max_resident_set_kb"] is not None:
        out["max_resident_set_mib"] = kb_to_mib(out["max_resident_set_kb"])
    return out


def parse_time_log(path: Path) -> dict[str, float | int | None]:
    """Parse a ``time_scale_rm_main.log`` file if it exists."""
    if not path.exists():
        return parse_time_text("")
    return parse_time_text(path.read_text(errors="ignore"))
