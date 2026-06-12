"""Parse SQUID job_metrics.json files for GMD2026 cases."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def parse_job_metrics(path: Path) -> dict[str, Any]:
    """Read ``job_metrics.json`` and compute core-hours when possible."""
    out: dict[str, Any] = {
        "wallclock_s": None,
        "mpi_ranks": None,
        "omp_threads": None,
        "node_count": None,
        "core_hours": None,
    }
    if not path.exists():
        return out

    try:
        raw = json.loads(path.read_text())
    except json.JSONDecodeError:
        return out

    for key in out:
        if key in raw:
            out[key] = raw[key]

    if out["core_hours"] is None:
        try:
            wallclock = float(out["wallclock_s"])
            mpi_ranks = float(out["mpi_ranks"])
            omp_threads = float(out["omp_threads"] or 1.0)
        except (TypeError, ValueError):
            return out
        out["core_hours"] = wallclock * mpi_ranks * omp_threads / 3600.0
    return out

