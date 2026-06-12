#!/usr/bin/env python3
"""Collect SQUID model-job status records into one failure summary."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
GMD_ROOT = SCRIPT_DIR.parent
STATUS_DIR = GMD_ROOT / "job_status"

COLUMNS = [
    "case_group",
    "case_name",
    "status",
    "exit_status",
    "expected_failure",
    "pbs_jobid",
    "wallclock_s",
    "mpi_ranks",
    "omp_threads",
    "node_count",
    "case_dir",
    "job_metrics",
    "note",
]


def model_case_dirs() -> list[Path]:
    """Return model case directories that should eventually report a status."""
    return sorted(
        run_conf.parent
        for run_conf in GMD_ROOT.glob("*/*/run.conf")
        if "postprocess" not in run_conf.parts and run_conf.parent.parent.name != "common"
    )


def load_status_files() -> dict[tuple[str, str], dict[str, Any]]:
    """Load per-job status JSON files written by squid_model_job.sh."""
    records: dict[tuple[str, str], dict[str, Any]] = {}
    if not STATUS_DIR.exists():
        return records
    for path in sorted(STATUS_DIR.glob("*.json")):
        try:
            record = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        group = str(record.get("case_group", ""))
        case = str(record.get("case_name", ""))
        if group and case:
            records[(group, case)] = record
    return records


def normalize_row(group: str, case: str, record: dict[str, Any] | None, case_dir: Path) -> dict[str, Any]:
    """Build one output row from a status record or a not-run placeholder."""
    if record is None:
        return {
            "case_group": group,
            "case_name": case,
            "status": "NOT_RUN",
            "exit_status": "NA",
            "expected_failure": "NA",
            "pbs_jobid": "NA",
            "wallclock_s": "NA",
            "mpi_ranks": "NA",
            "omp_threads": "NA",
            "node_count": "NA",
            "case_dir": str(case_dir),
            "job_metrics": str(case_dir / "job_metrics.json"),
            "note": "No job_status record was found.",
        }
    row = {column: record.get(column, "NA") for column in COLUMNS}
    row["note"] = ""
    if row["status"] == "FAIL":
        row["note"] = "Inspect LOG.pe*, *.e*, and time_scale_rm_main.log in the case directory."
    if row["status"] == "UNEXPECTED_PASS":
        row["note"] = "This case was expected to fail but exited with status 0."
    if row["status"] == "EXPECTED_FAILURE":
        row["note"] = "Intentional fail-fast test behaved as expected."
    return row


def write_outputs(rows: list[dict[str, Any]]) -> None:
    """Write CSV, JSON, and Markdown summaries."""
    STATUS_DIR.mkdir(exist_ok=True)
    csv_path = GMD_ROOT / "job_failure_summary.csv"
    json_path = GMD_ROOT / "job_failure_summary.json"
    md_path = GMD_ROOT / "job_failure_summary.md"

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(rows, indent=2))
    md_lines = [
        "| " + " | ".join(COLUMNS) + " |",
        "| " + " | ".join(["---"] * len(COLUMNS)) + " |",
    ]
    for row in rows:
        md_lines.append("| " + " | ".join(str(row.get(column, "NA")) for column in COLUMNS) + " |")
    md_path.write_text("\n".join(md_lines) + "\n")


def main() -> int:
    """Entry point."""
    records = load_status_files()
    rows: list[dict[str, Any]] = []
    for case_dir in model_case_dirs():
        group = case_dir.parent.name
        case = case_dir.name
        rows.append(normalize_row(group, case, records.get((group, case)), case_dir))

    write_outputs(rows)
    fail_count = sum(row["status"] in {"FAIL", "UNEXPECTED_PASS"} for row in rows)
    not_run_count = sum(row["status"] == "NOT_RUN" for row in rows)
    expected_failure_count = sum(row["status"] == "EXPECTED_FAILURE" for row in rows)
    print(f"FAIL: {fail_count} cases")
    print(f"EXPECTED_FAILURE: {expected_failure_count} cases")
    print(f"NOT_RUN: {not_run_count} cases")
    return 1 if fail_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
