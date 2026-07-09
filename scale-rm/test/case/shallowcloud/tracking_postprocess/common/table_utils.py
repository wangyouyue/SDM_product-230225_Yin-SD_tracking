"""Table writers and NA-safe helpers for tracking analysis outputs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

NA = "NA"


def is_missing(value: Any) -> bool:
    """Return True when a value should be serialized as a missing value."""
    if value is None or value == NA:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def safe_float(value: Any) -> float | None:
    """Convert a table value to float, returning None for NA-like values."""
    if is_missing(value):
        return None
    try:
        output = float(value)
    except (TypeError, ValueError):
        return None
    return output if math.isfinite(output) else None


def safe_ratio(numerator: Any, denominator: Any) -> float | None:
    """Compute a ratio without converting missing denominators into zero."""
    num = safe_float(numerator)
    den = safe_float(denominator)
    if num is None or den is None or den == 0.0:
        return None
    return num / den


def clean_value(value: Any) -> Any:
    """Prepare a value for CSV, Markdown, LaTeX, and JSON serialization."""
    if is_missing(value):
        return NA
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple, set)):
        return json.dumps([clean_value(item) for item in value], sort_keys=True)
    if isinstance(value, dict):
        return json.dumps({str(key): clean_value(val) for key, val in value.items()}, sort_keys=True)
    if isinstance(value, float):
        return value if math.isfinite(value) else NA
    return value


def clean_row(row: dict[str, Any], columns: Iterable[str]) -> dict[str, Any]:
    """Return a row with all requested columns and no Python None values."""
    return {column: clean_value(row.get(column)) for column in columns}


def infer_columns(rows: list[dict[str, Any]], preferred: Iterable[str] | None = None) -> list[str]:
    """Infer stable output columns from preferred columns and row keys."""
    columns: list[str] = []
    for column in preferred or []:
        if column not in columns:
            columns.append(column)
    for row in rows:
        for column in row:
            if column not in columns:
                columns.append(column)
    return columns


def latex_escape(value: Any) -> str:
    """Escape a compact LaTeX table cell."""
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for src, dst in replacements.items():
        text = text.replace(src, dst)
    return text


def write_table_bundle(
    rows: list[dict[str, Any]] | dict[str, Any],
    base_path: Path,
    columns: Iterable[str] | None = None,
) -> None:
    """Write a table as CSV, JSON, Markdown, and LaTeX."""
    row_list = [rows] if isinstance(rows, dict) else list(rows)
    fieldnames = list(columns) if columns is not None else infer_columns(row_list, None)
    cleaned = [clean_row(row, fieldnames) for row in row_list]

    base_path.parent.mkdir(parents=True, exist_ok=True)
    with base_path.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cleaned)

    base_path.with_suffix(".json").write_text(json.dumps(cleaned, indent=2, sort_keys=False) + "\n")

    md_lines = [
        "| " + " | ".join(fieldnames) + " |",
        "| " + " | ".join(["---"] * len(fieldnames)) + " |",
    ]
    for row in cleaned:
        md_lines.append("| " + " | ".join(str(row.get(column, NA)) for column in fieldnames) + " |")
    base_path.with_suffix(".md").write_text("\n".join(md_lines) + "\n")

    tex_lines = [
        "\\begin{tabular}{" + "l" * max(1, len(fieldnames)) + "}",
        "\\hline",
        " & ".join(latex_escape(column) for column in fieldnames) + r" \\",
        "\\hline",
    ]
    for row in cleaned:
        tex_lines.append(" & ".join(latex_escape(row.get(column, NA)) for column in fieldnames) + r" \\")
    tex_lines.extend(["\\hline", "\\end{tabular}"])
    base_path.with_suffix(".tex").write_text("\n".join(tex_lines) + "\n")

