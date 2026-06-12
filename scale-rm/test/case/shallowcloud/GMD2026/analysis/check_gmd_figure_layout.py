#!/usr/bin/env python3
"""Lightweight QA checks for regenerated GMD2026 candidate figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parent
import sys

sys.path.insert(0, str(ANALYSIS_DIR))

from common.table_utils import write_table_bundle  # noqa: E402


BAD_TERMS = [
    "3D representativeness proven",
    "representative cohort proven",
    "large-droplet formation mechanism",
    "coalescence episodes",
    "collision count",
    "full collision history",
    "TPHT measured full-BW reduction",
    "causal formation pathway",
]


def forbidden_hits(text: str) -> list[str]:
    """Return misleading terms while allowing explicit negative caveats."""
    lowered = text.lower()
    hits: list[str] = []
    for term in BAD_TERMS:
        term_lower = term.lower()
        if term_lower not in lowered:
            continue
        if f"not {term_lower}" in lowered:
            continue
        hits.append(term)
    return hits


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True, help="analysis_outputs_for_GMD directory.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned QA output only.")
    return parser.parse_args()


def read_manifest(outdir: Path) -> list[dict[str, str]]:
    """Read the generator manifest."""
    path = outdir / "tables" / "diagnostics" / "gmd_candidate_figure_manifest.csv"
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def add(rows: list[dict[str, Any]], stem: str, check: str, status: str, message: str) -> None:
    """Append one QA row."""
    rows.append({"figure_stem": stem, "check": check, "status": status, "message": message})


def check_files(outdir: Path, row: dict[str, str], rows: list[dict[str, Any]]) -> None:
    """Check expected figure files exist."""
    stem = row.get("figure_stem", "")
    category = row.get("category", "")
    base = outdir / "figures" / ("main_candidates" if category == "main_candidate" else "supplement_candidates")
    missing = [suffix for suffix in (".pdf", ".svg", ".png") if not (base / f"{stem}{suffix}").exists()]
    if missing:
        add(rows, stem, "saved_formats", "FAIL", f"Missing formats: {', '.join(missing)}")
    else:
        add(rows, stem, "saved_formats", "PASS", "pdf/svg/png present")


def check_manifest_row(row: dict[str, str], rows: list[dict[str, Any]]) -> None:
    """Check labels, units, wording, and candidate-specific caveats."""
    stem = row.get("figure_stem", "")
    text = " ".join(str(row.get(key, "")) for key in row)
    axis_labels = row.get("axis_labels", "")
    if not axis_labels or axis_labels == "NA":
        add(rows, stem, "axis_labels", "FAIL", "Axis labels are missing from manifest")
    elif stem != "main_candidate_framework_experiment_design" and "(" not in axis_labels:
        add(rows, stem, "axis_labels", "WARN", "No unit-like parentheses found in axis labels")
    else:
        add(rows, stem, "axis_labels", "PASS", "Axis labels recorded")

    panel_labels = row.get("panel_labels", "")
    if row.get("category") == "main_candidate" and not all(label in panel_labels for label in ("a", "b")):
        add(rows, stem, "panel_labels", "FAIL", "Main-candidate figure is missing panel labels")
    else:
        add(rows, stem, "panel_labels", "PASS", "Panel labels recorded")

    try:
        width_mm = float(row.get("width_mm", "0"))
    except ValueError:
        width_mm = 0.0
    if width_mm < 80.0:
        add(rows, stem, "figure_size", "WARN", "Figure width is below single-column target")
    else:
        add(rows, stem, "figure_size", "PASS", "Figure width is within expected range")

    try:
        tick_count = float(row.get("x_tick_count_max", "0") or 0)
    except ValueError:
        tick_count = 0.0
    try:
        tick_chars = float(row.get("x_tick_label_max_chars", "0") or 0)
    except ValueError:
        tick_chars = 0.0
    if tick_count > 8 or tick_chars > 20:
        add(rows, stem, "tick_labels", "WARN", "Potentially crowded categorical tick labels")
    else:
        add(rows, stem, "tick_labels", "PASS", "Tick-label density appears acceptable")

    bad_hits = forbidden_hits(text)
    if bad_hits:
        add(rows, stem, "misleading_wording", "FAIL", "Forbidden wording found: " + ", ".join(bad_hits))
    else:
        add(rows, stem, "misleading_wording", "PASS", "No forbidden wording in manifest")

    if stem == "main_candidate_TPHT_handoff_storage":
        marked = row.get("estimated_marked", "").lower()
        if "yes" in marked and "estimated" in marked:
            add(rows, stem, "estimated_quantity", "PASS", "Estimated full-BW output is marked")
        else:
            add(rows, stem, "estimated_quantity", "FAIL", "Estimated full-BW output is not clearly marked")

    if stem in {"main_candidate_TPHT_target_diagnostics", "main_candidate_TPHT_target_diagnostics_abc"}:
        included = row.get("unknown_category_included", "").lower()
        if included == "yes":
            add(rows, stem, "unknown_category", "PASS", "Unknown target category included")
        elif "grouped_with_radius_threshold" in included:
            add(rows, stem, "unknown_category", "PASS", "Output-level unknown category is grouped with the radius-threshold class")
        elif "source_na" in included:
            add(rows, stem, "unknown_category", "PASS", "Unknown category is source NA and is not plotted as a physical category")
        else:
            add(rows, stem, "unknown_category", "FAIL", "Unknown target category is missing")

    legend_status = row.get("legend_status", "").lower()
    if "cover" in legend_status:
        add(rows, stem, "legend", "WARN", "Legend may cover data")
    else:
        add(rows, stem, "legend", "PASS", "Legend status is acceptable")


def check_readme(outdir: Path, rows: list[dict[str, Any]]) -> None:
    """Check README caveats and numbering policy."""
    readme = outdir / "README_GMD_FIGURES.md"
    stem = "README_GMD_FIGURES"
    if not readme.exists():
        add(rows, stem, "readme", "FAIL", "README_GMD_FIGURES.md is missing")
        return
    text = readme.read_text()
    required = [
        "No final figure numbering has been assigned",
        "2D sampling-procedure verification",
        "not 3D representativeness",
        "controlled cold-start computational benchmark",
        "if_coal binary flag",
        "estimated full-BW output",
        "estimated storage reduction",
        "TPHT_ID_DIAG is absent or empty",
    ]
    missing = [phrase for phrase in required if phrase not in text]
    if missing:
        add(rows, stem, "readme_caveats", "FAIL", "Missing README caveats: " + "; ".join(missing))
    else:
        add(rows, stem, "readme_caveats", "PASS", "Required README caveats are present")
    bad_hits = forbidden_hits(text)
    if bad_hits:
        add(rows, stem, "readme_wording", "FAIL", "Forbidden wording found: " + ", ".join(bad_hits))
    else:
        add(rows, stem, "readme_wording", "PASS", "README wording caveats are acceptable")


def main() -> None:
    """CLI entry point."""
    args = parse_args()
    outdir = args.outdir.resolve()
    output_base = outdir / "tables" / "diagnostics" / "figure_layout_QA"
    if args.dry_run:
        print("[dry-run] check GMD figure layout")
        print(f"[dry-run] outdir={outdir}")
        print(f"[dry-run] would write {output_base}.csv and {output_base}.md")
        return
    manifest = read_manifest(outdir)
    qa_rows: list[dict[str, Any]] = []
    if not manifest:
        add(qa_rows, "manifest", "manifest", "FAIL", "gmd_candidate_figure_manifest.csv is missing or empty")
    for row in manifest:
        check_files(outdir, row, qa_rows)
        check_manifest_row(row, qa_rows)
    check_readme(outdir, qa_rows)
    write_table_bundle(qa_rows, output_base, ["figure_stem", "check", "status", "message"])


if __name__ == "__main__":
    main()
