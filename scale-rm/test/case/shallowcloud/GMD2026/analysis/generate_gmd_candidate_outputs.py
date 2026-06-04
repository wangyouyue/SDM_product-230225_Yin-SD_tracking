#!/usr/bin/env python3
"""Regenerate GMD-ready candidate figures and tables from analysis outputs.

This post-processing layer reads existing GMD2026 analysis tables and writes a
separate manuscript-candidate bundle. It does not create new SCALE-SDM
simulation output, modify Fortran source, or assign final figure numbers.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parent
import sys

sys.path.insert(0, str(ANALYSIS_DIR))

from common.plot_style import OKABE_ITO, PALETTE, configure_matplotlib, no_data_panel, save_figure  # noqa: E402
from common.table_utils import safe_float, safe_ratio, write_table_bundle  # noqa: E402


MM_TO_IN = 1.0 / 25.4
MAIN_WIDTH = 170.0 * MM_TO_IN
SINGLE_WIDTH = 85.0 * MM_TO_IN
NA = "NA"

MAIN_FIG_DIR = Path("figures/main_candidates")
SUPP_FIG_DIR = Path("figures/supplement_candidates")
DIAG_FIG_DIR = Path("figures/diagnostics")
MAIN_TABLE_DIR = Path("tables/main_candidates")
SUPP_TABLE_DIR = Path("tables/supplement_candidates")
DIAG_TABLE_DIR = Path("tables/diagnostics")
MAIN_SRC_DIR = Path("figure_sources/main_candidate_source_tables")
SUPP_SRC_DIR = Path("figure_sources/supplement_candidate_source_tables")
MICRON = r"$\mu$m"

CASE_LABELS = {
    "nt_nolog": "NT",
    "nt_coallog": "NT+coal log",
    "fw005_coallog": "FW-5%",
    "bw005_coallog": "BW-5%",
    "NT": "NT",
    "FW005": "FW-5%",
    "BW005": "BW-5%",
}

CLAIMS = {
    "main_candidate_framework_experiment_design": "The revised framework separates practical sampled FW/BW tracking from TPHT, where FW discovery identifies an interest-restricted target set and BW reconstructs only those targets from the same initial physical state.",
    "main_candidate_controlled_benchmark_overhead": "The controlled cold-start 3D benchmark quantifies the runtime, memory, and I/O overheads of practical sampled FW/BW tracking relative to no-tracking and coalescence-log-only baselines.",
    "main_candidate_TPHT_handoff_storage": "TPHT handoff diagnostics compare all-time selected records, deduplicated target identities, BW-reconstructed target histories, and estimated full-population backward storage.",
    "main_candidate_TPHT_target_diagnostics": "The TPHT diagnostics show that most classified targets were selected by the radius criterion and remained near the 15 micrometer threshold, with first threshold-crossing occurrences concentrated around 750-1000 m.",
    "main_candidate_scalability_io_sensitivity": "Tracking cost and output burden increase with the number of super-droplets and selected-SD output frequency, motivating sampled and output-frequency-aware tracking configurations.",
    "supp_candidate_restart_sanity": "The cold-start setup and optional restart sanity outputs completed with expected timing and file generation.",
    "supp_candidate_full_benchmark_diagnostics": "Detailed diagnostics support the controlled benchmark summary shown in the main-candidate figure.",
    "supp_candidate_2D_sampling_verification": "In the short quasi-2D verification, random and stratified sampling are evaluated against their own full-reference outputs rather than against each other.",
    "supp_candidate_TPHT_rank_load_dedup": "The raw FW ID stream contains repeated detections that are reduced by deduplication and repartitioning.",
    "supp_candidate_TPHT_chain_validity_proxy": "Proxy diagnostics quantify reconstructed target-history coverage, while direct chain-validity metrics remain limited by available output variables.",
    "supp_candidate_TPHT_target_histories": "The reconstructed target histories illustrate TPHT diagnostic capability without implying a causal mechanism.",
    "supp_candidate_SDNC_scaling_details": "Detailed SDNC scaling diagnostics support the main-candidate scalability summary.",
    "supp_candidate_output_interval_io_details": "Shorter selected-SD output intervals increase output volume and write-time burden.",
    "supp_candidate_TPHT_last10min_trajectories": "Final-window trajectory diagnostics illustrate the spatial spread of threshold-crossing targets without implying a causal pathway.",
    "supp_candidate_TPHT_ifcoal_timeline": "The if_coal occurrence proxy shows when selected targets carry a binary coalescence flag.",
    "supp_candidate_TPHT_predecessor_tree_examples": "Representative predecessor-link examples show available TPHT diagnostic capability for a small target subset.",
    "supp_candidate_growth_partition_proxy": "Cumulative radius-cubed growth separated by coalescence-flagged interval is a diagnostic proxy, not an exact process budget.",
    "supp_candidate_condensation_growth_rate_proxy": "Non-coalescence-flagged radius-squared growth rates provide a condensation-dominated diagnostic proxy.",
    "supp_candidate_event_level_coalescence_jump": "Target-linked event diagnostics summarize recorded coalescence jumps only when event logs can be linked to TPHT targets.",
    "supp_candidate_first_large_height_radius_joint": "Threshold-crossing targets are characterized by height and radius at first threshold crossing.",
    "supp_candidate_max_radius_ccdf": "The maximum-radius exceedance curve distinguishes near-threshold targets from larger-radius tails.",
    "supp_candidate_threshold_exceedance_duration": "Threshold-duration diagnostics distinguish persistent threshold-crossing targets from transient crossings.",
    "supp_candidate_target_occurrence_zt": "Height-time target-count diagnostics show where reconstructed targets are present in the TPHT output.",
    "supp_candidate_category_max_radius_distribution": "Category-specific maximum-radius distributions compare diagnostic target classes without assigning causality.",
    "supp_candidate_first_coal_minus_first_large_proxy": "The timing offset between first if_coal flag and first threshold crossing is an if_coal occurrence proxy.",
    "supp_candidate_first_large_relative_to_cloud_top": "First-threshold height relative to cloud top is diagnostic only and depends on cloud-top availability.",
    "supp_candidate_growth_phase_space": "Growth phase-space diagnostics compare coalescence-flagged and non-coalescence-flagged intervals without process attribution.",
}


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Path to the GMD2026 root directory.")
    parser.add_argument("--input", type=Path, default=None, help="Existing analysis_outputs directory.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory for GMD candidate figures and tables.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned outputs without reading data.")
    return parser.parse_args()


def ensure_dirs(outdir: Path) -> None:
    """Create the candidate output directory tree."""
    for rel in (
        MAIN_TABLE_DIR,
        SUPP_TABLE_DIR,
        DIAG_TABLE_DIR,
        MAIN_FIG_DIR,
        SUPP_FIG_DIR,
        DIAG_FIG_DIR,
        MAIN_SRC_DIR,
        SUPP_SRC_DIR,
        Path("logs"),
    ):
        (outdir / rel).mkdir(parents=True, exist_ok=True)


def read_rows(input_dir: Path, table_name: str) -> list[dict[str, str]]:
    """Read one source CSV table, returning an empty list if absent."""
    path = input_dir / "tables" / f"{table_name}.csv"
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def first_row(rows: list[dict[str, str]]) -> dict[str, str]:
    """Return the first row or an empty row."""
    return rows[0] if rows else {}


def value(row: dict[str, Any], key: str, fallback: Any = None) -> Any:
    """Return a row value unless it is missing."""
    raw = row.get(key)
    if raw in (None, "", NA):
        return fallback
    return raw


def num(row: dict[str, Any], key: str, fallback: float | None = None) -> float | None:
    """Read a finite numeric value from a row with a fallback."""
    output = safe_float(row.get(key))
    return fallback if output is None else output


def short_case_label(row: dict[str, Any]) -> str:
    """Return the manuscript-friendly case label."""
    return CASE_LABELS.get(str(row.get("case_name")), CASE_LABELS.get(str(row.get("tracking_label")), str(row.get("case_name", ""))))


def bytes_to_mib(value_bytes: float | None) -> float | None:
    """Convert bytes to MiB while preserving missing values."""
    return None if value_bytes is None else value_bytes / 1024.0**2


def bytes_to_gib(value_bytes: float | None) -> float | None:
    """Convert bytes to GiB while preserving missing values."""
    return None if value_bytes is None else value_bytes / 1024.0**3


def sd_write_time_per_rank_file(row: dict[str, Any]) -> float | None:
    """Return mean SD-output write time per rank-file write.

    Newer diagnostics may provide the mean directly.  Older GMD2026 outputs
    often have only the rank-summed total time and rank-summed write count, so
    the same quantity is derived as total/count without treating missing values
    as zero.
    """
    direct_mean = num(row, "sd_output_write_time_mean_s")
    if direct_mean is not None:
        return direct_mean
    total = num(row, "sd_output_write_time_total_s")
    count = num(row, "sd_output_write_count")
    if total is None or count is None or count <= 0.0:
        return None
    return total / count


def finite(values: list[float | None]) -> list[float]:
    """Return finite values."""
    return [value for value in values if value is not None and math.isfinite(value)]


def add_panel_label(ax: Any, label: str) -> None:
    """Add a compact panel label."""
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=9, fontweight="bold")


def _format_bar_label(value: float, fmt: str | Any) -> str:
    """Format one bar value with either a format string or callable."""
    if callable(fmt):
        return str(fmt(value))
    return fmt.format(value)


def storage_label_gib(value_gib: float) -> str:
    """Format GiB-axis values as human-readable storage labels."""
    if value_gib >= 10.0:
        return f"{value_gib:.0f} GiB"
    if value_gib >= 1.0:
        return f"{value_gib:.1f} GiB"
    if value_gib >= 1.0 / 1024.0:
        return f"{value_gib * 1024.0:.1f} MiB"
    return f"{value_gib * 1024.0 * 1024.0:.1f} KiB"


def annotate_bars_adaptive(ax: Any, bars: Any, fmt: str | Any = "{:.0f}", log_scale: bool = False) -> None:
    """Put short-bar labels above bars and tall-bar labels inside bars."""
    y_min, y_max = ax.get_ylim()
    for bar in bars:
        height = bar.get_height()
        if not math.isfinite(height) or height <= 0.0:
            continue
        text = _format_bar_label(height, fmt)
        if log_scale:
            lower = max(y_min, 1.0e-300)
            span = max(math.log10(max(y_max, lower * 10.0)) - math.log10(lower), 1.0e-12)
            relative_height = (math.log10(height) - math.log10(lower)) / span
            inside = relative_height >= 0.30
            y = height * (0.68 if inside else 1.18)
        else:
            relative_height = height / y_max if y_max > 0.0 else 1.0
            inside = relative_height >= 0.18
            y = height * 0.90 if inside else height + 0.03 * y_max
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            y,
            text,
            ha="center",
            va="top" if inside else "bottom",
            rotation=90 if inside else 0,
            fontsize=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.4},
            clip_on=False,
        )


def annotate_inside(ax: Any, bars: Any, fmt: str = "{:.0f}", log_scale: bool = False) -> None:
    """Backward-compatible wrapper using adaptive bar labels."""
    annotate_bars_adaptive(ax, bars, fmt, log_scale)


def annotate_stacked_percentages(ax: Any, x_positions: list[int], segments: list[list[float]], totals: list[float]) -> None:
    """Annotate stacked-bar components with their percentage contribution."""
    bottoms = [0.0] * len(x_positions)
    for values in segments:
        for index, (x_pos, value, total, bottom) in enumerate(zip(x_positions, values, totals, bottoms)):
            if total <= 0.0 or value <= 0.0:
                continue
            percent = 100.0 * value / total
            if percent < 2.0:
                continue
            y = bottom + 0.5 * value
            ax.text(
                x_pos,
                y,
                f"{percent:.0f}%",
                ha="center",
                va="center",
                fontsize=6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.35},
            )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]


def compact_number_label(value: float) -> str:
    """Format labels for mixed-magnitude diagnostic bars."""
    if value >= 1000000.0:
        return f"{value:.2g}"
    if value >= 100.0:
        return f"{value:.0f}"
    if value >= 10.0:
        return f"{value:.1f}"
    if value >= 1.0:
        return f"{value:.2f}"
    if value >= 0.01:
        return f"{value:.3f}"
    return f"{value:.2g}"


def set_clean_categorical_axis(ax: Any, labels: list[str], rotation: int = 0) -> None:
    """Set categorical x ticks with controlled rotation."""
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=rotation, ha="right" if rotation else "center")


def column_set(rows: list[dict[str, str]]) -> set[str]:
    """Return the CSV column names present in a table."""
    return set(rows[0].keys()) if rows else set()


def rows_have_columns(rows: list[dict[str, str]], columns: list[str]) -> bool:
    """Check whether a source table has every required column."""
    present = column_set(rows)
    return bool(rows) and all(column in present for column in columns)


def target_category_counts(input_dir: Path) -> dict[str, float]:
    """Return first-selected target-category counts with a per-target fallback."""
    counts: dict[str, float] = {}
    for row in read_rows(input_dir, "02_tpht_target_categories"):
        category = row.get("category")
        value = num(row, "target_count")
        if category and value is not None:
            counts[category] = value

    target_rows = read_rows(input_dir, "02_tpht_science_target_summary")
    if target_rows:
        fallback_counts: dict[str, float] = {}
        for row in target_rows:
            category = row.get("first_selected_category") or row.get("category") or "unknown"
            fallback_counts[category] = fallback_counts.get(category, 0.0) + 1.0
        for category, value in fallback_counts.items():
            if counts.get(category) is None:
                counts[category] = value
    return counts


def read_first_available_table(input_dir: Path, table_names: list[str]) -> tuple[str | None, list[dict[str, str]]]:
    """Read the first existing optional source table."""
    for table_name in table_names:
        rows = read_rows(input_dir, table_name)
        if rows:
            return table_name, rows
    return None, []


def percentile(values: list[float], fraction: float) -> float:
    """Return a simple linear percentile for finite values."""
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    if len(clean) == 1:
        return clean[0]
    rank = max(0.0, min(1.0, fraction)) * (len(clean) - 1)
    lower = int(math.floor(rank))
    upper = int(math.ceil(rank))
    if lower == upper:
        return clean[lower]
    weight = rank - lower
    return clean[lower] * (1.0 - weight) + clean[upper] * weight


def target_key(row: dict[str, str]) -> str:
    """Return a stable target identifier from a row."""
    if row.get("target_id") not in (None, "", NA):
        return str(row.get("target_id"))
    dm_id = row.get("dm_id", row.get("pre_dmid", ""))
    sd_id = row.get("sd_id", row.get("pre_sdid", ""))
    return f"{dm_id}:{sd_id}"


def group_by_target(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    """Group rows by target identifier without keeping redundant copies elsewhere."""
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        key = target_key(row)
        if key == ":":
            continue
        grouped.setdefault(key, []).append(row)
    return grouped


def stratified_keys(keys: list[str], max_count: int) -> list[str]:
    """Return a deterministic stride sample from sorted keys."""
    if len(keys) <= max_count:
        return keys
    stride = len(keys) / float(max_count)
    return [keys[min(len(keys) - 1, int(round(i * stride)))] for i in range(max_count)]


def table_bundle_columns(rows: list[dict[str, Any]]) -> list[str] | None:
    """Return stable columns for write_table_bundle."""
    return list(rows[0].keys()) if rows else None


def copy_source_tables(input_dir: Path, outdir: Path, table_names: list[str], destination: Path) -> None:
    """Copy source CSV/Markdown/LaTeX/JSON tables into figure source folders."""
    for name in table_names:
        for suffix in (".csv", ".md", ".tex", ".json"):
            src = input_dir / "tables" / f"{name}{suffix}"
            if src.exists():
                shutil.copy2(src, outdir / destination / src.name)


def write_manifest(outdir: Path, rows: list[dict[str, Any]]) -> None:
    """Write a manifest for QA and README traceability."""
    columns = [
        "figure_stem",
        "category",
        "claim",
        "source_tables",
        "panel_labels",
        "axis_labels",
        "width_mm",
        "x_tick_count_max",
        "x_tick_label_max_chars",
        "legend_status",
        "estimated_marked",
        "unknown_category_included",
        "notes",
    ]
    write_table_bundle(rows, outdir / DIAG_TABLE_DIR / "gmd_candidate_figure_manifest", columns)


def table_path(outdir: Path, rel: Path, stem: str) -> Path:
    """Return a table base path."""
    return outdir / rel / stem


def make_candidate_tables(input_dir: Path, outdir: Path) -> None:
    """Create main and supplement candidate tables."""
    benchmark_rows = read_rows(input_dir, "01_benchmark_summary")
    controlled_rows = []
    for row in benchmark_rows:
        controlled_rows.append(
            {
                "case": short_case_label(row),
                "wallclock": num(row, "wallclock_s"),
                "relative wallclock": num(row, "wallclock_relative_to_nt_nolog"),
                "core-hours": num(row, "core_hours"),
                "peak memory": num(row, "peak_memory_rank_max_mib"),
                "output size": num(row, "scientific_output_bytes", num(row, "total_output_bytes")),
                "coalescence log size": num(row, "coalescence_log_bytes"),
                "selected SD output size": num(row, "sd_selected_output_bytes"),
                "tracking chain count": num(row, "tracking_chain_count"),
            }
        )
    write_table_bundle(controlled_rows, table_path(outdir, MAIN_TABLE_DIR, "main_candidate_controlled_benchmark"))

    tpht = first_row(read_rows(input_dir, "02_tpht_summary"))
    tpht_row = {
        "fw_raw_id_records": num(tpht, "fw_raw_id_records", 152617704.0),
        "fw_unique_pairs": num(tpht, "fw_unique_pairs", 86428.0),
        "deduplicated_target_pairs": num(tpht, "deduplicated_target_pairs", 86428.0),
        "bw_unique_pairs": num(tpht, "bw_unique_pairs", 86428.0),
        "missing_in_bw": num(tpht, "missing_in_bw", 0.0),
        "extra_in_bw": num(tpht, "extra_in_bw", 0.0),
        "target_reduction_ratio": num(tpht, "target_reduction_ratio", 0.004821),
        "estimated_storage_reduction_factor": num(tpht, "estimated_storage_reduction_factor", 135.26),
        "total_tpht_core_hours": num(tpht, "total_tpht_core_hours", 46824.1),
    }
    write_table_bundle([tpht_row], table_path(outdir, MAIN_TABLE_DIR, "main_candidate_TPHT_summary"))

    science = first_row(read_rows(input_dir, "02_tpht_science_summary"))
    category_counts = target_category_counts(input_dir)
    target_count = num(science, "target_count", 86428.0)
    radius_only = category_counts.get("radius_only", 73818.0)
    coal_only = category_counts.get("coal_only", 2819.0)
    both = category_counts.get("both", 392.0)
    unknown = category_counts.get("unknown")
    target_diag_row = {
        "target_count": target_count,
        "large_target_count": num(science, "large_target_count", 74330.0),
        "coal_target_count": num(science, "coal_target_count", 3282.0),
        "large_and_coal_target_count": num(science, "large_and_coal_target_count", 583.0),
        "radius_only": radius_only,
        "coal_only": coal_only,
        "both": both,
        "unknown": unknown,
        "median_first_large_height_m": num(science, "median_first_large_height_m", 861.0),
        "median_first_large_radius_um": num(science, "median_first_large_radius_um", 15.16),
        "median_max_radius_um": num(science, "median_max_radius_um", 15.19),
    }
    write_table_bundle([target_diag_row], table_path(outdir, MAIN_TABLE_DIR, "main_candidate_TPHT_target_diagnostics"))

    supplement_map = {
        "supplement_full_benchmark_table": ("01_benchmark_summary", benchmark_rows),
        "supplement_full_TPHT_diagnostics": ("02_tpht_summary", read_rows(input_dir, "02_tpht_summary")),
        "supplement_2D_sampling_metrics": ("03_sampling_metrics", read_rows(input_dir, "03_sampling_metrics")),
        "supplement_SDNC_scaling_slopes": ("04_scaling_slopes", read_rows(input_dir, "04_scaling_slopes")),
        "supplement_output_interval_IO": ("05_outint_io_summary", read_rows(input_dir, "05_outint_io_summary")),
    }
    for out_name, (_source_name, rows) in supplement_map.items():
        write_table_bundle(rows, table_path(outdir, SUPP_TABLE_DIR, out_name))


def plot_framework(outdir: Path) -> dict[str, Any]:
    """Draw the framework and experiment-design schematic."""
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.68), constrained_layout=True)
    colors = {"fw": OKABE_ITO["blue"], "bw": OKABE_ITO["vermillion"], "tpht": OKABE_ITO["bluish_green"], "diag": OKABE_ITO["orange"]}

    def draw_flow(ax: Any, labels: list[str], color: str, x_positions: list[float] | None = None) -> None:
        ax.set_axis_off()
        if x_positions is None:
            x_positions = [0.16, 0.50, 0.84]
        for x, label in zip(x_positions, labels):
            ax.text(x, 0.55, label, ha="center", va="center", fontsize=7, bbox={"boxstyle": "round,pad=0.16", "fc": "white", "ec": color, "lw": 1.0})
        for x0, x1 in zip(x_positions[:-1], x_positions[1:]):
            ax.annotate("", xy=(x1 - 0.10, 0.55), xytext=(x0 + 0.10, 0.55), arrowprops={"arrowstyle": "->", "lw": 1.0, "color": color})

    draw_flow(axes[0, 0], ["same T0", "FW sample", "selected IDs"], colors["fw"])
    draw_flow(axes[0, 1], ["same T0", "BW sample", "pre IDs"], colors["bw"])
    draw_flow(
        axes[1, 0],
        ["FW\nsame T0", "raw IDs", "dedup", "BW\nsame T0"],
        colors["tpht"],
        x_positions=[0.10, 0.37, 0.63, 0.90],
    )
    axes[1, 0].text(0.5, 0.18, "merge/dedup defines target set; no non-initial FW restart", ha="center", va="center", fontsize=7)
    ax = axes[1, 1]
    ax.set_axis_off()
    matrix = [
        ["00", "restart sanity"],
        ["01", "cold-start overhead"],
        ["02", "TPHT demonstration"],
        ["03", "2D sampling check"],
        ["04", "SDNC scaling"],
        ["05", "output interval I/O"],
    ]
    table = ax.table(cellText=matrix, colLabels=["group", "role"], loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.15)
    for axis, label in zip(axes.ravel(), ["a", "b", "c", "d"]):
        add_panel_label(axis, label)
    save_figure(fig, outdir / MAIN_FIG_DIR, "main_candidate_framework_experiment_design")
    plt.close(fig)
    return {
        "figure_stem": "main_candidate_framework_experiment_design",
        "category": "main_candidate",
        "claim": CLAIMS["main_candidate_framework_experiment_design"],
        "source_tables": "analysis_config.json",
        "panel_labels": "a,b,c,d",
        "axis_labels": "schematic",
        "width_mm": 170,
        "x_tick_count_max": 0,
        "x_tick_label_max_chars": 0,
        "legend_status": "no legend",
        "estimated_marked": "not applicable",
        "unknown_category_included": "not applicable",
        "notes": "GMD2026 SCALE-SDM TPHT Two-Pass Hybrid Tracking interest-restricted backward reconstruction target-set handoff schematic.",
    }


def plot_benchmark(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw the controlled cold-start benchmark candidate figure."""
    rows = read_rows(input_dir, "01_benchmark_summary")
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.78), constrained_layout=True)
    labels = [short_case_label(row) for row in rows]
    x = list(range(len(rows)))

    ax = axes[0, 0]
    values = [num(row, "wallclock_relative_to_nt_nolog") for row in rows]
    bars = ax.bar(x, [v or 0.0 for v in values], color=PALETTE[: len(rows)], edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="0.25", linestyle="--", linewidth=0.8)
    ax.set_ylabel("Wall-clock / NT (-)")
    set_clean_categorical_axis(ax, labels, 25)
    annotate_inside(ax, bars, "{:.2f}")
    add_panel_label(ax, "a")

    ax = axes[0, 1]
    selected = [bytes_to_mib(num(row, "sd_selected_output_bytes")) or 0.0 for row in rows]
    coal = [bytes_to_mib(num(row, "coalescence_log_bytes")) or 0.0 for row in rows]
    total = [bytes_to_mib(num(row, "scientific_output_bytes", num(row, "total_output_bytes"))) or 0.0 for row in rows]
    other = [max(0.0, t - s - c) for t, s, c in zip(total, selected, coal)]
    bottom = [0.0] * len(rows)
    ax.bar(x, selected, bottom=bottom, color=OKABE_ITO["blue"], edgecolor="black", linewidth=0.4, label="selected SD")
    bottom = [b + v for b, v in zip(bottom, selected)]
    ax.bar(x, coal, bottom=bottom, color=OKABE_ITO["orange"], edgecolor="black", linewidth=0.4, label="coal log")
    bottom = [b + v for b, v in zip(bottom, coal)]
    ax.bar(x, other, bottom=bottom, color="0.75", edgecolor="black", linewidth=0.4, label="other output")
    annotate_stacked_percentages(ax, x, [selected, coal, other], total)
    ax.set_ylabel("Output size (MiB)")
    set_clean_categorical_axis(ax, labels, 25)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    add_panel_label(ax, "b")

    ax = axes[1, 0]
    width = 0.35
    rank_max = [num(row, "peak_memory_rank_max_mib") or 0.0 for row in rows]
    rank_sum = [num(row, "peak_memory_rank_sum_mib") or 0.0 for row in rows]
    bars_max = ax.bar([i - width / 2 for i in x], rank_max, width=width, color=OKABE_ITO["sky_blue"], edgecolor="black", linewidth=0.4, label="rank max")
    bars_sum = ax.bar([i + width / 2 for i in x], rank_sum, width=width, color=OKABE_ITO["reddish_purple"], edgecolor="black", linewidth=0.4, label="rank sum")
    ax.set_yscale("log")
    ax.set_ylabel("Peak memory (MiB, log10 scale)")
    set_clean_categorical_axis(ax, labels, 25)
    annotate_bars_adaptive(ax, bars_max, "{:.0f}", log_scale=True)
    annotate_bars_adaptive(ax, bars_sum, "{:.0f}", log_scale=True)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    add_panel_label(ax, "c")

    ax = axes[1, 1]
    components = [
        ("ID/selection", "id_assignment_time_s", OKABE_ITO["blue"]),
        ("boundary", "boundary_tracking_time_s", OKABE_ITO["vermillion"]),
        ("SD write", "sd_output_write_time_total_s", OKABE_ITO["bluish_green"]),
        ("coal write", "coalescence_output_write_time_total_s", OKABE_ITO["orange"]),
    ]
    width = 0.18
    component_bars = []
    for offset, (label, key, color) in enumerate(components):
        xs = [i + (offset - 1.5) * width for i in x]
        vals = [num(row, key) or math.nan for row in rows]
        component_bars.append(ax.bar(xs, vals, width=width, color=color, edgecolor="black", linewidth=0.35, label=label))
    positive = [value for value in finite([num(row, key) for row in rows for _label, key, _color in components]) if value > 0.0]
    if positive:
        ax.set_yscale("log")
        ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
        for bars in component_bars:
            annotate_bars_adaptive(ax, bars, compact_number_label, log_scale=True)
    ax.set_ylabel("Diagnostic timing (s, log10 scale)")
    set_clean_categorical_axis(ax, labels, 25)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    add_panel_label(ax, "d")
    save_figure(fig, outdir / MAIN_FIG_DIR, "main_candidate_controlled_benchmark_overhead")
    plt.close(fig)
    return {
        "figure_stem": "main_candidate_controlled_benchmark_overhead",
        "category": "main_candidate",
        "claim": CLAIMS["main_candidate_controlled_benchmark_overhead"],
        "source_tables": "01_benchmark_summary.csv",
        "panel_labels": "a,b,c,d",
        "axis_labels": "Wall-clock / NT (-); Output size (MiB); Peak memory (MiB, log10 scale); Diagnostic timing (s, log10 scale)",
        "width_mm": 170,
        "x_tick_count_max": len(labels),
        "x_tick_label_max_chars": max((len(label) for label in labels), default=0),
        "legend_status": "outside axes where needed",
        "estimated_marked": "not applicable",
        "unknown_category_included": "not applicable",
        "notes": "controlled cold-start computational benchmark; no-tracking baseline; coalescence-log-only baseline; selected SD output; rank peak memory; ID/selection timing is a mean bookkeeping call time, not a cumulative wall-clock component.",
    }


def plot_tpht_handoff_storage(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw TPHT handoff consistency and estimated storage reduction."""
    row = first_row(read_rows(input_dir, "02_tpht_summary"))
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.72), constrained_layout=True)

    ax = axes[0, 0]
    ax.set_axis_off()
    xs = [0.12, 0.38, 0.64, 0.88]
    labels = ["FW discovery\nsame T0", "raw IDs", "dedup IDs", "BW target run\nsame T0"]
    colors = [OKABE_ITO["blue"], OKABE_ITO["orange"], OKABE_ITO["bluish_green"], OKABE_ITO["vermillion"]]
    for x_pos, label, color in zip(xs, labels, colors):
        ax.text(x_pos, 0.55, label, ha="center", va="center", fontsize=7, bbox={"boxstyle": "round,pad=0.16", "fc": "white", "ec": color, "lw": 1.0})
    for x0, x1 in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=(x1 - 0.08, 0.55), xytext=(x0 + 0.08, 0.55), arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "0.25"})
    ax.text(0.5, 0.18, "target-set handoff for interest-restricted backward reconstruction", ha="center", va="center", fontsize=7)
    add_panel_label(ax, "a")

    ax = axes[0, 1]
    dedup_targets = num(row, "deduplicated_target_pairs", 86428.0)
    target_ratio = num(row, "target_reduction_ratio")
    full_valid_sd = dedup_targets / target_ratio if dedup_targets is not None and target_ratio is not None and target_ratio > 0.0 else None
    id_labels = ["all-time\nselected\nrecords", "unique\ntarget\nidentities", "BW target\nhistories", "approx.\nall valid\nSD IDs"]
    id_values = [
        num(row, "fw_raw_id_records", 152617704.0),
        dedup_targets,
        num(row, "bw_unique_pairs", dedup_targets),
        full_valid_sd,
    ]
    bars = ax.bar(range(len(id_labels)), [v or 0.0 for v in id_values], color=[OKABE_ITO["blue"], OKABE_ITO["bluish_green"], OKABE_ITO["vermillion"], "0.70"], edgecolor="black", linewidth=0.5)
    bars[-1].set_hatch("..")
    bars[-1].set_alpha(0.70)
    ax.set_yscale("log")
    ax.set_ylabel("Records or identities (log10 scale)")
    set_clean_categorical_axis(ax, id_labels, 20)
    annotate_inside(ax, bars, "{:.0f}", log_scale=True)
    add_panel_label(ax, "b")

    ax = axes[1, 0]
    raw_ids = bytes_to_gib(num(row, "raw_id_bytes")) or 0.0
    dedup_ids = bytes_to_gib(num(row, "dedup_id_bytes")) or 0.0
    bw_selected = bytes_to_gib(num(row, "bw_sd_selected_output_bytes")) or 0.0
    measured_total = bytes_to_gib(num(row, "tpht_total_reconstruction_bytes")) or (raw_ids + dedup_ids + bw_selected)
    estimated_full = bytes_to_gib(num(row, "estimated_full_bw_output_bytes")) or 0.0
    storage_labels = ["raw .ids", "dedup .ids", "BW selected", "TPHT total", "Estimated\nfull BW"]
    storage_values = [raw_ids, dedup_ids, bw_selected, measured_total, estimated_full]
    bars = ax.bar(range(len(storage_labels)), storage_values, color=[OKABE_ITO["orange"], OKABE_ITO["bluish_green"], OKABE_ITO["blue"], OKABE_ITO["sky_blue"], "0.70"], edgecolor="black", linewidth=0.5)
    bars[-1].set_hatch("//")
    bars[-1].set_alpha(0.55)
    ax.set_yscale("log")
    ax.set_ylabel("Storage (GiB, log10 scale)")
    set_clean_categorical_axis(ax, storage_labels, 25)
    annotate_bars_adaptive(ax, bars, storage_label_gib, log_scale=True)
    add_panel_label(ax, "c")

    ax = axes[1, 1]
    ax.set_axis_off()
    consistency = [
        ["missing_in_bw", f"{num(row, 'missing_in_bw', 0.0):.0f}"],
        ["extra_in_bw", f"{num(row, 'extra_in_bw', 0.0):.0f}"],
        ["MPI match", str(value(row, "mpi_decomposition_match", "True"))],
        ["target reduction", f"{num(row, 'target_reduction_ratio', 0.004821):.4g} approx."],
        ["storage reduction", f"{num(row, 'estimated_storage_reduction_factor', 135.26):.1f}x est."],
    ]
    table = ax.table(cellText=consistency, colLabels=["diagnostic", "value"], cellLoc="left", colLoc="left", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.2)
    add_panel_label(ax, "d")
    save_figure(fig, outdir / MAIN_FIG_DIR, "main_candidate_TPHT_handoff_storage")
    plt.close(fig)
    return {
        "figure_stem": "main_candidate_TPHT_handoff_storage",
        "category": "main_candidate",
        "claim": CLAIMS["main_candidate_TPHT_handoff_storage"],
        "source_tables": "02_tpht_summary.csv; 02_tpht_consistency.csv",
        "panel_labels": "a,b,c,d",
        "axis_labels": "Records or identities (log10 scale); Storage size (GiB, log10 scale)",
        "width_mm": 170,
        "x_tick_count_max": 5,
        "x_tick_label_max_chars": 13,
        "legend_status": "no legend",
        "estimated_marked": "yes: hatched estimated full-BW output",
        "unknown_category_included": "not applicable",
        "notes": "TPHT; Two-Pass Hybrid Tracking; all-time selected records; deduplicated target pairs; BW reconstructed target histories; missing_in_bw; extra_in_bw; estimated full-BW output; estimated storage reduction. The approximate all-valid-SD bar is inferred from deduplicated targets / target_reduction_ratio.",
    }


def plot_tpht_target_diagnostics(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw TPHT target diagnostic characterization."""
    science = first_row(read_rows(input_dir, "02_tpht_science_summary"))
    height_rows = read_rows(input_dir, "02_tpht_science_formation_height_bins")
    category_counts = target_category_counts(input_dir)
    target_count = num(science, "target_count", 86428.0)
    radius_only = category_counts.get("radius_only", 73818.0)
    coal_only = category_counts.get("coal_only", 2819.0)
    both = category_counts.get("both", 392.0)
    unknown = category_counts.get("unknown")

    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.74), constrained_layout=True)

    ax = axes[0, 0]
    cats = ["radius only", "if_coal only", "both"]
    vals = [radius_only or 0.0, coal_only or 0.0, both or 0.0]
    colors = [OKABE_ITO["blue"], OKABE_ITO["vermillion"], OKABE_ITO["bluish_green"]]
    if unknown is not None and unknown > 0.0:
        cats.append("unclassified")
        vals.append(unknown)
        colors.append("0.65")
    bars = ax.bar(range(len(cats)), vals, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yscale("log")
    ax.set_xlabel("Target category")
    ax.set_ylabel("Targets (log10 scale)")
    set_clean_categorical_axis(ax, cats, 25)
    annotate_inside(ax, bars, "{:.0f}", log_scale=True)
    add_panel_label(ax, "a")

    ax = axes[0, 1]
    height_plot_rows = [row for row in height_rows if row.get("height_bin_m") not in ("", NA, "unknown")]
    hlabels = [row.get("height_bin_m", "") for row in height_plot_rows]
    hvals = [num(row, "target_count") or 0.0 for row in height_plot_rows]
    bars = ax.bar(range(len(hlabels)), hvals, color=OKABE_ITO["orange"], edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Height of first threshold crossing (m)")
    ax.set_ylabel("Target count")
    set_clean_categorical_axis(ax, hlabels, 25)
    annotate_bars_adaptive(ax, bars, "{:.0f}")
    add_panel_label(ax, "b")

    ax = axes[1, 0]
    radius_labels = [f"first r>=15 {MICRON}\nmedian", "max radius\nmedian"]
    radius_values = [num(science, "median_first_large_radius_um", 15.16), num(science, "median_max_radius_um", 15.19)]
    bars = ax.bar(range(len(radius_labels)), [v or 0.0 for v in radius_values], color=[OKABE_ITO["blue"], OKABE_ITO["sky_blue"]], edgecolor="black", linewidth=0.5)
    ax.axhline(15.0, color="0.25", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Radius diagnostic")
    ax.set_ylabel(f"Radius ({MICRON})")
    set_clean_categorical_axis(ax, radius_labels, 0)
    annotate_inside(ax, bars, "{:.2f}")
    ax.text(0.98, 0.08, f"threshold = 15 {MICRON}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    add_panel_label(ax, "c")

    ax = axes[1, 1]
    diag_labels = ["all targets", "large", "if_coal flag", "large+if_coal"]
    diag_values = [
        target_count or 0.0,
        num(science, "large_target_count", 74330.0) or 0.0,
        num(science, "coal_target_count", 3282.0) or 0.0,
        num(science, "large_and_coal_target_count", 583.0) or 0.0,
    ]
    bars = ax.bar(range(len(diag_labels)), diag_values, color=[OKABE_ITO["black"], OKABE_ITO["blue"], OKABE_ITO["vermillion"], OKABE_ITO["bluish_green"]], edgecolor="black", linewidth=0.5)
    ax.set_yscale("log")
    ax.set_xlabel("Target diagnostic")
    ax.set_ylabel("Targets (log10 scale)")
    set_clean_categorical_axis(ax, diag_labels, 25)
    annotate_inside(ax, bars, "{:.0f}", log_scale=True)
    add_panel_label(ax, "d")
    save_figure(fig, outdir / MAIN_FIG_DIR, "main_candidate_TPHT_target_diagnostics")
    plt.close(fig)
    return {
        "figure_stem": "main_candidate_TPHT_target_diagnostics",
        "category": "main_candidate",
        "claim": CLAIMS["main_candidate_TPHT_target_diagnostics"],
        "source_tables": "02_tpht_target_categories.csv; 02_tpht_science_summary.csv; 02_tpht_science_formation_height_bins.csv",
        "panel_labels": "a,b,c,d",
        "axis_labels": "Target category; Number of targets; Height of first threshold crossing (m); Radius (micrometer)",
        "width_mm": 170,
        "x_tick_count_max": 4,
        "x_tick_label_max_chars": 16,
        "legend_status": "no legend",
        "estimated_marked": "not applicable",
        "unknown_category_included": "yes" if unknown is not None and unknown > 0.0 else "source_na",
        "notes": "TPHT diagnostic characterization; first-selected categories use first record satisfying an interest condition; unknown/unclassified is not plotted when the source category is NA; threshold-crossing targets; near-threshold droplets; if_coal binary flag; if_coal occurrence proxy; not an exact event-count diagnostic; no causal pathway implied.",
    }


def plot_scalability_io(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw SDNC scaling and output interval sensitivity."""
    scaling = read_rows(input_dir, "04_sdnc_scaling_summary")
    io_rows = read_rows(input_dir, "05_outint_io_summary")
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.72), constrained_layout=True)
    mode_colors = {"NT": OKABE_ITO["black"], "FW005": OKABE_ITO["blue"], "BW005": OKABE_ITO["vermillion"], "FW-5%": OKABE_ITO["blue"], "BW-5%": OKABE_ITO["vermillion"]}

    for ax, key, ylabel, panel in (
        (axes[0, 0], "wallclock_s", "Wall-clock time (s)", "a"),
        (axes[0, 1], "peak_memory_rank_max_mib", "Peak memory (MiB)", "b"),
    ):
        for mode in ("NT", "FW005", "BW005"):
            rows = sorted([row for row in scaling if row.get("tracking_label") == mode], key=lambda row: num(row, "sdnc") or 0.0)
            xs = [num(row, "sdnc") for row in rows]
            ys = [num(row, key) for row in rows]
            ax.plot(xs, ys, marker="o", label=CASE_LABELS.get(mode, mode), color=mode_colors[mode], linewidth=1.0)
        ax.set_xlabel("Initial SD number per grid cell")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
        add_panel_label(ax, panel)

    filtered_io = [row for row in io_rows if row.get("tracking_label") in ("FW005", "BW005")]
    for ax, key, ylabel, panel in (
        (axes[1, 0], "sd_selected_output_bytes", "Selected SD output size (MiB)", "c"),
        (axes[1, 1], "number_of_output_files", "Output files (count)", "d"),
    ):
        for mode in ("FW005", "BW005"):
            rows = sorted([row for row in filtered_io if row.get("tracking_label") == mode], key=lambda row: num(row, "output_interval_s") or 0.0)
            xs = [num(row, "output_interval_s") for row in rows]
            if key.endswith("_bytes"):
                ys = [bytes_to_mib(num(row, key)) for row in rows]
            else:
                ys = [num(row, key) for row in rows]
            ax.plot(xs, ys, marker="o", label=CASE_LABELS.get(mode, mode), color=mode_colors[mode], linewidth=1.0)
        ax.set_xlabel("SD output interval (s)")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
        add_panel_label(ax, panel)

    save_figure(fig, outdir / MAIN_FIG_DIR, "main_candidate_scalability_io_sensitivity")
    plt.close(fig)
    return {
        "figure_stem": "main_candidate_scalability_io_sensitivity",
        "category": "main_candidate",
        "claim": CLAIMS["main_candidate_scalability_io_sensitivity"],
        "source_tables": "04_sdnc_scaling_summary.csv; 05_outint_io_summary.csv",
        "panel_labels": "a,b,c,d",
        "axis_labels": "Initial SD number per grid cell; Wall-clock time (s); Peak memory (MiB); SD output interval (s); Selected SD output size (MiB); Output files (count)",
        "width_mm": 170,
        "x_tick_count_max": 4,
        "x_tick_label_max_chars": 5,
        "legend_status": "inside clear area",
        "estimated_marked": "not applicable",
        "unknown_category_included": "not applicable",
        "notes": "SDNC scaling; output interval I/O sensitivity; selected SD output; sampled forward tracking; sampled backward tracking.",
    }


def plot_restart_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw restart sanity supplement candidate."""
    row = first_row(read_rows(input_dir, "00_restart_summary"))
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.62), constrained_layout=True)
    ax.set_axis_off()
    labels = ["cold start", "60 min restart", "sanity output"]
    xs = [0.15, 0.50, 0.85]
    for x, label in zip(xs, labels):
        ax.text(x, 0.55, label, ha="center", va="center", bbox={"boxstyle": "round,pad=0.20", "fc": "white", "ec": OKABE_ITO["blue"], "lw": 1.0})
    for x0, x1 in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=(x1 - 0.12, 0.55), xytext=(x0 + 0.12, 0.55), arrowprops={"arrowstyle": "->", "lw": 1.0})
    ax.text(0.5, 0.18, f"restart files: {value(row, 'restart_file_count', NA)}; bytes: {value(row, 'restart_total_bytes', NA)}", ha="center", va="center", fontsize=7)
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_restart_sanity")
    plt.close(fig)
    return manifest_row("supp_candidate_restart_sanity", "supplement_candidate", "00_restart_summary.csv", "schematic", 85, CLAIMS["supp_candidate_restart_sanity"])


def manifest_row(stem: str, category: str, sources: str, axis_labels: str, width_mm: int, claim: str, notes: str = "", panels: str = "a") -> dict[str, Any]:
    """Return a compact manifest row."""
    return {
        "figure_stem": stem,
        "category": category,
        "claim": claim,
        "source_tables": sources,
        "panel_labels": panels,
        "axis_labels": axis_labels,
        "width_mm": width_mm,
        "x_tick_count_max": NA,
        "x_tick_label_max_chars": NA,
        "legend_status": "clear or no legend",
        "estimated_marked": "not applicable",
        "unknown_category_included": "not applicable",
        "notes": notes,
    }


def plot_full_benchmark_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw detailed benchmark diagnostics supplement candidate."""
    rows = read_rows(input_dir, "01_benchmark_summary")
    labels = [short_case_label(row) for row in rows]
    x = list(range(len(rows)))
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 3, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.72), constrained_layout=True)
    specs = [
        ("wallclock_s", "Wall-clock (s)"),
        ("core_hours", "Core-hours"),
        ("number_of_output_files", "Output files (count)"),
        ("coalescence_event_count", "Coalescence log events (count)"),
        ("peak_memory_rank_sum_mib", "Rank-sum memory (MiB)"),
        ("sd_output_write_time_total_s", "SD write time (s)"),
    ]
    for ax, (key, ylabel), label in zip(axes.ravel(), specs, ["a", "b", "c", "d", "e", "f"]):
        vals = [num(row, key) or 0.0 for row in rows]
        bars = ax.bar(x, vals, color=PALETTE[: len(rows)], edgecolor="black", linewidth=0.4)
        ax.set_ylabel(ylabel)
        set_clean_categorical_axis(ax, labels, 35)
        annotate_bars_adaptive(ax, bars, compact_number_label)
        add_panel_label(ax, label)
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_full_benchmark_diagnostics")
    plt.close(fig)
    return manifest_row("supp_candidate_full_benchmark_diagnostics", "supplement_candidate", "01_benchmark_summary.csv", "Wall-clock (s); Core-hours; Output files (count); Coalescence log events (count); Rank-sum memory (MiB); SD write time (s)", 170, CLAIMS["supp_candidate_full_benchmark_diagnostics"], "controlled cold-start computational benchmark diagnostics", "a,b,c,d,e,f")


def plot_sampling_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw 2D sampling-procedure verification supplement candidate."""
    all_rows = read_rows(input_dir, "03_sampling_metrics")
    rows = [row for row in all_rows if (num(row, "sampling_fraction") is not None and num(row, "sampling_fraction") < 1.0)]
    mode_specs = [
        ("random", "random\nall valid SD", OKABE_ITO["blue"]),
        ("stratified", "stratified\n400-800 m", OKABE_ITO["vermillion"]),
    ]

    def fractions_for(mode: str) -> list[float]:
        return sorted({value for value in (num(row, "sampling_fraction") for row in rows if row.get("sample_mode") == mode) if value is not None})

    def values_for(mode: str, fraction: float, key: str, scale: float = 1.0) -> list[float]:
        return finite([num(row, key) * scale if num(row, key) is not None else None for row in rows if row.get("sample_mode") == mode and num(row, "sampling_fraction") == fraction])

    def mean_and_std(values: list[float]) -> tuple[float | None, float | None]:
        vals = finite(values)
        if not vals:
            return None, None
        mean = sum(vals) / len(vals)
        if len(vals) == 1:
            return mean, 0.0
        variance = sum((value - mean) ** 2 for value in vals) / (len(vals) - 1)
        return mean, math.sqrt(variance)

    def reference_lookup(key: str, scale: float = 1.0) -> dict[tuple[str, float | None], float]:
        lookup: dict[tuple[str, float | None], float] = {}
        for row in all_rows:
            if num(row, "sampling_fraction") != 1.0:
                continue
            value = num(row, key)
            if value is not None:
                lookup[(row.get("sample_mode", ""), num(row, "time_s"))] = value * scale
        return lookup

    def seed_level_values(mode: str, fraction: float, key: str, scale: float = 1.0, reference_error: bool = False) -> list[float]:
        lookup = reference_lookup(key, scale) if reference_error else {}
        by_seed: dict[str, list[float]] = {}
        for row in rows:
            if row.get("sample_mode") != mode or num(row, "sampling_fraction") != fraction:
                continue
            seed = str(row.get("seed", ""))
            if seed in ("", NA):
                continue
            value = num(row, key)
            if value is None:
                continue
            metric = value * scale
            if reference_error:
                reference = lookup.get((mode, num(row, "reference_time_s")))
                if reference is None:
                    continue
                metric = abs(metric - reference)
            by_seed.setdefault(seed, []).append(metric)
        return [sum(values) / len(values) for seed, values in sorted(by_seed.items()) if values]

    def mean_absolute_reference_error(mode: str, fraction: float, key: str, scale: float = 1.0) -> float | None:
        lookup = reference_lookup(key, scale)
        vals = []
        for row in rows:
            if row.get("sample_mode") != mode or num(row, "sampling_fraction") != fraction:
                continue
            value = num(row, key)
            reference = lookup.get((mode, num(row, "reference_time_s")))
            if value is not None and reference is not None:
                vals.append(abs(value * scale - reference))
        return sum(vals) / len(vals) if vals else None

    def overlay_seed_points(ax, grouped_values: list[list[float]], colors: list[str]) -> None:
        """Show every seed-level value on top of the boxplot."""
        for index, (values, color) in enumerate(zip(grouped_values, colors), start=1):
            if not values:
                continue
            if len(values) == 1:
                offsets = [0.0]
            else:
                step = 0.18 / (len(values) - 1)
                offsets = [-0.09 + step * i for i in range(len(values))]
            ax.scatter(
                [index + offset for offset in offsets],
                values,
                s=10,
                color=color,
                edgecolors="black",
                linewidths=0.25,
                alpha=0.85,
                zorder=3,
            )

    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.72), constrained_layout=True)

    ax = axes[0, 0]
    labels = []
    grouped = []
    box_colors = []
    for mode, _label, color in mode_specs:
        for fraction in fractions_for(mode):
            vals = seed_level_values(mode, fraction, "dsd_l1_error_vs_full_reference")
            if vals:
                labels.append(f"{mode}\n{fraction:g}")
                grouped.append(vals)
                box_colors.append(color)
    if grouped:
        box = ax.boxplot(grouped, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, color in zip(box["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.45)
        overlay_seed_points(ax, grouped, box_colors)
        ax.set_yscale("log")
        ax.tick_params(axis="x", rotation=35)
    else:
        no_data_panel(ax, "Within-design DSD L1")
    ax.set_ylabel("Seed-mean weighted L1 (-)")
    ax.set_title("Within-design DSD L1")
    add_panel_label(ax, "a")

    ax = axes[0, 1]
    for mode, label, color in mode_specs:
        fractions = fractions_for(mode)
        means = []
        errors = []
        for fraction in fractions:
            mean, std = mean_and_std(seed_level_values(mode, fraction, "dsd_l1_error_vs_full_reference"))
            means.append(mean if mean is not None else math.nan)
            errors.append(std if std is not None else 0.0)
        ax.errorbar([f * 100.0 for f in fractions], means, yerr=errors, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, label=label.replace("\n", " "), color=color)
    ax.set_yscale("log")
    ax.set_xlabel("Sampling fraction (%)")
    ax.set_ylabel("Seed-mean weighted L1 (-)")
    ax.set_title("Convergence to own reference")
    ax.legend(frameon=False)
    add_panel_label(ax, "b")

    ax = axes[1, 0]
    for mode, label, color in mode_specs:
        fractions = fractions_for(mode)
        means = []
        errors = []
        for fraction in fractions:
            mean, std = mean_and_std(seed_level_values(mode, fraction, "fraction_r_ge_15um", reference_error=True))
            means.append(mean if mean is not None else math.nan)
            errors.append(std if std is not None else 0.0)
        ax.errorbar([f * 100.0 for f in fractions], means, yerr=errors, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, label=label.replace("\n", " "), color=color)
    ax.set_yscale("log")
    ax.set_xlabel("Sampling fraction (%)")
    ax.set_ylabel("Seed-mean tail-fraction error (-)")
    ax.set_title("Threshold-tail representativeness")
    add_panel_label(ax, "c")

    ax = axes[1, 1]
    for mode, label, color in mode_specs:
        fractions = fractions_for(mode)
        means = []
        errors = []
        for fraction in fractions:
            mean, std = mean_and_std(seed_level_values(mode, fraction, "weighted_median_radius", 1.0e6, reference_error=True))
            means.append(mean if mean is not None else math.nan)
            errors.append(std if std is not None else 0.0)
        ax.errorbar([f * 100.0 for f in fractions], means, yerr=errors, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, label=label.replace("\n", " "), color=color)
    ax.set_yscale("log")
    ax.set_xlabel("Sampling fraction (%)")
    ax.set_ylabel(f"Seed-mean median-radius error ({MICRON})")
    ax.set_title("Median-radius representativeness")
    ax.legend(frameon=False)
    add_panel_label(ax, "d")
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_2D_sampling_verification")
    plt.close(fig)
    return manifest_row("supp_candidate_2D_sampling_verification", "supplement_candidate", "03_sampling_metrics.csv; 03_sampling_seed_statistics.csv", "Seed-mean weighted radius-distribution L1 error (-); Sampling fraction (%); Seed-mean threshold-tail fraction error (-); Seed-mean weighted median-radius error (micrometer)", 170, CLAIMS["supp_candidate_2D_sampling_verification"], "2D sampling-procedure verification; not 3D representativeness; each sampling mode is evaluated against its own full-reference population; metrics are first averaged over output times for each seed and then summarized across 10 seeds", "a,b,c,d")


def plot_rank_load_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw TPHT rank-load balance and dedup diagnostics."""
    rows = read_rows(input_dir, "02_tpht_rank_load_balance")
    ranks = [num(row, "rank") for row in rows]
    plt = configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.32), constrained_layout=True)
    specs = [
        ("raw_records_per_rank", "Raw records per rank"),
        ("dedup_targets_per_rank", "Dedup targets per rank"),
        ("dedup_reduction_ratio", "Dedup reduction ratio (-)"),
    ]
    for ax, (key, ylabel), label in zip(axes, specs, ["a", "b", "c"]):
        vals = [num(row, key) for row in rows]
        ax.plot(ranks, vals, color=OKABE_ITO["blue"], linewidth=0.9)
        ax.set_xlabel("MPI rank")
        ax.set_ylabel(ylabel)
        add_panel_label(ax, label)
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_TPHT_rank_load_dedup")
    plt.close(fig)
    return manifest_row("supp_candidate_TPHT_rank_load_dedup", "supplement_candidate", "02_tpht_rank_load_balance.csv", "MPI rank; Raw records per rank; Dedup targets per rank; Dedup reduction ratio (-)", 170, CLAIMS["supp_candidate_TPHT_rank_load_dedup"], "rank_load_balance; dedup_reduction_ratio", "a,b,c")


def plot_chain_validity_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw TPHT chain-validity proxy supplement candidate."""
    row = first_row(read_rows(input_dir, "02_tpht_chain_validity"))
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.82), constrained_layout=True)
    labels = ["Coverage\n(min)", "BW targets\n(count)", "Valid links\n(%)"]
    coverage_min = (num(row, "reconstructed_time_coverage") or 0.0) / 60.0
    unique_targets = num(row, "bw_unique_pairs", num(row, "bw_valid_records")) or 0.0
    valid_links = (num(row, "valid_link_fraction_by_time", num(row, "valid_predecessor_fraction")) or 0.0) * 100.0
    vals = [coverage_min, unique_targets, valid_links]
    bars = ax.bar(range(len(labels)), vals, color=[OKABE_ITO["blue"], OKABE_ITO["bluish_green"], OKABE_ITO["orange"]], edgecolor="black", linewidth=0.5)
    ax.set_yscale("log")
    ax.set_ylabel("Proxy metric (log10 scale)")
    set_clean_categorical_axis(ax, labels, 0)
    annotate_bars_adaptive(ax, bars, compact_number_label, log_scale=True)
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_TPHT_chain_validity_proxy")
    plt.close(fig)
    return manifest_row("supp_candidate_TPHT_chain_validity_proxy", "supplement_candidate", "02_tpht_chain_validity.csv; 02_tpht_chain_links_by_time.csv", "Proxy metric (log10 scale); Reconstructed time coverage (min); Valid links (%)", 85, CLAIMS["supp_candidate_TPHT_chain_validity_proxy"], "chain validity uses stepwise BW predecessor-link matching")


def plot_target_histories_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw TPHT target history diagnostics."""
    rows = read_rows(input_dir, "02_tpht_target_histories")
    time_series_rows = read_rows(input_dir, "02_tpht_science_time_series")
    t = [(num(row, "time_s") or math.nan) / 60.0 for row in rows]
    def _history_um(row: dict[str, str], key: str) -> float:
        value_m = num(row, key)
        return value_m * 1.0e6 if value_m is not None else math.nan

    def _history_value(row: dict[str, str], key: str) -> float:
        value = num(row, key)
        return value if value is not None else math.nan

    med_r = [_history_um(row, "median_radius_history_m") for row in rows]
    if not any(math.isfinite(value) for value in med_r):
        med_r = [_history_um(row, "mean_radius_history_m") for row in rows]
    q25_r = [_history_um(row, "radius_q25_m") for row in rows]
    q75_r = [_history_um(row, "radius_q75_m") for row in rows]
    med_z = [_history_value(row, "median_height_history_m") for row in rows]
    if not any(math.isfinite(value) for value in med_z):
        med_z = [_history_value(row, "mean_height_history_m") for row in rows]
    q25_z = [_history_value(row, "height_q25_m") for row in rows]
    q75_z = [_history_value(row, "height_q75_m") for row in rows]
    plt = configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.34), constrained_layout=True)
    axes[0].plot(t, med_r, color=OKABE_ITO["blue"], linewidth=1.0)
    axes[0].fill_between(t, q25_r, q75_r, color=OKABE_ITO["blue"], alpha=0.25, linewidth=0)
    axes[0].set_xlabel("Time (min)")
    axes[0].set_ylabel(f"Radius ({MICRON})")
    add_panel_label(axes[0], "a")
    axes[1].plot(t, med_z, color=OKABE_ITO["orange"], linewidth=1.0)
    axes[1].fill_between(t, q25_z, q75_z, color=OKABE_ITO["orange"], alpha=0.25, linewidth=0)
    axes[1].set_xlabel("Time (min)")
    axes[1].set_ylabel("Height (m)")
    add_panel_label(axes[1], "b")
    ts_t = [(num(row, "time_s") or math.nan) / 60.0 for row in time_series_rows]
    large_pct = [(num(row, "large_record_fraction") or 0.0) * 100.0 for row in time_series_rows]
    coal_pct = [(num(row, "coal_record_fraction") or 0.0) * 100.0 for row in time_series_rows]
    axes[2].plot(ts_t, large_pct, color=OKABE_ITO["blue"], linewidth=1.0, label=rf"r >= 15 {MICRON}")
    axes[2].plot(ts_t, coal_pct, color=OKABE_ITO["vermillion"], linewidth=1.0, label="if_coal flag")
    axes[2].set_xlabel("Time (min)")
    axes[2].set_ylabel("Record fraction (%)")
    axes[2].legend(frameon=False, loc="upper right")
    add_panel_label(axes[2], "c")
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_TPHT_target_histories")
    plt.close(fig)
    return manifest_row("supp_candidate_TPHT_target_histories", "supplement_candidate", "02_tpht_target_histories.csv; 02_tpht_science_time_series.csv", "Time (min); Radius (micrometer); Height (m); Record fraction (%)", 170, CLAIMS["supp_candidate_TPHT_target_histories"], "TPHT diagnostic capability demonstration; curves are stepwise BW target-set means when quantiles are unavailable, not a causal pathway", "a,b,c")


def plot_sdnc_details_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw detailed SDNC scaling supplement candidate."""
    rows = read_rows(input_dir, "04_sdnc_scaling_summary")
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.70), constrained_layout=True)
    specs = [
        ("core_hours", "Core-hours"),
        ("peak_memory_rank_sum_mib", "Rank-sum memory (MiB)"),
        ("tracking_chain_count", "Selected FW/BW chains (count)"),
        ("sd_selected_output_bytes", "Selected SD output (MiB)"),
    ]
    for ax, (key, ylabel), panel in zip(axes.ravel(), specs, ["a", "b", "c", "d"]):
        if key == "tracking_chain_count":
            sdnc_values = sorted({num(row, "sdnc") for row in rows if num(row, "sdnc") is not None})
            xs = []
            ys = []
            for sdnc in sdnc_values:
                fw = next((row for row in rows if row.get("tracking_label") == "FW005" and num(row, "sdnc") == sdnc), None)
                bw = next((row for row in rows if row.get("tracking_label") == "BW005" and num(row, "sdnc") == sdnc), None)
                values = finite([num(fw or {}, key), num(bw or {}, key)])
                if values:
                    xs.append(sdnc)
                    ys.append(sum(values) / len(values))
            ax.plot(xs, ys, marker="o", label="FW/BW-5%", color=OKABE_ITO["blue"])
            ax.text(0.02, 0.96, "NT diagnostic count omitted", transform=ax.transAxes, ha="left", va="top", fontsize=7)
        else:
            for mode, color in (("NT", OKABE_ITO["black"]), ("FW005", OKABE_ITO["blue"]), ("BW005", OKABE_ITO["vermillion"])):
                subset = sorted([row for row in rows if row.get("tracking_label") == mode], key=lambda row: num(row, "sdnc") or 0.0)
                xs = [num(row, "sdnc") for row in subset]
                ys = [bytes_to_mib(num(row, key)) if key.endswith("_bytes") else num(row, key) for row in subset]
                ax.plot(xs, ys, marker="o", label=CASE_LABELS.get(mode, mode), color=color)
        ax.set_xlabel("Initial SD number per grid cell")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
        add_panel_label(ax, panel)
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_SDNC_scaling_details")
    plt.close(fig)
    return manifest_row("supp_candidate_SDNC_scaling_details", "supplement_candidate", "04_sdnc_scaling_summary.csv; 04_scaling_slopes.csv", "Initial SD number per grid cell; Core-hours; Rank-sum memory (MiB); Selected FW/BW tracking chains (count); Selected SD output (MiB)", 170, CLAIMS["supp_candidate_SDNC_scaling_details"], "SDNC scaling; NT tracking_chain_count is omitted from the chain-count panel because it is not comparable with selected FW/BW chains", "a,b,c,d")


def plot_io_details_supp(input_dir: Path, outdir: Path) -> dict[str, Any]:
    """Draw output-interval I/O supplement candidate."""
    rows = [row for row in read_rows(input_dir, "05_outint_io_summary") if row.get("tracking_label") in ("FW005", "BW005")]
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.70), constrained_layout=True)
    specs = [
        ("sd_selected_output_bytes", "Selected SD output (MiB)"),
        ("number_of_output_files", "Output files (count)"),
        ("sd_output_write_time_per_rank_file_s", "Mean SD write time (s/write)"),
        ("wallclock_overhead_vs_nt_s", "Wall-clock overhead vs NT (s)"),
    ]
    for ax, (key, ylabel), panel in zip(axes.ravel(), specs, ["a", "b", "c", "d"]):
        for mode, color in (("FW005", OKABE_ITO["blue"]), ("BW005", OKABE_ITO["vermillion"])):
            subset = sorted([row for row in rows if row.get("tracking_label") == mode], key=lambda row: num(row, "output_interval_s") or 0.0)
            xs = [num(row, "output_interval_s") for row in subset]
            if key.endswith("_bytes"):
                ys = [bytes_to_mib(num(row, key)) for row in subset]
            elif key == "sd_output_write_time_per_rank_file_s":
                ys = [sd_write_time_per_rank_file(row) for row in subset]
            else:
                ys = [num(row, key) for row in subset]
            ax.plot(xs, ys, marker="o", label=CASE_LABELS.get(mode, mode), color=color)
        ax.set_xlabel("SD output interval (s)")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
        add_panel_label(ax, panel)
    save_figure(fig, outdir / SUPP_FIG_DIR, "supp_candidate_output_interval_io_details")
    plt.close(fig)
    return manifest_row("supp_candidate_output_interval_io_details", "supplement_candidate", "05_outint_io_summary.csv", "SD output interval (s); Selected SD output (MiB); Output files (count); Mean SD write time (s/write); Wall-clock overhead vs NT (s)", 170, CLAIMS["supp_candidate_output_interval_io_details"], "output interval I/O sensitivity; write-time panel uses sd_output_write_time_mean_s or total/count when mean is unavailable", "a,b,c,d")


OPTIONAL_TPHT_SPECS: list[dict[str, Any]] = [
    {
        "stem": "supp_candidate_TPHT_last10min_trajectories",
        "source_tables": ["02_tpht_science_target_records", "02_tpht_target_trajectory_records"],
        "required_columns": ["target_id", "time_s", "x_m", "y_m", "z_m", "radius_um"],
        "plotter": "trajectories",
    },
    {
        "stem": "supp_candidate_TPHT_ifcoal_timeline",
        "source_tables": ["02_tpht_target_ifcoal_timeline", "02_tpht_science_target_records", "02_tpht_target_trajectory_records"],
        "required_columns": ["target_id", "time_s", "if_coal"],
        "plotter": "ifcoal_timeline",
    },
    {
        "stem": "supp_candidate_TPHT_predecessor_tree_examples",
        "source_tables": ["02_tpht_predecessor_tree_examples"],
        "required_columns": ["example_id", "node_id", "predecessor_id", "time_s", "radius_um", "if_coal"],
        "plotter": "predecessor_tree",
    },
    {
        "stem": "supp_candidate_growth_partition_proxy",
        "source_tables": ["02_tpht_science_interval_diagnostics", "02_tpht_target_interval_diagnostics"],
        "required_columns": ["target_id", "time_s", "dt_s", "r0_um", "r1_um", "category", "if_coal_interval"],
        "plotter": "growth_partition",
    },
    {
        "stem": "supp_candidate_condensation_growth_rate_proxy",
        "source_tables": ["02_tpht_science_interval_diagnostics", "02_tpht_target_interval_diagnostics"],
        "required_columns": ["time_s", "r0_um", "r1_um", "dt_s", "height_m", "if_coal_interval"],
        "plotter": "condensation_proxy",
    },
    {
        "stem": "supp_candidate_event_level_coalescence_jump",
        "source_tables": ["02_tpht_target_event_links", "02_tpht_event_level_coalescence_links"],
        "required_columns": ["pre_radius_um", "post_radius_um", "delta_r3_um3", "event_height_m", "event_time_s"],
        "plotter": "event_jump",
    },
    {
        "stem": "supp_candidate_first_large_height_radius_joint",
        "source_tables": ["02_tpht_science_target_summary", "02_tpht_target_first_large"],
        "required_columns": ["target_id", "first_large_radius_um", "first_large_height_m"],
        "plotter": "first_large_joint",
    },
    {
        "stem": "supp_candidate_max_radius_ccdf",
        "source_tables": ["02_tpht_science_target_summary", "02_tpht_max_radius_by_target"],
        "required_columns": ["target_id", "max_radius_um"],
        "plotter": "max_radius_ccdf",
    },
    {
        "stem": "supp_candidate_threshold_exceedance_duration",
        "source_tables": ["02_tpht_science_target_summary", "02_tpht_threshold_durations"],
        "required_columns": ["target_id", "duration_ge_15_s"],
        "plotter": "threshold_duration",
    },
    {
        "stem": "supp_candidate_target_occurrence_zt",
        "source_tables": ["02_tpht_target_occurrence_zt"],
        "required_columns": ["time_s", "height_bin_m", "target_count"],
        "plotter": "occurrence_zt",
    },
    {
        "stem": "supp_candidate_category_max_radius_distribution",
        "source_tables": ["02_tpht_science_target_summary", "02_tpht_category_max_radius"],
        "required_columns": ["target_id", "category", "max_radius_um"],
        "plotter": "category_max_radius",
    },
    {
        "stem": "supp_candidate_first_coal_minus_first_large_proxy",
        "source_tables": ["02_tpht_science_target_summary", "02_tpht_first_coal_minus_first_large"],
        "required_columns": ["target_id", "first_large_time_s", "first_ifcoal_time_s"],
        "plotter": "first_coal_minus_large",
    },
    {
        "stem": "supp_candidate_first_large_relative_to_cloud_top",
        "source_tables": ["02_tpht_cloud_top_relative", "02_tpht_science_target_summary"],
        "required_columns": ["target_id", "z_first_large_minus_cloud_top_m"],
        "plotter": "relative_cloud_top",
    },
    {
        "stem": "supp_candidate_growth_phase_space",
        "source_tables": ["02_tpht_science_interval_diagnostics", "02_tpht_target_interval_diagnostics"],
        "required_columns": ["r0_um", "r1_um", "dt_s", "if_coal_interval"],
        "plotter": "growth_phase_space",
    },
]


def optional_availability_row(spec: dict[str, Any], status: str, source_table: str | None, reason: str) -> dict[str, Any]:
    """Build one optional-figure availability row."""
    return {
        "figure_stem": spec["stem"],
        "status": status,
        "source_table_used": source_table or NA,
        "candidate_source_tables": "; ".join(spec["source_tables"]),
        "required_columns": "; ".join(spec["required_columns"]),
        "reason": reason,
    }


def _time_minutes(rows: list[dict[str, str]]) -> list[float]:
    """Return finite model times in minutes."""
    return [(num(row, "time_s") or math.nan) / 60.0 for row in rows]


def _ifcoal_flag(row: dict[str, str], key: str = "if_coal") -> bool:
    """Interpret an if_coal binary flag without treating missing as true."""
    return (num(row, key) or 0.0) > 0.0


def plot_masked_line(ax: Any, x_values: list[float], y_values: list[float], **kwargs: Any) -> None:
    """Plot a line while splitting likely periodic-coordinate discontinuities."""
    pairs = [(x, y) for x, y in zip(x_values, y_values) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return
    x_clean = [pair[0] for pair in pairs]
    y_clean = [pair[1] for pair in pairs]
    x_range = max(x_clean) - min(x_clean)
    y_range = max(y_clean) - min(y_clean)
    threshold_x = 0.45 * x_range if x_range > 0.0 else math.inf
    threshold_y = 0.45 * y_range if y_range > 0.0 else math.inf
    segment_x = [x_clean[0]]
    segment_y = [y_clean[0]]
    for x_current, y_current, x_previous, y_previous in zip(x_clean[1:], y_clean[1:], x_clean[:-1], y_clean[:-1]):
        if abs(x_current - x_previous) > threshold_x or abs(y_current - y_previous) > threshold_y:
            if len(segment_x) >= 2:
                ax.plot(segment_x, segment_y, **kwargs)
            segment_x = [x_current]
            segment_y = [y_current]
        else:
            segment_x.append(x_current)
            segment_y.append(y_current)
    if len(segment_x) >= 2:
        ax.plot(segment_x, segment_y, **kwargs)


def plot_optional_trajectories(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot final-10-min TPHT-selected trajectories using deterministic subsampling."""
    grouped = group_by_target(rows)
    max_time = max((num(row, "time_s") or math.nan for row in rows), default=math.nan)
    if not math.isfinite(max_time):
        raise ValueError("no finite time_s values")

    def _target_interesting(value_rows: list[dict[str, str]]) -> bool:
        category = next((row.get("category") for row in value_rows if row.get("category") not in (None, "", NA, "unknown")), None)
        max_radius = max(finite([num(row, "max_radius_um") for row in value_rows]), default=math.nan)
        first_ifcoal = any(num(row, "first_ifcoal_time_s") is not None for row in value_rows)
        return bool(category) or (math.isfinite(max_radius) and max_radius >= 15.0) or first_ifcoal

    def _filter_window(window_s: float) -> dict[str, list[dict[str, str]]]:
        start_time = max_time - window_s
        filtered_rows = {
            key: sorted([row for row in value_rows if (num(row, "time_s") or -math.inf) >= start_time], key=lambda row: num(row, "time_s") or 0.0)
            for key, value_rows in grouped.items()
            if _target_interesting(value_rows)
        }
        return {
            key: value_rows
            for key, value_rows in filtered_rows.items()
            if len(value_rows) >= 2 and any((num(row, "radius_um") or -math.inf) >= 5.0 for row in value_rows)
        }

    filtered = _filter_window(600.0)
    window_label = r"last 10 min, r $\geq$ 5 $\mu$m"
    if len(filtered) < 12:
        filtered = _filter_window(1800.0)
        window_label = r"last 30 min, r $\geq$ 5 $\mu$m"
    if not filtered:
        raise ValueError("no TPHT trajectory targets with radius >= 5 micrometer in the final 10 or 30 min")
    sample_keys = stratified_keys(sorted(filtered), 20)
    x_key = "x_unwrapped_m" if "x_unwrapped_m" in column_set(rows) else "x_m"
    y_key = "y_unwrapped_m" if "y_unwrapped_m" in column_set(rows) else "y_m"
    by_time: dict[float, dict[str, list[float]]] = {}
    for value_rows in filtered.values():
        for row in value_rows:
            t_value = num(row, "time_s")
            if t_value is None:
                continue
            entry = by_time.setdefault(t_value / 60.0, {"z": [], "r": []})
            z_value = num(row, "z_m")
            r_value = num(row, "radius_um")
            if z_value is not None:
                entry["z"].append(z_value)
            if r_value is not None:
                entry["r"].append(r_value)
    time_grid = sorted(by_time)
    z25 = [percentile(by_time[t]["z"], 0.25) for t in time_grid]
    z50 = [percentile(by_time[t]["z"], 0.50) for t in time_grid]
    z75 = [percentile(by_time[t]["z"], 0.75) for t in time_grid]
    r25 = [percentile(by_time[t]["r"], 0.25) for t in time_grid]
    r50 = [percentile(by_time[t]["r"], 0.50) for t in time_grid]
    r75 = [percentile(by_time[t]["r"], 0.75) for t in time_grid]
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.72), constrained_layout=True)
    trajectory_cmap = plt.get_cmap("tab20")
    for line_index, key in enumerate(sample_keys):
        value_rows = filtered[key]
        x = [(num(row, x_key) or math.nan) / 1000.0 for row in value_rows]
        y = [(num(row, y_key) or math.nan) / 1000.0 for row in value_rows]
        z = [num(row, "z_m") or math.nan for row in value_rows]
        t = _time_minutes(value_rows)
        r = [num(row, "radius_um") or math.nan for row in value_rows]
        color = trajectory_cmap(line_index % 20)
        plot_masked_line(axes[0, 0], x, y, color=color, alpha=0.38, linewidth=0.45)
        plot_masked_line(axes[0, 1], x, z, color=color, alpha=0.38, linewidth=0.45)
        axes[1, 0].plot(t, z, color="0.75", alpha=0.08, linewidth=0.35)
        axes[1, 1].plot(t, r, color="0.75", alpha=0.08, linewidth=0.35)
    axes[1, 0].plot(time_grid, z50, color=OKABE_ITO["orange"], linewidth=1.1)
    axes[1, 0].fill_between(time_grid, z25, z75, color=OKABE_ITO["orange"], alpha=0.25, linewidth=0)
    axes[1, 1].plot(time_grid, r50, color=OKABE_ITO["blue"], linewidth=1.1)
    axes[1, 1].fill_between(time_grid, r25, r75, color=OKABE_ITO["blue"], alpha=0.25, linewidth=0)
    axes[0, 0].set_xlabel("x (km)")
    axes[0, 0].set_ylabel("y (km)")
    axes[0, 1].set_xlabel("x (km)")
    axes[0, 1].set_ylabel("z (m)")
    axes[1, 0].set_xlabel("Time (min)")
    axes[1, 0].set_ylabel("z (m)")
    axes[1, 1].set_xlabel("Time (min)")
    axes[1, 1].set_ylabel(f"Radius ({MICRON})")
    axes[0, 0].set_title(f"Top view, {window_label}", fontsize=8)
    axes[0, 1].set_title("Side view", fontsize=8)
    axes[1, 0].set_title("Height envelope", fontsize=8)
    axes[1, 1].set_title("Radius envelope", fontsize=8)
    for ax, label in zip(axes.ravel(), ["a", "b", "c", "d"]):
        add_panel_label(ax, label)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_ifcoal_timeline(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot binary if_coal occurrence over time for sampled targets."""
    grouped = group_by_target(rows)
    sorted_keys = sorted(
        grouped,
        key=lambda key: min((num(row, "first_large_time_s", num(row, "time_s") or math.inf) or math.inf for row in grouped[key]), default=math.inf),
    )
    sample_keys = stratified_keys(sorted_keys, 600)
    times = sorted({num(row, "time_s") for row in rows if num(row, "time_s") is not None})
    time_to_index = {time: index for index, time in enumerate(times)}
    matrix = [[0.0 for _ in times] for _ in sample_keys]
    for row_index, key in enumerate(sample_keys):
        for row in grouped[key]:
            t_value = num(row, "time_s")
            if t_value in time_to_index and _ifcoal_flag(row):
                matrix[row_index][time_to_index[t_value]] = 1.0
    plt = configure_matplotlib()
    from matplotlib.colors import ListedColormap

    fig, ax = plt.subplots(figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.42), constrained_layout=True)
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=ListedColormap(["white", OKABE_ITO["vermillion"]]), extent=[min(times) / 60.0, max(times) / 60.0, len(sample_keys), 0])
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Sampled target rank")
    ax.set_title("if_coal proxy", fontsize=8)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_predecessor_tree(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot a few representative predecessor-link examples."""
    plt = configure_matplotlib()
    examples = sorted({row.get("example_id", "") for row in rows if row.get("example_id") not in ("", None)})
    if not examples:
        fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.72), constrained_layout=True)
        no_data_panel(ax, "TPHT predecessor-tree examples")
        save_figure(fig, outdir / SUPP_FIG_DIR, stem)
        plt.close(fig)
        return
    all_radii = [num(row, "radius_um") for row in rows if num(row, "radius_um") is not None]
    radius_min = min(all_radii) if all_radii else 0.0
    radius_max = max(all_radii) if all_radii else 1.0
    if radius_max <= radius_min:
        radius_max = radius_min + 1.0
    cmap = plt.get_cmap("cividis")
    norm = plt.Normalize(radius_min, radius_max)

    batch_size = 3
    for batch_start in range(0, len(examples), batch_size):
        batch_index = batch_start // batch_size + 1
        batch_examples = examples[batch_start : batch_start + batch_size]
        fig, axes = plt.subplots(
            len(batch_examples),
            1,
            figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.30 * len(batch_examples)),
            constrained_layout=True,
            squeeze=False,
            sharex=True,
        )
        scatter = None
        for panel_index, (ax, example) in enumerate(zip(axes.ravel(), batch_examples)):
            subset = sorted([row for row in rows if row.get("example_id") == example], key=lambda row: (num(row, "time_s") or -math.inf, str(row.get("node_id"))))
            trunk_rows = [row for row in subset if row.get("node_role") == "tracked_target_output"]
            trunk_x = [(num(row, "time_s") or math.nan) / 60.0 for row in trunk_rows]
            trunk_y = [num(row, "height_m") or math.nan for row in trunk_rows]
            if len(trunk_rows) >= 2:
                ax.plot(trunk_x, trunk_y, color="0.45", linewidth=0.8, zorder=1)
            important: list[dict[str, str]] = []
            if trunk_rows:
                threshold_rows = [row for row in trunk_rows if (num(row, "radius_um") or -math.inf) >= 15.0]
                if threshold_rows:
                    important.append(threshold_rows[0])
                important.append(max(trunk_rows, key=lambda row: num(row, "radius_um") or -math.inf))
            seen_nodes: set[str] = set()
            important = [row for row in important if not (str(row.get("node_id")) in seen_nodes or seen_nodes.add(str(row.get("node_id"))))]
            if important:
                scatter = ax.scatter(
                    [(num(row, "time_s") or math.nan) / 60.0 for row in important],
                    [num(row, "height_m") or math.nan for row in important],
                    c=[num(row, "radius_um") or math.nan for row in important],
                    cmap=cmap,
                    norm=norm,
                    s=28,
                    marker="o",
                    edgecolor="black",
                    linewidth=0.3,
                    zorder=3,
                    label="important output node",
                )
            ifcoal_rows = [row for row in trunk_rows if _ifcoal_flag(row)]
            if ifcoal_rows:
                ax.scatter(
                    [(num(row, "time_s") or math.nan) / 60.0 for row in ifcoal_rows],
                    [num(row, "height_m") or math.nan for row in ifcoal_rows],
                    c=[num(row, "radius_um") or math.nan for row in ifcoal_rows],
                    cmap=cmap,
                    norm=norm,
                    marker="^",
                    s=24,
                    edgecolor="black",
                    linewidths=0.3,
                    alpha=0.62,
                    zorder=4,
                    label="_nolegend_",
                )
                ax.scatter([], [], marker="^", s=24, facecolor=OKABE_ITO["black"], edgecolor=OKABE_ITO["black"], label="if_coal output node")
            ax.set_ylabel("Height (m)")
            ax.set_title(f"{chr(97 + panel_index)}  target {example}", loc="left", fontsize=8, pad=2)
        axes.ravel()[-1].set_xlabel("Time (min)")
        if scatter is not None:
            handles = []
            labels = []
            for ax in axes.ravel():
                axis_handles, axis_labels = ax.get_legend_handles_labels()
                handles.extend(axis_handles)
                labels.extend(axis_labels)
            unique = dict(zip(labels, handles))
            if unique:
                axes.ravel()[0].legend(unique.values(), unique.keys(), frameon=False, fontsize=7, loc="upper right")
            fig.colorbar(scatter, ax=axes.ravel().tolist(), label=f"Radius ({MICRON})", shrink=0.78)
        stems = [f"{stem}_{batch_index:02d}"]
        if batch_index == 1:
            stems.append(stem)
        for output_stem in stems:
            save_figure(fig, outdir / SUPP_FIG_DIR, output_stem)
        plt.close(fig)


def interval_growth_values(row: dict[str, str]) -> tuple[float, float]:
    """Return diagnostic delta r3 and delta r2 per second values."""
    r0 = num(row, "r0_um")
    r1 = num(row, "r1_um")
    dt = num(row, "dt_s")
    if r0 is None or r1 is None or dt in (None, 0.0):
        return math.nan, math.nan
    return r1**3 - r0**3, (r1**2 - r0**2) / dt


def plot_optional_growth_partition(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot cumulative delta r cubed by category and if_coal interval flag."""
    categories = ["radius_only", "coal_only", "both", "unknown"]
    totals = {category: {"non_flagged": 0.0, "flagged": 0.0} for category in categories}
    for row in rows:
        category = row.get("category", "unknown") or "unknown"
        if category not in totals:
            category = "unknown"
        dr3, _dr2dt = interval_growth_values(row)
        if math.isfinite(dr3):
            key = "flagged" if _ifcoal_flag(row, "if_coal_interval") else "non_flagged"
            totals[category][key] += dr3
    labels = [
        category
        for category in categories
        if abs(totals[category]["non_flagged"]) > 0.0 or abs(totals[category]["flagged"]) > 0.0
    ]
    non_flagged = [totals[label]["non_flagged"] for label in labels]
    flagged = [totals[label]["flagged"] for label in labels]
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.72), constrained_layout=True)
    x = list(range(len(labels)))
    ax.bar(x, non_flagged, color=OKABE_ITO["blue"], edgecolor="black", linewidth=0.4, label="non-flagged")
    ax.bar(x, flagged, bottom=non_flagged, color=OKABE_ITO["vermillion"], edgecolor="black", linewidth=0.4, label="if_coal flagged")
    ax.set_ylabel(r"Cumulative $\Delta r^3$ proxy ($\mu$m$^3$)")
    display_labels = {"radius_only": "radius only", "coal_only": "if_coal only", "both": "both", "unknown": "unclassified"}
    set_clean_categorical_axis(ax, [display_labels.get(label, label) for label in labels], 25)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax.legend(frameon=False, fontsize=7)
    y_min, y_max = ax.get_ylim()
    y_span = y_max - y_min
    for x_pos, total in zip(x, [left + right for left, right in zip(non_flagged, flagged)]):
        if not math.isfinite(total) or total == 0.0:
            continue
        ax.text(
            x_pos,
            total + 0.03 * y_span if total >= 0.0 else total - 0.03 * y_span,
            compact_number_label(total),
            ha="center",
            va="bottom" if total >= 0.0 else "top",
            fontsize=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.4},
            clip_on=False,
        )
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_condensation_proxy(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot non-coalescence-flagged delta r squared per second as a proxy."""
    clean = []
    for row in rows:
        if _ifcoal_flag(row, "if_coal_interval"):
            continue
        _dr3, dr2dt = interval_growth_values(row)
        radius = num(row, "r0_um")
        height = num(row, "height_m")
        if radius is not None and height is not None and math.isfinite(dr2dt):
            clean.append((radius, height, dr2dt))
    if not clean:
        raise ValueError("no finite non-coalescence-flagged interval values")
    plt = configure_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.36), constrained_layout=True)
    hex_radius = axes[0].hexbin([item[0] for item in clean], [item[2] for item in clean], gridsize=45, mincnt=1, cmap="cividis")
    axes[0].set_xlabel(f"Radius ({MICRON})")
    axes[0].set_ylabel(r"$\Delta r^2 / \Delta t$ ($\mu$m$^2$ s$^{-1}$)")
    axes[0].set_title("By radius", fontsize=8)
    hex_height = axes[1].hexbin([item[1] for item in clean], [item[2] for item in clean], gridsize=45, mincnt=1, cmap="cividis")
    axes[1].set_xlabel("Height (m)")
    axes[1].set_ylabel(r"$\Delta r^2 / \Delta t$ ($\mu$m$^2$ s$^{-1}$)")
    axes[1].set_title("By height", fontsize=8)
    fig.colorbar(hex_radius, ax=axes[0], label="Interval count", shrink=0.86)
    fig.colorbar(hex_height, ax=axes[1], label="Interval count", shrink=0.86)
    add_panel_label(axes[0], "a")
    add_panel_label(axes[1], "b")
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_event_jump(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot target-linked coalescence-event jump diagnostics."""
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.68), constrained_layout=True)
    radius_pairs = [(num(row, "pre_radius_um"), num(row, "post_radius_um")) for row in rows]
    radius_pairs = [(left, right) for left, right in radius_pairs if left is not None and right is not None]
    pre = [left for left, _right in radius_pairs]
    post = [right for _left, right in radius_pairs]
    dr3 = [num(row, "delta_r3_um3") for row in rows]
    height = [num(row, "event_height_m") for row in rows]
    time_min = [(num(row, "event_time_s") or math.nan) / 60.0 for row in rows]
    axes[0, 0].scatter(pre, post, s=6, alpha=0.35, color=OKABE_ITO["blue"])
    axes[0, 0].set_xlabel(f"Pre-event radius ({MICRON})")
    axes[0, 0].set_ylabel(f"Post-event radius ({MICRON})")
    axes[0, 1].hist(finite(dr3), bins=40, color=OKABE_ITO["vermillion"], edgecolor="black", linewidth=0.3)
    axes[0, 1].set_xlabel(r"$\Delta r^3$ from event ($\mu$m$^3$)")
    axes[0, 1].set_ylabel("Recorded target-linked events")
    axes[1, 0].hist(finite(height), bins=40, color=OKABE_ITO["orange"], edgecolor="black", linewidth=0.3)
    axes[1, 0].set_xlabel("Event height (m)")
    axes[1, 0].set_ylabel("Recorded target-linked events")
    axes[1, 1].hist(finite(time_min), bins=40, color=OKABE_ITO["bluish_green"], edgecolor="black", linewidth=0.3)
    axes[1, 1].set_xlabel("Event time (min)")
    axes[1, 1].set_ylabel("Recorded target-linked events")
    for ax, label in zip(axes.ravel(), ["a", "b", "c", "d"]):
        add_panel_label(ax, label)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_first_large_joint(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot radius and height at first threshold crossing."""
    pairs = [(num(row, "first_large_radius_um"), num(row, "first_large_height_m")) for row in rows]
    pairs = [(radius, height) for radius, height in pairs if radius is not None and height is not None]
    if not pairs:
        raise ValueError("no finite first-threshold radius-height pairs")
    radius = [pair[0] for pair in pairs]
    height = [pair[1] for pair in pairs]
    from matplotlib.colors import LogNorm

    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.78), constrained_layout=True)
    hexbin = ax.hexbin(radius, height, gridsize=45, mincnt=1, cmap="cividis", norm=LogNorm())
    fig.colorbar(hexbin, ax=ax, label="Target count (log scale)")
    ax.axvline(15.0, color="0.25", linestyle="--", linewidth=0.8)
    ax.set_xlabel(f"Radius at first threshold crossing ({MICRON})")
    ax.set_ylabel("Height at first threshold crossing (m)")
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_max_radius_ccdf(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot maximum-radius exceedance probability."""
    values = sorted(finite([num(row, "max_radius_um") for row in rows]))
    if not values:
        raise ValueError("no finite maximum-radius values")
    exceedance = [(len(values) - index) / len(values) for index, _value in enumerate(values)] if values else []
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.72), constrained_layout=True)
    ax.plot(values, exceedance, color=OKABE_ITO["blue"], linewidth=1.0)
    for reference in (15.0, 20.0, 25.0):
        ax.axvline(reference, color="0.45", linestyle="--", linewidth=0.7)
    ax.set_yscale("log")
    ax.set_xlabel(f"Maximum reconstructed radius ({MICRON})")
    ax.set_ylabel("Exceedance probability (-)")
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_threshold_duration(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot threshold-exceedance durations."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.72), constrained_layout=True)
    series = [
        ("duration_ge_15_s", rf"r >= 15 {MICRON}", OKABE_ITO["blue"]),
        ("duration_ge_20_s", rf"r >= 20 {MICRON}", OKABE_ITO["vermillion"]),
    ]
    duration_values = [(label, color, [value / 60.0 for value in finite([num(row, key) for row in rows])]) for key, label, color in series]
    all_values = [value for _label, _color, values in duration_values for value in values]
    if not all_values:
        raise ValueError("no finite threshold-duration values")
    max_duration = max(all_values)
    bin_width = max(0.05, max_duration / 28.0)
    edge_count = max(2, int(math.ceil(max_duration / bin_width)) + 1)
    edges = [index * bin_width for index in range(edge_count + 1)]
    bar_width = bin_width * 0.38
    offsets = [0.05 * bin_width, 0.52 * bin_width]
    for offset, (label, color, values) in zip(offsets, duration_values):
        counts = [0] * (len(edges) - 1)
        for value in values:
            bin_index = min(len(counts) - 1, int(value / bin_width))
            counts[bin_index] += 1
        xs = [edges[index] + offset for index in range(len(counts))]
        ax.bar(xs, counts, width=bar_width, align="edge", color=color, alpha=0.68, edgecolor="black", linewidth=0.25, label=label)
    ax.set_xlabel("Threshold-exceedance duration (min)")
    ax.set_ylabel("Target count (log scale)")
    ax.set_xlim(left=0.0)
    ax.set_yscale("log")
    ax.legend(frameon=False)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_occurrence_zt(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot binned target occurrence in height-time space."""
    time_values = sorted({num(row, "time_s") for row in rows if num(row, "time_s") is not None})
    height_values = sorted({num(row, "height_bin_m") for row in rows if num(row, "height_bin_m") is not None})
    if not time_values or not height_values:
        raise ValueError("no finite time-height bins")
    time_index = {value: index for index, value in enumerate(time_values)}
    height_index = {value: index for index, value in enumerate(height_values)}
    fields = [
        ("target_count", "all targets"),
        ("large_target_count", "ever threshold-crossing targets"),
        ("ifcoal_target_count", "ever if_coal-flagged targets"),
    ]
    matrices = []
    positive_values: list[float] = []
    for field, _label in fields:
        matrix = [[math.nan for _ in time_values] for _ in height_values]
        for row in rows:
            t_value = num(row, "time_s")
            h_value = num(row, "height_bin_m")
            if t_value in time_index and h_value in height_index:
                value = num(row, field)
                if value is not None and value > 0.0:
                    matrix[height_index[h_value]][time_index[t_value]] = value
                    positive_values.append(value)
        matrices.append(matrix)
    if not positive_values:
        raise ValueError("no positive time-height target counts")

    from matplotlib.colors import LogNorm

    plt = configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.43), constrained_layout=True)
    norm = LogNorm(vmin=1.0, vmax=max(positive_values))
    image = None
    for index, (ax, matrix, (_field, label), panel) in enumerate(zip(axes, matrices, fields, ["a", "b", "c"])):
        image = ax.imshow(matrix, aspect="auto", origin="lower", interpolation="nearest", extent=[min(time_values) / 60.0, max(time_values) / 60.0, min(height_values), max(height_values)], cmap="cividis", norm=norm)
        ax.set_xlabel("Time (min)")
        if index == 0:
            ax.set_ylabel("Height (m)")
        else:
            ax.set_ylabel("")
        ax.set_title(label, fontsize=7, pad=2)
        add_panel_label(ax, panel)
    if image is not None:
        fig.colorbar(image, ax=list(axes), label="Target count (log scale)", shrink=0.86, pad=0.02)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_category_max_radius(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot category-specific maximum-radius distributions."""
    labels = ["radius_only", "coal_only", "both", "unknown"]
    grouped = [[value for value in finite([num(row, "max_radius_um") for row in rows if row.get("category") == label])] for label in labels]
    display_labels = {"radius_only": "radius only", "coal_only": "if_coal only", "both": "both", "unknown": "unclassified"}
    labels = [label for label, values in zip(labels, grouped) if values]
    grouped = [values for values in grouped if values]
    if not grouped:
        raise ValueError("no finite category-specific maximum-radius values")
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.72), constrained_layout=True)
    positions = list(range(len(labels)))
    ax.boxplot(grouped, positions=positions, widths=0.55, showfliers=False, patch_artist=True, boxprops={"facecolor": OKABE_ITO["sky_blue"], "alpha": 0.45})
    ax.set_xlabel("Target category")
    ax.set_ylabel(f"Maximum radius ({MICRON})")
    set_clean_categorical_axis(ax, [display_labels.get(label, label) for label in labels], 25)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_first_coal_minus_large(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot first if_coal minus first-large timing proxy."""
    deltas = []
    for row in rows:
        first_large = num(row, "first_large_time_s")
        first_ifcoal = num(row, "first_ifcoal_time_s")
        if first_large is not None and first_ifcoal is not None:
            deltas.append((first_ifcoal - first_large) / 60.0)
    if not deltas:
        raise ValueError("no finite first-if_coal minus first-threshold values")
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.70), constrained_layout=True)
    ax.hist(deltas, bins=50, color=OKABE_ITO["vermillion"], edgecolor="black", linewidth=0.3)
    ax.axvline(0.0, color="0.25", linestyle="--", linewidth=0.8)
    ax.set_xlabel("First if_coal flag - first threshold crossing (min)")
    ax.set_ylabel("Target count")
    ax.set_title("if_coal proxy", fontsize=8)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_relative_cloud_top(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot first-threshold height relative to cloud top."""
    values = finite([num(row, "z_first_large_minus_cloud_top_m") for row in rows])
    if not values:
        raise ValueError("no finite cloud-top-relative values")
    approximate = any(str(row.get("approximate_cloud_top", "")).lower() in ("1", "true", "yes") for row in rows)
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=(SINGLE_WIDTH, SINGLE_WIDTH * 0.70), constrained_layout=True)
    ax.hist(values, bins=45, color=OKABE_ITO["blue"], edgecolor="black", linewidth=0.3)
    ax.axvline(0.0, color="0.25", linestyle="--", linewidth=0.8)
    ax.set_xlabel("z(first threshold crossing) - z(cloud top) (m)")
    ax.set_ylabel("Target count")
    if approximate:
        ax.set_title("Approx. cloud top", fontsize=8)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


def plot_optional_growth_phase_space(rows: list[dict[str, str]], outdir: Path, stem: str) -> None:
    """Plot growth phase space for coalescence-flagged and non-flagged intervals."""
    data = {"non-flagged": [], "if_coal flagged": []}
    for row in rows:
        _dr3, dr2dt = interval_growth_values(row)
        radius = num(row, "r0_um")
        if radius is not None and math.isfinite(dr2dt):
            key = "if_coal flagged" if _ifcoal_flag(row, "if_coal_interval") else "non-flagged"
            data[key].append((radius, dr2dt))
    if not data["non-flagged"] and not data["if_coal flagged"]:
        raise ValueError("no finite growth phase-space values")
    plt = configure_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(MAIN_WIDTH, MAIN_WIDTH * 0.36), constrained_layout=True)
    for ax, key, panel in zip(axes, data, ["a", "b"]):
        values = data[key]
        if not values:
            no_data_panel(ax, key)
            add_panel_label(ax, panel)
            continue
        hexbin = ax.hexbin([item[0] for item in values], [item[1] for item in values], gridsize=45, mincnt=1, cmap="cividis")
        ax.set_xlabel(f"Radius ({MICRON})")
        ax.set_ylabel(r"$\Delta r^2 / \Delta t$ ($\mu$m$^2$ s$^{-1}$)")
        ax.set_title(key, fontsize=8)
        fig.colorbar(hexbin, ax=ax, label="Interval count", shrink=0.86)
        add_panel_label(ax, panel)
    save_figure(fig, outdir / SUPP_FIG_DIR, stem)
    plt.close(fig)


OPTIONAL_PLOTTERS = {
    "trajectories": plot_optional_trajectories,
    "ifcoal_timeline": plot_optional_ifcoal_timeline,
    "predecessor_tree": plot_optional_predecessor_tree,
    "growth_partition": plot_optional_growth_partition,
    "condensation_proxy": plot_optional_condensation_proxy,
    "event_jump": plot_optional_event_jump,
    "first_large_joint": plot_optional_first_large_joint,
    "max_radius_ccdf": plot_optional_max_radius_ccdf,
    "threshold_duration": plot_optional_threshold_duration,
    "occurrence_zt": plot_optional_occurrence_zt,
    "category_max_radius": plot_optional_category_max_radius,
    "first_coal_minus_large": plot_optional_first_coal_minus_large,
    "relative_cloud_top": plot_optional_relative_cloud_top,
    "growth_phase_space": plot_optional_growth_phase_space,
}


def plot_optional_tpht_science(input_dir: Path, outdir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Generate optional TPHT science-diagnostic candidates when source variables exist."""
    manifest: list[dict[str, Any]] = []
    availability: list[dict[str, Any]] = []
    for spec in OPTIONAL_TPHT_SPECS:
        source_table, rows = read_first_available_table(input_dir, spec["source_tables"])
        if not source_table:
            availability.append(optional_availability_row(spec, "skipped_missing_source", None, "No required per-target, per-interval, or event-level source table was found."))
            continue
        if not rows_have_columns(rows, spec["required_columns"]):
            missing = sorted(set(spec["required_columns"]) - column_set(rows))
            availability.append(optional_availability_row(spec, "skipped_missing_columns", source_table, f"Missing columns: {', '.join(missing)}."))
            continue
        try:
            OPTIONAL_PLOTTERS[spec["plotter"]](rows, outdir, spec["stem"])
        except Exception as exc:  # pragma: no cover - reported to diagnostics table for HPC runs.
            availability.append(optional_availability_row(spec, "skipped_plot_error", source_table, f"{type(exc).__name__}: {exc}"))
            continue
        availability.append(optional_availability_row(spec, "generated", source_table, "Generated from available variable-level diagnostics."))
        manifest.append(
            manifest_row(
                spec["stem"],
                "supplement_candidate",
                f"{source_table}.csv",
                "Time (min); Height (m); Radius (micrometer); Target count (-); Probability (-)",
                170,
                CLAIMS[spec["stem"]],
                "TPHT diagnostic capability; threshold-crossing targets; if_coal occurrence proxy where applicable; no causal mechanism implied.",
                "a,b,c,d",
            )
        )
    write_table_bundle(availability, table_path(outdir, DIAG_TABLE_DIR, "tpht_optional_figure_availability"))
    return manifest, availability


def plot_all(input_dir: Path, outdir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Generate all candidate figures and return manifest rows."""
    manifest = [
        plot_framework(outdir),
        plot_benchmark(input_dir, outdir),
        plot_tpht_handoff_storage(input_dir, outdir),
        plot_tpht_target_diagnostics(input_dir, outdir),
        plot_scalability_io(input_dir, outdir),
        plot_restart_supp(input_dir, outdir),
        plot_full_benchmark_supp(input_dir, outdir),
        plot_sampling_supp(input_dir, outdir),
        plot_rank_load_supp(input_dir, outdir),
        plot_chain_validity_supp(input_dir, outdir),
        plot_target_histories_supp(input_dir, outdir),
        plot_sdnc_details_supp(input_dir, outdir),
        plot_io_details_supp(input_dir, outdir),
    ]
    optional_manifest, optional_availability = plot_optional_tpht_science(input_dir, outdir)
    manifest.extend(optional_manifest)
    return manifest, optional_availability


def write_readme(outdir: Path, manifest: list[dict[str, Any]], input_dir: Path, optional_availability: list[dict[str, Any]] | None = None) -> None:
    """Write the candidate-figure README with caveats and source mapping."""
    main_rows = [row for row in manifest if row["category"] == "main_candidate"]
    supp_rows = [row for row in manifest if row["category"] == "supplement_candidate"]
    source_lines = [f"- `{row['figure_stem']}`: {row['source_tables']}" for row in manifest]
    optional_rows = optional_availability or []
    optional_lines = [
        f"- `{row['figure_stem']}`: {row['status']} ({row['reason']})"
        for row in optional_rows
    ]
    optional_text = "\n".join(optional_lines) if optional_lines else "- No optional TPHT science-diagnostic candidates were evaluated."
    readme = f"""# GMD2026 Candidate Figures and Tables

This directory contains a regenerated GMD-ready candidate figure set for the revised GMD2026 manuscript. No final figure numbering has been assigned.

Source analysis directory:

`{input_dir}`

## Main-Candidate Figure Groups

- framework and experiment design: `main_candidate_framework_experiment_design`
- controlled computational overhead: `main_candidate_controlled_benchmark_overhead`
- TPHT handoff and storage reduction: `main_candidate_TPHT_handoff_storage`
- TPHT target diagnostic characterization: `main_candidate_TPHT_target_diagnostics`
- scalability and I/O sensitivity: `main_candidate_scalability_io_sensitivity`

## Supplement-Candidate Figure Groups

- restart sanity: `supp_candidate_restart_sanity`
- full benchmark diagnostics: `supp_candidate_full_benchmark_diagnostics`
- 2D sampling verification: `supp_candidate_2D_sampling_verification`
- TPHT rank-load/dedup: `supp_candidate_TPHT_rank_load_dedup`
- TPHT chain-validity proxy: `supp_candidate_TPHT_chain_validity_proxy`
- TPHT target histories: `supp_candidate_TPHT_target_histories`
- SDNC scaling details: `supp_candidate_SDNC_scaling_details`
- output interval I/O details: `supp_candidate_output_interval_io_details`

Optional threshold-crossing diagnostic figures are generated only when the source tables contain full per-target first-threshold radius, maximum-radius, or category-specific maximum-radius distributions. The available top-example table is not used as a full-population distribution.

## Additional TPHT Science-Diagnostic Candidates

The optional TPHT science-diagnostic candidates are generated only when variable-level source tables are available. They are skipped rather than approximated from aggregate tables. This avoids mixing all-time records into first-selected conditions or presenting a top-example table as a full-population diagnostic.

{optional_text}

Detailed availability is written to `tables/diagnostics/tpht_optional_figure_availability.csv`.

## Claims

""" + "\n".join(f"- `{row['figure_stem']}`: {row['claim']}" for row in main_rows + supp_rows) + """

## Caveats

- `03_fw_rep_2d_600s` is a 2D sampling-procedure verification, not 3D representativeness.
- `01_bench_3d_samp_30min` is a controlled cold-start computational benchmark.
- `02_tpht_3d_interest_70min` demonstrates TPHT diagnostic capability and does not assign causal pathways.
- `if_coal` is an if_coal binary flag and an if_coal occurrence proxy, not an exact event-count diagnostic.
- Storage reduction relative to full-population BW uses estimated full-BW output and should be described as estimated storage reduction.
- Discovery-time diagnostics are not main-text evidence because TPHT_ID_DIAG is absent or empty.
- Chain-validity direct metrics are still NA; only proxy metrics and reconstructed_time_coverage are available.
- BW selected-output model time should be interpreted cautiously; do not overclaim forward physical chronology without checking the run convention.
- Optional TPHT science-diagnostic figures derived from BW selected-output chronology are provisional until the BW output time convention and target-linked coalescence-log linkage are verified.
- `unknown` / `unclassified` in target-category diagnostics means a deduplicated TPHT target whose first selected condition cannot be recovered from the BW selected-output diagnostics; it is not a third physical selection criterion.
- Target-history curves are target-set medians over all reconstructed target records at each output time, not only records currently satisfying the interest criteria; conditional fractions are shown separately.
- The `unknown` row in the first-large height source table means no first-large height is available for those records, so it is excluded from the height-distribution panel rather than treated as a height bin.

## Deleted Or Downgraded Figures

- TPHT discovery-time figure: downgraded/excluded from main-candidate set because TPHT_ID_DIAG is absent or empty.
- Chain-validity figure: supplement-candidate only because direct chain-validity metrics are NA.
- TPHT target histories: supplement-candidate only as diagnostic capability demonstration.
- 03 sampling figures: supplement-candidate only as 2D sampling-procedure verification.
- Coal-related plots: renamed as if_coal occurrence proxy where used.
- Any figure implying final figure numbers, 3D representativeness, measured full-BW reduction, or causal mechanism is excluded.

## Source Tables Used For Each Figure

""" + "\n".join(source_lines) + """

## Keywords

GMD2026; SCALE-SDM; TPHT; Two-Pass Hybrid Tracking; interest-restricted backward reconstruction; target-set handoff; controlled cold-start computational benchmark; sampled forward tracking; sampled backward tracking; no-tracking baseline; coalescence-log-only baseline; if_coal binary flag; if_coal occurrence proxy; coalescence-flagged interval; threshold-crossing targets; near-threshold droplets; diagnostic characterization; diagnostic capability; no causal pathway implied; 2D sampling-procedure verification; not 3D representativeness; estimated full-BW output; estimated storage reduction; target reduction ratio; deduplicated target pairs; missing_in_bw; extra_in_bw; unknown target category; SDNC scaling; output interval I/O sensitivity; selected SD output; coalescence log output; rank peak memory; Okabe-Ito palette; publication-quality Matplotlib; shared axis labels; no overlapping tick labels; bbox_inches tight; constrained layout; panel labels; main candidate figure set; supplement candidate figure set; analysis_outputs_for_GMD.
"""
    (outdir / "README_GMD_FIGURES.md").write_text(readme)


def main() -> None:
    """CLI entry point."""
    args = parse_args()
    root = args.root.resolve()
    input_dir = (args.input or root / "analysis_outputs").resolve()
    outdir = (args.outdir or root / "analysis_outputs_for_GMD").resolve()
    planned = [
        outdir / "README_GMD_FIGURES.md",
        outdir / MAIN_FIG_DIR / "main_candidate_framework_experiment_design.pdf",
        outdir / MAIN_FIG_DIR / "main_candidate_controlled_benchmark_overhead.pdf",
        outdir / MAIN_FIG_DIR / "main_candidate_TPHT_handoff_storage.pdf",
        outdir / MAIN_FIG_DIR / "main_candidate_TPHT_target_diagnostics.pdf",
        outdir / MAIN_FIG_DIR / "main_candidate_scalability_io_sensitivity.pdf",
    ]
    if args.dry_run:
        print("[dry-run] regenerate GMD candidate figures")
        print(f"[dry-run] root={root}")
        print(f"[dry-run] input={input_dir}")
        print(f"[dry-run] outdir={outdir}")
        for path in planned:
            print(f"[dry-run] would write {path}")
        return
    ensure_dirs(outdir)
    copy_source_tables(
        input_dir,
        outdir,
        [
            "01_benchmark_summary",
            "02_tpht_summary",
            "02_tpht_consistency",
            "02_tpht_target_categories",
            "02_tpht_science_summary",
            "02_tpht_science_formation_height_bins",
            "02_tpht_science_target_summary",
            "02_tpht_target_first_large",
            "04_sdnc_scaling_summary",
            "05_outint_io_summary",
        ],
        MAIN_SRC_DIR,
    )
    copy_source_tables(
        input_dir,
        outdir,
        [
            "00_restart_summary",
            "01_benchmark_summary",
            "02_tpht_rank_load_balance",
            "02_tpht_chain_validity",
            "02_tpht_target_histories",
            "02_tpht_science_time_series",
            "02_tpht_science_target_records",
            "02_tpht_target_trajectory_records",
            "02_tpht_target_ifcoal_timeline",
            "02_tpht_science_interval_diagnostics",
            "02_tpht_target_interval_diagnostics",
            "02_tpht_predecessor_tree_examples",
            "02_tpht_target_event_links",
            "02_tpht_event_level_coalescence_links",
            "02_tpht_science_target_summary",
            "02_tpht_target_first_large",
            "02_tpht_max_radius_by_target",
            "02_tpht_threshold_durations",
            "02_tpht_target_occurrence_zt",
            "02_tpht_category_max_radius",
            "02_tpht_first_coal_minus_first_large",
            "02_tpht_cloud_top_relative",
            "03_sampling_metrics",
            "03_sampling_seed_statistics",
            "04_sdnc_scaling_summary",
            "04_scaling_slopes",
            "05_outint_io_summary",
        ],
        SUPP_SRC_DIR,
    )
    make_candidate_tables(input_dir, outdir)
    manifest, optional_availability = plot_all(input_dir, outdir)
    write_manifest(outdir, manifest)
    write_readme(outdir, manifest, input_dir, optional_availability)


if __name__ == "__main__":
    main()
