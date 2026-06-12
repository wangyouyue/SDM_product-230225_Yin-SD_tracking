#!/usr/bin/env python3
"""Plot TPHT raw-output sanity checks for mode-setting validation.

This GMD2026 diagnostic script is intentionally separate from the manuscript
candidate figures.  It reads already generated TPHT diagnostic tables and makes
foreground-safe plots for checking whether the SCALE-SDM TPHT selected output,
predecessor links, and reconstructed histories look internally consistent.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir, warn_or_raise  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, label_bars, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float, write_table_bundle  # noqa: E402


FIGURE_OUTPUTS = [
    "figures/02_tpht_mode_sanity_overview.{pdf,svg,png}",
    "figures/02_tpht_mode_sanity_time_order.{pdf,svg,png}",
    "figures/02_tpht_mode_sanity_sample_chains.{pdf,svg,png}",
    "figures/02_tpht_mode_sanity_chain_summary.{pdf,svg,png}",
]
TABLE_OUTPUTS = [
    "tables/02_tpht_mode_sanity_findings.{csv,md,tex,json}",
]

FINDING_COLUMNS = [
    "check_name",
    "source_table",
    "value",
    "unit",
    "interpretation",
    "warnings",
]


def _number(row: dict[str, Any], key: str) -> float | None:
    """Return a finite numeric value from a table row."""
    return safe_float(row.get(key))


def _time_min(row: dict[str, Any]) -> float | None:
    """Return model/output time in minutes from a row."""
    value = _number(row, "time_min")
    if value is not None:
        return value
    seconds = _number(row, "time_s")
    return seconds / 60.0 if seconds is not None else None


def _sorted_points(rows: list[dict[str, Any]], value_key: str) -> tuple[list[float], list[float]]:
    """Return finite x-y points sorted by model/output time."""
    points = []
    for row in rows:
        time_min = _time_min(row)
        value = _number(row, value_key)
        if time_min is not None and value is not None:
            points.append((time_min, value))
    points.sort(key=lambda item: item[0])
    return [point[0] for point in points], [point[1] for point in points]


def _format_value(value: float | None, precision: int = 4) -> str:
    """Format a scalar finding while preserving missing values."""
    if value is None or not math.isfinite(value):
        return "NA"
    if abs(value) >= 1.0e5 or (0.0 < abs(value) < 1.0e-3):
        return f"{value:.{precision}e}"
    return f"{value:.{precision}g}"


def _warn_missing(path: Path, warnings: list[str], strict: bool) -> None:
    """Record that a source table required for a sanity plot is missing."""
    warn_or_raise(warnings, f"Missing source table: {path}", strict=strict)


def _read_source_table(outdir: Path, table_name: str, warnings: list[str], strict: bool) -> list[dict[str, str]]:
    """Read one source CSV table and warn if it is unavailable."""
    path = outdir / "tables" / table_name
    rows = read_table(path)
    if not rows:
        _warn_missing(path, warnings, strict)
    return rows


def _plot_overview(outdir: Path, rows: list[dict[str, str]]) -> None:
    """Plot count and bulk property histories from TPHT selected-output tables."""
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=figure_size("double", 0.70), sharex=True, constrained_layout=True)
    panels = [
        ("record_count", "Selected-output records (count)", "Selected records", PALETTE[0]),
        ("mean_radius_um", r"Mean radius ($\mu$m)", "Mean radius", PALETTE[1]),
        ("large_record_fraction", "Record fraction (%)", r"r $\geq$ 15 $\mu$m", PALETTE[2]),
        ("coal_record_fraction", "Record fraction (%)", "if_coal occurrence proxy", PALETTE[3]),
    ]
    for ax, (key, ylabel, title, color) in zip(axes.flat, panels):
        xs, ys = _sorted_points(rows, key)
        if key.endswith("_fraction"):
            ys = [value * 100.0 for value in ys]
        if not xs or not has_any_number(ys):
            no_data_panel(ax, title)
            continue
        ax.plot(xs, ys, color=color, marker="o", markersize=1.8, linewidth=1.0)
        if key == "mean_radius_um":
            ax.axhline(15.0, color="0.35", linestyle="--", linewidth=0.8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    for ax in axes[1]:
        ax.set_xlabel("Model/output time (min)")
    save_figure(fig, outdir / "figures", "02_tpht_mode_sanity_overview")
    plt.close(fig)


def _direction_rows(rows: list[dict[str, str]], direction: str) -> list[dict[str, str]]:
    """Return rows for one time-order assumption."""
    return [row for row in rows if str(row.get("direction", "")) == direction]


def _plot_time_order(outdir: Path, rows: list[dict[str, str]]) -> None:
    """Plot raw predecessor-link consistency for ascending and descending time."""
    plt = configure_matplotlib()
    fig, axes = plt.subplots(3, 1, figsize=figure_size("double", 0.82), sharex=True, constrained_layout=True)
    direction_styles = [
        ("ascending", PALETTE[0], "ascending model time"),
        ("descending", PALETTE[1], "descending model time"),
    ]
    for direction, color, label in direction_styles:
        direction_table = _direction_rows(rows, direction)
        xs, fractions = _sorted_points(direction_table, "valid_link_fraction")
        _, unmatched = _sorted_points(direction_table, "unmatched_records")
        _, counts = _sorted_points(direction_table, "record_count")
        if xs and has_any_number(fractions):
            axes[0].plot(xs, [value * 100.0 for value in fractions], color=color, label=label, linewidth=1.0)
        if xs and has_any_number(unmatched):
            axes[1].plot(xs, unmatched, color=color, label=label, linewidth=1.0)
        if xs and has_any_number(counts):
            axes[2].plot(xs, counts, color=color, label=label, linewidth=1.0)
    if rows:
        axes[0].set_ylabel("Valid links (%)")
        axes[1].set_ylabel("Unmatched records (count)")
        axes[2].set_ylabel("Selected records (count)")
        axes[2].set_xlabel("Model/output time (min)")
        for ax in axes:
            ax.legend(loc="best", frameon=False)
    else:
        for ax in axes:
            no_data_panel(ax, "Raw time-order check")
    save_figure(fig, outdir / "figures", "02_tpht_mode_sanity_time_order")
    plt.close(fig)


def _chain_groups(rows: list[dict[str, str]], max_chains: int) -> list[tuple[str, list[dict[str, str]]]]:
    """Group raw chain-sample rows by target ID."""
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        target_id = str(row.get("target_id", ""))
        if target_id in {"", "NA", "None"}:
            continue
        grouped.setdefault(target_id, []).append(row)
    output = []
    for target_id, group_rows in sorted(grouped.items())[:max(1, max_chains)]:
        group_rows.sort(key=lambda row: _time_min(row) if _time_min(row) is not None else math.inf)
        output.append((target_id, group_rows))
    return output


def _plot_sample_chains(outdir: Path, rows: list[dict[str, str]], max_chains: int) -> None:
    """Plot raw radius and height histories for a small set of traced chains."""
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 1, figsize=figure_size("double", 0.70), sharex=True, constrained_layout=True)
    groups = _chain_groups(rows, max_chains)
    if not groups:
        for ax in axes:
            no_data_panel(ax, "Raw chain samples")
    else:
        cmap = plt.get_cmap("tab20")
        for index, (target_id, group_rows) in enumerate(groups):
            xs, radius = _sorted_points(group_rows, "radius_um")
            _, height = _sorted_points(group_rows, "height_m")
            if not xs:
                continue
            label = target_id.replace("chain_", "c")
            color = cmap(index % 20)
            if has_any_number(radius):
                axes[0].plot(xs, radius, color=color, linewidth=0.9, label=label)
            if has_any_number(height):
                axes[1].plot(xs, height, color=color, linewidth=0.9, label=label)
        axes[0].axhline(5.0, color="0.5", linewidth=0.7, linestyle=":")
        axes[0].axhline(15.0, color="0.25", linewidth=0.8, linestyle="--")
        axes[0].set_ylabel(r"Radius ($\mu$m)")
        axes[1].set_ylabel("Height (m)")
        axes[1].set_xlabel("Model/output time (min)")
        axes[0].legend(loc="upper right", bbox_to_anchor=(1.18, 1.02), frameon=False, ncol=1)
    save_figure(fig, outdir / "figures", "02_tpht_mode_sanity_sample_chains")
    plt.close(fig)


def _plot_chain_summary(outdir: Path, rows: list[dict[str, str]], max_chains: int) -> None:
    """Plot compact summaries of raw traced chain samples."""
    plot_rows = rows[: max(1, max_chains)]
    plt = configure_matplotlib()
    fig, axes = plt.subplots(3, 1, figsize=figure_size("double", 0.90), sharex=True, constrained_layout=True)
    if not plot_rows:
        for ax in axes:
            no_data_panel(ax, "Chain summary")
    else:
        labels = [str(row.get("target_id", "")).replace("chain_", "c") for row in plot_rows]
        x_positions = list(range(len(plot_rows)))
        series = [
            ("peak_time_s", "Peak time (min)", PALETTE[0], 1.0 / 60.0),
            ("peak_radius_um", r"Peak radius ($\mu$m)", PALETTE[1], 1.0),
            ("last_radius_um", r"Last radius ($\mu$m)", PALETTE[2], 1.0),
        ]
        for ax, (key, ylabel, color, scale) in zip(axes, series):
            values = []
            for row in plot_rows:
                value = _number(row, key)
                values.append(value * scale if value is not None else 0.0)
            bars = ax.bar(x_positions, values, color=color, edgecolor="black", linewidth=0.5)
            positive = [value for value in values if value > 0.0]
            ax.set_ylim(0.0, max(positive) * 1.18 if positive else 1.0)
            ax.set_ylabel(ylabel)
            label_bars(ax, bars)
        axes[-1].set_xticks(x_positions, labels, rotation=40, ha="right")
        axes[-1].set_xlabel("Sampled chain ID")
    save_figure(fig, outdir / "figures", "02_tpht_mode_sanity_chain_summary")
    plt.close(fig)


def _append_finding(
    rows: list[dict[str, Any]],
    check_name: str,
    source_table: str,
    value: Any,
    unit: str,
    interpretation: str,
    warnings: str | None = None,
) -> None:
    """Append one compact finding row."""
    rows.append(
        {
            "check_name": check_name,
            "source_table": source_table,
            "value": value if value not in (None, "") else "NA",
            "unit": unit,
            "interpretation": interpretation,
            "warnings": warnings,
        }
    )


def _findings_from_tables(
    science_rows: list[dict[str, str]],
    order_summary_rows: list[dict[str, str]],
    chain_summary_rows: list[dict[str, str]],
    warnings: list[str],
) -> list[dict[str, Any]]:
    """Build a compact sanity findings table from available diagnostics."""
    rows: list[dict[str, Any]] = []
    if science_rows:
        times = [_time_min(row) for row in science_rows]
        times = [value for value in times if value is not None]
        record_counts = [_number(row, "record_count") for row in science_rows]
        record_counts = [value for value in record_counts if value is not None]
        _append_finding(rows, "science_time_group_count", "02_tpht_science_time_series.csv", len(times), "count", "Number of output times represented in the science time-series table.")
        _append_finding(rows, "selected_record_count_min", "02_tpht_science_time_series.csv", _format_value(min(record_counts) if record_counts else None), "records", "Minimum selected-output records per time.")
        _append_finding(rows, "selected_record_count_max", "02_tpht_science_time_series.csv", _format_value(max(record_counts) if record_counts else None), "records", "Maximum selected-output records per time.")
    if order_summary_rows:
        row = order_summary_rows[0]
        _append_finding(rows, "direction_inference", "02_tpht_raw_time_order_check.csv", row.get("direction_inference", "NA"), "-", "Preferred raw selected-output time ordering from predecessor-link consistency.")
        _append_finding(rows, "ascending_mean_valid_link_fraction", "02_tpht_raw_time_order_check.csv", row.get("ascending_mean_valid_link_fraction", "NA"), "fraction", "Mean adjacent-output predecessor-link consistency in ascending model time.")
        _append_finding(rows, "descending_mean_valid_link_fraction", "02_tpht_raw_time_order_check.csv", row.get("descending_mean_valid_link_fraction", "NA"), "fraction", "Mean adjacent-output predecessor-link consistency in descending model time.")
        _append_finding(rows, "ascending_unmatched_total", "02_tpht_raw_time_order_check.csv", row.get("ascending_unmatched_total", "NA"), "records", "Unmatched selected-output records in the ascending-time link test.")
    if chain_summary_rows:
        peak_times = [_number(row, "peak_time_s") for row in chain_summary_rows]
        peak_times = [value for value in peak_times if value is not None]
        first_ge_15 = [_number(row, "first_ge_15um_time_s") for row in chain_summary_rows]
        first_ge_15 = [value for value in first_ge_15 if value is not None]
        last_radius = [_number(row, "last_radius_um") for row in chain_summary_rows]
        last_radius = [value for value in last_radius if value is not None]
        last_ge_5 = sum(1 for row in chain_summary_rows if str(row.get("last_radius_ge_5um", "")).lower() == "true")
        last_ge_15 = sum(1 for row in chain_summary_rows if str(row.get("last_radius_ge_15um", "")).lower() == "true")
        _append_finding(rows, "raw_sample_chain_count", "02_tpht_raw_chain_sample_summary.csv", len(chain_summary_rows), "chains", "Number of raw chains summarized.")
        _append_finding(rows, "raw_sample_peak_time_range", "02_tpht_raw_chain_sample_summary.csv", f"{_format_value(min(peak_times) / 60.0 if peak_times else None)}-{_format_value(max(peak_times) / 60.0 if peak_times else None)}", "min", "Range of peak-radius times in sampled raw chains.")
        _append_finding(rows, "raw_sample_first_ge_15um_time_range", "02_tpht_raw_chain_sample_summary.csv", f"{_format_value(min(first_ge_15) / 60.0 if first_ge_15 else None)}-{_format_value(max(first_ge_15) / 60.0 if first_ge_15 else None)}", "min", "Range of first r >= 15 micrometer times in sampled raw chains.")
        _append_finding(rows, "raw_sample_last_radius_max", "02_tpht_raw_chain_sample_summary.csv", _format_value(max(last_radius) if last_radius else None), "micrometer", "Largest final-time radius among sampled raw chains.")
        _append_finding(rows, "raw_sample_last_ge_5um_count", "02_tpht_raw_chain_sample_summary.csv", last_ge_5, "chains", "Sampled chains with final-time radius >= 5 micrometer.")
        _append_finding(rows, "raw_sample_last_ge_15um_count", "02_tpht_raw_chain_sample_summary.csv", last_ge_15, "chains", "Sampled chains with final-time radius >= 15 micrometer.")
    if warnings:
        _append_finding(rows, "source_warnings", "multiple", len(warnings), "warnings", "Input tables missing or incomplete; inspect warnings column.", "; ".join(warnings))
    return rows


def plot_and_write(outdir: Path, strict: bool, max_chains: int) -> None:
    """Create TPHT mode sanity figures and a compact findings table."""
    warnings: list[str] = []
    science_rows = _read_source_table(outdir, "02_tpht_science_time_series.csv", warnings, strict)
    order_rows = _read_source_table(outdir, "02_tpht_raw_time_order_by_level.csv", warnings, strict)
    order_summary_rows = _read_source_table(outdir, "02_tpht_raw_time_order_check.csv", warnings, strict)
    chain_rows = _read_source_table(outdir, "02_tpht_raw_chain_samples.csv", warnings, strict)
    chain_summary_rows = _read_source_table(outdir, "02_tpht_raw_chain_sample_summary.csv", warnings, strict)

    _plot_overview(outdir, science_rows)
    _plot_time_order(outdir, order_rows)
    _plot_sample_chains(outdir, chain_rows, max_chains)
    _plot_chain_summary(outdir, chain_summary_rows, max_chains)

    findings = _findings_from_tables(science_rows, order_summary_rows, chain_summary_rows, warnings)
    write_table_bundle(findings, outdir / "tables" / "02_tpht_mode_sanity_findings", FINDING_COLUMNS)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    parser.add_argument("--max-chains", type=int, default=12, help="Maximum raw chain samples to draw.")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, FIGURE_OUTPUTS + TABLE_OUTPUTS)
        return
    ensure_output_dirs(outdir)
    plot_and_write(outdir, strict=bool(args.strict), max_chains=max(1, args.max_chains))


if __name__ == "__main__":
    main()
