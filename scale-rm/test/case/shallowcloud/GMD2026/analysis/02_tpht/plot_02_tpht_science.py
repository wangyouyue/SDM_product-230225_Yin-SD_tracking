#!/usr/bin/env python3
"""Plot TPHT threshold-crossing and history diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, label_bars, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402


PATHWAY_LABELS = {
    "coal_before_large": "coal -> large",
    "large_before_coal": "large -> coal",
    "simultaneous_large_and_coal": "same output",
    "large_without_coal": "large, no coal",
    "coal_without_large": "coal, not large",
    "no_large_or_coal": "neither diagnosed",
    "timing_unknown": "timing unknown",
}


def _bar_plot(rows: list[dict], label_key: str, value_key: str, ylabel: str, title: str, stem: str, outdir: Path, log_y: bool = False) -> None:
    """Draw a compact absolute-count bar plot with internal labels."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.74))
    labels = [str(row.get(label_key, "")) for row in rows]
    values = [safe_float(row.get(value_key)) for row in rows]
    if not rows or not has_any_number(values):
        no_data_panel(ax, title)
    else:
        x_positions = list(range(len(rows)))
        bars = ax.bar(
            x_positions,
            [value or 0.0 for value in values],
            color=[PALETTE[index % len(PALETTE)] for index in x_positions],
            edgecolor="black",
            linewidth=0.5,
        )
        positive = [value for value in values if value is not None and value > 0.0]
        if log_y and positive:
            ax.set_yscale("log")
            ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
            label_bars(ax, bars, log_y=True)
        else:
            ax.set_ylim(0.0, max(positive) * 1.08 if positive else 1.0)
            label_bars(ax, bars)
        ax.set_xticks(x_positions, labels, rotation=32, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    save_figure(fig, outdir / "figures", stem)
    plt.close(fig)


def _plot_formation_pathways(outdir: Path) -> None:
    """Plot the diagnostic order of threshold crossing and coalescence flags."""
    rows = read_table(outdir / "tables" / "02_tpht_science_pathways.csv")
    plot_rows = []
    for row in rows:
        pathway = str(row.get("pathway", ""))
        count = safe_float(row.get("target_count"))
        if pathway in {"", "NA", "None"} or count is None or count <= 0.0:
            continue
        plot_rows.append({"label": PATHWAY_LABELS.get(pathway, pathway), "target_count": count})
    _bar_plot(plot_rows, "label", "target_count", "Targets (count)", "TPHT diagnostic timing class", "02_tpht_science_formation_pathways", outdir, log_y=True)


def _plot_formation_height(outdir: Path) -> None:
    """Plot where targets first reached the radius threshold."""
    rows = read_table(outdir / "tables" / "02_tpht_science_formation_height_bins.csv")
    plot_rows = []
    for row in rows:
        label = str(row.get("height_bin_m", ""))
        count = safe_float(row.get("target_count"))
        if label in {"unknown", "", "NA", "None"} or count is None or count <= 0.0:
            continue
        plot_rows.append({"height_bin_m": label, "target_count": count})
    _bar_plot(plot_rows, "height_bin_m", "target_count", "Targets (count)", "First r >= 15 micrometer height", "02_tpht_science_formation_height", outdir)


def _plot_coalescence_counts(outdir: Path) -> None:
    """Plot a proxy for how often selected targets carried coalescence flags."""
    rows = read_table(outdir / "tables" / "02_tpht_science_coalescence_counts.csv")
    plot_rows = []
    for row in rows:
        label = str(row.get("coal_episode_count_proxy", ""))
        count = safe_float(row.get("target_count"))
        if label in {"", "NA", "None"} or count is None or count <= 0.0:
            continue
        plot_rows.append({"coal_episode_count_proxy": label, "target_count": count})
    _bar_plot(plot_rows, "coal_episode_count_proxy", "target_count", "Targets (count)", "if_coal occurrence proxy", "02_tpht_science_coalescence_counts", outdir, log_y=True)


def _plot_history(outdir: Path) -> None:
    """Plot mean reconstructed radius and height histories when time is available."""
    rows = read_table(outdir / "tables" / "02_tpht_science_time_series.csv")
    times = [(safe_float(row.get("time_s")) / 60.0) if safe_float(row.get("time_s")) is not None else None for row in rows]
    radius = [safe_float(row.get("mean_radius_um")) for row in rows]
    height = [safe_float(row.get("mean_height_m")) for row in rows]
    large_fraction = [safe_float(row.get("large_record_fraction")) for row in rows]
    coal_fraction = [safe_float(row.get("coal_record_fraction")) for row in rows]

    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=figure_size("double", 0.68), sharex=True)
    flat_axes = [axes[0][0], axes[0][1], axes[1][0], axes[1][1]]
    series = [
        (radius, r"Mean radius ($\mu$m)", "All-record mean radius", PALETTE[0]),
        (height, "Mean height (m)", "All-record mean height", PALETTE[1]),
        (large_fraction, "Fraction (-)", r"Records with r $\geq$ 15 $\mu$m", PALETTE[2]),
        (coal_fraction, "Fraction (-)", "Records with if_coal flag", PALETTE[3]),
    ]
    valid_times = [time for time in times if time is not None]
    if not has_any_number(valid_times):
        for ax in flat_axes:
            no_data_panel(ax, "TPHT history", "No time coordinate available")
    else:
        for ax, (values, ylabel, title, color) in zip(flat_axes, series):
            points = sorted((time, value) for time, value in zip(times, values) if time is not None and value is not None)
            if not points:
                no_data_panel(ax, title)
                continue
            xs, ys = zip(*points)
            ax.plot(xs, ys, color=color, marker="o", markersize=2.0)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
        for ax in axes[1]:
            ax.set_xlabel("TPHT output time (min)")
    save_figure(fig, outdir / "figures", "02_tpht_science_history")
    plt.close(fig)


def _plot_interest_condition_fallback(outdir: Path) -> None:
    """Plot first-selected interest conditions, excluding unknown diagnostics."""
    rows = read_table(outdir / "tables" / "02_tpht_target_categories.csv")
    keep = {"radius_only": "r >= 15 micrometer", "coal_only": "if_coal only", "both": "both"}
    plot_rows = []
    for row in rows:
        category = str(row.get("category", ""))
        count = safe_float(row.get("target_count"))
        if category in keep and count is not None and count > 0.0:
            plot_rows.append({"label": keep[category], "target_count": count})
    _bar_plot(plot_rows, "label", "target_count", "Targets (count)", "First-selected interest condition", "02_tpht_science_condition_composition", outdir, log_y=True)


def plot(outdir: Path) -> None:
    """Create TPHT science-oriented diagnostic figures."""
    _plot_formation_pathways(outdir)
    _plot_formation_height(outdir)
    _plot_coalescence_counts(outdir)
    _plot_history(outdir)
    _plot_interest_condition_fallback(outdir)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(
            Path(__file__).name,
            root,
            outdir,
            [
                "figures/02_tpht_science_formation_pathways.{pdf,svg,png}",
                "figures/02_tpht_science_formation_height.{pdf,svg,png}",
                "figures/02_tpht_science_coalescence_counts.{pdf,svg,png}",
                "figures/02_tpht_science_history.{pdf,svg,png}",
                "figures/02_tpht_science_condition_composition.{pdf,svg,png}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
