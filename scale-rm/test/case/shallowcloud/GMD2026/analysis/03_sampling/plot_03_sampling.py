#!/usr/bin/env python3
"""Plot short 2D sampling-procedure verification diagnostics."""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402


def _partial_rows(rows: list[dict]) -> list[dict]:
    """Return sampled cases only."""
    return [row for row in rows if safe_float(row.get("sampling_fraction")) is not None and (safe_float(row.get("sampling_fraction")) or 0.0) < 1.0]


def _boxplot(outdir: Path, rows: list[dict]) -> None:
    """Plot DSD L1 error distributions by mode and fraction."""
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in _partial_rows(rows):
        value = safe_float(row.get("dsd_l1_error_vs_full_reference"))
        fraction = safe_float(row.get("sampling_fraction"))
        if value is not None and fraction is not None:
            grouped[f"{row.get('sample_mode')}\nf={fraction:g}"].append(value)
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("double", 0.38))
    if not grouped:
        no_data_panel(ax, "Sampling error boxplot")
    else:
        labels = list(grouped)
        ax.boxplot([grouped[label] for label in labels], tick_labels=labels, patch_artist=True, boxprops={"facecolor": PALETTE[0], "alpha": 0.45})
        ax.set_ylabel("DSD L1 error (-)")
        ax.set_title("Sampling error")
        ax.tick_params(axis="x", rotation=35)
    save_figure(fig, outdir / "figures", "03_sampling_error_boxplot")
    plt.close(fig)


def _error_vs_fraction(outdir: Path, rows: list[dict]) -> None:
    """Plot mean DSD L1 error versus sampling fraction."""
    grouped: dict[str, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in _partial_rows(rows):
        value = safe_float(row.get("dsd_l1_error_vs_full_reference"))
        fraction = safe_float(row.get("sampling_fraction"))
        if value is not None and fraction is not None:
            grouped[str(row.get("sample_mode"))][fraction].append(value)
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    if not grouped:
        no_data_panel(ax, "Error vs fraction")
    else:
        for index, (mode, by_fraction) in enumerate(sorted(grouped.items())):
            xs = sorted(by_fraction)
            ys = [sum(by_fraction[x]) / len(by_fraction[x]) for x in xs]
            ax.plot(xs, ys, marker="o", color=PALETTE[index], label=mode)
        ax.set_xlabel("Sampling fraction (-)")
        ax.set_ylabel("Mean DSD L1 error (-)")
        ax.set_title("Error vs fraction")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "03_sampling_error_vs_fraction")
    plt.close(fig)


def _threshold_error(outdir: Path, rows: list[dict]) -> None:
    """Plot radius-threshold fractions by sampling fraction and mode."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    fractions = [safe_float(row.get("sampling_fraction")) for row in _partial_rows(rows)]
    values = [safe_float(row.get("fraction_r_ge_15um")) for row in _partial_rows(rows)]
    if not has_any_number(fractions) or not has_any_number(values):
        no_data_panel(ax, "Radius-threshold diagnostics")
    else:
        for index, mode in enumerate(("random", "stratified")):
            mode_rows = [row for row in _partial_rows(rows) if row.get("sample_mode") == mode]
            by_fraction: dict[float, list[float]] = defaultdict(list)
            for row in mode_rows:
                fraction = safe_float(row.get("sampling_fraction"))
                value = safe_float(row.get("fraction_r_ge_15um"))
                if fraction is not None and value is not None:
                    by_fraction[fraction].append(value)
            xs = sorted(by_fraction)
            ys = [sum(by_fraction[x]) / len(by_fraction[x]) for x in xs]
            if xs:
                ax.plot(xs, ys, marker="o", color=PALETTE[index], label=mode)
        ax.set_xlabel("Sampling fraction (-)")
        ax.set_ylabel("Fraction r >= 15 um (-)")
        ax.set_title("Radius-threshold diagnostic")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "03_sampling_radius_threshold_error")
    plt.close(fig)


def _dsd_example(outdir: Path, rows: list[dict]) -> None:
    """Plot a compact DSD proxy using median radius versus fraction."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    values = [safe_float(row.get("weighted_median_radius")) for row in rows]
    if not has_any_number(values):
        no_data_panel(ax, "DSD example")
    else:
        for index, mode in enumerate(("random", "stratified")):
            by_fraction: dict[float, list[float]] = defaultdict(list)
            for row in rows:
                if row.get("sample_mode") != mode:
                    continue
                fraction = safe_float(row.get("sampling_fraction"))
                value = safe_float(row.get("weighted_median_radius"))
                if fraction is not None and value is not None:
                    by_fraction[fraction].append(value)
            xs = sorted(by_fraction)
            ys = [sum(by_fraction[x]) / len(by_fraction[x]) for x in xs]
            if xs:
                ax.plot(xs, ys, marker="o", color=PALETTE[index], label=mode)
        ax.set_xlabel("Sampling fraction (-)")
        ax.set_ylabel("Median radius (m)")
        ax.set_title("DSD proxy example")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "03_sampling_dsd_example")
    plt.close(fig)


def plot(outdir: Path) -> None:
    """Create all group-03 figures."""
    rows = read_table(outdir / "tables" / "03_sampling_metrics.csv")
    _boxplot(outdir, rows)
    _error_vs_fraction(outdir, rows)
    _threshold_error(outdir, rows)
    _dsd_example(outdir, rows)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/03_sampling_*.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
