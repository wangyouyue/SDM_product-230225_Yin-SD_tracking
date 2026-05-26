#!/usr/bin/env python3
"""Plot GMD2026 SDNC scaling diagnostics."""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402


def _plot_series(ax, rows: list[dict], metric: str, label: str, color: str) -> bool:
    """Plot one tracking-label series and return whether anything was drawn."""
    label_rows = sorted([row for row in rows if row.get("tracking_label") == label], key=lambda row: safe_float(row.get("sdnc")) or -1)
    points = [(safe_float(row.get("sdnc")), safe_float(row.get(metric))) for row in label_rows]
    points = [(x_value, y_value) for x_value, y_value in points if x_value is not None and y_value is not None]
    if not points:
        return False
    xs, ys = zip(*points)
    ax.plot(xs, ys, marker="o", color=color, label=label)
    return True


def _merged_fw_bw_points(rows: list[dict], metric: str, tolerance: float = 0.01) -> tuple[list[float], list[float]] | None:
    """Return averaged FW/BW points when both modes are effectively identical."""
    by_label: dict[str, dict[float, float]] = defaultdict(dict)
    for row in rows:
        label = row.get("tracking_label")
        if label not in {"FW005", "BW005"}:
            continue
        sdnc = safe_float(row.get("sdnc"))
        value = safe_float(row.get(metric))
        if sdnc is None or value is None:
            continue
        by_label[label][sdnc] = value

    common_sdnc = sorted(set(by_label["FW005"]) & set(by_label["BW005"]))
    if not common_sdnc:
        return None

    xs: list[float] = []
    ys: list[float] = []
    for sdnc in common_sdnc:
        fw_value = by_label["FW005"][sdnc]
        bw_value = by_label["BW005"][sdnc]
        scale = max(abs(fw_value), abs(bw_value), 1.0)
        if abs(fw_value - bw_value) / scale > tolerance:
            return None
        xs.append(sdnc)
        ys.append((fw_value + bw_value) / 2.0)
    return xs, ys


def _line_metric(
    rows: list[dict],
    metric: str,
    ylabel: str,
    title: str,
    stem: str,
    outdir: Path,
    log_y: bool = False,
    merge_fw_bw: bool = False,
) -> None:
    """Plot a metric against SDNC for NT/FW/BW."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    all_values = [safe_float(row.get(metric)) for row in rows]
    if not has_any_number(all_values):
        no_data_panel(ax, title)
    else:
        plotted = False
        plotted = _plot_series(ax, rows, metric, "NT", PALETTE[0]) or plotted
        merged_points = _merged_fw_bw_points(rows, metric) if merge_fw_bw else None
        if merged_points is not None:
            xs, ys = merged_points
            ax.plot(xs, ys, marker="o", color=PALETTE[2], label="BW005/FW005")
            plotted = True
        else:
            plotted = _plot_series(ax, rows, metric, "FW005", PALETTE[1]) or plotted
            plotted = _plot_series(ax, rows, metric, "BW005", PALETTE[2]) or plotted
        if not plotted:
            no_data_panel(ax, title)
            save_figure(fig, outdir / "figures", stem)
            plt.close(fig)
            return
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel(r"SDNC (# grid-cell$^{-1}$)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", stem)
    plt.close(fig)


def _overhead_ratio(rows: list[dict], outdir: Path) -> None:
    """Plot FW/BW overhead ratios against SDNC."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    values = [safe_float(row.get("fw_overhead_ratio")) for row in rows] + [safe_float(row.get("bw_overhead_ratio")) for row in rows]
    if not has_any_number(values):
        no_data_panel(ax, "Tracking overhead ratio")
    else:
        by_sdnc: dict[float, dict[str, float]] = defaultdict(dict)
        for row in rows:
            sdnc = safe_float(row.get("sdnc"))
            if sdnc is None:
                continue
            if safe_float(row.get("fw_overhead_ratio")) is not None:
                by_sdnc[sdnc]["FW005"] = safe_float(row.get("fw_overhead_ratio")) or 0.0
            if safe_float(row.get("bw_overhead_ratio")) is not None:
                by_sdnc[sdnc]["BW005"] = safe_float(row.get("bw_overhead_ratio")) or 0.0
        xs = sorted(by_sdnc)
        for index, label in enumerate(("FW005", "BW005")):
            ys = [by_sdnc[x].get(label, 0.0) for x in xs]
            ax.plot(xs, ys, marker="o", color=PALETTE[index], label=label)
        ax.axhline(1.0, color="0.3", linestyle="--", linewidth=0.8)
        ax.set_xlabel(r"SDNC (# grid-cell$^{-1}$)")
        ax.set_ylabel("Wallclock ratio vs NT (-)")
        ax.set_title("Tracking overhead ratio")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "04_tracking_overhead_ratio")
    plt.close(fig)


def plot(outdir: Path) -> None:
    """Create all group-04 figures."""
    rows = read_table(outdir / "tables" / "04_sdnc_scaling_summary.csv")
    _line_metric(rows, "wallclock_s", "Wallclock (s)", "Wallclock vs SDNC", "04_wallclock_vs_sdnc", outdir)
    _line_metric(rows, "core_hours", "Core-hours (h)", "Core-hours vs SDNC", "04_core_hours_vs_sdnc", outdir)
    _line_metric(rows, "peak_memory_rank_max_mib", "Peak rank memory (MiB)", "Memory vs SDNC", "04_memory_vs_sdnc", outdir)
    _overhead_ratio(rows, outdir)
    _line_metric(rows, "tracking_chain_count", "Tracking chains (count)", "Chain count vs SDNC", "04_chain_count_vs_sdnc", outdir, log_y=True, merge_fw_bw=True)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/04_*.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
