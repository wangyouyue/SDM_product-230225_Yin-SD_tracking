#!/usr/bin/env python3
"""Plot group-01 controlled 3D sampled benchmark diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, label_stacked_percentages, no_data_panel, read_table, save_figure, simple_bar  # noqa: E402
from common.table_utils import safe_float  # noqa: E402


def _labels(rows: list[dict]) -> list[str]:
    """Return compact case labels."""
    return [str(row.get("case_name", "")).replace("_coallog", "").replace("_", "\n") for row in rows]


def _stacked_output(rows: list[dict], outdir: Path) -> None:
    """Plot scientific output byte components using MiB units."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.78))
    labels = _labels(rows)
    to_mib = 1.0 / 1024.0**2
    selected = [(safe_float(row.get("sd_selected_output_bytes")) or 0.0) * to_mib for row in rows]
    sd_all = [(safe_float(row.get("sd_all_output_bytes")) or 0.0) * to_mib for row in rows]
    coal = [(safe_float(row.get("coalescence_log_bytes")) or 0.0) * to_mib for row in rows]
    history = [(safe_float(row.get("history_output_bytes")) or 0.0) * to_mib for row in rows]
    ids = [(safe_float(row.get("tpht_id_bytes")) or 0.0) * to_mib for row in rows]
    total = [(safe_float(row.get("scientific_output_bytes")) or safe_float(row.get("total_output_bytes")) or 0.0) * to_mib for row in rows]
    other = [max(0.0, t - s - a - c - h - i) for t, s, a, c, h, i in zip(total, selected, sd_all, coal, history, ids)]
    if not has_any_number(total):
        no_data_panel(ax, "Output size")
    else:
        bottoms = [0.0] * len(rows)
        x_positions = list(range(len(rows)))
        for values, label, color in (
            (selected, "selected SD", PALETTE[0]),
            (sd_all, "all SD", PALETTE[2]),
            (coal, "coalescence log", PALETTE[1]),
            (history, "history.pe*", PALETTE[3]),
            (ids, "TPHT ids", PALETTE[4]),
            (other, "other scientific", "0.75"),
        ):
            if any(value > 0.0 for value in values):
                ax.bar(x_positions, values, bottom=bottoms, label=label, color=color, edgecolor="black", linewidth=0.4)
                label_stacked_percentages(ax, x_positions, values, bottoms, total, min_fraction=0.05)
                bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        ax.set_xticks(x_positions, labels, rotation=35, ha="right")
        ax.set_ylabel("Scientific output size (MiB)")
        ax.set_title("Output size")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "01_output_size")
    plt.close(fig)


def _runtime_components(rows: list[dict], outdir: Path) -> None:
    """Plot measured runtime components and remaining model runtime."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("double", 0.45))
    labels = _labels(rows)
    keys = [
        ("id_assignment_time_s", "ID assignment", PALETTE[0]),
        ("boundary_tracking_time_s", "boundary tracking", PALETTE[1]),
        ("sd_output_write_time_total_s", "SD write", PALETTE[2]),
        ("coalescence_output_write_time_total_s", "coal write", PALETTE[3]),
    ]
    bottoms = [0.0] * len(rows)
    if not has_any_number([safe_float(row.get("wallclock_s")) for row in rows]):
        no_data_panel(ax, "Runtime components")
    else:
        x_positions = list(range(len(rows)))
        wall = [safe_float(row.get("wallclock_s")) or 0.0 for row in rows]
        for key, label, color in keys:
            values = [safe_float(row.get(key)) or 0.0 for row in rows]
            ax.bar(x_positions, values, bottom=bottoms, label=label, color=color, edgecolor="black", linewidth=0.4)
            label_stacked_percentages(ax, x_positions, values, bottoms, wall, min_fraction=0.03)
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        other = [max(0.0, wall_value - bottom) for wall_value, bottom in zip(wall, bottoms)]
        ax.bar(x_positions, other, bottom=bottoms, label="other model runtime", color="0.82", edgecolor="black", linewidth=0.4)
        label_stacked_percentages(ax, x_positions, other, bottoms, wall, min_fraction=0.05)
        ax.set_xticks(x_positions, labels, rotation=35, ha="right")
        ax.set_ylabel("Wallclock component (s)")
        ax.set_title("Runtime components")
        ax.legend(frameon=False, ncol=3)
    save_figure(fig, outdir / "figures", "01_runtime_components")
    plt.close(fig)


def plot(outdir: Path) -> None:
    """Create all group-01 figures."""
    rows = read_table(outdir / "tables" / "01_benchmark_summary.csv")
    simple_bar(rows, "case_name", "wallclock_relative_to_nt_nolog", "Relative wallclock (-)", "Wallclock relative to nt_nolog", outdir / "figures", "01_wallclock_relative", reference=1.0, annotate=True)
    simple_bar(rows, "case_name", "core_hours", "Core-hours (h)", "Core-hours", outdir / "figures", "01_core_hours", annotate=True)
    _stacked_output(rows, outdir)
    simple_bar(rows, "case_name", "peak_memory_rank_max_mib", "Peak rank memory (MiB)", "Peak memory", outdir / "figures", "01_peak_memory", annotate=True)
    _runtime_components(rows, outdir)
    simple_bar(rows, "case_name", "tracking_chain_count", "Tracking chains (count)", "Tracking chain count", outdir / "figures", "01_tracking_chain_count", annotate=True, log_y=True)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/01_*.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
