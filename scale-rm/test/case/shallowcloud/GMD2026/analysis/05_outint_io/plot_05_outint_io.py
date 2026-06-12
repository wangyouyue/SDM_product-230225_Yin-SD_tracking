#!/usr/bin/env python3
"""Plot GMD2026 selected-output interval I/O sensitivity."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402


def _line_metric(rows: list[dict], metric: str, ylabel: str, title: str, stem: str, outdir: Path, scale: float = 1.0) -> None:
    """Plot an interval-sensitivity metric for FW005 and BW005."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    values = [(safe_float(row.get(metric)) * scale) if safe_float(row.get(metric)) is not None else None for row in rows]
    if not has_any_number(values):
        no_data_panel(ax, title)
    else:
        plotted = False
        for index, label in enumerate(("FW005", "BW005")):
            label_rows = sorted(
                [row for row in rows if row.get("tracking_label") == label],
                key=lambda row: safe_float(row.get("output_interval_s")) or -1,
            )
            xs = [safe_float(row.get("output_interval_s")) for row in label_rows]
            ys = [(safe_float(row.get(metric)) * scale) if safe_float(row.get(metric)) is not None else None for row in label_rows]
            if has_any_number(xs) and has_any_number(ys):
                ax.plot([x or 0.0 for x in xs], [y or 0.0 for y in ys], marker="o", color=PALETTE[index], label=label)
                plotted = True
        if not plotted:
            no_data_panel(ax, title)
            save_figure(fig, outdir / "figures", stem)
            plt.close(fig)
            return
        ax.set_xlabel("Output interval (s)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", stem)
    plt.close(fig)


def plot(outdir: Path) -> None:
    """Create all group-05 figures."""
    rows = read_table(outdir / "tables" / "05_outint_io_summary.csv")
    _line_metric(rows, "sd_selected_output_bytes", "Selected output (MiB)", "Selected output size vs interval", "05_output_size_vs_interval", outdir, scale=1.0 / 1024.0**2)
    _line_metric(rows, "sd_output_write_time_total_s", "Write time total (s)", "Write time vs interval", "05_write_time_vs_interval", outdir)
    _line_metric(rows, "sd_selected_file_count", "Selected SD file count", "Selected SD file count vs interval", "05_file_count_vs_interval", outdir)
    _line_metric(rows, "wallclock_overhead_vs_nt_s", "Wallclock overhead (s)", "Wallclock overhead vs interval", "05_wallclock_overhead_vs_interval", outdir)
    _line_metric(rows, "output_bytes_per_tracked_chain", "Bytes per tracked chain", "Output efficiency", "05_output_efficiency", outdir)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/05_*.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
