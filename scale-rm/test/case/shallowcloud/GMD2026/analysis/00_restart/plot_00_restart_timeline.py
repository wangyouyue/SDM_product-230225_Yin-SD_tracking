#!/usr/bin/env python3
"""Plot the GMD2026 restart sanity timeline schematic."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import configure_matplotlib, figure_size, save_figure  # noqa: E402


def plot(outdir: Path) -> None:
    """Draw cold-start to restart-output schematic."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.42))
    ax.annotate("", xy=(1.0, 0.5), xytext=(0.08, 0.5), arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "black"})
    ax.scatter([0.08, 1.0], [0.5, 0.5], s=28, color=["#0072B2", "#D55E00"], zorder=3)
    ax.text(0.08, 0.63, "cold start", ha="center", va="bottom")
    ax.text(1.0, 0.63, "60 min restart", ha="center", va="bottom")
    ax.text(0.54, 0.36, "Group 00 restart sanity only", ha="center", va="top", fontsize=7)
    ax.set_xlim(0, 1.08)
    ax.set_ylim(0, 1)
    ax.axis("off")
    save_figure(fig, outdir / "figures", "00_restart_timeline")
    plt.close(fig)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/00_restart_timeline.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()

