"""Publication-quality Matplotlib style for GMD2026 figures."""

from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

from .table_utils import read_csv_rows, safe_float

OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}
PALETTE = [
    OKABE_ITO["blue"],
    OKABE_ITO["vermillion"],
    OKABE_ITO["bluish_green"],
    OKABE_ITO["orange"],
    OKABE_ITO["reddish_purple"],
    OKABE_ITO["sky_blue"],
    OKABE_ITO["yellow"],
    OKABE_ITO["black"],
]


def configure_matplotlib() -> Any:
    """Configure Matplotlib with a compact publication style."""
    mpl_config = Path(tempfile.gettempdir()) / "gmd2026_matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
    xdg_cache = Path(tempfile.gettempdir()) / "gmd2026_xdg_cache"
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )
    return plt


def figure_size(width: str = "single", aspect: float = 0.68) -> tuple[float, float]:
    """Return figure size in inches for single- or double-column GMD figures."""
    width_mm = 85.0 if width == "single" else 170.0
    width_in = width_mm / 25.4
    return width_in, width_in * aspect


def save_figure(fig: Any, outdir: Path, stem: str, dpi: int = 600) -> None:
    """Save one figure as PDF, SVG, and high-resolution PNG."""
    outdir.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in ((".pdf", {}), (".svg", {}), ((".png"), {"dpi": max(600, dpi)})):
        fig.savefig(outdir / f"{stem}{suffix}", **kwargs)


def no_data_panel(ax: Any, title: str, message: str = "No available data") -> None:
    """Draw a minimal placeholder panel when all metrics are NA."""
    ax.set_title(title)
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def numeric_series(rows: list[dict[str, Any]], key: str) -> list[float | None]:
    """Return a numeric series with None for NA."""
    return [safe_float(row.get(key)) for row in rows]


def read_table(path: Path) -> list[dict[str, str]]:
    """Read a CSV table for plotting."""
    return read_csv_rows(path)


def has_any_number(values: Iterable[float | None]) -> bool:
    """Return True when a plotted series contains at least one finite number."""
    return any(value is not None and math.isfinite(value) for value in values)


def _format_bar_value(value: float) -> str:
    """Format absolute bar values compactly without hiding their scale."""
    absolute = abs(value)
    if absolute >= 1.0e6 or (0.0 < absolute < 1.0e-2):
        return f"{value:.2e}"
    if absolute >= 100.0:
        return f"{value:.0f}"
    if absolute >= 10.0:
        return f"{value:.1f}"
    return f"{value:.2g}"


def label_bars(ax: Any, bars: Any, fmt: str | None = None, log_y: bool = False) -> None:
    """Add compact absolute-value labels inside bars when explicitly requested."""
    for bar in bars:
        height = bar.get_height()
        if not math.isfinite(height) or height <= 0.0:
            continue
        text = fmt.format(height) if fmt is not None else _format_bar_value(height)
        y_position = height * (0.72 if log_y else 0.94)
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            y_position,
            text,
            ha="center",
            va="top",
            fontsize=5.5,
            rotation=90,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.5},
        )


def label_stacked_percentages(
    ax: Any,
    x_positions: list[float],
    values: list[float],
    bottoms: list[float],
    totals: list[float],
    min_fraction: float = 0.06,
) -> None:
    """Label sizable stacked-bar components with their percentage contribution."""
    for x_position, value, bottom, total in zip(x_positions, values, bottoms, totals):
        if not (math.isfinite(value) and math.isfinite(total)) or value <= 0.0 or total <= 0.0:
            continue
        fraction = value / total
        if fraction < min_fraction:
            continue
        ax.text(
            x_position,
            bottom + value / 2.0,
            f"{fraction * 100.0:.0f}%",
            ha="center",
            va="center",
            fontsize=5.5,
            color="black",
        )


def simple_bar(
    rows: list[dict[str, Any]],
    label_key: str,
    value_key: str,
    ylabel: str,
    title: str,
    outdir: Path,
    stem: str,
    reference: float | None = None,
    annotate: bool = False,
    log_y: bool = False,
) -> None:
    """Create a compact bar plot."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.78))
    labels = [str(row.get(label_key, "")) for row in rows]
    values = numeric_series(rows, value_key)
    if not has_any_number(values):
        no_data_panel(ax, title)
    else:
        plotted = [value if value is not None else 0.0 for value in values]
        colors = [PALETTE[index % len(PALETTE)] for index in range(len(plotted))]
        bars = ax.bar(range(len(plotted)), plotted, color=colors, edgecolor="black", linewidth=0.5)
        positive = [value for value in plotted if value > 0.0]
        if log_y and positive:
            ax.set_yscale("log")
            ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
        elif plotted:
            ax.set_ylim(0.0, max(plotted) * 1.08 if max(plotted) > 0.0 else 1.0)
        if annotate:
            label_bars(ax, bars, log_y=log_y)
        ax.set_xticks(range(len(labels)), labels, rotation=35, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if reference is not None:
            ax.axhline(reference, color="0.3", linewidth=0.8, linestyle="--")
    save_figure(fig, outdir, stem)
    plt.close(fig)
