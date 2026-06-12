#!/usr/bin/env python3
"""Plot TPHT target-linked coalescence participant radius pairs.

This GMD2026 diagnostic uses existing analysis source tables, not raw model
output.  It is similar in spirit to the older scatter_sd_r1_vs_sd_r2.py plot,
but it is restricted to coalescence records represented in the TPHT
target-linked event table.  The filled density is a log-binned frequency map,
which avoids the honeycomb appearance introduced by Matplotlib hexbin.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.plot_style import OKABE_ITO, configure_matplotlib, figure_size, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402

EVENT_TABLE = "02_tpht_target_event_links.csv"
SINGLE_FIGURE_STEM = "supp_candidate_TPHT_coalescence_radius_pairs_all_scale"
COMPARISON_FIGURE_STEM = "supp_candidate_Hall_kernel_vs_TPHT_coalescence_radius_pairs"

# Hall kernel table copied from the reference halls_kernel.py script so the
# comparison panel can be redrawn with the same font and axis style as the TPHT
# diagnostic panel.
HALL_ECOLL = np.array(
    [
        0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010,
        0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0030, 0.0030, 0.0030, 0.0040, 0.0050,
        0.0050, 0.0050, 0.0100, 0.1000, 0.0500, 0.2000, 0.5000, 0.7700, 0.8700, 0.9700,
        0.0070, 0.0070, 0.0070, 0.0080, 0.0090, 0.0100, 0.0100, 0.0700, 0.4000, 0.4300,
        0.5800, 0.7900, 0.9300, 0.9600, 1.0000, 0.0090, 0.0090, 0.0090, 0.0120, 0.0150,
        0.0100, 0.0200, 0.2800, 0.6000, 0.6400, 0.7500, 0.9100, 0.9700, 0.9800, 1.0000,
        0.0140, 0.0140, 0.0140, 0.0150, 0.0160, 0.0300, 0.0600, 0.5000, 0.7000, 0.7700,
        0.8400, 0.9500, 0.9700, 1.0000, 1.0000, 0.0170, 0.0170, 0.0170, 0.0200, 0.0220,
        0.0600, 0.1000, 0.6200, 0.7800, 0.8400, 0.8800, 0.9500, 1.0000, 1.0000, 1.0000,
        0.0300, 0.0300, 0.0240, 0.0220, 0.0320, 0.0620, 0.2000, 0.6800, 0.8300, 0.8700,
        0.9000, 0.9500, 1.0000, 1.0000, 1.0000, 0.0250, 0.0250, 0.0250, 0.0360, 0.0430,
        0.1300, 0.2700, 0.7400, 0.8600, 0.8900, 0.9200, 1.0000, 1.0000, 1.0000, 1.0000,
        0.0270, 0.0270, 0.0270, 0.0400, 0.0520, 0.2000, 0.4000, 0.7800, 0.8800, 0.9000,
        0.9400, 1.0000, 1.0000, 1.0000, 1.0000, 0.0300, 0.0300, 0.0300, 0.0470, 0.0640,
        0.2500, 0.5000, 0.8000, 0.9000, 0.9100, 0.9500, 1.0000, 1.0000, 1.0000, 1.0000,
        0.0400, 0.0400, 0.0330, 0.0370, 0.0680, 0.2400, 0.5500, 0.8000, 0.9000, 0.9100,
        0.9500, 1.0000, 1.0000, 1.0000, 1.0000, 0.0350, 0.0350, 0.0350, 0.0550, 0.0790,
        0.2900, 0.5800, 0.8000, 0.9000, 0.9100, 0.9500, 1.0000, 1.0000, 1.0000, 1.0000,
        0.0370, 0.0370, 0.0370, 0.0620, 0.0820, 0.2900, 0.5900, 0.7800, 0.9000, 0.9100,
        0.9500, 1.0000, 1.0000, 1.0000, 1.0000, 0.0370, 0.0370, 0.0370, 0.0600, 0.0800,
        0.2900, 0.5800, 0.7700, 0.8900, 0.9100, 0.9500, 1.0000, 1.0000, 1.0000, 1.0000,
        0.0370, 0.0370, 0.0370, 0.0410, 0.0750, 0.2500, 0.5400, 0.7600, 0.8800, 0.9200,
        0.9500, 1.0000, 1.0000, 1.0000, 1.0000, 0.0370, 0.0370, 0.0370, 0.0520, 0.0670,
        0.2500, 0.5100, 0.7700, 0.8800, 0.9300, 0.9700, 1.0000, 1.0000, 1.0000, 1.0000,
        0.0370, 0.0370, 0.0370, 0.0470, 0.0570, 0.2500, 0.4900, 0.7700, 0.8900, 0.9500,
        1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0360, 0.0360, 0.0360, 0.0420, 0.0480,
        0.2300, 0.4700, 0.7800, 0.9200, 1.0000, 1.0200, 1.0200, 1.0200, 1.0200, 1.0200,
        0.0400, 0.0400, 0.0350, 0.0330, 0.0400, 0.1120, 0.4500, 0.7900, 1.0100, 1.0300,
        1.0400, 1.0400, 1.0400, 1.0400, 1.0400, 0.0330, 0.0330, 0.0330, 0.0330, 0.0330,
        0.1190, 0.4700, 0.9500, 1.3000, 1.7000, 2.3000, 2.3000, 2.3000, 2.3000, 2.3000,
        0.0270, 0.0270, 0.0270, 0.0270, 0.0270, 0.1250, 0.5200, 1.4000, 2.3000, 3.0000,
        4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
    ]
).reshape(21, 15).T
HALL_R0COL = np.array([6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 100.0, 150.0, 200.0, 300.0])
HALL_RATCOL = np.array([0.00, 0.050, 0.10, 0.150, 0.20, 0.250, 0.30, 0.350, 0.40, 0.450, 0.50, 0.550, 0.60, 0.650, 0.70, 0.750, 0.80, 0.850, 0.90, 0.950, 1.00])


def _event_table_candidates(gmd_outdir: Path) -> list[Path]:
    """Return likely locations for the target-linked event table."""
    return [
        gmd_outdir / "figure_sources" / "supplement_candidate_source_tables" / EVENT_TABLE,
        gmd_outdir / "tables" / EVENT_TABLE,
        gmd_outdir / EVENT_TABLE,
    ]


def _resolve_event_table(gmd_outdir: Path, explicit: Path | None) -> Path:
    """Find the existing event-link CSV table."""
    if explicit is not None:
        return explicit
    for candidate in _event_table_candidates(gmd_outdir):
        if candidate.exists():
            return candidate
    return _event_table_candidates(gmd_outdir)[0]


def _iter_paired_events(path: Path):
    """Yield adjacent participant-1/participant-2 coalescence records."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        previous: dict[str, str] | None = None
        for row in reader:
            if previous is None:
                previous = row
                continue
            if (
                previous.get("participant") == "1"
                and row.get("participant") == "2"
                and previous.get("source_file") == row.get("source_file")
                and previous.get("event_time_s") == row.get("event_time_s")
            ):
                yield previous, row
                previous = None
            else:
                previous = row


def _event_weight(first: dict[str, str], second: dict[str, str]) -> tuple[float, bool]:
    """Return event weight and whether multiplicity was available."""
    num_col = safe_float(first.get("num_col"))
    if num_col is None or num_col <= 0.0:
        num_col = safe_float(second.get("num_col"))
    if num_col is None or num_col <= 0.0:
        num_col = 1.0

    n1 = safe_float(first.get("sd_n1") or first.get("multiplicity_1") or first.get("sd_n"))
    n2 = safe_float(second.get("sd_n2") or second.get("multiplicity_2") or second.get("sd_n"))
    if n1 is None or n2 is None or n1 <= 0.0 or n2 <= 0.0:
        return num_col, False
    return num_col * min(n1, n2), True


def _read_radius_pairs(
    path: Path,
    radius_min_um: float | None,
    radius_max_um: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Read larger/smaller pre-coalescence radii for paired TPHT target-linked events."""
    larger: list[float] = []
    smaller: list[float] = []
    weights: list[float] = []
    paired_count = 0
    invalid_count = 0
    outside_count = 0
    multiplicity_weight_count = 0
    max_radius = None
    min_radius = None

    for first, second in _iter_paired_events(path):
        paired_count += 1
        radius_1 = safe_float(first.get("pre_radius_um"))
        radius_2 = safe_float(second.get("pre_radius_um"))
        if radius_1 is None or radius_2 is None or radius_1 <= 0.0 or radius_2 <= 0.0:
            invalid_count += 1
            continue

        big = max(radius_1, radius_2)
        small = min(radius_1, radius_2)
        max_radius = big if max_radius is None else max(max_radius, big)
        min_radius = small if min_radius is None else min(min_radius, small)
        if radius_min_um is not None and small < radius_min_um:
            outside_count += 1
            continue
        if radius_max_um is not None and big > radius_max_um:
            outside_count += 1
            continue

        weight, used_multiplicity = _event_weight(first, second)
        if used_multiplicity:
            multiplicity_weight_count += 1
        larger.append(big)
        smaller.append(small)
        weights.append(weight)

    if multiplicity_weight_count:
        weighting = "num_col * min(sd_n1, sd_n2) where multiplicity columns were available"
    else:
        weighting = "num_col only; multiplicity columns were not available in this source table"

    summary = {
        "event_table": str(path),
        "paired_event_count": paired_count,
        "plotted_pair_count": len(larger),
        "invalid_pair_count": invalid_count,
        "outside_radius_limit_count": outside_count,
        "radius_min_um": radius_min_um if radius_min_um is not None else "data_min_positive",
        "radius_max_um": radius_max_um if radius_max_um is not None else "data_max_positive",
        "minimum_smaller_radius_um": min_radius,
        "maximum_larger_radius_um": max_radius,
        "weighting": weighting,
        "multiplicity_weighted_pair_count": multiplicity_weight_count,
        "frequency_normalization": "histogram values are divided by total plotted weight",
        "scope": "coalescence records represented in the TPHT target-linked event table with both participants present",
    }
    return np.asarray(larger), np.asarray(smaller), np.asarray(weights), summary


def _write_summary(outdir: Path, summary: dict[str, Any]) -> None:
    """Write compact machine-readable summaries for figure provenance."""
    table_dir = outdir / "tables" / "diagnostics"
    table_dir.mkdir(parents=True, exist_ok=True)
    for stem in (SINGLE_FIGURE_STEM, COMPARISON_FIGURE_STEM):
        csv_path = table_dir / f"{stem}_summary.csv"
        json_path = table_dir / f"{stem}_summary.json"
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary))
            writer.writeheader()
            writer.writerow(summary)
        with json_path.open("w") as handle:
            json.dump(summary, handle, indent=2)


def _log_limits(values: np.ndarray, radius_min_um: float | None, radius_max_um: float | None) -> tuple[float, float]:
    """Choose positive log-axis limits from requested limits and data."""
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        raise RuntimeError("no positive radius values are available for log-scale plotting")
    lower = radius_min_um if radius_min_um is not None else float(positive.min())
    upper = radius_max_um if radius_max_um is not None else float(positive.max())
    lower = max(lower, np.nextafter(0.0, 1.0))
    if upper <= lower:
        upper = lower * 10.0
    # Use a small data-relative margin instead of full decade rounding.  This
    # keeps the all-scale plot faithful while avoiding the empty area caused by
    # expanding 0.007--247 um to 0.001--1000 um.
    return lower / 1.12, upper * 1.12


def _hall_collection_efficiency(larger_radius_um: float, ratio: float) -> float:
    """Interpolate the Hall collection-efficiency table from halls_kernel.py."""
    irr = int(np.searchsorted(HALL_R0COL, larger_radius_um, side="left")) + 1
    iqq = int(np.searchsorted(HALL_RATCOL[1:], ratio, side="left")) + 2
    iqq = max(2, min(iqq, len(HALL_RATCOL)))

    q = (ratio - HALL_RATCOL[iqq - 2]) / (HALL_RATCOL[iqq - 1] - HALL_RATCOL[iqq - 2])
    if irr >= 16:
        efficiency = (1.0 - q) * HALL_ECOLL[14, iqq - 2] + q * HALL_ECOLL[14, iqq - 1]
        return min(efficiency, 1.0)
    if 2 <= irr < 16:
        p = (larger_radius_um - HALL_R0COL[irr - 2]) / (HALL_R0COL[irr - 1] - HALL_R0COL[irr - 2])
        return (
            (1.0 - p) * (1.0 - q) * HALL_ECOLL[irr - 2, iqq - 2]
            + p * (1.0 - q) * HALL_ECOLL[irr - 1, iqq - 2]
            + (1.0 - p) * q * HALL_ECOLL[irr - 2, iqq - 1]
            + p * q * HALL_ECOLL[irr - 1, iqq - 1]
        )
    return (1.0 - q) * HALL_ECOLL[0, iqq - 2] + q * HALL_ECOLL[0, iqq - 1]


def _hall_kernel_grid(
    radius_limits: tuple[float, float],
    grid_size: int = 280,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the Hall kernel field on the original symmetric R_j/R_k axes."""
    radius_j_values = np.logspace(math.log10(radius_limits[0]), math.log10(radius_limits[1]), grid_size)
    radius_k_values = np.logspace(math.log10(radius_limits[0]), math.log10(radius_limits[1]), grid_size)
    radius_j, radius_k = np.meshgrid(radius_j_values, radius_k_values)
    larger = np.maximum(radius_j, radius_k)
    smaller = np.minimum(radius_j, radius_k)
    ratio = np.clip(smaller / larger, HALL_RATCOL[1] * 0.2, 1.0)

    vector_efficiency = np.vectorize(_hall_collection_efficiency, otypes=[float])
    efficiency = vector_efficiency(larger, ratio)
    delta_v = 1.19e8 * ((larger * 1e-4) ** 2 - (smaller * 1e-4) ** 2)
    sweep_area = np.pi * ((larger + smaller) * 1e-4) ** 2
    kernel = efficiency * delta_v * sweep_area
    kernel = np.ma.masked_invalid(kernel)
    kernel = np.ma.masked_less_equal(kernel, 0.0)
    return radius_j, radius_k, kernel


def _histogram_frequency(
    larger: np.ndarray,
    smaller: np.ndarray,
    weights: np.ndarray,
    bins: int,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a normalized log-binned 2-D histogram."""
    x_edges = np.logspace(math.log10(x_limits[0]), math.log10(x_limits[1]), bins + 1)
    y_edges = np.logspace(math.log10(y_limits[0]), math.log10(y_limits[1]), bins + 1)
    hist, x_edges, y_edges = np.histogram2d(larger, smaller, bins=(x_edges, y_edges), weights=weights)
    total = float(np.nansum(hist))
    if total > 0.0:
        hist = hist / total
    hist = np.ma.masked_less_equal(hist.T, 0.0)
    return hist, x_edges, y_edges


def _frequency_label(summary: dict[str, Any]) -> str:
    """Return a colorbar label that matches the available weighting evidence."""
    if int(summary.get("multiplicity_weighted_pair_count") or 0) > 0:
        return "Multiplicity-weighted frequency"
    return "Record-weighted frequency"


def _plot_radius_pair_panel(
    fig: Any,
    ax: Any,
    larger: np.ndarray,
    smaller: np.ndarray,
    weights: np.ndarray,
    summary: dict[str, Any],
    bins: int,
    radius_min_um: float | None,
    radius_max_um: float | None,
    *,
    add_panel_label: str | None = None,
    add_colorbar: bool = True,
    radius_limits: tuple[float, float] | None = None,
) -> Any:
    """Plot one log-binned TPHT coalescence radius-pair panel."""
    from matplotlib.colors import LogNorm

    if radius_limits is None:
        radius_limits = _log_limits(np.concatenate([larger, smaller]), radius_min_um, radius_max_um)
    x_limits = radius_limits
    y_limits = radius_limits
    hist, x_edges, y_edges = _histogram_frequency(larger, smaller, weights, bins, x_limits, y_limits)
    positive_hist = hist.compressed()
    norm = LogNorm(vmin=float(positive_hist.min()), vmax=float(positive_hist.max())) if positive_hist.size else None

    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        hist,
        cmap="plasma",
        norm=norm,
        shading="auto",
        rasterized=True,
    )
    line_min = max(x_limits[0], y_limits[0])
    line_max = min(x_limits[1], y_limits[1])
    if line_min < line_max:
        ax.plot([line_min, line_max], [line_min, line_max], color=OKABE_ITO["black"], linestyle="--", linewidth=0.9)
    ax.text(0.86, 0.91, "1:1", transform=ax.transAxes, ha="center", va="center", fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, which="major", color="0.82", linewidth=0.45)
    ax.grid(False, which="minor")
    ax.tick_params(top=True, right=True)
    ax.set_xlabel("Larger super-droplet radius\nbefore coalescence ($\\mu$m)")
    ax.set_ylabel("Smaller super-droplet radius\nbefore coalescence ($\\mu$m)")
    if add_panel_label:
        ax.text(
            -0.10,
            1.04,
            add_panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            fontweight="bold",
            clip_on=False,
        )

    info = (
        f"paired records: {summary['paired_event_count']:,}\n"
        f"plotted pairs: {summary['plotted_pair_count']:,}\n"
        f"{_frequency_label(summary).lower()}"
    )
    ax.text(
        0.04,
        0.94,
        info,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.75", "alpha": 0.90},
    )
    if add_colorbar:
        cbar = fig.colorbar(mesh, ax=ax, pad=0.02, shrink=0.68)
        cbar.set_label(_frequency_label(summary))
    return mesh


def _plot_hall_kernel_panel(
    ax: Any,
    radius_limits: tuple[float, float],
    *,
    add_panel_label: str | None = None,
) -> None:
    """Plot the Hall kernel panel on the TPHT radius range."""
    from matplotlib.colors import LogNorm

    radius_j, radius_k, kernel = _hall_kernel_grid(radius_limits)
    positive = kernel.compressed()
    if positive.size == 0:
        raise RuntimeError("Hall kernel grid did not contain positive values")
    low_level = max(10.0 ** math.floor(math.log10(float(positive.min()))), 1e-14)
    high_level = 10.0 ** math.ceil(math.log10(float(positive.max())))
    mesh = ax.pcolormesh(
        radius_j,
        radius_k,
        kernel,
        shading="auto",
        cmap="Blues",
        norm=LogNorm(vmin=low_level, vmax=high_level),
        rasterized=True,
    )
    levels = np.logspace(math.log10(low_level), math.log10(high_level), 9)
    ax.contour(radius_j, radius_k, kernel, levels=levels, colors="#313695", linewidths=0.55, alpha=0.55)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*radius_limits)
    ax.set_ylim(*radius_limits)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, which="major", color="0.82", linewidth=0.45)
    ax.grid(False, which="minor")
    ax.tick_params(top=True, right=True)
    ax.set_xlabel("Droplet radius $R_j$ ($\\mu$m)")
    ax.set_ylabel("Droplet radius $R_k$ ($\\mu$m)")
    ax.set_facecolor("#f8f9fa")
    cbar = ax.figure.colorbar(mesh, ax=ax, pad=0.02, shrink=0.68)
    cbar.set_label("Hall kernel (cm$^3$ s$^{-1}$)")
    if add_panel_label:
        ax.text(
            -0.10,
            1.04,
            add_panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            fontweight="bold",
            clip_on=False,
        )


def plot(
    gmd_outdir: Path,
    outdir: Path,
    event_table: Path | None,
    radius_min_um: float | None,
    radius_max_um: float | None,
    bins: int,
    hall_image: Path | None,
) -> None:
    """Create all-scale TPHT radius-pair figures."""
    resolved_table = _resolve_event_table(gmd_outdir, event_table)
    if not resolved_table.exists():
        raise FileNotFoundError(f"event table not found: {resolved_table}")

    larger, smaller, weights, summary = _read_radius_pairs(resolved_table, radius_min_um, radius_max_um)
    if larger.size == 0:
        raise RuntimeError("no valid target-linked coalescence radius pairs were available for plotting")
    summary["hall_kernel_reference_image"] = str(hall_image) if hall_image is not None else "NA"
    summary["plot_type"] = "double-log pcolormesh from log-spaced 2-D histogram, not hexbin"
    summary["hall_kernel_panel"] = "redrawn from halls_kernel.py table and formula with matched Matplotlib style"

    plt = configure_matplotlib()
    figure_dir = outdir / "figures" / "supplement_candidates"
    radius_limits = _log_limits(np.concatenate([larger, smaller]), radius_min_um, radius_max_um)
    summary["x_axis_radius_um"] = f"{radius_limits[0]:.6g}--{radius_limits[1]:.6g}"
    summary["y_axis_radius_um"] = f"{radius_limits[0]:.6g}--{radius_limits[1]:.6g}"

    fig, ax = plt.subplots(figsize=figure_size("single", 1.06), constrained_layout=True)
    _plot_radius_pair_panel(
        fig,
        ax,
        larger,
        smaller,
        weights,
        summary,
        bins,
        radius_min_um,
        radius_max_um,
        radius_limits=radius_limits,
    )
    save_figure(fig, figure_dir, SINGLE_FIGURE_STEM)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=figure_size("double", 0.56), constrained_layout=True)
    _plot_hall_kernel_panel(axes[0], radius_limits, add_panel_label="a")
    mesh = _plot_radius_pair_panel(
        fig,
        axes[1],
        larger,
        smaller,
        weights,
        summary,
        bins,
        radius_min_um,
        radius_max_um,
        add_panel_label="b",
        add_colorbar=False,
        radius_limits=radius_limits,
    )
    cbar = fig.colorbar(mesh, ax=axes[1], pad=0.02, shrink=0.68)
    cbar.set_label(_frequency_label(summary))
    save_figure(fig, figure_dir, COMPARISON_FIGURE_STEM)
    plt.close(fig)

    _write_summary(outdir, summary)
    print(f"wrote {figure_dir / (SINGLE_FIGURE_STEM + '.png')}")
    print(f"wrote {figure_dir / (COMPARISON_FIGURE_STEM + '.png')}")
    print(json.dumps(summary, indent=2))


def _optional_float(value: str) -> float | None:
    """Parse an optional float CLI value."""
    if value.upper() in {"NA", "NONE", "AUTO"}:
        return None
    return float(value)


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmd-outdir", type=Path, required=True, help="analysis_outputs_for_GMD directory with figure_sources tables")
    parser.add_argument("--outdir", type=Path, required=True, help="output directory for generated figure and summary")
    parser.add_argument("--event-table", type=Path, default=None, help="optional explicit 02_tpht_target_event_links.csv path")
    parser.add_argument("--radius-min-um", type=_optional_float, default=None, help="minimum positive radius; default uses data minimum")
    parser.add_argument("--radius-max-um", type=_optional_float, default=None, help="maximum radius; default uses data maximum")
    parser.add_argument("--bins", type=int, default=90, help="number of log-spaced bins on each axis")
    parser.add_argument("--hall-image", type=Path, default=None, help="optional hall_kernel_contour.png path recorded for provenance only")
    parser.add_argument("--dry-run", action="store_true", help="show resolved inputs without writing figures")
    args = parser.parse_args()

    resolved_table = _resolve_event_table(args.gmd_outdir, args.event_table)
    if args.dry_run:
        print(f"[dry-run] event_table={resolved_table}")
        print(f"[dry-run] outdir={args.outdir}")
        print(f"[dry-run] radius_min_um={args.radius_min_um}")
        print(f"[dry-run] radius_max_um={args.radius_max_um}")
        print(f"[dry-run] bins={args.bins}")
        print(f"[dry-run] hall_image={args.hall_image}")
        return
    plot(
        args.gmd_outdir,
        args.outdir,
        args.event_table,
        args.radius_min_um,
        args.radius_max_um,
        args.bins,
        args.hall_image,
    )


if __name__ == "__main__":
    main()
