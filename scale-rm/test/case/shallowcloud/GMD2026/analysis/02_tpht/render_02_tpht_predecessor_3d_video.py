#!/usr/bin/env python3
"""Render a GMD2026 TPHT 3-D QC/QR video with selected target trajectories.

The default target selection reproduces the three target IDs represented by
``supp_candidate_TPHT_predecessor_tree_examples_final10_07``: rows 19--21 of
the final-window selected-target table when three targets are drawn per figure.
Target markers are colored by radius using a cyan-magenta scale that is visually
distinct from the QC greyscale and QR yellow-red scale.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "gmd2026_matplotlib"))

try:  # Optional for --dry-run on machines without the rendering stack.
    import imageio.v2 as imageio
except Exception:  # pragma: no cover
    imageio = None

try:  # Optional for --dry-run on machines without Matplotlib.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
except Exception:  # pragma: no cover
    matplotlib = None
    plt = None
    ScalarMappable = None
    LinearSegmentedColormap = None
    LogNorm = None
    Normalize = None

try:  # Optional for --dry-run before the compact export exists.
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None

try:  # Optional, used only for visualization smoothing.
    from scipy.ndimage import gaussian_filter
except Exception:  # pragma: no cover
    gaussian_filter = None


EXPORT_FILE = "tpht_bw_qc_qr_3d_fullgrid.nc"
TRAJECTORY_TABLE = "02_tpht_target_trajectory_records.csv"
SELECTED_TABLE = "supp_candidate_TPHT_predecessor_tree_examples_final10_selected_targets.csv"
MARKERS = ("o", "^", "D")
LINE_COLORS = ("#00A6D6", "#CC79A7", "#009E73")


@dataclass(frozen=True)
class TargetTrajectory:
    """One selected TPHT target trajectory."""

    target_id: str
    marker: str
    line_color: str
    time_s: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    radius_um: np.ndarray


@dataclass(frozen=True)
class ScanStats:
    """Robust color scale statistics for one 3-D field export."""

    qc_vmax_gkg: float
    qr_vmax_gkg: float
    radius_vmin_um: float
    radius_vmax_um: float


def parse_args() -> argparse.Namespace:
    """Return command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Path to the GMD2026 directory.")
    parser.add_argument("--analysis-outdir", type=Path, default=None, help="Existing analysis_outputs directory with TPHT science tables.")
    parser.add_argument("--gmd-outdir", type=Path, default=None, help="analysis_outputs_for_GMD directory.")
    parser.add_argument("--export-file", type=Path, default=None, help="Compact QC/QR export NetCDF.")
    parser.add_argument("--figure-index", type=int, default=7, help="Final-window predecessor-tree figure index.")
    parser.add_argument("--targets-per-figure", type=int, default=3, help="Targets per predecessor-tree figure.")
    parser.add_argument("--target-ids", nargs="*", default=None, help="Override target IDs.")
    parser.add_argument("--qc-threshold-gkg", type=float, default=0.001, help="Visual QC threshold in g/kg.")
    parser.add_argument("--qr-threshold-gkg", type=float, default=0.0001, help="Visual QR threshold in g/kg.")
    parser.add_argument("--max-qc-points", type=int, default=80_000, help="Maximum plotted QC points per frame.")
    parser.add_argument("--max-qr-points", type=int, default=80_000, help="Maximum plotted QR points per frame.")
    parser.add_argument("--smooth-qc-sigma-cells", type=float, default=1.25, help="Gaussian smoothing sigma for QC visualization in grid-cell units.")
    parser.add_argument("--smooth-qr-sigma-cells", type=float, default=0.0, help="Gaussian smoothing sigma for QR visualization in grid-cell units.")
    parser.add_argument("--qc-alpha", type=float, default=0.09, help="QC point-cloud alpha; lower values are more transparent.")
    parser.add_argument("--qr-alpha", type=float, default=0.70, help="QR point-cloud alpha.")
    parser.add_argument("--qc-marker-size", type=float, default=5.0, help="QC square marker size.")
    parser.add_argument("--qr-marker-size", type=float, default=11.0, help="QR square marker size.")
    parser.add_argument("--color-percentile", type=float, default=99.5, help="Robust QC/QR color upper percentile.")
    parser.add_argument("--fps", type=int, default=8, help="Video frame rate.")
    parser.add_argument("--dpi", type=int, default=160, help="Frame DPI.")
    parser.add_argument("--fig-width", type=float, default=16.0)
    parser.add_argument("--fig-height", type=float, default=9.6)
    parser.add_argument("--start-time-s", type=float, default=0.0, help="First model time to render.")
    parser.add_argument("--trajectory-start-time-s", type=float, default=None, help="First model time shown for target trajectories; defaults to --start-time-s.")
    parser.add_argument("--periodic-break-fraction", type=float, default=0.5, help="Break trajectory lines when modulo x/y jumps exceed this fraction of the domain.")
    parser.add_argument("--write-frame-pdfs", action="store_true", help="Write every rendered frame as a PDF.")
    parser.add_argument("--frame-pdf-dir", type=Path, default=None, help="Directory for per-frame PDF output.")
    parser.add_argument("--frames-only", action="store_true", help="Write per-frame outputs and summary without writing an MP4.")
    parser.add_argument("--times-min", nargs="+", type=float, default=None, help="Render only the nearest frames to these model times in minutes.")
    parser.add_argument("--max-frames", type=int, default=None, help="Debugging limit on rendered frames.")
    parser.add_argument("--output-stem", default="tpht_predecessor_final10_07_3d_0-4200s")
    parser.add_argument("--dry-run", action="store_true", help="Print planned inputs/outputs without rendering.")
    return parser.parse_args()


def configure_fonts() -> str:
    """Use the same robust font approach as the UAEREP Matplotlib renderer."""
    if matplotlib is None or plt is None:
        raise RuntimeError("Matplotlib is required for rendering")
    arial = Path("/Library/Fonts/Microsoft/Arial.ttf")
    if arial.exists():
        matplotlib.font_manager.fontManager.addfont(str(arial))
        plt.rcParams["font.family"] = "Arial"
        return str(arial)
    plt.rcParams["font.family"] = "DejaVu Sans"
    return "DejaVu Sans fallback"


def safe_float(value: Any) -> float | None:
    """Convert a CSV/NetCDF scalar to finite float or None."""
    if value in (None, "", "NA"):
        return None
    try:
        output = float(value)
    except (TypeError, ValueError):
        return None
    return output if math.isfinite(output) else None


def open_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV table."""
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def default_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """Resolve analysis, GMD-candidate, and export paths."""
    root = args.root.resolve()
    analysis_outdir = (args.analysis_outdir or root / "analysis_outputs").resolve()
    gmd_outdir = (args.gmd_outdir or root / "analysis_outputs_for_GMD").resolve()
    export_file = (args.export_file or gmd_outdir / "video_sources" / "tpht_3d" / EXPORT_FILE).resolve()
    return analysis_outdir, gmd_outdir, export_file


def selected_table_candidates(gmd_outdir: Path) -> list[Path]:
    """Return likely selected-target table locations."""
    return [
        gmd_outdir / "tables" / "diagnostics" / SELECTED_TABLE,
        gmd_outdir / "tables" / "supplement_candidates" / SELECTED_TABLE,
        gmd_outdir / "figure_sources" / "supplement_candidate_source_tables" / SELECTED_TABLE,
    ]


def trajectory_table_candidates(analysis_outdir: Path, gmd_outdir: Path) -> list[Path]:
    """Return likely target-trajectory table locations."""
    return [
        analysis_outdir / "tables" / TRAJECTORY_TABLE,
        gmd_outdir / "figure_sources" / "supplement_candidate_source_tables" / TRAJECTORY_TABLE,
    ]


def first_existing(paths: list[Path], label: str) -> Path:
    """Return the first existing path or raise a clear error."""
    for path in paths:
        if path.exists():
            return path
    joined = "\n  ".join(str(path) for path in paths)
    raise FileNotFoundError(f"No {label} found. Checked:\n  {joined}")


def target_ids_from_selected_table(path: Path, figure_index: int, targets_per_figure: int) -> list[str]:
    """Read the target IDs corresponding to one predecessor-tree figure index."""
    rows = [row for row in open_csv(path) if row.get("target_id") not in (None, "", "NA")]
    rows.sort(key=lambda row: safe_float(row.get("rank")) or math.inf)
    start = max(0, (figure_index - 1) * targets_per_figure)
    selected = rows[start : start + targets_per_figure]
    if len(selected) < targets_per_figure:
        raise RuntimeError(f"Selected-target table has only {len(rows)} rows; cannot read figure index {figure_index}")
    return [row["target_id"] for row in selected]


def read_trajectories(path: Path, target_ids: list[str], domain_x_m: float, domain_y_m: float) -> list[TargetTrajectory]:
    """Read selected TPHT trajectories for the full 70 min interval."""
    wanted = set(target_ids)
    grouped: dict[str, list[dict[str, str]]] = {target_id: [] for target_id in target_ids}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            target_id = row.get("target_id")
            if target_id in wanted:
                grouped[target_id].append(row)
    trajectories: list[TargetTrajectory] = []
    for index, target_id in enumerate(target_ids):
        rows = sorted(grouped.get(target_id, []), key=lambda row: safe_float(row.get("time_s")) or -math.inf)
        if not rows:
            raise RuntimeError(f"No trajectory rows found for target {target_id}")
        time_s = np.asarray([safe_float(row.get("time_s")) for row in rows], dtype=float)
        x_m = np.asarray([safe_float(row.get("x_m")) for row in rows], dtype=float)
        y_m = np.asarray([safe_float(row.get("y_m")) for row in rows], dtype=float)
        z_m = np.asarray([safe_float(row.get("z_m")) for row in rows], dtype=float)
        radius_um = np.asarray([safe_float(row.get("radius_um")) for row in rows], dtype=float)
        x_m = np.mod(x_m, domain_x_m) if domain_x_m > 0 else x_m
        y_m = np.mod(y_m, domain_y_m) if domain_y_m > 0 else y_m
        finite = np.isfinite(time_s) & np.isfinite(x_m) & np.isfinite(y_m) & np.isfinite(z_m) & np.isfinite(radius_um)
        if not np.count_nonzero(finite):
            raise RuntimeError(f"Trajectory rows for {target_id} contain no finite x/y/z/radius values")
        trajectories.append(
            TargetTrajectory(
                target_id=target_id,
                marker=MARKERS[index % len(MARKERS)],
                line_color=LINE_COLORS[index % len(LINE_COLORS)],
                time_s=time_s[finite],
                x_m=x_m[finite],
                y_m=y_m[finite],
                z_m=z_m[finite],
                radius_um=radius_um[finite],
            )
        )
    return trajectories


def deterministic_subset(indices: np.ndarray, max_points: int) -> np.ndarray:
    """Return a deterministic pseudo-random subset of flat point indices."""
    if indices.size <= max_points:
        return indices
    rng = np.random.default_rng(20260612)
    selected = rng.choice(indices.size, size=max_points, replace=False)
    return np.sort(indices[selected])


def threshold_points(values_gkg: np.ndarray, threshold_gkg: float, max_points: int, x_km: np.ndarray, y_km: np.ndarray, z_km: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Return thresholded point cloud coordinates and values."""
    flat = values_gkg.ravel()
    indices = np.flatnonzero(np.isfinite(flat) & (flat >= threshold_gkg))
    full_count = int(indices.size)
    indices = deterministic_subset(indices, max_points)
    if indices.size == 0:
        return np.empty(0), np.empty(0), np.empty(0), np.empty(0), full_count
    k, j, i = np.unravel_index(indices, values_gkg.shape)
    return x_km[i], y_km[j], z_km[k], flat[indices], full_count


def sample_threshold_values(values_gkg: np.ndarray, threshold_gkg: float, max_samples: int) -> np.ndarray:
    """Sample thresholded values for robust color scales."""
    flat = values_gkg.ravel()
    indices = np.flatnonzero(np.isfinite(flat) & (flat >= threshold_gkg))
    indices = deterministic_subset(indices, max_samples)
    if indices.size == 0:
        return np.empty(0, dtype=np.float32)
    return np.asarray(flat[indices], dtype=np.float32)


def robust_color_limit(samples: list[np.ndarray], threshold_gkg: float, actual_max_gkg: float, percentile: float) -> float:
    """Return a robust upper colorbar limit."""
    if actual_max_gkg <= threshold_gkg:
        return threshold_gkg * 1.01
    valid = [sample[np.isfinite(sample) & (sample >= threshold_gkg)] for sample in samples if sample.size]
    if not valid:
        return actual_max_gkg
    merged = np.concatenate(valid)
    if merged.size == 0:
        return actual_max_gkg
    robust = float(np.nanpercentile(merged, percentile))
    robust = min(actual_max_gkg, robust)
    return max(threshold_gkg * 1.01, robust)


def smooth_visual_field(values_gkg: np.ndarray, sigma_cells: float, label: str) -> np.ndarray:
    """Apply optional Gaussian smoothing for visualization only."""
    if sigma_cells <= 0.0:
        return values_gkg
    if gaussian_filter is None:
        raise RuntimeError(f"scipy is required when --smooth-{label}-sigma-cells is greater than zero")
    smoothed = gaussian_filter(values_gkg, sigma=float(sigma_cells), mode="nearest")
    return np.asarray(smoothed, dtype=np.float32)


def light_greys_cmap() -> LinearSegmentedColormap:
    """Return a light greyscale colormap that keeps QC visually secondary."""
    if plt is None or LinearSegmentedColormap is None:
        raise RuntimeError("Matplotlib is required for the QC colormap")
    base = plt.get_cmap("Greys")
    return LinearSegmentedColormap.from_list("GMD2026LightGreys", base(np.linspace(0.04, 0.52, 256)))


def scan_export(path: Path, trajectories: list[TargetTrajectory], args: argparse.Namespace) -> ScanStats:
    """Scan the export for robust QC/QR and radius color limits."""
    qc_samples: list[np.ndarray] = []
    qr_samples: list[np.ndarray] = []
    max_qc = 0.0
    max_qr = 0.0
    with Dataset(path, "r") as dataset:
        times = np.asarray(dataset.variables["time"][:], dtype=float)
        frame_indices = frame_indices_for_args(times, args)
        for time_index in frame_indices:
            qc_gkg = np.asarray(dataset.variables["QC_hyd"][time_index], dtype=np.float32) * 1000.0
            qr_gkg = np.asarray(dataset.variables["QR_hyd"][time_index], dtype=np.float32) * 1000.0
            qc_plot_gkg = smooth_visual_field(qc_gkg, args.smooth_qc_sigma_cells, "qc")
            qr_plot_gkg = smooth_visual_field(qr_gkg, args.smooth_qr_sigma_cells, "qr")
            if np.isfinite(qc_plot_gkg).any():
                max_qc = max(max_qc, float(np.nanmax(qc_plot_gkg)))
            if np.isfinite(qr_plot_gkg).any():
                max_qr = max(max_qr, float(np.nanmax(qr_plot_gkg)))
            qc_samples.append(sample_threshold_values(qc_plot_gkg, args.qc_threshold_gkg, 50_000))
            qr_samples.append(sample_threshold_values(qr_plot_gkg, args.qr_threshold_gkg, 50_000))
    trajectory_start = effective_trajectory_start_time(args)
    radius_values = []
    for item in trajectories:
        mask = np.isfinite(item.radius_um) & (item.time_s >= trajectory_start - 1.0e-9)
        if np.count_nonzero(mask):
            radius_values.append(item.radius_um[mask])
    if not radius_values:
        radius_values = [item.radius_um[np.isfinite(item.radius_um)] for item in trajectories]
    radii = np.concatenate(radius_values)
    radius_vmin = max(0.0, float(np.nanmin(radii))) if radii.size else 0.0
    radius_vmax = float(np.nanmax(radii)) if radii.size else 1.0
    if radius_vmax <= radius_vmin:
        radius_vmax = radius_vmin + 1.0
    return ScanStats(
        qc_vmax_gkg=robust_color_limit(qc_samples, args.qc_threshold_gkg, max_qc, args.color_percentile),
        qr_vmax_gkg=robust_color_limit(qr_samples, args.qr_threshold_gkg, max_qr, args.color_percentile),
        radius_vmin_um=radius_vmin,
        radius_vmax_um=radius_vmax,
    )


def format_log_ticks(max_value: float, threshold: float) -> list[float]:
    """Format compact log colorbar ticks."""
    ticks = [threshold]
    for value in (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
        if threshold < value < max_value:
            ticks.append(value)
    if max_value > threshold:
        ticks.append(max_value)
    return ticks


def draw_reference_box(axis: Any, xlim: tuple[float, float], ylim: tuple[float, float], zlim: tuple[float, float]) -> None:
    """Draw a simple 3-D domain box."""
    xmin, xmax = xlim
    ymin, ymax = ylim
    zmin, zmax = zlim
    edge_color = "0.62"
    for z in (zmin, zmax):
        axis.plot([xmin, xmax], [ymin, ymin], [z, z], color=edge_color, lw=0.8)
        axis.plot([xmin, xmax], [ymax, ymax], [z, z], color=edge_color, lw=0.8)
        axis.plot([xmin, xmin], [ymin, ymax], [z, z], color=edge_color, lw=0.8)
        axis.plot([xmax, xmax], [ymin, ymax], [z, z], color=edge_color, lw=0.8)
    for x in (xmin, xmax):
        for y in (ymin, ymax):
            axis.plot([x, x], [y, y], [zmin, zmax], color=edge_color, lw=0.8)


def effective_trajectory_start_time(args: argparse.Namespace) -> float:
    """Return the first time used for trajectory display."""
    if args.trajectory_start_time_s is None:
        return float(args.start_time_s)
    return float(args.trajectory_start_time_s)


def frame_indices_for_args(times: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    """Return frame indices after applying start-time and debug limits."""
    times_min = getattr(args, "times_min", None)
    if times_min:
        requested = [int(np.argmin(np.abs(times - float(time_min) * 60.0))) for time_min in times_min]
        frame_indices = np.asarray(list(dict.fromkeys(requested)), dtype=int)
    else:
        frame_indices = np.flatnonzero(times >= float(args.start_time_s) - 1.0e-9)
    if args.max_frames is not None:
        frame_indices = frame_indices[: max(0, args.max_frames)]
    if frame_indices.size == 0:
        raise RuntimeError(f"No export frames found for the requested time selection")
    return frame_indices.astype(int)


def current_trajectory_slice(trajectory: TargetTrajectory, time_s: float, start_time_s: float) -> np.ndarray:
    """Return indices up to the current frame time."""
    if time_s < start_time_s - 1.0e-9:
        return np.empty(0, dtype=int)
    indices = np.flatnonzero((trajectory.time_s >= start_time_s - 1.0e-9) & (trajectory.time_s <= time_s + 1.0e-9))
    if indices.size == 0:
        return np.empty(0, dtype=int)
    return indices


def plot_periodic_segments(
    axis: Any,
    trajectory: TargetTrajectory,
    indices: np.ndarray,
    domain_x_m: float,
    domain_y_m: float,
    break_fraction: float,
) -> None:
    """Plot trajectory line segments without connecting across periodic jumps."""
    if indices.size < 2:
        return
    x_m = trajectory.x_m[indices]
    y_m = trajectory.y_m[indices]
    z_m = trajectory.z_m[indices]
    x_jump = np.abs(np.diff(x_m)) > max(0.0, break_fraction) * domain_x_m
    y_jump = np.abs(np.diff(y_m)) > max(0.0, break_fraction) * domain_y_m
    break_after = np.flatnonzero(x_jump | y_jump)
    segment_starts = np.concatenate(([0], break_after + 1))
    segment_ends = np.concatenate((break_after + 1, [indices.size]))
    for start, end in zip(segment_starts, segment_ends):
        if end - start < 2:
            continue
        axis.plot(
            x_m[start:end] / 1000.0,
            y_m[start:end] / 1000.0,
            z_m[start:end] / 1000.0,
            color=trajectory.line_color,
            lw=1.0,
            alpha=0.50,
        )


def render_frame(
    path: Path,
    time_index: int,
    trajectories: list[TargetTrajectory],
    qc_norm: LogNorm,
    qr_norm: LogNorm,
    radius_norm: Normalize,
    stats: ScanStats,
    args: argparse.Namespace,
    domain_x_m: float,
    domain_y_m: float,
    pdf_path: Path | None = None,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Render one video frame."""
    with Dataset(path, "r") as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        x_km = np.asarray(dataset.variables["x"][:], dtype=float) / 1000.0
        y_km = np.asarray(dataset.variables["y"][:], dtype=float) / 1000.0
        z_km = np.asarray(dataset.variables["z"][:], dtype=float) / 1000.0
        qc_gkg = np.asarray(dataset.variables["QC_hyd"][time_index], dtype=np.float32) * 1000.0
        qr_gkg = np.asarray(dataset.variables["QR_hyd"][time_index], dtype=np.float32) * 1000.0

    qc_plot_gkg = smooth_visual_field(qc_gkg, args.smooth_qc_sigma_cells, "qc")
    qr_plot_gkg = smooth_visual_field(qr_gkg, args.smooth_qr_sigma_cells, "qr")
    qcx, qcy, qcz, qcv, qc_count = threshold_points(qc_plot_gkg, args.qc_threshold_gkg, args.max_qc_points, x_km, y_km, z_km)
    qrx, qry, qrz, qrv, qr_count = threshold_points(qr_plot_gkg, args.qr_threshold_gkg, args.max_qr_points, x_km, y_km, z_km)

    fig = plt.figure(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    axis = fig.add_axes([0.14, 0.16, 0.72, 0.74], projection="3d")
    axis.view_init(elev=18.0, azim=-58.0)
    xlim = (float(np.nanmin(x_km)), float(np.nanmax(x_km)))
    ylim = (float(np.nanmin(y_km)), float(np.nanmax(y_km)))
    zlim = (0.0, float(np.nanmax(z_km)))
    axis.set_xlim(*xlim)
    axis.set_ylim(*ylim)
    axis.set_zlim(*zlim)
    axis.set_box_aspect((max(1.0e-9, xlim[1] - xlim[0]), max(1.0e-9, ylim[1] - ylim[0]), max(1.0e-9, zlim[1] - zlim[0])))
    axis.set_xlabel("x (km)", labelpad=8)
    axis.set_ylabel("y (km)", labelpad=8)
    axis.set_zlabel("z AGL (km)", labelpad=8)
    axis.grid(True, color="0.85", linewidth=0.6)
    axis.tick_params(axis="both", which="major", labelsize=8)
    draw_reference_box(axis, xlim, ylim, zlim)

    qc_cmap = light_greys_cmap()
    if qcx.size:
        axis.scatter(qcx, qcy, qcz, c=qcv, cmap=qc_cmap, norm=qc_norm, s=args.qc_marker_size, marker="s", alpha=args.qc_alpha, linewidths=0, rasterized=True)
    if qrx.size:
        axis.scatter(qrx, qry, qrz, c=qrv, cmap="YlOrRd", norm=qr_norm, s=args.qr_marker_size, marker="s", alpha=args.qr_alpha, linewidths=0, rasterized=True)

    radius_cmap = plt.get_cmap("cool")
    legend_handles = []
    legend_labels = []
    frame_time = float(time[time_index])
    trajectory_start = effective_trajectory_start_time(args)
    for trajectory_index, trajectory in enumerate(trajectories):
        target_label = f"Target {chr(ord('a') + trajectory_index)}"
        indices = current_trajectory_slice(trajectory, frame_time, trajectory_start)
        if indices.size == 0:
            handle = plt.Line2D([0], [0], marker=trajectory.marker, color=trajectory.line_color, linestyle="", markerfacecolor=trajectory.line_color, markeredgewidth=0, markersize=8)
            legend_handles.append(handle)
            legend_labels.append(target_label)
            continue
        plot_periodic_segments(axis, trajectory, indices, domain_x_m, domain_y_m, args.periodic_break_fraction)
        axis.scatter(
            trajectory.x_m[indices] / 1000.0,
            trajectory.y_m[indices] / 1000.0,
            trajectory.z_m[indices] / 1000.0,
            c=trajectory.radius_um[indices],
            cmap=radius_cmap,
            norm=radius_norm,
            s=20,
            marker=trajectory.marker,
            alpha=0.72,
            linewidths=0,
            depthshade=False,
        )
        current = indices[-1]
        axis.scatter(
            [trajectory.x_m[current] / 1000.0],
            [trajectory.y_m[current] / 1000.0],
            [trajectory.z_m[current] / 1000.0],
            c=[trajectory.radius_um[current]],
            cmap=radius_cmap,
            norm=radius_norm,
            s=90,
            marker=trajectory.marker,
            linewidths=0,
            depthshade=False,
        )
        handle = plt.Line2D([0], [0], marker=trajectory.marker, color=trajectory.line_color, linestyle="", markerfacecolor=trajectory.line_color, markeredgewidth=0, markersize=8)
        legend_handles.append(handle)
        legend_labels.append(target_label)

    axis.legend(legend_handles, legend_labels, frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0.0, 1.03), ncol=1)
    fig.suptitle(f"TPHT backward-reconstruction diagnostics, t = {frame_time / 60.0:.1f} min", fontsize=18, y=0.965)
    if trajectory_start > float(np.nanmin(time)) + 1.0e-9:
        trajectory_note = f"target histories traced from {trajectory_start / 60.0:.1f} min"
    else:
        trajectory_note = "target histories traced from 0.0 min"
    fig.text(0.5, 0.050, f"QC/QR Eulerian fields with TPHT selected-target trajectories; marker colour = target radius; {trajectory_note}", ha="center", va="center", fontsize=12)

    cax_qc = fig.add_axes([0.105, 0.31, 0.018, 0.40])
    cbar_qc = fig.colorbar(ScalarMappable(norm=qc_norm, cmap=qc_cmap), cax=cax_qc)
    cbar_qc.set_label(r"$Q_C$ (g kg$^{-1}$)", fontsize=13, labelpad=14)
    cbar_qc.ax.yaxis.set_label_position("left")
    cbar_qc.set_ticks(format_log_ticks(stats.qc_vmax_gkg, args.qc_threshold_gkg))
    cbar_qc.ax.set_yticklabels([f"{tick:g}" for tick in cbar_qc.get_ticks()])
    cbar_qc.ax.tick_params(labelsize=9)

    cax_qr = fig.add_axes([0.875, 0.31, 0.018, 0.40])
    cbar_qr = fig.colorbar(ScalarMappable(norm=qr_norm, cmap="YlOrRd"), cax=cax_qr)
    cbar_qr.set_label(r"$Q_R$ (g kg$^{-1}$)", fontsize=13, labelpad=14)
    cbar_qr.set_ticks(format_log_ticks(stats.qr_vmax_gkg, args.qr_threshold_gkg))
    cbar_qr.ax.set_yticklabels([f"{tick:g}" for tick in cbar_qr.get_ticks()])
    cbar_qr.ax.tick_params(labelsize=9)

    cax_radius = fig.add_axes([0.34, 0.082, 0.32, 0.018])
    cbar_radius = fig.colorbar(ScalarMappable(norm=radius_norm, cmap=radius_cmap), cax=cax_radius, orientation="horizontal")
    cbar_radius.set_label("Target radius (µm)", fontsize=12, labelpad=4)
    cbar_radius.ax.xaxis.set_label_position("top")
    cbar_radius.ax.tick_params(labelsize=9)

    if pdf_path is not None:
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(pdf_path, dpi=args.dpi, bbox_inches="tight")
    fig.canvas.draw()
    image = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
    plt.close(fig)
    return image, {
        "time_s": frame_time,
        "max_QC_gkg": float(np.nanmax(qc_plot_gkg)) if np.isfinite(qc_plot_gkg).any() else math.nan,
        "max_QR_gkg": float(np.nanmax(qr_plot_gkg)) if np.isfinite(qr_plot_gkg).any() else math.nan,
        "QC_points": qc_count,
        "QR_points": qr_count,
    }


def write_summary(output_dir: Path, video_path: Path | None, key_frame_dir: Path, source_export: Path, target_ids: list[str], rows: list[dict[str, float | int]], args: argparse.Namespace, font_source: str, frame_pdf_dir: Path | None) -> None:
    """Write a compact Markdown summary for traceability."""
    lines = [
        "# TPHT Backward 3-D Trajectory Video",
        "",
        f"- video: `{video_path}`" if video_path is not None else "- video: not written (`--frames-only`)",
        f"- key frames: `{key_frame_dir}`",
        f"- per-frame PDFs: `{frame_pdf_dir}`" if frame_pdf_dir is not None else "- per-frame PDFs: not requested",
        f"- source export: `{source_export}`",
        f"- target IDs: {', '.join(f'`{target_id}`' for target_id in target_ids)}",
        f"- selected predecessor-tree source: final10 figure index {args.figure_index}, {args.targets_per_figure} targets per figure",
        f"- rendered start time: {args.start_time_s:g} s ({args.start_time_s / 60.0:g} min)",
        f"- trajectory start time: {effective_trajectory_start_time(args):g} s ({effective_trajectory_start_time(args) / 60.0:g} min)",
        f"- periodic boundary line break fraction: {args.periodic_break_fraction:g}",
        f"- fps: {args.fps}",
        f"- QC threshold: {args.qc_threshold_gkg:g} g/kg",
        f"- QR threshold: {args.qr_threshold_gkg:g} g/kg",
        f"- QC smoothing: Gaussian sigma {args.smooth_qc_sigma_cells:g} grid cells for visualization only",
        f"- QR smoothing: Gaussian sigma {args.smooth_qr_sigma_cells:g} grid cells for visualization only",
        f"- QC alpha: {args.qc_alpha:g}",
        f"- QR alpha: {args.qr_alpha:g}",
        f"- QC square marker size: {args.qc_marker_size:g}",
        f"- QR square marker size: {args.qr_marker_size:g}",
        f"- requested model times: {', '.join(f'{time_min:g} min' for time_min in args.times_min)}" if args.times_min else "- requested model times: all selected frames",
        f"- particle colormap: `cool`; marker colour represents target radius in µm",
        f"- font family: `{plt.rcParams['font.family'][0]}` from `{font_source}`",
        "",
        "## Per-frame Diagnostics",
        "",
        "| time_s | time_min | max_QC_gkg | max_QR_gkg | QC points | QR points |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        time_s = float(row["time_s"])
        lines.append(
            f"| {time_s:.0f} | {time_s / 60.0:.2f} | {float(row['max_QC_gkg']):.6e} | "
            f"{float(row['max_QR_gkg']):.6e} | {int(row['QC_points'])} | {int(row['QR_points'])} |"
        )
    (output_dir / "tpht_predecessor_3d_video_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_video(args: argparse.Namespace) -> Path:
    """Render the complete 70 min TPHT 3-D diagnostic video."""
    analysis_outdir, gmd_outdir, export_file = default_paths(args)
    trajectory_table = first_existing(trajectory_table_candidates(analysis_outdir, gmd_outdir), "TPHT trajectory table")
    selected_table: Path | None = None
    if args.target_ids:
        target_ids = list(args.target_ids)
    else:
        selected_table = first_existing(selected_table_candidates(gmd_outdir), "final-window selected-target table")
        target_ids = target_ids_from_selected_table(selected_table, args.figure_index, args.targets_per_figure)
    output_dir = gmd_outdir / "videos" / "tpht_3d"
    video_path = None if args.frames_only else output_dir / f"{args.output_stem}.mp4"
    frame_pdf_dir = args.frame_pdf_dir or output_dir / "pdf_frames" / args.output_stem
    frame_pdf_dir = frame_pdf_dir.resolve() if args.write_frame_pdfs else None
    if args.dry_run:
        print("[dry-run] GMD2026 TPHT 3-D video render")
        print(f"[dry-run] export_file={export_file}")
        print(f"[dry-run] selected_table={selected_table if selected_table is not None else 'not required when --target-ids is set'}")
        print(f"[dry-run] trajectory_table={trajectory_table}")
        print(f"[dry-run] target_ids={','.join(target_ids)}")
        print(f"[dry-run] start_time_s={args.start_time_s:g}")
        print(f"[dry-run] trajectory_start_time_s={effective_trajectory_start_time(args):g}")
        print(f"[dry-run] times_min={','.join(str(item) for item in args.times_min) if args.times_min else 'all'}")
        print(f"[dry-run] output_video={video_path if video_path is not None else 'not written (--frames-only)'}")
        print(f"[dry-run] frame_pdf_dir={frame_pdf_dir if frame_pdf_dir is not None else 'not requested'}")
        return video_path or output_dir

    if Dataset is None:
        raise RuntimeError("netCDF4 is required for rendering. Activate the GMD2026 sdm_env environment on SQUID.")
    if imageio is None and not args.frames_only:
        raise RuntimeError("imageio is required for rendering. Activate the GMD2026 sdm_env environment on SQUID.")
    if matplotlib is None or plt is None or LogNorm is None or Normalize is None:
        raise RuntimeError("Matplotlib is required for rendering. Activate the GMD2026 sdm_env environment on SQUID.")

    with Dataset(export_file, "r") as dataset:
        x = np.asarray(dataset.variables["x"][:], dtype=float)
        y = np.asarray(dataset.variables["y"][:], dtype=float)
        times = np.asarray(dataset.variables["time"][:], dtype=float)
        domain_x_m = float(np.nanmax(x) + 0.5 * getattr(dataset, "dx", 0.0))
        domain_y_m = float(np.nanmax(y) + 0.5 * getattr(dataset, "dy", 0.0))
    trajectories = read_trajectories(trajectory_table, target_ids, domain_x_m, domain_y_m)

    key_frame_dir = output_dir / "key_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    key_frame_dir.mkdir(parents=True, exist_ok=True)
    font_source = configure_fonts()
    stats = scan_export(export_file, trajectories, args)
    qc_norm = LogNorm(vmin=args.qc_threshold_gkg, vmax=stats.qc_vmax_gkg)
    qr_norm = LogNorm(vmin=args.qr_threshold_gkg, vmax=stats.qr_vmax_gkg)
    radius_norm = Normalize(vmin=stats.radius_vmin_um, vmax=stats.radius_vmax_um)

    frame_indices = frame_indices_for_args(times, args)
    rows: list[dict[str, float | int]] = []
    writer_context = imageio.get_writer(video_path, fps=args.fps, codec="libx264", quality=8, macro_block_size=16) if video_path is not None else None
    if writer_context is not None:
        writer_context.__enter__()
    try:
        for count, time_index in enumerate(frame_indices, start=1):
            frame_time = float(times[int(time_index)])
            pdf_path = frame_pdf_dir / f"{args.output_stem}_t{int(round(frame_time)):05d}s.pdf" if frame_pdf_dir is not None else None
            image, row = render_frame(export_file, int(time_index), trajectories, qc_norm, qr_norm, radius_norm, stats, args, domain_x_m, domain_y_m, pdf_path)
            if writer_context is not None:
                writer_context.append_data(image)
            rows.append(row)
            if count == 1 or count % 50 == 0 or count == len(frame_indices):
                print(f"[render_3d] rendered {count}/{len(frame_indices)} frames", flush=True)
    finally:
        if writer_context is not None:
            writer_context.__exit__(None, None, None)

    key_times = [float(times[frame_indices[0]]), float(times[frame_indices[len(frame_indices) // 2]]), float(times[frame_indices[-1]])]
    for key_time in key_times:
        nearest_index = int(np.argmin(np.abs(times - key_time)))
        image, _ = render_frame(export_file, nearest_index, trajectories, qc_norm, qr_norm, radius_norm, stats, args, domain_x_m, domain_y_m)
        imageio.imwrite(key_frame_dir / f"{args.output_stem}_t{int(round(times[nearest_index])):05d}s.png", image)

    write_summary(output_dir, video_path, key_frame_dir, export_file, target_ids, rows, args, font_source, frame_pdf_dir)
    if video_path is not None:
        print(f"[render_3d] wrote {video_path}", flush=True)
    if frame_pdf_dir is not None:
        print(f"[render_3d] wrote per-frame PDFs in {frame_pdf_dir}", flush=True)
    return video_path or output_dir


def main() -> None:
    """CLI entry point."""
    render_video(parse_args())


if __name__ == "__main__":
    main()
