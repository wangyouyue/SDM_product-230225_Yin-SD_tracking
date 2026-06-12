#!/usr/bin/env python3
"""Render a few filled-voxel TPHT QC/QR preview frames.

This script is intentionally for visual inspection only.  It renders each
thresholded Eulerian grid cell as a filled 3-D cuboid so that the user can
compare this style with the lighter cell-center marker video renderer.
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np

from render_02_tpht_predecessor_3d_video import (
    EXPORT_FILE,
    LINE_COLORS,
    MARKERS,
    ScanStats,
    configure_fonts,
    default_paths,
    effective_trajectory_start_time,
    first_existing,
    format_log_ticks,
    frame_indices_for_args,
    light_greys_cmap,
    plot_periodic_segments,
    read_trajectories,
    scan_export,
    smooth_visual_field,
    trajectory_table_candidates,
)

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import LogNorm, Normalize
except Exception:  # pragma: no cover
    matplotlib = None
    plt = None
    ScalarMappable = None
    LogNorm = None
    Normalize = None


def parse_args() -> argparse.Namespace:
    """Return command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Path to the GMD2026 directory.")
    parser.add_argument("--analysis-outdir", type=Path, default=None, help="Existing analysis_outputs directory with TPHT science tables.")
    parser.add_argument("--gmd-outdir", type=Path, default=None, help="analysis_outputs_for_GMD directory.")
    parser.add_argument("--export-file", type=Path, default=None, help="Compact QC/QR export NetCDF.")
    parser.add_argument("--target-ids", nargs="*", default=["chain_02008800", "chain_00803600", "chain_03182000"], help="Target IDs to overlay.")
    parser.add_argument("--times-min", nargs="+", type=float, default=[60.0, 65.0, 70.0], help="Model times in minutes to render.")
    parser.add_argument("--trajectory-start-time-s", type=float, default=3600.0, help="First model time shown for target trajectories.")
    parser.add_argument("--periodic-break-fraction", type=float, default=0.5, help="Break trajectory lines when modulo x/y jumps exceed this fraction of the domain.")
    parser.add_argument("--qc-threshold-gkg", type=float, default=0.001, help="Visual QC threshold in g/kg.")
    parser.add_argument("--qr-threshold-gkg", type=float, default=0.0001, help="Visual QR threshold in g/kg.")
    parser.add_argument("--smooth-qc-sigma-cells", type=float, default=1.25, help="Gaussian smoothing sigma for QC visualization in grid-cell units.")
    parser.add_argument("--smooth-qr-sigma-cells", type=float, default=0.0, help="Gaussian smoothing sigma for QR visualization in grid-cell units.")
    parser.add_argument("--qc-alpha", type=float, default=0.18, help="QC cuboid alpha.")
    parser.add_argument("--qr-alpha", type=float, default=0.70, help="QR cuboid alpha.")
    parser.add_argument("--color-percentile", type=float, default=99.5, help="Robust QC/QR color upper percentile.")
    parser.add_argument("--fig-width", type=float, default=16.0)
    parser.add_argument("--fig-height", type=float, default=9.6)
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory for preview PNG/PDF files.")
    parser.add_argument("--output-stem", default="tpht_predecessor_final10_07_3d_voxel_preview")
    parser.add_argument("--dry-run", action="store_true", help="Print planned output paths without rendering.")
    return parser.parse_args()


def axis_edges(centers: np.ndarray, fallback_spacing: float) -> np.ndarray:
    """Return grid-cell edges from cell centers."""
    centers = np.asarray(centers, dtype=float)
    if centers.size == 1:
        half = 0.5 * fallback_spacing
        return np.asarray([centers[0] - half, centers[0] + half], dtype=float)
    middle = 0.5 * (centers[:-1] + centers[1:])
    first = centers[0] - (middle[0] - centers[0])
    last = centers[-1] + (centers[-1] - middle[-1])
    return np.concatenate(([first], middle, [last]))


def rgba_voxel_colors(values: np.ndarray, mask: np.ndarray, norm: LogNorm, cmap: object, alpha: float) -> np.ndarray:
    """Create a float RGBA color array for Matplotlib voxels."""
    colors = np.zeros(mask.shape + (4,), dtype=float)
    colors[mask] = cmap(norm(values[mask]))
    colors[..., 3] = np.where(mask, alpha, 0.0)
    return colors


def current_indices(time_s: np.ndarray, frame_time: float, start_time_s: float) -> np.ndarray:
    """Return trajectory rows from start_time_s through frame_time."""
    indices = np.flatnonzero((time_s >= start_time_s - 1.0e-9) & (time_s <= frame_time + 1.0e-9))
    if indices.size:
        return indices
    future = np.flatnonzero(time_s >= start_time_s - 1.0e-9)
    candidates = future if future.size else np.arange(time_s.size)
    nearest = int(candidates[np.argmin(np.abs(time_s[candidates] - frame_time))])
    return np.asarray([nearest], dtype=int)


def draw_voxel_frame(
    export_file: Path,
    time_index: int,
    trajectories: list[object],
    stats: ScanStats,
    args: argparse.Namespace,
    output_dir: Path,
    domain_x_m: float,
    domain_y_m: float,
) -> Path:
    """Render one filled-voxel preview frame."""
    with Dataset(export_file, "r") as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        x_km = np.asarray(dataset.variables["x"][:], dtype=float) / 1000.0
        y_km = np.asarray(dataset.variables["y"][:], dtype=float) / 1000.0
        z_km = np.asarray(dataset.variables["z"][:], dtype=float) / 1000.0
        qc_gkg = np.asarray(dataset.variables["QC_hyd"][time_index], dtype=np.float32) * 1000.0
        qr_gkg = np.asarray(dataset.variables["QR_hyd"][time_index], dtype=np.float32) * 1000.0

    qc_plot_gkg = smooth_visual_field(qc_gkg, args.smooth_qc_sigma_cells, "qc")
    qr_plot_gkg = smooth_visual_field(qr_gkg, args.smooth_qr_sigma_cells, "qr")
    qc_mask_zyx = np.isfinite(qc_plot_gkg) & (qc_plot_gkg >= args.qc_threshold_gkg)
    qr_mask_zyx = np.isfinite(qr_plot_gkg) & (qr_plot_gkg >= args.qr_threshold_gkg)

    # Axes3D.voxels expects arrays ordered as x, y, z.
    qc_values_xyz = np.transpose(qc_plot_gkg, (2, 1, 0))
    qr_values_xyz = np.transpose(qr_plot_gkg, (2, 1, 0))
    qc_mask_xyz = np.transpose(qc_mask_zyx, (2, 1, 0))
    qr_mask_xyz = np.transpose(qr_mask_zyx, (2, 1, 0))
    x_edges = axis_edges(x_km, 0.05)
    y_edges = axis_edges(y_km, 0.05)
    z_edges = axis_edges(z_km, 0.005)
    x_grid = x_edges[:, None, None]
    y_grid = y_edges[None, :, None]
    z_grid = z_edges[None, None, :]

    qc_cmap = light_greys_cmap()
    qr_cmap = plt.get_cmap("YlOrRd")
    radius_cmap = plt.get_cmap("cool")
    qc_norm = LogNorm(vmin=args.qc_threshold_gkg, vmax=stats.qc_vmax_gkg)
    qr_norm = LogNorm(vmin=args.qr_threshold_gkg, vmax=stats.qr_vmax_gkg)
    radius_norm = Normalize(vmin=stats.radius_vmin_um, vmax=stats.radius_vmax_um)

    fig = plt.figure(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    axis = fig.add_axes([0.14, 0.16, 0.72, 0.74], projection="3d")
    axis.view_init(elev=18.0, azim=-58.0)
    axis.set_xlim(float(x_edges[0]), float(x_edges[-1]))
    axis.set_ylim(float(y_edges[0]), float(y_edges[-1]))
    axis.set_zlim(0.0, float(z_edges[-1]))
    axis.set_box_aspect((x_edges[-1] - x_edges[0], y_edges[-1] - y_edges[0], z_edges[-1] - z_edges[0]))
    axis.set_xlabel("x (km)", labelpad=8)
    axis.set_ylabel("y (km)", labelpad=8)
    axis.set_zlabel("z AGL (km)", labelpad=8)
    axis.tick_params(axis="both", which="major", labelsize=8)
    axis.grid(True, color="0.86", linewidth=0.6)

    if np.count_nonzero(qc_mask_xyz):
        axis.voxels(
            x_grid,
            y_grid,
            z_grid,
            qc_mask_xyz,
            facecolors=rgba_voxel_colors(qc_values_xyz, qc_mask_xyz, qc_norm, qc_cmap, args.qc_alpha),
            edgecolor="none",
            shade=False,
        )
    if np.count_nonzero(qr_mask_xyz):
        axis.voxels(
            x_grid,
            y_grid,
            z_grid,
            qr_mask_xyz,
            facecolors=rgba_voxel_colors(qr_values_xyz, qr_mask_xyz, qr_norm, qr_cmap, args.qr_alpha),
            edgecolor="none",
            shade=False,
        )

    frame_time = float(time[time_index])
    trajectory_start = effective_trajectory_start_time(args)
    handles = []
    labels = []
    for trajectory_index, trajectory in enumerate(trajectories):
        indices = current_indices(trajectory.time_s, frame_time, trajectory_start)
        plot_periodic_segments(axis, trajectory, indices, domain_x_m, domain_y_m, args.periodic_break_fraction)
        axis.scatter(
            trajectory.x_m[indices] / 1000.0,
            trajectory.y_m[indices] / 1000.0,
            trajectory.z_m[indices] / 1000.0,
            c=trajectory.radius_um[indices],
            cmap=radius_cmap,
            norm=radius_norm,
            s=24,
            marker=MARKERS[trajectory_index % len(MARKERS)],
            alpha=0.82,
            linewidths=0,
            depthshade=False,
        )
        handles.append(
            plt.Line2D(
                [0],
                [0],
                marker=MARKERS[trajectory_index % len(MARKERS)],
                color=LINE_COLORS[trajectory_index % len(LINE_COLORS)],
                linestyle="",
                markerfacecolor=LINE_COLORS[trajectory_index % len(LINE_COLORS)],
                markeredgewidth=0,
                markersize=8,
            )
        )
        labels.append(trajectory.target_id)

    axis.legend(handles, labels, frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0.0, 1.03), ncol=1)
    fig.suptitle(f"GMD2026 TPHT filled-voxel preview, t = {frame_time / 60.0:.1f} min", fontsize=18, y=0.965)
    fig.text(0.5, 0.050, "Filled cuboids show thresholded QC/QR grid cells; marker color = target radius", ha="center", va="center", fontsize=12)

    cax_qc = fig.add_axes([0.105, 0.31, 0.018, 0.40])
    cbar_qc = fig.colorbar(ScalarMappable(norm=qc_norm, cmap=qc_cmap), cax=cax_qc)
    cbar_qc.set_label("QC (g/kg)", fontsize=13, labelpad=14)
    cbar_qc.ax.yaxis.set_label_position("left")
    cbar_qc.set_ticks(format_log_ticks(stats.qc_vmax_gkg, args.qc_threshold_gkg))
    cbar_qc.ax.set_yticklabels([f"{tick:g}" for tick in cbar_qc.get_ticks()])
    cbar_qc.ax.tick_params(labelsize=9)

    cax_qr = fig.add_axes([0.875, 0.31, 0.018, 0.40])
    cbar_qr = fig.colorbar(ScalarMappable(norm=qr_norm, cmap=qr_cmap), cax=cax_qr)
    cbar_qr.set_label("QR (g/kg)", fontsize=13, labelpad=14)
    cbar_qr.set_ticks(format_log_ticks(stats.qr_vmax_gkg, args.qr_threshold_gkg))
    cbar_qr.ax.set_yticklabels([f"{tick:g}" for tick in cbar_qr.get_ticks()])
    cbar_qr.ax.tick_params(labelsize=9)

    cax_radius = fig.add_axes([0.34, 0.082, 0.32, 0.018])
    cbar_radius = fig.colorbar(ScalarMappable(norm=radius_norm, cmap=radius_cmap), cax=cax_radius, orientation="horizontal")
    cbar_radius.set_label("Target radius (µm)", fontsize=12, labelpad=4)
    cbar_radius.ax.xaxis.set_label_position("top")
    cbar_radius.ax.tick_params(labelsize=9)

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / f"{args.output_stem}_t{int(round(frame_time)):05d}s.png"
    pdf_path = output_dir / f"{args.output_stem}_t{int(round(frame_time)):05d}s.pdf"
    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    return png_path


def main() -> None:
    """CLI entry point."""
    args = parse_args()
    if Dataset is None or matplotlib is None or plt is None or LogNorm is None or Normalize is None:
        raise RuntimeError("netCDF4 and Matplotlib are required for filled-voxel preview rendering")
    analysis_outdir, gmd_outdir, export_file = default_paths(args)
    output_dir = args.output_dir or gmd_outdir / "videos" / "tpht_3d" / "voxel_previews"
    output_dir = output_dir.resolve()
    trajectory_table = first_existing(trajectory_table_candidates(analysis_outdir, gmd_outdir), "TPHT trajectory table")
    if args.dry_run:
        print("[dry-run] GMD2026 TPHT filled-voxel preview")
        print(f"[dry-run] export_file={export_file}")
        print(f"[dry-run] trajectory_table={trajectory_table}")
        print(f"[dry-run] target_ids={','.join(args.target_ids)}")
        print(f"[dry-run] output_dir={output_dir}")
        print(f"[dry-run] times_min={','.join(str(item) for item in args.times_min)}")
        return

    configure_fonts()
    with Dataset(export_file, "r") as dataset:
        x_m = np.asarray(dataset.variables["x"][:], dtype=float)
        y_m = np.asarray(dataset.variables["y"][:], dtype=float)
        times_s = np.asarray(dataset.variables["time"][:], dtype=float)
        domain_x_m = float(np.nanmax(x_m) + 0.5 * getattr(dataset, "dx", 0.0))
        domain_y_m = float(np.nanmax(y_m) + 0.5 * getattr(dataset, "dy", 0.0))
    args.start_time_s = min(args.times_min) * 60.0
    args.max_frames = None
    frame_indices_for_args(times_s, args)
    trajectories = read_trajectories(trajectory_table, args.target_ids, domain_x_m, domain_y_m)
    stats = scan_export(export_file, trajectories, args)
    written = []
    for time_min in args.times_min:
        nearest = int(np.argmin(np.abs(times_s - time_min * 60.0)))
        written.append(draw_voxel_frame(export_file, nearest, trajectories, stats, args, output_dir, domain_x_m, domain_y_m))
        print(f"[voxel_preview] wrote {written[-1]}", flush=True)


if __name__ == "__main__":
    tempfile.tempdir = tempfile.gettempdir()
    main()
