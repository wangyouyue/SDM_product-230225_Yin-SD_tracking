#!/usr/bin/env python3
"""Plot identity-safe SCALE-SDM trajectories produced by ``sd_output.py``."""

import argparse
import math
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sdm_tracking_matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "sdm_tracking_xdg_cache"))
    import matplotlib

    if "--show" not in sys.argv and not os.environ.get("MPLBACKEND"):
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LogNorm, Normalize
    from netCDF4 import Dataset
except ImportError as exc:  # pragma: no cover - exercised by deployment checks
    raise SystemExit(
        "traj_plot.py requires numpy, netCDF4, and matplotlib. Activate the "
        "project Python environment before running it."
    ) from exc


def parse_args(argv=None):
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="tracking_trajectories.nc")
    parser.add_argument("--output-base", default="particle_trajectories")
    parser.add_argument("--target", action="append", default=[], metavar="DM_ID:SD_ID")
    parser.add_argument("--trajectory-index", action="append", default=[], type=int)
    parser.add_argument("--number", type=int, default=10, help="number of automatic trajectories")
    parser.add_argument("--initial-radius-max-um", type=float, default=5.0)
    parser.add_argument("--max-radius-min-um", type=float, default=15.0)
    parser.add_argument("--x-period-m", type=float, default=None, help="periodic x length; file metadata is used by default")
    parser.add_argument(
        "--full-domain",
        action="store_true",
        help="show x=0..period and z from zero instead of zooming to selected trajectories",
    )
    parser.add_argument(
        "--relative-position",
        action="store_true",
        help="plot displacement from each trajectory's first record instead of absolute position",
    )
    parser.add_argument("--z-min-m", type=float, default=None, help="optional lower height limit")
    parser.add_argument("--z-max-m", type=float, default=None, help="optional upper height limit")
    parser.add_argument("--reverse-time", action="store_true", help="reverse start/end direction for BW presentation")
    parser.add_argument(
        "--style",
        choices=("gmd", "diagnostic"),
        default="gmd",
        help="publication or larger-screen diagnostic styling",
    )
    parser.add_argument(
        "--figure-width",
        choices=("single", "double"),
        default="double",
        help="GMD single-column (85 mm) or double-column (170 mm) width",
    )
    parser.add_argument("--title", default=None, help="optional figure title; omitted by default in GMD style")
    parser.add_argument("--panel-label", default=None, help="optional panel label, for example '(a)'")
    parser.add_argument("--show", action="store_true")
    return parser.parse_args(argv)


def _parse_pair(text):
    fields = text.replace(",", ":").split(":")
    if len(fields) != 2:
        raise ValueError("invalid target pair: {0}".format(text))
    return int(fields[0]), int(fields[1])


def _read_rows(path):
    """Read tidy trajectory rows and file metadata."""
    with Dataset(path, "r") as handle:
        required = ("trajectory_index", "time_s", "dm_id", "sd_id", "sd_x", "sd_z", "sd_r")
        missing = [name for name in required if name not in handle.variables]
        if missing:
            raise RuntimeError("trajectory file is missing variables: {0}".format(", ".join(missing)))
        arrays = {}
        for name in required:
            variable = handle.variables[name]
            fill_value = -999 if np.issubdtype(variable.dtype, np.integer) else np.nan
            arrays[name] = np.ma.filled(variable[:], fill_value)
        metadata = {
            "tracking_mode": getattr(handle, "tracking_mode", "unknown"),
            "x_period_m": float(getattr(handle, "x_period_m", np.nan)),
        }
    rows = defaultdict(list)
    for index in range(len(arrays["time_s"])):
        key = int(arrays["trajectory_index"][index])
        rows[key].append(
            {
                "time_s": float(arrays["time_s"][index]),
                "dm_id": int(arrays["dm_id"][index]),
                "sd_id": int(arrays["sd_id"][index]),
                "x_m": float(arrays["sd_x"][index]),
                "z_m": float(arrays["sd_z"][index]),
                "radius_um": float(arrays["sd_r"][index]) * 1.0e6,
            }
        )
    for trajectory_rows in rows.values():
        trajectory_rows.sort(key=lambda row: row["time_s"])
    return rows, metadata


def _choose_trajectories(rows, args):
    """Choose exact targets or a deterministic growth-oriented sample."""
    pair_to_index = {}
    for index, trajectory_rows in rows.items():
        if trajectory_rows:
            pair_to_index[(trajectory_rows[0]["dm_id"], trajectory_rows[0]["sd_id"])] = index

    chosen = []
    for target in args.target:
        pair = _parse_pair(target)
        if pair not in pair_to_index:
            raise RuntimeError("target {0}:{1} is absent from trajectory file".format(*pair))
        chosen.append(pair_to_index[pair])
    chosen.extend(args.trajectory_index)
    if chosen:
        return list(dict.fromkeys(chosen))

    candidates = []
    for index, trajectory_rows in rows.items():
        radii = [row["radius_um"] for row in trajectory_rows if math.isfinite(row["radius_um"])]
        if not radii:
            continue
        initial_radius = radii[0]
        maximum_radius = max(radii)
        if initial_radius < args.initial_radius_max_um and maximum_radius >= args.max_radius_min_um:
            candidates.append((index, maximum_radius, initial_radius))
    if not candidates:
        print(
            "WARNING: no trajectories satisfy the radius-growth filter; plotting the largest available targets",
            file=sys.stderr,
        )
        for index, trajectory_rows in rows.items():
            radii = [row["radius_um"] for row in trajectory_rows if math.isfinite(row["radius_um"])]
            if radii:
                candidates.append((index, max(radii), radii[0]))
    candidates.sort(key=lambda item: (-item[1], item[0]))
    return [item[0] for item in candidates[: max(0, args.number)]]


def _gap_limit(rows):
    """Infer a conservative gap threshold from actual output times."""
    times = sorted(set(row["time_s"] for trajectory_rows in rows.values() for row in trajectory_rows))
    differences = [b - a for a, b in zip(times[:-1], times[1:]) if b > a]
    return 1.5 * float(np.median(differences)) if differences else math.inf


def _line_segments(trajectory_rows, x_period, gap_limit, relative_position=False):
    """Build only physically adjacent segments, breaking gaps and boundary wraps."""
    valid_rows = [
        row
        for row in trajectory_rows
        if all(math.isfinite(row[name]) for name in ("time_s", "x_m", "z_m", "radius_um"))
        and row["radius_um"] > 0.0
    ]
    if x_period and math.isfinite(x_period) and x_period > 0.0:
        for row in valid_rows:
            row["plot_x_m"] = row["x_m"] % x_period
    else:
        for row in valid_rows:
            row["plot_x_m"] = row["x_m"]
    if valid_rows and relative_position:
        valid_rows[0]["display_x_m"] = 0.0
        valid_rows[0]["display_z_m"] = 0.0
        for previous, current in zip(valid_rows[:-1], valid_rows[1:]):
            delta_x = current["plot_x_m"] - previous["plot_x_m"]
            if x_period and math.isfinite(x_period) and x_period > 0.0:
                if delta_x > 0.5 * x_period:
                    delta_x -= x_period
                elif delta_x < -0.5 * x_period:
                    delta_x += x_period
            current["display_x_m"] = previous["display_x_m"] + delta_x
            current["display_z_m"] = current["z_m"] - valid_rows[0]["z_m"]
    else:
        for row in valid_rows:
            row["display_x_m"] = row["plot_x_m"]
            row["display_z_m"] = row["z_m"]
    segments = []
    colors = []
    for first, second in zip(valid_rows[:-1], valid_rows[1:]):
        if abs(second["time_s"] - first["time_s"]) > gap_limit:
            continue
        if x_period and math.isfinite(x_period) and x_period > 0.0:
            if abs(second["plot_x_m"] - first["plot_x_m"]) > 0.5 * x_period:
                continue
        segments.append(
            [
                [first["display_x_m"], first["display_z_m"]],
                [second["display_x_m"], second["display_z_m"]],
            ]
        )
        colors.append(0.5 * (first["radius_um"] + second["radius_um"]))
    return valid_rows, segments, colors


def _configure_style(style):
    """Apply the repository's GMD2026 style contract or a screen style."""
    if style == "gmd":
        parameters = {
            "font.family": "DejaVu Sans",
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
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    else:
        parameters = {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    plt.rcParams.update(parameters)


def _figure_size(style, width):
    """Return a manuscript column size or the legacy screen-oriented size."""
    if style == "diagnostic":
        return 7.2, 4.4
    width_mm = 85.0 if width == "single" else 170.0
    width_in = width_mm / 25.4
    return width_in, width_in * 0.58


def main(argv=None):
    """Create PDF, SVG, and high-resolution PNG trajectory figures."""
    args = parse_args(argv)
    if args.number < 1 and not args.target and not args.trajectory_index:
        raise SystemExit("--number must be positive when no exact target is supplied")
    if args.full_domain and args.relative_position:
        raise SystemExit("--full-domain and --relative-position are mutually exclusive")
    input_path = Path(args.input).resolve()
    try:
        rows, metadata = _read_rows(input_path)
        chosen = _choose_trajectories(rows, args)
    except (OSError, RuntimeError, ValueError) as exc:
        raise SystemExit("ERROR: {0}".format(exc))
    if not chosen:
        raise SystemExit("ERROR: no trajectories are available for plotting")

    x_period = args.x_period_m if args.x_period_m is not None else metadata["x_period_m"]
    gap_limit = _gap_limit(rows)
    plotted = []
    all_radii = []
    for index in chosen:
        if index not in rows:
            raise SystemExit("ERROR: trajectory index {0} is absent".format(index))
        trajectory_rows = list(rows[index])
        if args.reverse_time:
            trajectory_rows.reverse()
        valid_rows, segments, colors = _line_segments(
            trajectory_rows, x_period, gap_limit, relative_position=args.relative_position
        )
        if valid_rows:
            plotted.append((index, valid_rows, segments, colors))
            all_radii.extend(row["radius_um"] for row in valid_rows)
    if not plotted:
        raise SystemExit("ERROR: selected trajectories contain no finite positive-radius records")

    positive_radii = np.asarray([value for value in all_radii if value > 0.0], dtype=float)
    vmin = float(np.min(positive_radii))
    vmax = float(np.max(positive_radii))
    norm = LogNorm(vmin=vmin, vmax=vmax) if vmax > vmin else Normalize(vmin=0.9 * vmin, vmax=1.1 * vmax)
    cmap = plt.get_cmap("viridis")

    _configure_style(args.style)
    figure, axis = plt.subplots(
        figsize=_figure_size(args.style, args.figure_width), constrained_layout=True
    )
    point_size = 5 if args.style == "gmd" else 8
    start_size = 20 if args.style == "gmd" else 28
    end_size = 24 if args.style == "gmd" else 34
    for plot_number, (index, trajectory_rows, segments, colors) in enumerate(plotted):
        if segments:
            collection = LineCollection(segments, cmap=cmap, norm=norm, linewidths=1.0, alpha=0.9)
            collection.set_array(np.asarray(colors, dtype=float))
            axis.add_collection(collection)
        x_values = [row["display_x_m"] for row in trajectory_rows]
        z_values = [row["display_z_m"] for row in trajectory_rows]
        radius_values = [row["radius_um"] for row in trajectory_rows]
        axis.scatter(
            x_values,
            z_values,
            c=radius_values,
            cmap=cmap,
            norm=norm,
            s=point_size,
            linewidths=0.2,
            edgecolors="white",
        )
        axis.scatter(
            x_values[0], z_values[0], s=start_size, marker="o", facecolor="#0072B2", edgecolor="white", linewidth=0.5,
            label="Start" if plot_number == 0 else None, zorder=4,
        )
        axis.scatter(
            x_values[-1], z_values[-1], s=end_size, marker="x", color="#D55E00", linewidth=1.2,
            label="End" if plot_number == 0 else None, zorder=5,
        )

    axis.autoscale()
    axis.margins(x=0.04, y=0.06)
    if args.full_domain and x_period and math.isfinite(x_period) and x_period > 0.0:
        axis.set_xlim(0.0, x_period)
    if args.full_domain and args.z_min_m is None:
        axis.set_ylim(bottom=0.0)
    if args.z_min_m is not None or args.z_max_m is not None:
        current_bottom, current_top = axis.get_ylim()
        axis.set_ylim(
            args.z_min_m if args.z_min_m is not None else current_bottom,
            args.z_max_m if args.z_max_m is not None else current_top,
        )
    if args.relative_position:
        axis.set_xlabel("Horizontal displacement, Δx (m)")
        axis.set_ylabel("Vertical displacement, Δz (m)")
    else:
        axis.set_xlabel("X position (m)")
        axis.set_ylabel("Z position (m)")
    automatic_title = None
    if args.style == "diagnostic":
        if len(plotted) == 1:
            first = plotted[0][1][0]
            automatic_title = "SD trajectory ({0}:{1})".format(first["dm_id"], first["sd_id"])
        else:
            automatic_title = "{0} {1}-tracking trajectories".format(
                len(plotted), metadata["tracking_mode"]
            )
    if args.title is not None or automatic_title is not None:
        axis.set_title(args.title if args.title is not None else automatic_title)
    if args.panel_label:
        axis.text(
            -0.10,
            1.02,
            args.panel_label,
            transform=axis.transAxes,
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    axis.legend(frameon=False, loc="best")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = figure.colorbar(scalar, ax=axis, pad=0.02)
    colorbar.set_label("Radius (µm)")
    colorbar.ax.tick_params(direction="in", width=0.8, labelsize=7 if args.style == "gmd" else 8)

    output_base = Path(args.output_base).resolve()
    output_base.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix, kwargs in ((".pdf", {}), (".svg", {}), (".png", {"dpi": 600})):
        output = output_base.with_suffix(suffix)
        figure.savefig(output, bbox_inches="tight", **kwargs)
        outputs.append(output)
    if args.show:
        plt.show()
    plt.close(figure)
    for output in outputs:
        print("wrote: {0}".format(output))


if __name__ == "__main__":
    main()
