#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PALETTE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#7F7F7F",
    "light_gray": "#D9D9D9",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-style TPHT figures from analyze_tpht_tracks.py outputs. "
            "The script writes both PDF and PNG versions for manuscript and quick-view use."
        )
    )
    parser.add_argument(
        "--analysis-dir",
        required=True,
        help="Directory containing tpht_* CSV/JSON outputs from analyze_tpht_tracks.py.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="Directory where figures will be written. Defaults to <analysis-dir>/figures.",
    )
    return parser.parse_args()


def set_publication_style():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.minor.width": 0.6,
            "ytick.minor.width": 0.6,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.minor.size": 2.0,
            "ytick.minor.size": 2.0,
            "lines.linewidth": 1.6,
            "lines.markersize": 4.0,
            "figure.dpi": 180,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": "#EAEAEA",
            "grid.linewidth": 0.6,
        }
    )


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def read_csv(analysis_dir: str, name: str) -> pd.DataFrame:
    path = os.path.join(analysis_dir, name)
    if not os.path.exists(path):
        raise SystemExit(f"Required analysis file was not found: {path}")
    return pd.read_csv(path)


def add_panel_label(ax, label: str):
    ax.text(
        -0.12,
        1.03,
        label,
        transform=ax.transAxes,
        fontsize=9.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=PALETTE["black"],
    )


def save_figure(fig: plt.Figure, output_dir: str, stem: str):
    pdf_path = os.path.join(output_dir, f"{stem}.pdf")
    png_path = os.path.join(output_dir, f"{stem}.png")
    fig.savefig(pdf_path)
    fig.savefig(png_path)
    plt.close(fig)
    return pdf_path, png_path


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(False, index=frame.index)
    return frame[column].astype(str).str.lower().isin(["true", "1", "t", "yes"])


def build_pair_key(frame: pd.DataFrame) -> pd.Series:
    return frame["pre_dmid"].astype(int).astype(str) + ":" + frame["pre_sdid"].astype(int).astype(str)


def compute_collision_cumulative(events: pd.DataFrame):
    tracked_events = events[events["tracked_any"]].copy()
    tracked_events["pair1_key"] = (
        tracked_events["pre_dmid1"].astype(int).astype(str)
        + ":"
        + tracked_events["pre_sdid1"].astype(int).astype(str)
    )
    tracked_events["pair2_key"] = (
        tracked_events["pre_dmid2"].astype(int).astype(str)
        + ":"
        + tracked_events["pre_sdid2"].astype(int).astype(str)
    )

    seen = set()
    rows = []
    for time_label, group in tracked_events.groupby("time_label", sort=True):
        for _, row in group.iterrows():
            if bool(row["pair1_tracked"]):
                seen.add(row["pair1_key"])
            if bool(row["pair2_tracked"]):
                seen.add(row["pair2_key"])
        rows.append({"time_label": time_label, "cumulative_tracked_pairs": len(seen)})
    return pd.DataFrame(rows)


def prepare_inputs(analysis_dir: str):
    time_summary = read_csv(analysis_dir, "tpht_time_summary.csv").sort_values("time_seconds_from_start")
    pair_timeseries = read_csv(analysis_dir, "tpht_pair_timeseries.csv").sort_values(
        ["pre_dmid", "pre_sdid", "time_seconds_from_start"]
    )
    pair_overview = read_csv(analysis_dir, "tpht_pair_overview.csv")
    collision_time = read_csv(analysis_dir, "tpht_collision_time_summary.csv")
    collision_events = read_csv(analysis_dir, "tpht_collision_events.csv")

    pair_overview["pair_key"] = build_pair_key(pair_overview)
    pair_timeseries["pair_key"] = build_pair_key(pair_timeseries)

    collision_events["pair1_tracked"] = bool_series(collision_events, "pair1_tracked")
    collision_events["pair2_tracked"] = bool_series(collision_events, "pair2_tracked")
    collision_events["tracked_any"] = bool_series(collision_events, "tracked_any")
    collision_events["tracked_both"] = bool_series(collision_events, "tracked_both")

    pair_overview["collided"] = (pair_overview["event_count"] > 0) | (
        pair_overview["if_coal_time_count"] > 0
    )
    pair_overview["growth_ratio"] = pair_overview["last_sd_r"] / np.maximum(
        pair_overview["first_sd_r"], 1.0e-20
    )

    pair_meta = pair_overview[
        [
            "pair_key",
            "pre_dmid",
            "pre_sdid",
            "collided",
            "event_count",
            "num_col_sum",
            "partner_pair_count",
            "first_sd_r",
            "last_sd_r",
            "delta_sd_r",
            "first_sd_z",
            "last_sd_z",
            "delta_sd_z",
            "max_radius",
            "growth_ratio",
        ]
    ]
    pair_timeseries = pair_timeseries.merge(pair_meta, on=["pair_key", "pre_dmid", "pre_sdid"], how="left")

    return time_summary, pair_timeseries, pair_overview, collision_time, collision_events


def plot_overview_timeseries(output_dir: str, time_summary: pd.DataFrame, collision_time: pd.DataFrame):
    merged = time_summary.merge(
        collision_time[
            [
                "time_label",
                "event_count_tracked_any",
                "num_col_sum_tracked_any",
                "tracked_pair_count_any",
            ]
        ],
        on="time_label",
        how="left",
    ).fillna(0)

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.8), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.plot(
        merged["time_seconds_from_start"],
        merged["tracked_count"],
        color=PALETTE["black"],
        label="Tracked pairs",
    )
    ax.plot(
        merged["time_seconds_from_start"],
        merged["nonlocal_pairs"],
        color=PALETTE["orange"],
        linestyle="--",
        label="Cross-rank pairs",
    )
    ax.set_ylabel("Pair count")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, axis="y")
    add_panel_label(ax, "A")

    ax = axes[1]
    ax.plot(
        merged["time_seconds_from_start"],
        merged["if_coal_count"],
        color=PALETTE["blue"],
        label="Tracked pairs with if_coal > 0",
    )
    ax.bar(
        merged["time_seconds_from_start"],
        merged["event_count_tracked_any"],
        width=3.6,
        color=PALETTE["vermillion"],
        alpha=0.35,
        label="Collision events involving tracked pairs",
    )
    ax.set_ylabel("Collision activity")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, axis="y")
    add_panel_label(ax, "B")

    ax = axes[2]
    x = merged["time_seconds_from_start"].to_numpy()
    y_mean = merged["sd_r_mean"].to_numpy()
    y_min = merged["sd_r_min"].to_numpy()
    y_max = merged["sd_r_max"].to_numpy()
    ax.fill_between(x, y_min, y_max, color=PALETTE["sky"], alpha=0.22, label="Min-max envelope")
    ax.plot(x, y_mean, color=PALETTE["green"], label="Mean radius")
    ax.set_yscale("log")
    ax.set_ylabel("Radius [m]")
    ax.set_xlabel("Time from BW start [s]")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, axis="y", which="both")
    add_panel_label(ax, "C")

    return save_figure(fig, output_dir, "fig01_tpht_overview_timeseries")


def split_wrapped_segments(group: pd.DataFrame, domain_x: float) -> Iterable[pd.DataFrame]:
    if len(group) <= 1:
        yield group
        return

    split_index = [0]
    x = group["sd_x"].to_numpy()
    for idx in range(1, len(group)):
        if abs(x[idx] - x[idx - 1]) > 0.5 * domain_x:
            split_index.append(idx)
    split_index.append(len(group))

    for start, end in zip(split_index[:-1], split_index[1:]):
        yield group.iloc[start:end]


def plot_trajectory_panels(output_dir: str, time_summary: pd.DataFrame, pair_timeseries: pd.DataFrame, pair_overview: pd.DataFrame):
    collided = pair_overview[pair_overview["collided"]].sort_values(
        ["num_col_sum", "max_radius"], ascending=[False, False]
    )
    large = pair_overview.sort_values(["max_radius", "delta_sd_r"], ascending=[False, False])
    selected_keys = list(dict.fromkeys(list(collided["pair_key"].head(12)) + list(large["pair_key"].head(12))))
    sample = pair_timeseries[pair_timeseries["pair_key"].isin(selected_keys)].copy()
    sample.sort_values(["pair_key", "time_seconds_from_start"], inplace=True)

    domain_x = float(time_summary["sd_x_max"].max() - time_summary["sd_x_min"].min())
    if domain_x <= 0.0:
        domain_x = float(pair_timeseries["sd_x"].max() - pair_timeseries["sd_x"].min())

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7), constrained_layout=True)

    for pair_key, group in sample.groupby("pair_key", sort=False):
        collided_flag = bool(group["collided"].iloc[0])
        color = PALETTE["vermillion"] if collided_flag else PALETTE["blue"]
        alpha = 0.95 if collided_flag else 0.55

        for segment in split_wrapped_segments(group, domain_x):
            axes[0].plot(segment["sd_x"], segment["sd_z"], color=color, alpha=alpha, linewidth=1.2)
        axes[1].plot(
            group["time_seconds_from_start"],
            group["sd_z"],
            color=color,
            alpha=alpha,
            linewidth=1.2,
        )
        axes[0].scatter(group["sd_x"].iloc[0], group["sd_z"].iloc[0], s=14, color=color, marker="o", zorder=3)
        axes[0].scatter(group["sd_x"].iloc[-1], group["sd_z"].iloc[-1], s=18, color=color, marker="^", zorder=3)

    axes[0].set_xlabel("x [m]")
    axes[0].set_ylabel("z [m]")
    axes[0].set_title("Representative TPHT trajectories")
    axes[0].grid(True)
    add_panel_label(axes[0], "A")

    axes[1].set_xlabel("Time from BW start [s]")
    axes[1].set_ylabel("z [m]")
    axes[1].set_title("Vertical evolution of the same pairs")
    axes[1].grid(True)
    add_panel_label(axes[1], "B")

    legend_handles = [
        mpl.lines.Line2D([], [], color=PALETTE["vermillion"], linewidth=1.8, label="Collision-involved pairs"),
        mpl.lines.Line2D([], [], color=PALETTE["blue"], linewidth=1.8, label="Non-colliding large pairs"),
    ]
    axes[1].legend(handles=legend_handles, frameon=False, loc="best")

    return save_figure(fig, output_dir, "fig02_tpht_trajectory_panels")


def plot_radius_phase_space(output_dir: str, pair_overview: pd.DataFrame):
    collided = pair_overview[pair_overview["collided"]]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), constrained_layout=True)

    ax = axes[0]
    hb = ax.hexbin(
        pair_overview["first_sd_r"],
        pair_overview["last_sd_r"],
        gridsize=55,
        xscale="log",
        yscale="log",
        mincnt=1,
        cmap="cividis",
    )
    ax.plot(
        [pair_overview["first_sd_r"].min(), pair_overview["last_sd_r"].max()],
        [pair_overview["first_sd_r"].min(), pair_overview["last_sd_r"].max()],
        linestyle="--",
        linewidth=1.0,
        color=PALETTE["gray"],
    )
    if not collided.empty:
        ax.scatter(
            collided["first_sd_r"],
            collided["last_sd_r"],
            s=12,
            color=PALETTE["vermillion"],
            alpha=0.75,
            edgecolors="white",
            linewidths=0.25,
            label="Collision-involved pairs",
        )
        ax.legend(frameon=False, loc="upper left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Initial radius [m]")
    ax.set_ylabel("Final radius [m]")
    ax.set_title("Radius growth phase space")
    add_panel_label(ax, "A")
    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label("Pair count per bin")

    ax = axes[1]
    growth_non = pair_overview.loc[~pair_overview["collided"], "growth_ratio"]
    growth_col = pair_overview.loc[pair_overview["collided"], "growth_ratio"]
    bins = np.geomspace(
        max(pair_overview["growth_ratio"].min(), 1.0e-2),
        max(pair_overview["growth_ratio"].max(), 1.0e0),
        40,
    )
    ax.hist(
        growth_non,
        bins=bins,
        color=PALETTE["blue"],
        alpha=0.45,
        label="Non-colliding",
    )
    if len(growth_col) > 0:
        ax.hist(
            growth_col,
            bins=bins,
            color=PALETTE["vermillion"],
            alpha=0.55,
            label="Collision-involved",
        )
    ax.set_xscale("log")
    ax.set_xlabel("Final / initial radius")
    ax.set_ylabel("Pair count")
    ax.set_title("Radius-amplification distribution")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(True, axis="y")
    add_panel_label(ax, "B")

    return save_figure(fig, output_dir, "fig03_tpht_radius_phase_space")


def plot_collision_timeline(
    output_dir: str,
    collision_time: pd.DataFrame,
    collision_events: pd.DataFrame,
):
    cumulative = compute_collision_cumulative(collision_events)
    merged = collision_time.merge(cumulative, on="time_label", how="left").fillna(0)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.bar(
        merged["time_seconds_from_start"],
        merged["event_count_tracked_any"],
        width=3.6,
        color=PALETTE["vermillion"],
        alpha=0.45,
        label="Events involving tracked pairs",
    )
    ax.plot(
        merged["time_seconds_from_start"],
        merged["num_col_sum_tracked_any"],
        color=PALETTE["purple"],
        linewidth=1.6,
        marker="o",
        label=r"$\Sigma$ num_col (tracked-any)",
    )
    ax.set_ylabel("Event intensity")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, axis="y")
    add_panel_label(ax, "A")

    ax = axes[1]
    ax.plot(
        merged["time_seconds_from_start"],
        merged["cumulative_tracked_pairs"],
        color=PALETTE["black"],
        linewidth=1.8,
        marker="o",
        label="Cumulative collision-involved tracked pairs",
    )
    ax.plot(
        merged["time_seconds_from_start"],
        merged["tracked_pair_count_any"],
        color=PALETTE["orange"],
        linewidth=1.4,
        linestyle="--",
        marker="s",
        label="Tracked pairs seen in events at this output time",
    )
    ax.set_xlabel("Time from BW start [s]")
    ax.set_ylabel("Tracked pair count")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(True, axis="y")
    add_panel_label(ax, "B")

    return save_figure(fig, output_dir, "fig04_tpht_collision_timeline")


def quantile_summary(frame: pd.DataFrame, value_col: str) -> pd.DataFrame:
    grouped = frame.groupby(["time_seconds_from_start", "collided"])[value_col]
    summary = grouped.quantile([0.25, 0.5, 0.75]).unstack()
    summary = summary.rename(columns={0.25: "q25", 0.5: "q50", 0.75: "q75"}).reset_index()
    return summary


def plot_group_evolution(output_dir: str, pair_timeseries: pd.DataFrame):
    radius_summary = quantile_summary(pair_timeseries, "sd_r")
    height_summary = quantile_summary(pair_timeseries, "sd_z")

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.3), sharex=True, constrained_layout=True)

    for collided_flag, color, label in [
        (False, PALETTE["blue"], "Non-colliding pairs"),
        (True, PALETTE["vermillion"], "Collision-involved pairs"),
    ]:
        subset = height_summary[height_summary["collided"] == collided_flag]
        if subset.empty:
            continue
        axes[0].fill_between(
            subset["time_seconds_from_start"],
            subset["q25"],
            subset["q75"],
            color=color,
            alpha=0.16,
        )
        axes[0].plot(
            subset["time_seconds_from_start"],
            subset["q50"],
            color=color,
            label=label,
        )

    axes[0].set_ylabel("z [m]")
    axes[0].set_title("Vertical-position evolution")
    axes[0].legend(frameon=False, loc="best")
    axes[0].grid(True, axis="y")
    add_panel_label(axes[0], "A")

    for collided_flag, color in [
        (False, PALETTE["blue"]),
        (True, PALETTE["vermillion"]),
    ]:
        subset = radius_summary[radius_summary["collided"] == collided_flag]
        if subset.empty:
            continue
        axes[1].fill_between(
            subset["time_seconds_from_start"],
            subset["q25"],
            subset["q75"],
            color=color,
            alpha=0.16,
        )
        axes[1].plot(subset["time_seconds_from_start"], subset["q50"], color=color)

    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time from BW start [s]")
    axes[1].set_ylabel("Radius [m]")
    axes[1].set_title("Radius evolution")
    axes[1].grid(True, axis="y", which="both")
    add_panel_label(axes[1], "B")

    return save_figure(fig, output_dir, "fig05_tpht_group_evolution")


def main():
    args = parse_args()
    analysis_dir = os.path.abspath(args.analysis_dir)
    output_dir = os.path.abspath(args.output_dir or os.path.join(analysis_dir, "figures"))

    ensure_dir(output_dir)
    set_publication_style()

    (
        time_summary,
        pair_timeseries,
        pair_overview,
        collision_time,
        collision_events,
    ) = prepare_inputs(analysis_dir)

    written = []
    written.extend(plot_overview_timeseries(output_dir, time_summary, collision_time))
    written.extend(plot_trajectory_panels(output_dir, time_summary, pair_timeseries, pair_overview))
    written.extend(plot_radius_phase_space(output_dir, pair_overview))
    written.extend(plot_collision_timeline(output_dir, collision_time, collision_events))
    written.extend(plot_group_evolution(output_dir, pair_timeseries))

    print(f"figure_count={len(written) // 2}")
    print(f"output_dir={output_dir}")
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
