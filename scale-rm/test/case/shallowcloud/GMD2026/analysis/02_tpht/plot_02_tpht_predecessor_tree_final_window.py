#!/usr/bin/env python3
"""Plot final-window TPHT predecessor-tree diagnostic examples.

This script uses existing GMD2026 analysis source tables, not raw simulation
output.  It selects threshold-crossing targets that remain diagnostically
interesting during the final analysis window and draws one compact target
history per figure.  Coalescence markers are target-linked event diagnostics,
not a complete causal formation tree.  Triangle colors represent the
coalescing partner radius from the target-linked event table.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.plot_style import OKABE_ITO, configure_matplotlib, figure_size, save_figure  # noqa: E402
from common.table_utils import safe_float, write_table_bundle  # noqa: E402


TARGET_TABLE = "02_tpht_science_target_summary.csv"
TRAJECTORY_TABLE = "02_tpht_target_trajectory_records.csv"
EVENT_TABLE = "02_tpht_target_event_links.csv"
LINK_TABLE = "02_tpht_chain_links_by_time.csv"

SELECTED_COLUMNS = [
    "rank",
    "target_id",
    "selection_window_start_s",
    "selection_window_end_s",
    "category",
    "window_record_count",
    "window_large_record_count",
    "window_ifcoal_record_count",
    "window_max_radius_um",
    "window_max_height_m",
    "event_row_count_window",
    "event_row_count_total",
    "num_col_sum_window",
    "num_col_sum_total",
    "max_event_delta_r3_um3",
    "first_large_time_s",
    "first_ifcoal_time_s",
    "warnings",
]


def _open_rows(path: Path):
    """Yield CSV rows from an existing analysis table."""
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle)


def _resolve_table_dir(analysis_outdir: Path, gmd_outdir: Path) -> Path:
    """Find a table directory after analysis_outputs has been archived or removed."""
    candidates = [
        analysis_outdir / "tables",
        analysis_outdir,
        gmd_outdir / "figure_sources" / "supplement_candidate_source_tables",
        gmd_outdir / "figure_sources" / "main_candidate_source_tables",
    ]
    for candidate in candidates:
        if (candidate / TRAJECTORY_TABLE).exists() and (candidate / EVENT_TABLE).exists():
            return candidate
    return analysis_outdir / "tables"


def _max_output_time_s(tables_dir: Path) -> float:
    """Read the maximum TPHT output time from the compact chain-link table."""
    path = tables_dir / LINK_TABLE
    maximum = None
    if path.exists():
        for row in _open_rows(path):
            time_s = safe_float(row.get("time_s"))
            if time_s is not None:
                maximum = time_s if maximum is None else max(maximum, time_s)
    if maximum is None:
        maximum = 4200.0
    return maximum


def _scan_trajectory_candidates(path: Path, start_s: float, end_s: float, min_radius_um: float) -> dict[str, dict[str, Any]]:
    """Find sampled targets with final-window large-radius or if_coal records."""
    candidates: dict[str, dict[str, Any]] = {}
    for row in _open_rows(path):
        target_id = row.get("target_id")
        time_s = safe_float(row.get("time_s"))
        if not target_id or target_id == "NA" or time_s is None or time_s < start_s or time_s > end_s:
            continue
        radius_um = safe_float(row.get("radius_um"))
        height_m = safe_float(row.get("z_m"))
        if_coal = safe_float(row.get("if_coal")) or 0.0
        is_large = radius_um is not None and radius_um >= min_radius_um
        is_ifcoal = if_coal > 0.0
        if not (is_large or is_ifcoal):
            continue
        state = candidates.setdefault(
            target_id,
            {
                "target_id": target_id,
                "category": row.get("category"),
                "window_record_count": 0,
                "window_large_record_count": 0,
                "window_ifcoal_record_count": 0,
                "window_max_radius_um": None,
                "window_max_height_m": None,
                "first_large_time_s": safe_float(row.get("first_large_time_s")),
                "first_ifcoal_time_s": safe_float(row.get("first_ifcoal_time_s")),
            },
        )
        state["window_record_count"] += 1
        if is_large:
            state["window_large_record_count"] += 1
        if is_ifcoal:
            state["window_ifcoal_record_count"] += 1
        if radius_um is not None:
            current = state.get("window_max_radius_um")
            state["window_max_radius_um"] = radius_um if current is None else max(current, radius_um)
        if height_m is not None:
            current_h = state.get("window_max_height_m")
            state["window_max_height_m"] = height_m if current_h is None else max(current_h, height_m)
    return candidates


def _scan_events(path: Path, candidates: dict[str, dict[str, Any]], start_s: float, end_s: float) -> None:
    """Attach target-linked event counts to selected candidates."""
    if not candidates:
        return
    candidate_ids = set(candidates)
    for row in _open_rows(path):
        target_id = row.get("target_id")
        if target_id not in candidate_ids:
            continue
        state = candidates[target_id]
        event_time_s = safe_float(row.get("event_time_s"))
        delta_r3 = safe_float(row.get("delta_r3_um3"))
        num_col = safe_float(row.get("num_col")) or 0.0
        state["event_row_count_total"] = int(state.get("event_row_count_total") or 0) + 1
        state["num_col_sum_total"] = float(state.get("num_col_sum_total") or 0.0) + num_col
        if event_time_s is not None and start_s <= event_time_s <= end_s:
            state["event_row_count_window"] = int(state.get("event_row_count_window") or 0) + 1
            state["num_col_sum_window"] = float(state.get("num_col_sum_window") or 0.0) + num_col
        if delta_r3 is not None:
            current = state.get("max_event_delta_r3_um3")
            state["max_event_delta_r3_um3"] = delta_r3 if current is None else max(current, delta_r3)


def _select_targets(candidates: dict[str, dict[str, Any]], count: int) -> list[dict[str, Any]]:
    """Choose targets that best show final-window threshold/coalescence diagnostics."""
    ranked = sorted(
        candidates.values(),
        key=lambda row: (
            int(row.get("event_row_count_window") or 0),
            int(row.get("window_ifcoal_record_count") or 0),
            safe_float(row.get("window_max_radius_um")) or -math.inf,
            int(row.get("event_row_count_total") or 0),
            safe_float(row.get("max_event_delta_r3_um3")) or -math.inf,
        ),
        reverse=True,
    )
    return ranked[:count]


def _collect_rows(path: Path, target_ids: set[str]) -> dict[str, list[dict[str, Any]]]:
    """Collect rows for a small set of selected targets."""
    output = {target_id: [] for target_id in target_ids}
    for row in _open_rows(path):
        target_id = row.get("target_id")
        if target_id in output:
            output[target_id].append(row)
    return output


def _event_text(row: dict[str, Any], key: str) -> str:
    """Return a stable text value for event-row pairing."""
    value = row.get(key)
    if value in (None, "", "NA"):
        return ""
    return str(value).strip()


def _participant_text(row: dict[str, Any]) -> str:
    """Normalize the participant number for event-row pairing."""
    value = safe_float(row.get("participant"))
    if value is None:
        return _event_text(row, "participant")
    return str(int(round(value)))


def _event_row_key(row: dict[str, Any]) -> tuple[str, str, str, str, str, str]:
    """Return a key for one participant row in the event-link table."""
    return (
        _event_text(row, "source_file"),
        _event_text(row, "event_time_s"),
        _event_text(row, "num_col"),
        _participant_text(row),
        _event_text(row, "pre_radius_um"),
        _event_text(row, "partner_radius_um"),
    )


def _counterpart_event_key(row: dict[str, Any]) -> tuple[str, str, str, str, str, str] | None:
    """Return the expected key for the other participant row in the same event."""
    participant = _participant_text(row)
    if participant == "1":
        counterpart = "2"
    elif participant == "2":
        counterpart = "1"
    else:
        return None
    pre_radius = _event_text(row, "pre_radius_um")
    partner_radius = _event_text(row, "partner_radius_um")
    if not pre_radius or not partner_radius:
        return None
    return (
        _event_text(row, "source_file"),
        _event_text(row, "event_time_s"),
        _event_text(row, "num_col"),
        counterpart,
        partner_radius,
        pre_radius,
    )


def _attach_partner_heights(path: Path, events: dict[str, list[dict[str, Any]]]) -> None:
    """Attach formal partner heights when the matching participant row exists.

    SD_coal_output_NetCDF records do not carry participant heights.  The
    target-linked event table, however, contains one row per linked
    participant.  For plotting only, we recover the partner's selected-output
    height by matching the opposite participant row from the same event.
    """
    wanted: dict[tuple[str, str, str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for rows in events.values():
        for row in rows:
            key = _counterpart_event_key(row)
            if key is not None:
                wanted[key].append(row)
    remaining = sum(len(rows) for rows in wanted.values())
    if remaining == 0:
        return

    for row in _open_rows(path):
        key = _event_row_key(row)
        matches = wanted.get(key)
        if not matches:
            continue
        target_row = matches.pop(0)
        target_row["partner_height_m"] = row.get("event_height_m")
        target_row["partner_target_id"] = row.get("target_id")
        remaining -= 1
        if not matches:
            wanted.pop(key, None)
        if remaining <= 0:
            break


def _validate_partner_radius(events: dict[str, list[dict[str, Any]]]) -> None:
    """Require partner-radius values before drawing event triangles.

    Older source tables only contained the target participant radius.  Using
    those rows would silently color triangles by the wrong particle, so fail
    with an actionable message instead.
    """
    missing = 0
    total = 0
    for rows in events.values():
        for row in rows:
            if safe_float(row.get("event_time_s")) is None or safe_float(row.get("event_height_m")) is None:
                continue
            total += 1
            if "partner_radius_um" not in row or safe_float(row.get("partner_radius_um")) is None:
                missing += 1
    if total and missing:
        raise ValueError(
            "02_tpht_target_event_links.csv lacks finite partner_radius_um values "
            f"for {missing}/{total} plotted coalescence records. Rerun "
            "analysis/02_tpht/analyze_02_tpht_science.py with the updated scripts "
            "before regenerating predecessor-tree figures."
        )


def _finite(values: list[float | None]) -> list[float]:
    """Return finite values."""
    return [value for value in values if value is not None and math.isfinite(value)]


def _first_condition_label(category: Any) -> str:
    """Return an explicit first-selected interest-condition label."""
    labels = {
        "radius_only": "r>=15 µm",
        "coal_only": "if_coal",
        "both": "r>=15 µm and if_coal",
        "unknown": "unknown",
    }
    return labels.get(str(category), str(category))


def _highlighted_title(ax: Any, selected: dict[str, Any], panel_label: str, coalescence_record_count: int) -> None:
    """Add a two-line panel title with colored diagnostic values."""
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea, VPacker

    condition = _first_condition_label(selected.get("category"))
    max_radius = safe_float(selected.get("window_max_radius_um"))
    max_radius_text = "NA" if max_radius is None else f"{max_radius:.2f}"
    coalescence_record_text = f"{coalescence_record_count:d}"

    def text_area(text: str, color: str = OKABE_ITO["black"], weight: str = "normal") -> TextArea:
        return TextArea(text, textprops={"fontsize": 8.5, "color": color, "fontweight": weight})

    condition_color = {
        "radius_only": OKABE_ITO["blue"],
        "coal_only": OKABE_ITO["vermillion"],
        "both": OKABE_ITO["bluish_green"],
        "unknown": "0.35",
    }.get(str(selected.get("category")), "0.35")
    line_one = HPacker(
        children=[
            text_area(f"{panel_label} {selected['target_id']}; first condition: "),
            text_area(condition, condition_color, "bold"),
        ],
        align="baseline",
        pad=0,
        sep=0,
    )
    line_two = HPacker(
        children=[
            text_area("window max r="),
            text_area(max_radius_text, OKABE_ITO["vermillion"], "bold"),
            text_area(" µm; coalescence records="),
            text_area(coalescence_record_text, OKABE_ITO["bluish_green"], "bold"),
        ],
        align="baseline",
        pad=0,
        sep=0,
    )
    title_box = VPacker(children=[line_one, line_two], align="left", pad=0, sep=2)
    anchored_title = AnchoredOffsetbox(
        loc="lower left",
        child=title_box,
        bbox_to_anchor=(0.0, 1.02),
        bbox_transform=ax.transAxes,
        frameon=False,
        borderpad=0.0,
        pad=0.0,
    )
    ax.add_artist(anchored_title)
    # Reserve space for the anchored colored title under constrained layout.
    ax.set_title(" \n ", loc="left", fontsize=8.5)


def _plot_colored_line(ax: Any, plt: Any, trajectory: list[dict[str, Any]], radius_range: tuple[float, float]) -> Any:
    """Plot target height history colored by radius."""
    import numpy as np
    from matplotlib.collections import LineCollection

    points = []
    radii = []
    for row in trajectory:
        time_s = safe_float(row.get("time_s"))
        height_m = safe_float(row.get("z_m"))
        radius_um = safe_float(row.get("radius_um"))
        if time_s is None or height_m is None or radius_um is None:
            continue
        points.append((time_s / 60.0, height_m))
        radii.append(radius_um)
    if len(points) < 2:
        return None
    point_array = np.asarray(points, dtype=float)
    segments = np.stack([point_array[:-1], point_array[1:]], axis=1)
    segment_radii = np.asarray([(left + right) * 0.5 for left, right in zip(radii[:-1], radii[1:])], dtype=float)
    collection = LineCollection(segments, cmap=plt.get_cmap("cividis"), linewidth=1.0)
    collection.set_array(segment_radii)
    collection.set_clim(radius_range[0], radius_range[1])
    ax.add_collection(collection)
    ax.autoscale_view()
    return collection


def _plot_target_panel(
    ax: Any,
    plt: Any,
    selected: dict[str, Any],
    trajectory: list[dict[str, Any]],
    events: list[dict[str, Any]],
    panel_label: str,
    start_s: float,
    end_s: float,
    radius_range: tuple[float, float],
) -> Any:
    """Draw one target-history and all linked coalescence-event markers."""
    import numpy as np
    from matplotlib.colors import Normalize

    trajectory = sorted(trajectory, key=lambda row: safe_float(row.get("time_s")) or -math.inf)
    collection = _plot_colored_line(ax, plt, trajectory, radius_range)
    ax.axvspan(start_s / 60.0, end_s / 60.0, color="0.92", zorder=-5)
    radius_cmap = plt.get_cmap("cividis")
    radius_norm = Normalize(vmin=radius_range[0], vmax=radius_range[1])

    def radius_color(radius_um: float | None) -> Any:
        """Map a finite particle radius to the same color scale as the trajectory."""
        if radius_um is None or not math.isfinite(radius_um):
            return "0.55"
        return radius_cmap(radius_norm(np.clip(radius_um, radius_range[0], radius_range[1])))

    times = [safe_float(row.get("time_s")) for row in trajectory]
    heights = [safe_float(row.get("z_m")) for row in trajectory]
    radii = [safe_float(row.get("radius_um")) for row in trajectory]
    finite_times = _finite(times)
    finite_heights = _finite(heights)
    finite_radii = _finite(radii)
    if finite_times and finite_heights:
        ax.set_xlim(min(finite_times) / 60.0, max(finite_times) / 60.0)
        height_pad = max(20.0, (max(finite_heights) - min(finite_heights)) * 0.08)
        ax.set_ylim(min(finite_heights) - height_pad, max(finite_heights) + height_pad)

    event_rows = [
        row
        for row in events
        if safe_float(row.get("event_time_s")) is not None and safe_float(row.get("event_height_m")) is not None
    ]
    partner_heights = _finite([safe_float(row.get("partner_height_m")) for row in event_rows])
    if partner_heights:
        lower, upper = ax.get_ylim()
        height_pad = max(20.0, (max(finite_heights + partner_heights) - min(finite_heights + partner_heights)) * 0.08)
        ax.set_ylim(
            min(lower, min(partner_heights) - height_pad),
            max(upper, max(partner_heights) + height_pad),
        )

    large_rows = [row for row in trajectory if (safe_float(row.get("radius_um")) or -math.inf) >= 15.0]
    first_large_row = min(large_rows, key=lambda row: safe_float(row.get("time_s")) or math.inf) if large_rows else None
    max_radius_row = max(trajectory, key=lambda row: safe_float(row.get("radius_um")) or -math.inf) if trajectory else None
    final_row = max(trajectory, key=lambda row: safe_float(row.get("time_s")) or -math.inf) if trajectory else None
    key_markers = [
        (first_large_row, "o", "first r>=15 µm"),
        (max_radius_row, "*", "maximum radius output"),
        (final_row, "s", "final output"),
    ]
    for row, marker, label in key_markers:
        if row is None:
            continue
        time_s = safe_float(row.get("time_s"))
        height_m = safe_float(row.get("z_m"))
        radius_um = safe_float(row.get("radius_um"))
        if time_s is None or height_m is None:
            continue
        ax.scatter(
            [time_s / 60.0],
            [height_m],
            marker=marker,
            s=48 if marker == "*" else 32,
            facecolor=radius_color(radius_um),
            edgecolors=OKABE_ITO["black"],
            linewidths=0.30,
            zorder=6,
            label=label,
        )

    if event_rows:
        marker_x: list[float] = []
        marker_y: list[float] = []
        marker_colors: list[Any] = []
        used_offset_marker = False
        y_lower, y_upper = ax.get_ylim()
        x_lower, x_upper = ax.get_xlim()
        y_span = max(y_upper - y_lower, 1.0)
        x_span = max(x_upper - x_lower, 1.0e-6)
        for event_index, row in enumerate(event_rows):
            event_time_min = (safe_float(row.get("event_time_s")) or 0.0) / 60.0
            event_height = safe_float(row.get("event_height_m")) or math.nan
            partner_height = safe_float(row.get("partner_height_m"))
            if partner_height is None:
                used_offset_marker = True
                sign = 1.0 if event_index % 2 == 0 else -1.0
                candidate_y = event_height + sign * 0.11 * y_span
                if candidate_y > y_upper - 0.035 * y_span:
                    candidate_y = event_height - 0.11 * y_span
                if candidate_y < y_lower + 0.035 * y_span:
                    candidate_y = event_height + 0.11 * y_span
                candidate_y = min(max(candidate_y, y_lower + 0.03 * y_span), y_upper - 0.03 * y_span)
                sign_x = 1.0 if (event_index // 2) % 2 == 0 else -1.0
                candidate_x = min(max(event_time_min + sign_x * 0.012 * x_span, x_lower), x_upper)
                ax.plot(
                    [event_time_min, candidate_x],
                    [event_height, candidate_y],
                    color="0.45",
                    linestyle=(0, (2.0, 2.0)),
                    linewidth=0.45,
                    alpha=0.55,
                    zorder=5,
                )
            else:
                candidate_y = partner_height
                if abs(partner_height - event_height) < max(8.0, 0.015 * y_span):
                    used_offset_marker = True
                    sign_x = 1.0 if event_index % 2 == 0 else -1.0
                    candidate_x = min(max(event_time_min + sign_x * 0.018 * x_span, x_lower), x_upper)
                    ax.plot(
                        [event_time_min, candidate_x],
                        [event_height, candidate_y],
                        color="0.45",
                        linestyle=(0, (2.0, 2.0)),
                        linewidth=0.45,
                        alpha=0.55,
                        zorder=5,
                    )
                else:
                    candidate_x = event_time_min
            marker_x.append(candidate_x)
            marker_y.append(candidate_y)
            marker_colors.append(radius_color(safe_float(row.get("partner_radius_um"))))
        ax.scatter(
            marker_x,
            marker_y,
            marker="^",
            s=14,
            c=marker_colors,
            edgecolors=OKABE_ITO["black"],
            linewidths=0.25,
            alpha=0.52,
            zorder=7,
            label="_nolegend_",
        )
        ax.scatter([], [], marker="^", s=24, facecolor=OKABE_ITO["black"], edgecolors=OKABE_ITO["black"], linewidths=0.25, label="coalescence partner record")
        if used_offset_marker:
            ax.plot([], [], color="0.45", linestyle=(0, (2.0, 2.0)), linewidth=0.6, label="dashed link to event point")

    _highlighted_title(ax, selected, panel_label, len(event_rows))
    ax.set_ylabel("Height (m)")
    ax.tick_params(direction="in")
    return collection


def _plot_target_batch(
    selected_rows: list[dict[str, Any]],
    trajectories: dict[str, list[dict[str, Any]]],
    events: dict[str, list[dict[str, Any]]],
    output_dir: Path,
    figure_rank: int,
    start_s: float,
    end_s: float,
) -> None:
    """Draw one three-panel predecessor-tree diagnostic figure."""
    plt = configure_matplotlib()
    fig = plt.figure(figsize=figure_size("double", 1.02), constrained_layout=True)
    grid = fig.add_gridspec(4, 1, height_ratios=[0.34, 1.0, 1.0, 1.0])
    legend_ax = fig.add_subplot(grid[0])
    legend_ax.axis("off")
    axes = []
    for index in range(3):
        axes.append(fig.add_subplot(grid[index + 1], sharex=axes[0] if axes else None))
    collections = []
    all_radii = [
        radius
        for selected in selected_rows
        for row in trajectories.get(selected["target_id"], [])
        for radius in [safe_float(row.get("radius_um"))]
        if radius is not None
    ]
    radius_min = min(all_radii) if all_radii else 0.0
    radius_max = max(all_radii) if all_radii and max(all_radii) > radius_min else radius_min + 1.0
    radius_range = (radius_min, radius_max)
    for panel_index, (ax, selected) in enumerate(zip(axes, selected_rows)):
        panel_label = f"({chr(97 + panel_index)})"
        collection = _plot_target_panel(
            ax,
            plt,
            selected,
            trajectories.get(selected["target_id"], []),
            events.get(selected["target_id"], []),
            panel_label,
            start_s,
            end_s,
            radius_range,
        )
        if collection is not None:
            collections.append(collection)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("Time (min)")
    handles: list[Any] = []
    labels: list[str] = []
    for ax in axes:
        axis_handles, axis_labels = ax.get_legend_handles_labels()
        for handle, label in zip(axis_handles, axis_labels):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        legend_ax.legend(
            handles,
            labels,
            frameon=False,
            fontsize=7,
            loc="center",
            ncol=3,
        )
    if collections:
        cbar = fig.colorbar(collections[0], ax=axes, pad=0.015, shrink=0.92)
        cbar.set_label(r"Particle radius ($\mu$m)")
    stem = f"supp_candidate_TPHT_predecessor_tree_examples_final10_{figure_rank:02d}"
    save_figure(fig, output_dir, stem)
    save_figure(fig, output_dir, f"supp_candidate_TPHT_predecessor_tree_examples_{figure_rank:02d}")
    plt.close(fig)


def _write_selected_table(rows: list[dict[str, Any]], tables_dir: Path, start_s: float, end_s: float) -> None:
    """Write the selected final-window target list for traceability."""
    output_rows = []
    for rank, row in enumerate(rows, start=1):
        output = {column: row.get(column) for column in SELECTED_COLUMNS}
        output["rank"] = rank
        output["selection_window_start_s"] = start_s
        output["selection_window_end_s"] = end_s
        output["warnings"] = "selected from sampled trajectory source table; coalescence markers use partner_radius_um and partner selected-output height when the paired participant row is available"
        output_rows.append(output)
    write_table_bundle(output_rows, tables_dir / "supp_candidate_TPHT_predecessor_tree_examples_final10_selected_targets", SELECTED_COLUMNS)


def _read_selected_table(path: Path, required_count: int) -> list[dict[str, Any]]:
    """Read a previously selected final-window target list."""
    if not path.exists():
        return []
    rows = [row for row in _open_rows(path) if row.get("target_id") not in (None, "", "NA")]
    rows.sort(key=lambda row: safe_float(row.get("rank")) or math.inf)
    return rows[:required_count] if len(rows) >= required_count else []


def build_parser() -> argparse.ArgumentParser:
    """Return the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-outdir", type=Path, required=True, help="Existing GMD2026 analysis_outputs directory.")
    parser.add_argument("--gmd-outdir", type=Path, required=True, help="analysis_outputs_for_GMD directory to receive candidate figures.")
    parser.add_argument("--window-minutes", type=float, default=10.0, help="Final-window length in minutes.")
    parser.add_argument("--max-figures", type=int, default=10, help="Number of target figures to generate.")
    parser.add_argument("--targets-per-figure", type=int, default=3, help="Number of target panels per figure.")
    parser.add_argument("--min-radius-um", type=float, default=15.0, help="Radius threshold for a final-window target.")
    parser.add_argument("--rescan", action="store_true", help="Ignore an existing selected-target table and rescan source tables.")
    parser.add_argument("--dry-run", action="store_true", help="Report inputs and outputs without reading large source tables.")
    return parser


def main() -> None:
    """CLI entry point."""
    args = build_parser().parse_args()
    analysis_tables = _resolve_table_dir(args.analysis_outdir, args.gmd_outdir)
    gmd_figures = args.gmd_outdir / "figures" / "supplement_candidates"
    gmd_tables = args.gmd_outdir / "tables" / "diagnostics"
    trajectory_path = analysis_tables / TRAJECTORY_TABLE
    event_path = analysis_tables / EVENT_TABLE
    end_s = _max_output_time_s(analysis_tables)
    start_s = max(0.0, end_s - args.window_minutes * 60.0)

    target_count = max(1, args.max_figures) * max(1, args.targets_per_figure)
    expected_outputs = [gmd_figures / f"supp_candidate_TPHT_predecessor_tree_examples_final10_{index:02d}.png" for index in range(1, args.max_figures + 1)]
    if args.dry_run:
        print(f"[dry-run] final window: {start_s:g}-{end_s:g} s")
        print(f"[dry-run] trajectory source: {trajectory_path}")
        print(f"[dry-run] event source: {event_path}")
        for output in expected_outputs:
            print(f"[dry-run] would write {output.with_suffix('.pdf')}")
            print(f"[dry-run] would write {output.with_suffix('.svg')}")
            print(f"[dry-run] would write {output}")
        return

    if not trajectory_path.exists():
        raise FileNotFoundError(f"missing source table: {trajectory_path}")
    if not event_path.exists():
        raise FileNotFoundError(f"missing source table: {event_path}")

    selected_table = gmd_tables / "supp_candidate_TPHT_predecessor_tree_examples_final10_selected_targets.csv"
    gmd_figures.mkdir(parents=True, exist_ok=True)
    gmd_tables.mkdir(parents=True, exist_ok=True)
    selected = [] if args.rescan else _read_selected_table(selected_table, target_count)
    if not selected:
        candidates = _scan_trajectory_candidates(trajectory_path, start_s, end_s, args.min_radius_um)
        _scan_events(event_path, candidates, start_s, end_s)
        selected = _select_targets(candidates, target_count)
        if len(selected) < target_count:
            raise RuntimeError(f"only {len(selected)} final-window candidates found; expected {target_count}")
        _write_selected_table(selected, gmd_tables, start_s, end_s)
    selected_ids = {row["target_id"] for row in selected}
    trajectories = _collect_rows(trajectory_path, selected_ids)
    events = _collect_rows(event_path, selected_ids)
    _validate_partner_radius(events)
    _attach_partner_heights(event_path, events)
    batch_size = max(1, args.targets_per_figure)
    for figure_rank, start_index in enumerate(range(0, min(len(selected), target_count), batch_size), start=1):
        if figure_rank > args.max_figures:
            break
        _plot_target_batch(selected[start_index : start_index + batch_size], trajectories, events, gmd_figures, figure_rank, start_s, end_s)


if __name__ == "__main__":
    main()
