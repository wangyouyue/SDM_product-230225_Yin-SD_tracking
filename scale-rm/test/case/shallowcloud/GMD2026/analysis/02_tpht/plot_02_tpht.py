#!/usr/bin/env python3
"""Plot GMD2026 TPHT workflow, cost, storage, IDs, and target diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, label_bars, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402


MIB = 1024.0**2


def _one_row(outdir: Path, name: str) -> dict:
    """Read the first row from a table."""
    rows = read_table(outdir / "tables" / name)
    return rows[0] if rows else {}


def _workflow(outdir: Path) -> None:
    """Plot the TPHT workflow schematic."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("double", 0.30))
    labels = ["FW discovery", "raw IDs", "merge/dedup", "BW reconstruction"]
    x_positions = [0.08, 0.36, 0.64, 0.92]
    for index, (x_pos, label) in enumerate(zip(x_positions, labels)):
        ax.text(x_pos, 0.55, label, ha="center", va="center", bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": PALETTE[index], "lw": 1.0})
    for left, right in zip(x_positions[:-1], x_positions[1:]):
        ax.annotate("", xy=(right - 0.09, 0.55), xytext=(left + 0.09, 0.55), arrowprops={"arrowstyle": "->", "lw": 1.0})
    ax.text(0.5, 0.18, "Two-Pass Hybrid Tracking: interest-restricted backward reconstruction", ha="center", fontsize=8)
    ax.axis("off")
    save_figure(fig, outdir / "figures", "02_tpht_workflow")
    plt.close(fig)


def _cost(outdir: Path, row: dict) -> None:
    """Plot TPHT cost as core-hours only."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    core_values = [safe_float(row.get("fw_core_hours")), safe_float(row.get("bw_core_hours")), safe_float(row.get("total_tpht_core_hours"))]
    if not has_any_number(core_values):
        no_data_panel(ax, "Core-hour cost")
    else:
        bars = ax.bar(range(3), [value or 0.0 for value in core_values], color=PALETTE[:3], edgecolor="black", linewidth=0.5)
        ax.set_ylim(0.0, max(value or 0.0 for value in core_values) * 1.08)
        label_bars(ax, bars)
        ax.set_xticks(range(3), ["FW", "BW", "total"], rotation=25, ha="right")
        ax.set_ylabel("Core-hours (h)")
        ax.set_title("TPHT core-hour cost")
    save_figure(fig, outdir / "figures", "02_tpht_cost")
    plt.close(fig)


def _storage(outdir: Path, row: dict) -> None:
    """Plot TPHT storage components."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    labels = ["raw .ids", "dedup .ids", "BW selected", "FW coal", "BW coal"]
    keys = ["raw_id_bytes", "dedup_id_bytes", "bw_sd_selected_output_bytes", "coalescence_log_bytes_fw", "coalescence_log_bytes_bw"]
    values = [(safe_float(row.get(key)) / MIB) if safe_float(row.get(key)) is not None else None for key in keys]
    if not has_any_number(values):
        no_data_panel(ax, "TPHT storage")
    else:
        bars = ax.bar(range(len(keys)), [value or 0.0 for value in values], color=PALETTE[: len(keys)], edgecolor="black", linewidth=0.5)
        positive = [value for value in values if value is not None and value > 0.0]
        if positive and max(positive) / min(positive) > 50.0:
            ax.set_yscale("log")
            ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
            label_bars(ax, bars, log_y=True)
        else:
            ax.set_ylim(0.0, max(positive) * 1.08 if positive else 1.0)
            label_bars(ax, bars)
        ax.set_xticks(range(len(keys)), labels, rotation=35, ha="right")
        ax.set_ylabel("Storage (MiB)")
        ax.set_title("TPHT storage")
    save_figure(fig, outdir / "figures", "02_tpht_storage")
    plt.close(fig)


def _id_counts(outdir: Path, row: dict) -> None:
    """Plot raw, unique, deduplicated, and BW ID counts."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    labels = ["raw records", "FW unique", "dedup targets", "BW unique"]
    keys = ["fw_raw_id_records", "fw_unique_pairs", "deduplicated_target_pairs", "bw_unique_pairs"]
    values = [safe_float(row.get(key)) for key in keys]
    if not has_any_number(values):
        no_data_panel(ax, "TPHT ID counts")
    else:
        bars = ax.bar(range(len(keys)), [value or 0.0 for value in values], color=PALETTE[: len(keys)], edgecolor="black", linewidth=0.5)
        positive = [value for value in values if value is not None and value > 0.0]
        if positive and max(positive) / min(positive) > 50.0:
            ax.set_yscale("log")
            ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
            label_bars(ax, bars, log_y=True)
        else:
            ax.set_ylim(0.0, max(positive) * 1.08 if positive else 1.0)
            label_bars(ax, bars)
        ax.set_xticks(range(len(keys)), labels, rotation=35, ha="right")
        ax.set_ylabel("ID count")
        ax.set_title("TPHT ID counts")
    save_figure(fig, outdir / "figures", "02_tpht_id_counts")
    plt.close(fig)


def _target_categories(outdir: Path) -> None:
    """Plot TPHT target category counts."""
    rows = [
        row
        for row in read_table(outdir / "tables" / "02_tpht_target_categories.csv")
        if (safe_float(row.get("target_count")) or 0.0) > 0.0
    ]
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.70))
    values = [safe_float(row.get("target_count")) for row in rows]
    if not has_any_number(values):
        no_data_panel(ax, "Target categories")
    else:
        bars = ax.bar(range(len(rows)), [value or 0.0 for value in values], color=[PALETTE[index % len(PALETTE)] for index in range(len(rows))], edgecolor="black", linewidth=0.5)
        positive = [value for value in values if value is not None and value > 0.0]
        if positive and max(positive) / min(positive) > 50.0:
            ax.set_yscale("log")
            ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
            label_bars(ax, bars, log_y=True)
        else:
            ax.set_ylim(0.0, max(positive) * 1.08 if positive else 1.0)
            label_bars(ax, bars)
        ax.set_xticks(range(len(rows)), [row.get("category", "") for row in rows], rotation=30, ha="right")
        ax.set_ylabel("Targets (count)")
        ax.set_title("First-selected target conditions")
    save_figure(fig, outdir / "figures", "02_tpht_target_categories")
    plt.close(fig)


def _discovery(outdir: Path) -> None:
    """Plot target discovery through time with first-hour cap shading."""
    rows = read_table(outdir / "tables" / "02_tpht_discovery_time.csv")
    times = [(safe_float(row.get("time_s")) / 60.0) if safe_float(row.get("time_s")) is not None else None for row in rows]
    new_targets = [safe_float(row.get("new_targets_per_time")) for row in rows]
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    if not has_any_number(times) or not has_any_number(new_targets):
        no_data_panel(ax, "Discovery time")
    else:
        ax.axvspan(0.0, 60.0, color="0.9", zorder=0, label="first-hour cap region")
        ax.plot([time or 0.0 for time in times], [value or 0.0 for value in new_targets], marker="o", color=PALETTE[0])
        ax.set_xlabel("Model time (min)")
        ax.set_ylabel("New targets (count)")
        ax.set_title("Target discovery time")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "02_tpht_discovery_time")
    plt.close(fig)


def _rank_load(outdir: Path) -> None:
    """Plot raw and deduplicated target load by rank."""
    rows = read_table(outdir / "tables" / "02_tpht_rank_load_balance.csv")
    ranks = [safe_float(row.get("rank")) for row in rows]
    raw = [safe_float(row.get("raw_records_per_rank")) for row in rows]
    dedup = [safe_float(row.get("dedup_targets_per_rank")) for row in rows]
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("double", 0.35))
    if not has_any_number(ranks) or not (has_any_number(raw) or has_any_number(dedup)):
        no_data_panel(ax, "Rank load balance")
    else:
        x_positions = list(range(len(rows)))
        ax.plot(x_positions, [value or 0.0 for value in raw], label="raw", color=PALETTE[0], marker="o", markersize=2)
        ax.plot(x_positions, [value or 0.0 for value in dedup], label="dedup", color=PALETTE[1], marker="s", markersize=2)
        positive = [value for value in raw + dedup if value is not None and value > 0.0]
        if positive and max(positive) / min(positive) > 50.0:
            ax.set_yscale("log")
        ax.set_xlabel("Rank index")
        ax.set_ylabel("Records or targets (count)")
        ax.set_title("Rank load balance")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "02_tpht_rank_load_balance")
    plt.close(fig)


def _chain_validity(outdir: Path) -> None:
    """Plot chain-validity fractions."""
    row = _one_row(outdir, "02_tpht_chain_validity.csv")
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.65))
    labels = ["valid predecessor", "invalid predecessor"]
    values = [safe_float(row.get("valid_predecessor_fraction")), safe_float(row.get("invalid_predecessor_fraction"))]
    if not has_any_number(values):
        count_labels = ["unique BW targets"]
        count_values = [safe_float(row.get("bw_unique_pairs", row.get("bw_valid_records")))]
        if not has_any_number(count_values):
            no_data_panel(ax, "Chain validity", "No predecessor-fraction diagnostics")
        else:
            bars = ax.bar(range(len(count_values)), [value or 0.0 for value in count_values], color=PALETTE[:1], edgecolor="black", linewidth=0.5)
            positive = [value for value in count_values if value is not None and value > 0.0]
            ax.set_yscale("log")
            ax.set_ylim(min(positive) / 3.0, max(positive) * 3.0)
            label_bars(ax, bars, log_y=True)
            ax.set_xticks(range(len(count_values)), count_labels, rotation=25, ha="right")
            ax.set_ylabel("Count")
            ax.set_title("BW reconstruction ID coverage")
    else:
        bars = ax.bar(range(2), [value or 0.0 for value in values], color=PALETTE[:2], edgecolor="black", linewidth=0.5)
        label_bars(ax, bars)
        ax.set_xticks(range(2), labels, rotation=25, ha="right")
        ax.set_ylabel("Fraction (-)")
        ax.set_ylim(0, 1)
        ax.set_title("Chain validity")
    save_figure(fig, outdir / "figures", "02_tpht_chain_validity")
    plt.close(fig)


def _target_histories(outdir: Path) -> None:
    """Plot median radius and height histories when available."""
    rows = read_table(outdir / "tables" / "02_tpht_target_histories.csv")
    times = [(safe_float(row.get("time_s")) / 60.0) if safe_float(row.get("time_s")) is not None else None for row in rows]
    radius = [safe_float(row.get("median_radius_history_m")) for row in rows]
    if not has_any_number(radius):
        radius = [safe_float(row.get("mean_radius_history_m")) for row in rows]
    height = [safe_float(row.get("median_height_history_m")) for row in rows]
    if not has_any_number(height):
        height = [safe_float(row.get("mean_height_history_m")) for row in rows]
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    if not has_any_number(times) or not (has_any_number(radius) or has_any_number(height)):
        row = rows[0] if rows else {}
        radius_median = safe_float(row.get("median_radius_history_m"))
        radius_q25 = safe_float(row.get("radius_q25_m"))
        radius_q75 = safe_float(row.get("radius_q75_m"))
        height_median = safe_float(row.get("median_height_history_m"))
        height_q25 = safe_float(row.get("height_q25_m"))
        height_q75 = safe_float(row.get("height_q75_m"))
        aggregate_values = [radius_median, radius_q25, radius_q75, height_median, height_q25, height_q75]
        if not has_any_number(aggregate_values):
            no_data_panel(ax, "Target histories", "No time-resolved target-history diagnostics")
        else:
            radius_um = (radius_median or 0.0) * 1.0e6
            radius_err_low = max(0.0, radius_um - (radius_q25 or radius_median or 0.0) * 1.0e6)
            radius_err_high = max(0.0, (radius_q75 or radius_median or 0.0) * 1.0e6 - radius_um)
            height_value = height_median or 0.0
            height_err_low = max(0.0, height_value - (height_q25 or height_median or 0.0))
            height_err_high = max(0.0, (height_q75 or height_median or 0.0) - height_value)
            ax.errorbar([0], [radius_um], yerr=[[radius_err_low], [radius_err_high]], fmt="o", color=PALETTE[0], capsize=3)
            ax.set_xticks([0], ["radius"])
            ax.set_ylabel("Radius (um)")
            ax.set_title("Aggregate target-history radius")
            twin = ax.twinx()
            twin.errorbar([0.18], [height_value], yerr=[[height_err_low], [height_err_high]], fmt="s", color=PALETTE[1], capsize=3)
            twin.set_ylabel("Height (m)")
            twin.tick_params(direction="in")
    else:
        if has_any_number(radius):
            ax.plot([time or 0.0 for time in times], [value or 0.0 for value in radius], color=PALETTE[0], label="radius (m)")
        if has_any_number(height):
            ax.plot([time or 0.0 for time in times], [value or 0.0 for value in height], color=PALETTE[1], label="height (m)")
        ax.set_xlabel("Model time (min)")
        ax.set_ylabel("Diagnostic value (m)")
        ax.set_title("Stepwise BW target histories")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "02_tpht_target_histories")
    plt.close(fig)


def plot(outdir: Path) -> None:
    """Create all TPHT figures."""
    row = _one_row(outdir, "02_tpht_summary.csv")
    _workflow(outdir)
    _cost(outdir, row)
    _storage(outdir, row)
    _id_counts(outdir, row)
    _target_categories(outdir)
    _discovery(outdir)
    _rank_load(outdir)
    _chain_validity(outdir)
    _target_histories(outdir)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/02_tpht_*.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
