#!/usr/bin/env python3
"""Plot compact TPHT predecessor-tree examples from BW selected output."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

from common.parse_logs import collect_case_metrics  # noqa: E402
from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.plot_style import OKABE_ITO, configure_matplotlib, figure_size, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import read_csv_rows, safe_float, write_table_bundle  # noqa: E402
from tpht_chain_utils import TREE_COLUMNS, build_predecessor_tree_examples  # noqa: E402


FIGURE_STEM = "02_tpht_predecessor_tree_examples"
TABLE_STEM = "02_tpht_predecessor_tree_examples"


def _target_sort_key(row: dict[str, str]) -> tuple[float, float, float, float]:
    """Rank targets by available coalescence-linked and size diagnostics."""
    return (
        safe_float(row.get("coal_event_count")) or 0.0,
        safe_float(row.get("coal_episode_count_proxy")) or 0.0,
        safe_float(row.get("max_radius_um")) or 0.0,
        safe_float(row.get("record_count")) or 0.0,
    )


def choose_target_ids(outdir: Path, max_examples: int) -> tuple[set[str], list[str]]:
    """Choose a few target chains for tree examples from cached TPHT summaries."""
    rows = read_csv_rows(outdir / "tables" / "02_tpht_science_target_summary.csv")
    warnings: list[str] = []
    candidates = [row for row in rows if row.get("target_id") not in (None, "", "NA")]
    if not candidates:
        return set(), ["02_tpht_science_target_summary.csv is missing or has no target_id values"]
    ranked = sorted(candidates, key=_target_sort_key, reverse=True)
    selected = {str(row["target_id"]) for row in ranked[: max(1, max_examples)]}
    if all((safe_float(row.get("coal_event_count")) or 0.0) <= 0.0 for row in ranked[: max(1, max_examples)]):
        warnings.append("selected predecessor-tree examples have no linked coalescence events; plot will show target histories only")
    return selected, warnings


def load_or_generate_tree_rows(root: Path, outdir: Path, options: dict[str, Any], strict: bool, max_examples: int) -> list[dict[str, Any]]:
    """Read cached tree rows or generate a compact source table from BW output."""
    table_path = outdir / "tables" / f"{TABLE_STEM}.csv"
    cached_rows = read_table(table_path)
    cached_examples = {row.get("example_id") for row in cached_rows if row.get("node_id") not in (None, "", "NA")}
    cache_has_partner_radius = all("partner_radius_um" in row for row in cached_rows) if cached_rows else False
    if cached_rows and len(cached_examples) >= max_examples and cache_has_partner_radius:
        return cached_rows

    config = load_config()
    bw_case = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}["bw_reconstruction"]
    target_ids, target_warnings = choose_target_ids(outdir, max_examples)
    _bw_row, bw_files, _records, metric_warnings = collect_case_metrics(root, bw_case, options)
    rows, tree_warnings = build_predecessor_tree_examples(
        bw_files.get("sd_selected_output_bytes", []),
        bw_files.get("coalescence_log_bytes", []),
        target_ids,
        options,
    )
    warnings = target_warnings + metric_warnings + tree_warnings
    if not rows:
        rows = [{column: None for column in TREE_COLUMNS}]
        rows[0]["warnings"] = warning_text(warnings)
    else:
        if warnings:
            rows[0]["warnings"] = warning_text(warnings)
    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle(rows, outdir / "tables" / TABLE_STEM, TREE_COLUMNS)
    return rows


def _group_rows(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Group valid tree rows by example target."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        example_id = row.get("example_id")
        node_id = row.get("node_id")
        if example_id in (None, "", "NA") or node_id in (None, "", "NA"):
            continue
        grouped.setdefault(str(example_id), []).append(row)
    return grouped


def _finite(values: list[float | None]) -> list[float]:
    """Return finite values."""
    return [value for value in values if value is not None and math.isfinite(value)]


def _nan_if_none(value: float | None) -> float:
    """Convert optional finite values to NaN for Matplotlib."""
    return value if value is not None and math.isfinite(value) else math.nan


def _plot_tree_batch(plt: Any, grouped: dict[str, list[dict[str, Any]]], example_ids: list[str], outdir: Path, stems: list[str]) -> None:
    """Draw one batch of three predecessor-tree proxy examples."""
    fig, axes = plt.subplots(len(example_ids), 1, figsize=figure_size("double", 0.34 * len(example_ids)), constrained_layout=True, sharex=True)
    if len(example_ids) == 1:
        axes = [axes]
    all_radii = _finite([safe_float(row.get("radius_um")) for example_rows in grouped.values() for row in example_rows])
    radius_min = min(all_radii) if all_radii else 0.0
    radius_max = max(all_radii) if all_radii else 1.0
    if radius_max <= radius_min:
        radius_max = radius_min + 1.0
    cmap = plt.get_cmap("cividis")
    norm = plt.Normalize(radius_min, radius_max)

    scatter = None
    for panel_index, (ax, example_id) in enumerate(zip(axes, example_ids)):
        example_rows = sorted(grouped[example_id], key=lambda row: (safe_float(row.get("time_s")) or -math.inf, str(row.get("node_id"))))
        node_by_id = {str(row.get("node_id")): row for row in example_rows if row.get("node_id") not in (None, "", "NA")}
        trunk_rows = [row for row in example_rows if row.get("node_role") == "tracked_target_output"]
        trunk_rows.sort(key=lambda row: safe_float(row.get("time_s")) or -math.inf)
        trunk_x = [_nan_if_none(safe_float(row.get("time_s"))) / 60.0 for row in trunk_rows]
        trunk_y = [_nan_if_none(safe_float(row.get("height_m"))) for row in trunk_rows]
        if len(trunk_rows) >= 2:
            ax.plot(trunk_x, trunk_y, color="0.45", linewidth=0.8, zorder=1)
        important_rows: list[dict[str, Any]] = []
        if trunk_rows:
            threshold_rows = [row for row in trunk_rows if (safe_float(row.get("radius_um")) or -math.inf) >= 15.0]
            if threshold_rows:
                important_rows.append(threshold_rows[0])
            important_rows.append(max(trunk_rows, key=lambda row: safe_float(row.get("radius_um")) or -math.inf))
            important_rows.extend([row for row in trunk_rows if safe_float(row.get("if_coal")) and (safe_float(row.get("if_coal")) or 0.0) > 0.0])
        seen_nodes: set[str] = set()
        important_rows = [row for row in important_rows if not (str(row.get("node_id")) in seen_nodes or seen_nodes.add(str(row.get("node_id"))))]
        if important_rows:
            scatter = ax.scatter(
                [_nan_if_none(safe_float(row.get("time_s"))) / 60.0 for row in important_rows],
                [_nan_if_none(safe_float(row.get("height_m"))) for row in important_rows],
                c=[_nan_if_none(safe_float(row.get("radius_um"))) for row in important_rows],
                cmap=cmap,
                norm=norm,
                s=28,
                marker="o",
                edgecolors="none",
                linewidths=0.0,
                zorder=3,
                label="important output node",
            )
        ifcoal_rows = [row for row in trunk_rows if safe_float(row.get("if_coal")) and (safe_float(row.get("if_coal")) or 0.0) > 0.0]
        if ifcoal_rows:
            ax.scatter(
                [_nan_if_none(safe_float(row.get("time_s"))) / 60.0 for row in ifcoal_rows],
                [_nan_if_none(safe_float(row.get("height_m"))) for row in ifcoal_rows],
                marker="x",
                color=OKABE_ITO["vermillion"],
                s=28,
                linewidths=0.9,
                zorder=4,
                label="x: if_coal flag",
            )

        branch_rows = [row for row in example_rows if row.get("node_role") == "coalescence_partner_proxy"]
        for row in branch_rows:
            predecessor = node_by_id.get(str(row.get("predecessor_id")))
            x0 = safe_float(row.get("time_s"))
            y0 = safe_float(row.get("height_m"))
            radius = safe_float(row.get("radius_um"))
            if predecessor is None or x0 is None or y0 is None:
                continue
            x1 = safe_float(predecessor.get("time_s"))
            y1 = safe_float(predecessor.get("height_m"))
            if x1 is None or y1 is None:
                continue
            color = cmap(norm(radius if radius is not None else radius_min))
            ax.plot([x0 / 60.0, x1 / 60.0], [y0, y1], color=color, linewidth=0.8, alpha=0.75, zorder=0)
            ax.scatter(
                [x0 / 60.0],
                [y0],
                marker="^",
                s=34,
                color=color,
                edgecolors="none",
                linewidths=0.0,
                zorder=4,
                label="triangle: linked partner radius",
            )

        ax.set_title(f"{chr(97 + panel_index)}  target {example_id}", loc="left", fontsize=9, pad=2)
        ax.set_ylabel("Height (m)")
        ax.tick_params(direction="in")
    axes[-1].set_xlabel("Time (min)")
    handles: list[Any] = []
    labels: list[str] = []
    for ax in axes:
        axis_handles, axis_labels = ax.get_legend_handles_labels()
        handles.extend(axis_handles)
        labels.extend(axis_labels)
    unique = dict(zip(labels, handles))
    if unique:
        axes[0].legend(unique.values(), unique.keys(), frameon=False, fontsize=7, loc="upper right")
    if scatter is not None:
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(mappable, ax=axes, location="right", shrink=0.88, pad=0.02)
        cbar.set_label(r"Radius ($\mu$m)")
    for stem in stems:
        save_figure(fig, outdir / "figures", stem)
    plt.close(fig)


def plot_tree_rows(rows: list[dict[str, Any]], outdir: Path) -> None:
    """Draw predecessor-tree proxy examples in batches of three targets."""
    grouped = _group_rows(rows)
    plt = configure_matplotlib()
    if not grouped:
        fig, ax = plt.subplots(figsize=figure_size("single", 0.75), constrained_layout=True)
        no_data_panel(ax, "TPHT predecessor-tree examples")
        save_figure(fig, outdir / "figures", FIGURE_STEM)
        plt.close(fig)
        return

    example_ids = list(grouped)
    batch_size = 3
    for batch_start in range(0, len(example_ids), batch_size):
        batch_index = batch_start // batch_size + 1
        batch_examples = example_ids[batch_start : batch_start + batch_size]
        stems = [f"{FIGURE_STEM}_{batch_index:02d}"]
        if batch_index == 1:
            stems.append(FIGURE_STEM)
        _plot_tree_batch(plt, grouped, batch_examples, outdir, stems)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    parser.add_argument("--max-examples", type=int, default=12, help="Maximum number of target chains to plot.")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(
            Path(__file__).name,
            root,
            outdir,
            [
                "tables/02_tpht_predecessor_tree_examples.{csv,md,tex,json}",
                "figures/02_tpht_predecessor_tree_examples.{pdf,svg,png}",
                "figures/02_tpht_predecessor_tree_examples_*.{pdf,svg,png}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    rows = load_or_generate_tree_rows(root, outdir, analysis_options(args), args.strict, args.max_examples)
    plot_tree_rows(rows, outdir)


if __name__ == "__main__":
    main()
