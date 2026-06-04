#!/usr/bin/env python3
"""Plot short 2D sampling-procedure verification diagnostics."""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float  # noqa: E402

MODE_ORDER = ("random", "stratified")
MODE_LABELS = {
    "random": "random\nall valid SD",
    "stratified": "stratified\n400-800 m",
}
MODE_COLORS = {
    "random": PALETTE[0],
    "stratified": PALETTE[1],
}


def _partial_rows(rows: list[dict]) -> list[dict]:
    """Return sampled cases only."""
    return [row for row in rows if safe_float(row.get("sampling_fraction")) is not None and (safe_float(row.get("sampling_fraction")) or 0.0) < 1.0]


def _mode_label(mode: str) -> str:
    """Return a label that makes each sampling design's reference population explicit."""
    return MODE_LABELS.get(mode, mode)


def _finite_values(rows: list[dict], key: str, scale: float = 1.0) -> list[float]:
    """Return finite numeric values for one table column."""
    values: list[float] = []
    for row in rows:
        value = safe_float(row.get(key))
        if value is not None:
            values.append(value * scale)
    return values


def _fractions(rows: list[dict], mode: str, include_full: bool = False) -> list[float]:
    """Return sorted sampling fractions for one sampling design."""
    output = set()
    for row in rows:
        if row.get("sample_mode") != mode:
            continue
        fraction = safe_float(row.get("sampling_fraction"))
        if fraction is None:
            continue
        if include_full or fraction < 1.0:
            output.add(fraction)
    return sorted(output)


def _mean_for_fraction(rows: list[dict], mode: str, fraction: float, key: str, scale: float = 1.0) -> float | None:
    """Return a mean metric for a sampling design and fraction."""
    values = _finite_values(
        [row for row in rows if row.get("sample_mode") == mode and safe_float(row.get("sampling_fraction")) == fraction],
        key,
        scale,
    )
    return sum(values) / len(values) if values else None


def _mean_and_std(values: list[float]) -> tuple[float | None, float | None]:
    """Return mean and seed-to-seed sample standard deviation."""
    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return None, None
    mean = sum(finite_values) / len(finite_values)
    if len(finite_values) == 1:
        return mean, 0.0
    variance = sum((value - mean) ** 2 for value in finite_values) / (len(finite_values) - 1)
    return mean, math.sqrt(variance)


def _reference_lookup(rows: list[dict], key: str, scale: float = 1.0) -> dict[tuple[str, float | None], float]:
    """Return full-reference values keyed by sampling design and model time."""
    lookup: dict[tuple[str, float | None], float] = {}
    for row in rows:
        if safe_float(row.get("sampling_fraction")) != 1.0:
            continue
        value = safe_float(row.get(key))
        if value is None:
            continue
        lookup[(str(row.get("sample_mode")), safe_float(row.get("time_s")))] = value * scale
    return lookup


def _seed_level_values(
    rows: list[dict],
    mode: str,
    fraction: float,
    key: str,
    scale: float = 1.0,
    reference_error: bool = False,
) -> list[float]:
    """Return one value per seed after averaging over output times."""
    reference = _reference_lookup(rows, key, scale) if reference_error else {}
    by_seed: dict[str, list[float]] = defaultdict(list)
    for row in _partial_rows(rows):
        if row.get("sample_mode") != mode or safe_float(row.get("sampling_fraction")) != fraction:
            continue
        seed = str(row.get("seed", ""))
        if seed in ("", "NA"):
            continue
        value = safe_float(row.get(key))
        if value is None:
            continue
        metric = value * scale
        if reference_error:
            ref_value = reference.get((mode, safe_float(row.get("reference_time_s"))))
            if ref_value is None:
                continue
            metric = abs(metric - ref_value)
        by_seed[seed].append(metric)
    seed_values = [sum(values) / len(values) for seed, values in sorted(by_seed.items()) if values]
    return seed_values


def _mean_absolute_reference_error(rows: list[dict], mode: str, fraction: float, key: str, scale: float = 1.0) -> float | None:
    """Return mean absolute error against the same-design, same-time full reference."""
    reference = _reference_lookup(rows, key, scale)
    errors: list[float] = []
    for row in _partial_rows(rows):
        if row.get("sample_mode") != mode or safe_float(row.get("sampling_fraction")) != fraction:
            continue
        value = safe_float(row.get(key))
        ref_value = reference.get((mode, safe_float(row.get("reference_time_s"))))
        if value is None or ref_value is None:
            continue
        errors.append(abs(value * scale - ref_value))
    return sum(errors) / len(errors) if errors else None


def _full_reference_mean(rows: list[dict], mode: str, key: str, scale: float = 1.0) -> float | None:
    """Return the full-reference mean for one sampling design."""
    values = _finite_values(
        [row for row in rows if row.get("sample_mode") == mode and safe_float(row.get("sampling_fraction")) == 1.0],
        key,
        scale,
    )
    return sum(values) / len(values) if values else None


def _add_panel_label(ax, label: str) -> None:
    """Add a compact panel label."""
    ax.text(-0.08, 1.04, f"({label})", transform=ax.transAxes, ha="left", va="bottom", fontsize=9, fontweight="bold")


def _overlay_seed_points(ax, grouped: list[list[float]], colors: list[str]) -> None:
    """Overlay all seed-level values on a boxplot with deterministic jitter."""
    for index, (values, color) in enumerate(zip(grouped, colors), start=1):
        if not values:
            continue
        if len(values) == 1:
            offsets = [0.0]
        else:
            step = 0.18 / (len(values) - 1)
            offsets = [-0.09 + step * i for i in range(len(values))]
        ax.scatter(
            [index + offset for offset in offsets],
            values,
            s=10,
            color=color,
            edgecolors="black",
            linewidths=0.25,
            alpha=0.85,
            zorder=3,
        )


def _boxplot(outdir: Path, rows: list[dict]) -> None:
    """Plot multiplicity-weighted radius-distribution L1 errors."""
    grouped: dict[str, list[float]] = defaultdict(list)
    colors_by_label: dict[str, str] = {}
    for mode in MODE_ORDER:
        for fraction in _fractions(rows, mode):
            values = _seed_level_values(rows, mode, fraction, "dsd_l1_error_vs_full_reference")
            if values:
                label = f"{_mode_label(mode)}\nf={fraction:g}"
                grouped[label].extend(values)
                colors_by_label[label] = MODE_COLORS[mode]
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("double", 0.38))
    if not grouped:
        no_data_panel(ax, "Sampling error boxplot")
    else:
        labels = list(grouped)
        values = [grouped[label] for label in labels]
        colors = [colors_by_label[label] for label in labels]
        box = ax.boxplot(values, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, label in zip(box["boxes"], labels):
            patch.set_facecolor(colors_by_label[label])
            patch.set_alpha(0.45)
        _overlay_seed_points(ax, values, colors)
        ax.set_yscale("log")
        ax.set_ylabel("Seed-mean weighted L1 (-)")
        ax.set_title("Within-design sampling error")
        ax.tick_params(axis="x", rotation=35)
    save_figure(fig, outdir / "figures", "03_sampling_error_boxplot")
    plt.close(fig)


def _error_vs_fraction(outdir: Path, rows: list[dict]) -> None:
    """Plot mean multiplicity-weighted L1 error versus sampling fraction."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    if not _partial_rows(rows):
        no_data_panel(ax, "Error vs fraction")
    else:
        for mode in MODE_ORDER:
            xs = _fractions(rows, mode)
            stats = [_mean_and_std(_seed_level_values(rows, mode, x, "dsd_l1_error_vs_full_reference")) for x in xs]
            ys = [mean if mean is not None else math.nan for mean, _std in stats]
            yerr = [std if std is not None else 0.0 for _mean, std in stats]
            ax.errorbar(xs, ys, yerr=yerr, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, color=MODE_COLORS[mode], label=_mode_label(mode).replace("\n", " "))
        ax.set_yscale("log")
        ax.set_xlabel("Sampling fraction (-)")
        ax.set_ylabel("Seed-mean weighted L1 (-)")
        ax.set_title("Within-design convergence")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "03_sampling_error_vs_fraction")
    plt.close(fig)


def _threshold_error(outdir: Path, rows: list[dict]) -> None:
    """Plot radius-threshold fraction errors by sampling fraction and mode."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    fractions = [safe_float(row.get("sampling_fraction")) for row in _partial_rows(rows)]
    values = [safe_float(row.get("fraction_r_ge_15um")) for row in _partial_rows(rows)]
    if not has_any_number(fractions) or not has_any_number(values):
        no_data_panel(ax, "Radius-threshold diagnostics")
    else:
        for mode in MODE_ORDER:
            xs = _fractions(rows, mode)
            stats = [_mean_and_std(_seed_level_values(rows, mode, x, "fraction_r_ge_15um", reference_error=True)) for x in xs]
            ys = [mean if mean is not None else math.nan for mean, _std in stats]
            yerr = [std if std is not None else 0.0 for _mean, std in stats]
            if xs:
                ax.errorbar(xs, ys, yerr=yerr, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, color=MODE_COLORS[mode], label=_mode_label(mode).replace("\n", " "))
        ax.set_yscale("log")
        ax.set_xlabel("Sampling fraction (-)")
        ax.set_ylabel("Seed-mean tail-fraction error (-)")
        ax.set_title("Threshold-tail representativeness")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "03_sampling_radius_threshold_error")
    plt.close(fig)


def _dsd_example(outdir: Path, rows: list[dict]) -> None:
    """Plot a compact weighted-radius proxy using median radius versus fraction."""
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("single", 0.72))
    values = [safe_float(row.get("weighted_median_radius")) for row in rows]
    if not has_any_number(values):
        no_data_panel(ax, "Weighted median radius")
    else:
        for mode in MODE_ORDER:
            xs = _fractions(rows, mode)
            if xs:
                stats = [_mean_and_std(_seed_level_values(rows, mode, x, "weighted_median_radius", 1.0e6, reference_error=True)) for x in xs]
                ys = [mean if mean is not None else math.nan for mean, _std in stats]
                yerr = [std if std is not None else 0.0 for _mean, std in stats]
                ax.errorbar(xs, ys, yerr=yerr, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, color=MODE_COLORS[mode], label=_mode_label(mode).replace("\n", " "))
        ax.set_yscale("log")
        ax.set_xlabel("Sampling fraction (-)")
        ax.set_ylabel("Seed-mean median-radius error (µm)")
        ax.set_title("Median-radius representativeness")
        ax.legend(frameon=False)
    save_figure(fig, outdir / "figures", "03_sampling_dsd_example")
    plt.close(fig)


def _combined_verification(outdir: Path, rows: list[dict]) -> None:
    """Plot a four-panel 2D sampling-procedure verification summary."""
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=figure_size("double", 0.72), constrained_layout=True)

    ax = axes[0, 0]
    labels: list[str] = []
    grouped: list[list[float]] = []
    colors: list[str] = []
    for mode in MODE_ORDER:
        for fraction in _fractions(rows, mode):
            values = _seed_level_values(rows, mode, fraction, "dsd_l1_error_vs_full_reference")
            if values:
                labels.append(f"{mode}\n{fraction:g}")
                grouped.append(values)
                colors.append(MODE_COLORS[mode])
    if grouped:
        box = ax.boxplot(grouped, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.45)
        _overlay_seed_points(ax, grouped, colors)
        ax.set_yscale("log")
        ax.set_ylabel("Seed-mean weighted L1 (-)")
        ax.tick_params(axis="x", rotation=35)
    else:
        no_data_panel(ax, "L1 error")
    ax.set_title("Within-design DSD L1")
    _add_panel_label(ax, "a")

    ax = axes[0, 1]
    for mode in MODE_ORDER:
        xs = _fractions(rows, mode)
        stats = [_mean_and_std(_seed_level_values(rows, mode, fraction, "dsd_l1_error_vs_full_reference")) for fraction in xs]
        ys = [mean if mean is not None else math.nan for mean, _std in stats]
        yerr = [std if std is not None else 0.0 for _mean, std in stats]
        ax.errorbar(xs, ys, yerr=yerr, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, color=MODE_COLORS[mode], label=_mode_label(mode).replace("\n", " "))
    ax.set_yscale("log")
    ax.set_xlabel("Sampling fraction (-)")
    ax.set_ylabel("Seed-mean weighted L1 (-)")
    ax.set_title("Convergence to own reference")
    ax.legend(frameon=False)
    _add_panel_label(ax, "b")

    ax = axes[1, 0]
    for mode in MODE_ORDER:
        xs = _fractions(rows, mode)
        stats = [_mean_and_std(_seed_level_values(rows, mode, fraction, "fraction_r_ge_15um", reference_error=True)) for fraction in xs]
        ys = [mean if mean is not None else math.nan for mean, _std in stats]
        yerr = [std if std is not None else 0.0 for _mean, std in stats]
        ax.errorbar(xs, ys, yerr=yerr, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, color=MODE_COLORS[mode], label=_mode_label(mode).replace("\n", " "))
    ax.set_yscale("log")
    ax.set_xlabel("Sampling fraction (-)")
    ax.set_ylabel("Seed-mean tail-fraction error (-)")
    ax.set_title("Threshold-tail representativeness")
    _add_panel_label(ax, "c")

    ax = axes[1, 1]
    for mode in MODE_ORDER:
        xs = _fractions(rows, mode)
        stats = [_mean_and_std(_seed_level_values(rows, mode, fraction, "weighted_median_radius", 1.0e6, reference_error=True)) for fraction in xs]
        ys = [mean if mean is not None else math.nan for mean, _std in stats]
        yerr = [std if std is not None else 0.0 for _mean, std in stats]
        ax.errorbar(xs, ys, yerr=yerr, marker="o", markersize=3.2, capsize=2.5, elinewidth=0.8, capthick=0.8, color=MODE_COLORS[mode], label=_mode_label(mode).replace("\n", " "))
    ax.set_yscale("log")
    ax.set_xlabel("Sampling fraction (-)")
    ax.set_ylabel("Seed-mean median-radius error (µm)")
    ax.set_title("Median-radius representativeness")
    _add_panel_label(ax, "d")

    save_figure(fig, outdir / "figures", "03_sampling_procedure_verification")
    plt.close(fig)


def plot(outdir: Path) -> None:
    """Create all group-03 figures."""
    rows = read_table(outdir / "tables" / "03_sampling_metrics.csv")
    _boxplot(outdir, rows)
    _error_vs_fraction(outdir, rows)
    _threshold_error(outdir, rows)
    _dsd_example(outdir, rows)
    _combined_verification(outdir, rows)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["figures/03_sampling_*.{pdf,svg,png}"])
        return
    ensure_output_dirs(outdir)
    plot(outdir)


if __name__ == "__main__":
    main()
