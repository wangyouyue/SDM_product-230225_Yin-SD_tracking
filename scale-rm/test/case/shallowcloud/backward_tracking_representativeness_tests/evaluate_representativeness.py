from __future__ import annotations

"""Evaluate SD backward-tracking representativeness and export metrics/figures."""

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class CaseMeta:
    case_name: str
    mode: str
    fraction: float
    seed_scale: int
    tracking_height_min: float
    tracking_height_max: float
    tracking_radius_min: float
    tracking_nz_bin: int
    tracking_nr_bin: int


@dataclass
class DistributionPack:
    z_edges: np.ndarray
    r_edges: np.ndarray
    zr_pdf: np.ndarray
    z_pdf: np.ndarray
    r_pdf: np.ndarray


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for metric and figure outputs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "representativeness_metrics.csv")
    parser.add_argument("--figure-dir", type=Path, default=Path(__file__).resolve().parent / "figures")
    parser.add_argument("--max-chain-samples", type=int, default=3000)
    return parser.parse_args()


def _extract_value(text: str, key: str, as_int: bool = False) -> float | int:
    match = re.search(rf"^\s*{re.escape(key)}\s*=\s*([^,\n]+)", text, flags=re.MULTILINE)
    if not match:
        return 0 if as_int else 0.0
    raw = match.group(1).strip().strip('"').strip("'")
    lower = raw.lower()
    if as_int:
        return int(float(lower.replace("d", "e")))
    return float(lower.replace("d", "e"))


def parse_case_meta(case_dir: Path) -> CaseMeta:
    """Read one case's run.conf and extract sampling metadata."""
    run_conf = (case_dir / "run.conf").read_text()
    mode_match = re.search(r'^\s*tracking_selection_mode\s*=\s*"([^"]+)"', run_conf, flags=re.MULTILINE)
    mode = mode_match.group(1).strip() if mode_match else "unknown"
    seed = int(_extract_value(run_conf, "RANDOM_SEED_SCALE", as_int=True))
    return CaseMeta(
        case_name=case_dir.name,
        mode=mode,
        fraction=float(_extract_value(run_conf, "tracking_fraction")),
        seed_scale=seed,
        tracking_height_min=float(_extract_value(run_conf, "tracking_height_min")),
        tracking_height_max=float(_extract_value(run_conf, "tracking_height_max")),
        tracking_radius_min=float(_extract_value(run_conf, "tracking_radius_min")),
        tracking_nz_bin=int(_extract_value(run_conf, "tracking_nz_bin", as_int=True)),
        tracking_nr_bin=int(_extract_value(run_conf, "tracking_nr_bin", as_int=True)),
    )


def collect_case_dirs(root: Path) -> list[Path]:
    """Collect all case directories under root."""
    return sorted([p for p in root.iterdir() if p.is_dir() and (p / "run.conf").exists()])


def parse_sd_filename(file_name: str) -> tuple[str, int] | None:
    match = re.search(r"_(\d{8}-\d{6}\.\d{3})\.pe(\d{6})$", file_name)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def collect_sd_files(case_dir: Path) -> dict[str, dict[int, Path]]:
    """Group SD NetCDF files by simulation time and MPI domain."""
    selected = sorted(case_dir.glob("SD_selected_NetCDF_*"))
    target = selected if selected else sorted(case_dir.glob("SD_all_NetCDF_*"))
    grouped: dict[str, dict[int, Path]] = {}
    for path in target:
        parsed = parse_sd_filename(path.name)
        if not parsed:
            continue
        time_key, domain = parsed
        grouped.setdefault(time_key, {})[domain] = path
    return grouped


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Compute weighted quantile for SD properties."""
    if values.size == 0:
        return float("nan")
    sorter = np.argsort(values)
    values_sorted = values[sorter]
    weights_sorted = weights[sorter]
    cdf = np.cumsum(weights_sorted)
    if cdf[-1] <= 0:
        return float("nan")
    threshold = q * cdf[-1]
    idx = int(np.searchsorted(cdf, threshold, side="left"))
    idx = max(0, min(idx, values_sorted.size - 1))
    return float(values_sorted[idx])


def _ks_distance(p: np.ndarray, q: np.ndarray) -> float:
    cdf_p = np.cumsum(p)
    cdf_q = np.cumsum(q)
    return float(np.max(np.abs(cdf_p - cdf_q)))


def _emd_1d(p: np.ndarray, q: np.ndarray, edges: np.ndarray) -> float:
    cdf_diff = np.abs(np.cumsum(p) - np.cumsum(q))
    widths = np.diff(edges)
    return float(np.sum(cdf_diff * widths))


def _safe_normalize(hist: np.ndarray) -> np.ndarray:
    total = np.sum(hist)
    if total <= 0:
        return np.zeros_like(hist, dtype=np.float64)
    return hist.astype(np.float64) / float(total)


def _read_var(ds: Dataset, name: str) -> np.ndarray | None:
    if name in ds.variables:
        return np.array(ds.variables[name][:])
    return None


def build_distribution(meta: CaseMeta, z_values: np.ndarray, r_values: np.ndarray, weights: np.ndarray) -> DistributionPack:
    """Build weighted height-radius distribution from SD samples."""
    z_min = meta.tracking_height_min
    z_max = meta.tracking_height_max
    z_n = max(1, int(meta.tracking_nz_bin))
    r_n = max(1, int(meta.tracking_nr_bin))
    r_positive = r_values[r_values > 0]
    r_min_data = float(np.min(r_positive)) if r_positive.size > 0 else 1.0e-9
    r_max_data = float(np.max(r_values)) if r_values.size > 0 else max(r_min_data * 10.0, 1.0e-9)
    r_min = max(meta.tracking_radius_min, r_min_data, 1.0e-9)
    if r_max_data <= r_min:
        r_max_data = r_min * 10.0
    z_edges = np.linspace(z_min, z_max, z_n + 1)
    r_edges = np.geomspace(r_min, r_max_data, r_n + 1)
    hist, _, _ = np.histogram2d(z_values, r_values, bins=[z_edges, r_edges], weights=weights)
    zr_pdf = _safe_normalize(hist)
    z_pdf = _safe_normalize(np.sum(hist, axis=1))
    r_pdf = _safe_normalize(np.sum(hist, axis=0))
    return DistributionPack(z_edges=z_edges, r_edges=r_edges, zr_pdf=zr_pdf, z_pdf=z_pdf, r_pdf=r_pdf)


def compute_chain_metrics(grouped_files: dict[str, dict[int, Path]], max_chain_samples: int, key_hmin: float, key_hmax: float) -> dict[str, float]:
    """Estimate trajectory-chain representativeness metrics."""
    if not grouped_files:
        return {
            "chain_length_mean": float("nan"),
            "chain_length_p90": float("nan"),
            "chain_length_p99": float("nan"),
            "coal_chain_fraction": float("nan"),
            "key_height_coverage_fraction": float("nan"),
        }
    time_keys = sorted(grouped_files.keys(), reverse=True)
    latest = time_keys[0]
    latest_nodes: list[tuple[int, int, float]] = []
    node_maps: list[dict[int, dict[str, np.ndarray]]] = []
    for t_idx, t_key in enumerate(time_keys):
        domain_map: dict[int, dict[str, np.ndarray]] = {}
        for domain, file_path in grouped_files[t_key].items():
            with Dataset(file_path, "r") as ds:
                pre_sdid = _read_var(ds, "pre_sdid")
                pre_dmid = _read_var(ds, "pre_dmid")
                if_coal = _read_var(ds, "if_coal")
                sd_z = _read_var(ds, "sd_z")
                sd_n = _read_var(ds, "sd_n")
                if pre_sdid is None or pre_dmid is None or if_coal is None or sd_z is None:
                    continue
                domain_map[domain] = {
                    "pre_sdid": np.asarray(pre_sdid).astype(np.int64),
                    "pre_dmid": np.asarray(pre_dmid).astype(np.int64),
                    "if_coal": np.asarray(if_coal).astype(np.int64),
                    "sd_z": np.asarray(sd_z).astype(np.float64),
                }
                if t_idx == 0:
                    weights = np.asarray(sd_n).astype(np.float64) if sd_n is not None else np.ones_like(pre_sdid, dtype=np.float64)
                    for idx in range(pre_sdid.shape[0]):
                        latest_nodes.append((domain, idx + 1, max(weights[idx], 1.0)))
        node_maps.append(domain_map)
    if not latest_nodes:
        return {
            "chain_length_mean": float("nan"),
            "chain_length_p90": float("nan"),
            "chain_length_p99": float("nan"),
            "coal_chain_fraction": float("nan"),
            "key_height_coverage_fraction": float("nan"),
        }
    choose_n = min(max_chain_samples, len(latest_nodes))
    weights = np.array([n[2] for n in latest_nodes], dtype=np.float64)
    probs = weights / np.sum(weights)
    rng = np.random.default_rng(42)
    selected_idx = rng.choice(len(latest_nodes), size=choose_n, replace=False, p=probs)
    lengths: list[int] = []
    coal_flags: list[int] = []
    coverage_flags: list[int] = []
    for idx in selected_idx:
        cur_domain, cur_sdid, _ = latest_nodes[int(idx)]
        chain_len = 0
        has_coal = 0
        covered = 0
        for t_idx, domain_map in enumerate(node_maps):
            if cur_domain not in domain_map:
                break
            rec = domain_map[cur_domain]
            arr_idx = int(cur_sdid) - 1
            if arr_idx < 0 or arr_idx >= rec["pre_sdid"].shape[0]:
                break
            chain_len += 1
            z_val = rec["sd_z"][arr_idx]
            if key_hmin <= z_val <= key_hmax:
                covered = 1
            if rec["if_coal"][arr_idx] > 0:
                has_coal = 1
            next_domain = int(rec["pre_dmid"][arr_idx])
            next_sdid = int(rec["pre_sdid"][arr_idx])
            if next_domain <= 0 or next_sdid <= 0:
                break
            cur_domain = next_domain
            cur_sdid = next_sdid
            if t_idx == len(node_maps) - 1:
                break
        lengths.append(chain_len)
        coal_flags.append(has_coal)
        coverage_flags.append(covered)
    arr_len = np.array(lengths, dtype=np.float64)
    return {
        "chain_length_mean": float(np.mean(arr_len)),
        "chain_length_p90": float(np.percentile(arr_len, 90)),
        "chain_length_p99": float(np.percentile(arr_len, 99)),
        "coal_chain_fraction": float(np.mean(np.array(coal_flags, dtype=np.float64))),
        "key_height_coverage_fraction": float(np.mean(np.array(coverage_flags, dtype=np.float64))),
    }


def compute_case_metrics(case_dir: Path, max_chain_samples: int) -> tuple[CaseMeta, dict[str, float], DistributionPack]:
    """Compute all scalar metrics for one test case."""
    meta = parse_case_meta(case_dir)
    grouped = collect_sd_files(case_dir)
    if not grouped:
        metrics = {
            "files_count": 0.0,
            "particle_count": 0.0,
            "weighted_nd_sum": float("nan"),
            "qc_proxy_mass": float("nan"),
            "radius_m0": float("nan"),
            "radius_m1": float("nan"),
            "radius_m2": float("nan"),
            "radius_m3": float("nan"),
            "radius_p90": float("nan"),
            "radius_p99": float("nan"),
            "aerosol_weighted_sum": float("nan"),
            "aerosol_positive_weighted_frac": float("nan"),
            "chain_length_mean": float("nan"),
            "chain_length_p90": float("nan"),
            "chain_length_p99": float("nan"),
            "coal_chain_fraction": float("nan"),
            "key_height_coverage_fraction": float("nan"),
        }
        empty = np.array([], dtype=np.float64)
        z_edges = np.linspace(meta.tracking_height_min, meta.tracking_height_max, max(meta.tracking_nz_bin, 1) + 1)
        r_edges = np.geomspace(1.0e-9, 1.0e-8, max(meta.tracking_nr_bin, 1) + 1)
        pack = DistributionPack(z_edges, r_edges, np.zeros((len(z_edges) - 1, len(r_edges) - 1)), empty, empty)
        return meta, metrics, pack
    latest_time = sorted(grouped.keys())[-1]
    files = list(grouped[latest_time].values())
    sd_r_list: list[np.ndarray] = []
    sd_n_list: list[np.ndarray] = []
    sd_z_list: list[np.ndarray] = []
    sd_asl_sum_list: list[np.ndarray] = []
    for file_path in files:
        with Dataset(file_path, "r") as ds:
            sd_r = np.array(ds.variables["sd_r"][:], dtype=np.float64)
            sd_n = np.array(ds.variables["sd_n"][:], dtype=np.float64)
            sd_z = np.array(ds.variables["sd_z"][:], dtype=np.float64)
            sd_r_list.append(sd_r)
            sd_n_list.append(sd_n)
            sd_z_list.append(sd_z)
            if "sd_asl" in ds.variables:
                sd_asl = np.array(ds.variables["sd_asl"][:], dtype=np.float64)
                if sd_asl.ndim == 2:
                    sd_asl_sum = np.sum(sd_asl, axis=0)
                else:
                    sd_asl_sum = sd_asl
                sd_asl_sum_list.append(sd_asl_sum)
            else:
                sd_asl_sum_list.append(np.zeros_like(sd_r))
    sd_r_all = np.concatenate(sd_r_list) if sd_r_list else np.array([], dtype=np.float64)
    sd_n_all = np.concatenate(sd_n_list) if sd_n_list else np.array([], dtype=np.float64)
    sd_z_all = np.concatenate(sd_z_list) if sd_z_list else np.array([], dtype=np.float64)
    sd_asl_all = np.concatenate(sd_asl_sum_list) if sd_asl_sum_list else np.array([], dtype=np.float64)
    weight_sum = float(np.sum(sd_n_all)) if sd_n_all.size > 0 else float("nan")
    radius_moments = {}
    for order in range(4):
        if sd_r_all.size == 0 or sd_n_all.size == 0 or np.sum(sd_n_all) <= 0:
            radius_moments[f"radius_m{order}"] = float("nan")
        else:
            radius_moments[f"radius_m{order}"] = float(np.sum(sd_n_all * np.power(sd_r_all, order)) / np.sum(sd_n_all))
    if sd_asl_all.size > 0 and sd_n_all.size > 0 and np.sum(sd_n_all) > 0:
        aerosol_weighted_sum = float(np.sum(sd_asl_all * sd_n_all))
        aerosol_positive_weighted_frac = float(np.sum(sd_n_all[sd_asl_all > 0]) / np.sum(sd_n_all))
    else:
        aerosol_weighted_sum = float("nan")
        aerosol_positive_weighted_frac = float("nan")
    if sd_r_all.size > 0 and sd_n_all.size > 0:
        p90 = weighted_quantile(sd_r_all, sd_n_all, 0.90)
        p99 = weighted_quantile(sd_r_all, sd_n_all, 0.99)
    else:
        p90 = float("nan")
        p99 = float("nan")
    rho_w = 1000.0
    qc_proxy_mass = float(np.sum(sd_n_all * (4.0 / 3.0) * math.pi * rho_w * np.power(sd_r_all, 3))) if sd_r_all.size > 0 else float("nan")
    dist = build_distribution(meta, sd_z_all, sd_r_all, sd_n_all if sd_n_all.size > 0 else np.ones_like(sd_r_all))
    chain = compute_chain_metrics(grouped, max_chain_samples, meta.tracking_height_min, meta.tracking_height_max)
    metrics = {
        "files_count": float(len(files)),
        "particle_count": float(sd_r_all.size),
        "weighted_nd_sum": weight_sum,
        "qc_proxy_mass": qc_proxy_mass,
        "radius_p90": p90,
        "radius_p99": p99,
        "aerosol_weighted_sum": aerosol_weighted_sum,
        "aerosol_positive_weighted_frac": aerosol_positive_weighted_frac,
    }
    metrics.update(radius_moments)
    metrics.update(chain)
    return meta, metrics, dist


def compare_to_reference(case_dist: DistributionPack, ref_dist: DistributionPack) -> dict[str, float]:
    """Compare sampled distributions against mode-matched reference."""
    eps = 1.0e-16
    n_z = min(case_dist.zr_pdf.shape[0], ref_dist.zr_pdf.shape[0])
    n_r = min(case_dist.zr_pdf.shape[1], ref_dist.zr_pdf.shape[1])
    case_zr = case_dist.zr_pdf[:n_z, :n_r]
    ref_zr = ref_dist.zr_pdf[:n_z, :n_r]
    if case_zr.size == 0 or ref_zr.size == 0:
        return {
            "zr_rel_error_mean": float("nan"),
            "zr_l1_distance": float("nan"),
            "z_ks_distance": float("nan"),
            "r_ks_distance": float("nan"),
            "z_emd_distance": float("nan"),
            "r_emd_distance": float("nan"),
        }
    zr_rel_error_mean = float(np.mean(np.abs(case_zr - ref_zr) / (ref_zr + eps)))
    zr_l1_distance = float(np.sum(np.abs(case_zr - ref_zr)))
    case_z = _safe_normalize(np.sum(case_zr, axis=1))
    ref_z = _safe_normalize(np.sum(ref_zr, axis=1))
    case_r = _safe_normalize(np.sum(case_zr, axis=0))
    ref_r = _safe_normalize(np.sum(ref_zr, axis=0))
    z_edges = case_dist.z_edges[: n_z + 1]
    r_edges = case_dist.r_edges[: n_r + 1]
    return {
        "zr_rel_error_mean": zr_rel_error_mean,
        "zr_l1_distance": zr_l1_distance,
        "z_ks_distance": _ks_distance(case_z, ref_z),
        "r_ks_distance": _ks_distance(case_r, ref_r),
        "z_emd_distance": _emd_1d(case_z, ref_z, z_edges),
        "r_emd_distance": _emd_1d(case_r, ref_r, r_edges),
    }


def main() -> None:
    args = parse_args()
    case_dirs = collect_case_dirs(args.root)
    all_results: list[dict[str, Any]] = []
    dist_by_case: dict[str, DistributionPack] = {}
    meta_by_case: dict[str, CaseMeta] = {}
    for case_dir in case_dirs:
        meta, metrics, dist = compute_case_metrics(case_dir, args.max_chain_samples)
        meta_by_case[meta.case_name] = meta
        dist_by_case[meta.case_name] = dist
        row: dict[str, Any] = {
            "case_name": meta.case_name,
            "mode": meta.mode,
            "tracking_fraction": meta.fraction,
            "random_seed_scale": meta.seed_scale,
        }
        row.update(metrics)
        all_results.append(row)
    reference_by_mode: dict[str, str] = {}
    for case_name, meta in meta_by_case.items():
        if abs(meta.fraction - 1.0) < 1.0e-12 and meta.mode not in reference_by_mode:
            reference_by_mode[meta.mode] = case_name
    for row in all_results:
        ref_name = reference_by_mode.get(str(row["mode"]), "")
        row["reference_case"] = ref_name
        if not ref_name or row["case_name"] == ref_name:
            row["zr_rel_error_mean"] = float("nan")
            row["zr_l1_distance"] = float("nan")
            row["z_ks_distance"] = float("nan")
            row["r_ks_distance"] = float("nan")
            row["z_emd_distance"] = float("nan")
            row["r_emd_distance"] = float("nan")
            continue
        comp = compare_to_reference(dist_by_case[row["case_name"]], dist_by_case[ref_name])
        row.update(comp)
    fieldnames = [
        "case_name",
        "mode",
        "tracking_fraction",
        "random_seed_scale",
        "reference_case",
        "files_count",
        "particle_count",
        "weighted_nd_sum",
        "qc_proxy_mass",
        "radius_m0",
        "radius_m1",
        "radius_m2",
        "radius_m3",
        "radius_p90",
        "radius_p99",
        "aerosol_weighted_sum",
        "aerosol_positive_weighted_frac",
        "chain_length_mean",
        "chain_length_p90",
        "chain_length_p99",
        "coal_chain_fraction",
        "key_height_coverage_fraction",
        "zr_rel_error_mean",
        "zr_l1_distance",
        "z_ks_distance",
        "r_ks_distance",
        "z_emd_distance",
        "r_emd_distance",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(all_results, key=lambda x: str(x["case_name"])):
            writer.writerow(row)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    _plot_summary_figures(all_results, args.figure_dir)
    print(str(args.output))


def _apply_journal_style() -> None:
    """Apply clean journal-like plotting style."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "savefig.dpi": 300,
        }
    )


def _group_rows(rows: list[dict[str, Any]], mode: str, fraction: float) -> list[dict[str, Any]]:
    grouped: list[dict[str, Any]] = []
    for row in rows:
        if str(row["mode"]) != mode:
            continue
        if abs(float(row["tracking_fraction"]) - fraction) > 1.0e-12:
            continue
        grouped.append(row)
    return grouped


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _extract_series(rows: list[dict[str, Any]], key: str) -> list[float]:
    out: list[float] = []
    for row in rows:
        value = _safe_float(row.get(key))
        if np.isfinite(value):
            out.append(value)
    return out


def _boxplot_metric(rows: list[dict[str, Any]], metric_key: str, y_label: str, out_path: Path) -> None:
    fractions = [0.05, 0.10, 0.20]
    modes = ["stratified", "random"]
    colors = {"stratified": "#1f78b4", "random": "#d95f02"}
    labels: list[str] = []
    data: list[list[float]] = []
    box_colors: list[str] = []
    for fraction in fractions:
        for mode in modes:
            subset = _group_rows(rows, mode, fraction)
            series = _extract_series(subset, metric_key)
            if not series:
                continue
            labels.append(f"{mode[0].upper()}-{int(fraction * 100)}%")
            data.append(series)
            box_colors.append(colors[mode])
    if not data:
        return
    fig, ax = plt.subplots(figsize=(6.6, 2.6))
    bp = ax.boxplot(data, patch_artist=True, widths=0.65, showfliers=False)
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor(color)
        patch.set_linewidth(0.8)
    for item in ["whiskers", "caps", "medians"]:
        for artist in bp[item]:
            artist.set_linewidth(0.8)
            artist.set_color("#333333")
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel(y_label)
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def _plot_reference_scatter(rows: list[dict[str, Any]], x_key: str, y_key: str, out_path: Path, x_label: str, y_label: str) -> None:
    fractions = [0.05, 0.10, 0.20]
    mode_markers = {"stratified": "o", "random": "s"}
    mode_colors = {"stratified": "#1f78b4", "random": "#d95f02"}
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    for mode in ["stratified", "random"]:
        for fraction in fractions:
            subset = _group_rows(rows, mode, fraction)
            xvals = _extract_series(subset, x_key)
            yvals = _extract_series(subset, y_key)
            n = min(len(xvals), len(yvals))
            if n == 0:
                continue
            ax.scatter(
                xvals[:n],
                yvals[:n],
                s=18,
                marker=mode_markers[mode],
                c=mode_colors[mode],
                alpha=0.75,
                edgecolors="none",
                label=f"{mode}-{int(fraction*100)}%",
            )
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(color="#e6e6e6", linewidth=0.6)
    handles, labels = ax.get_legend_handles_labels()
    uniq: dict[str, Any] = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    ax.legend(uniq.values(), uniq.keys(), frameon=False, ncol=2, handletextpad=0.3, columnspacing=0.8)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def _plot_summary_figures(rows: list[dict[str, Any]], figure_dir: Path) -> None:
    """Create publication-style PDF figures from metrics."""
    _apply_journal_style()
    _boxplot_metric(rows, "radius_p99", "Weighted radius P99 (m)", figure_dir / "radius_p99_by_group.pdf")
    _boxplot_metric(rows, "qc_proxy_mass", "qc proxy mass (kg)", figure_dir / "qc_proxy_mass_by_group.pdf")
    _boxplot_metric(rows, "coal_chain_fraction", "Coal-chain fraction", figure_dir / "coal_chain_fraction_by_group.pdf")
    _plot_reference_scatter(
        rows,
        "zr_l1_distance",
        "r_emd_distance",
        figure_dir / "distribution_distance_scatter.pdf",
        "2D z-r L1 distance to reference",
        "Radius EMD to reference (m)",
    )


if __name__ == "__main__":
    main()
