#!/usr/bin/env python3
"""Summarize Cold TPHT event-stream sizes and trigger counts for I/O calibration."""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


INT_RE = re.compile(r"-?\d+")
FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")


STREAM_GLOBS = {
    "collision": "SD_event_collision_NetCDF_*",
    "singleproc": "SD_event_singleproc_NetCDF_*",
    "diag": "SD_event_diag_NetCDF_*",
    "all": "SD_all_NetCDF_*",
    "selected": "SD_selected_NetCDF_*",
    "lifecycle": "SD_lifecycle_NetCDF_*",
}


def ncdump_payload(path: Path, variable: str) -> str:
    result = subprocess.run(
        ["ncdump", "-v", variable, str(path)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        return ""
    data = result.stdout.split("data:", 1)[-1]
    return data.split("=", 1)[-1]


def ncdump_int_var(path: Path, variable: str) -> list[int]:
    data = ncdump_payload(path, variable)
    if not data:
        return []
    return [int(value) for value in INT_RE.findall(data)]


def ncdump_float_var(path: Path, variable: str) -> list[float]:
    data = ncdump_payload(path, variable)
    if not data:
        return []
    return [float(value.replace("D", "E").replace("d", "e")) for value in FLOAT_RE.findall(data)]


def rank_from_path(path: Path) -> str:
    match = re.search(r"\.pe(\d+)", path.name)
    return match.group(1) if match else "unknown"


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lo = int(position)
    hi = min(lo + 1, len(ordered) - 1)
    weight = position - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def summarize_distribution(values: list[float]) -> dict:
    qs = {
        "p50": 0.50,
        "p75": 0.75,
        "p90": 0.90,
        "p95": 0.95,
        "p99": 0.99,
        "p99.9": 0.999,
    }
    result = {"count": len(values)}
    for label, q in qs.items():
        value = percentile(values, q)
        if value is not None:
            result[label] = value
    return result


def add_relative_metric(
    distributions: defaultdict[str, list[float]],
    denominator_floor_stats: dict[str, dict[str, int]],
    key: str,
    numerator: float,
    denominator: float,
    denominator_floor: float,
) -> None:
    abs_numerator = abs(numerator)
    abs_denominator = abs(denominator)
    stats = denominator_floor_stats.setdefault(key, {"count": 0, "denominator_below_floor_count": 0})
    stats["count"] += 1
    if abs_denominator < denominator_floor:
        stats["denominator_below_floor_count"] += 1

    distributions[f"{key}_raw_relative"].append(abs_numerator / max(abs_denominator, sys.float_info.min))
    distributions[f"{key}_absolute"].append(abs_numerator)
    if abs_denominator >= denominator_floor:
        distributions[f"{key}_denominator_filtered_relative"].append(abs_numerator / abs_denominator)


def summarize_case(case_dir: Path, denominator_floor: float) -> dict:
    output: dict = {
        "case_dir": str(case_dir),
        "streams": {},
        "trigger_code_counts": {},
        "trigger_level_counts": {},
        "trigger_pair_counts": {},
        "sd_event_mask_bit_counts": {},
        "sd_event_sig_mask_bit_counts": {},
        "sd_diag_mask_bit_counts": {},
        "denominator_floor": denominator_floor,
        "denominator_floor_stats": {},
        "relative_mass_change_distributions": {},
        "denominator_filtered_relative_mass_change_distributions": {},
        "absolute_mass_change_distributions": {},
    }
    trigger_counts: Counter[int] = Counter()
    trigger_level_counts: Counter[int] = Counter()
    trigger_pair_counts: Counter[tuple[int, int]] = Counter()
    event_bit_counts: Counter[int] = Counter()
    event_sig_bit_counts: Counter[int] = Counter()
    diag_bit_counts: Counter[int] = Counter()
    distributions: defaultdict[str, list[float]] = defaultdict(list)
    denominator_floor_stats: dict[str, dict[str, int]] = {}

    for stream, pattern in STREAM_GLOBS.items():
        paths = sorted(case_dir.glob(pattern))
        record_counts_by_rank: Counter[str] = Counter()
        output["streams"][stream] = {
            "files": len(paths),
            "bytes": sum(path.stat().st_size for path in paths if path.is_file()),
            "max_file_bytes": max((path.stat().st_size for path in paths if path.is_file()), default=0),
        }
        if stream in {"collision", "singleproc", "diag"}:
            for path in paths:
                triggers = ncdump_int_var(path, "trigger_code")
                trigger_counts.update(triggers)
                levels = ncdump_int_var(path, "trigger_level")
                if levels:
                    trigger_level_counts.update(levels)
                    trigger_pair_counts.update(zip(triggers, levels))
                record_counts_by_rank[rank_from_path(path)] += len(triggers)
                if stream == "collision":
                    mass1 = ncdump_float_var(path, "hydro_mass1_pre")
                    mass2 = ncdump_float_var(path, "hydro_mass2_pre")
                    for code, m1, m2 in zip(triggers, mass1, mass2):
                        add_relative_metric(
                            distributions,
                            denominator_floor_stats,
                            f"collision_pair_mass_ratio_trigger_{code}",
                            min(abs(m1), abs(m2)),
                            max(abs(m1), abs(m2)),
                            denominator_floor,
                        )
                elif stream == "singleproc":
                    liq_pre = ncdump_float_var(path, "liq_mass_pre")
                    liq_post = ncdump_float_var(path, "liq_mass_post")
                    ice_pre = ncdump_float_var(path, "ice_mass_pre")
                    ice_post = ncdump_float_var(path, "ice_mass_post")
                    for code, lp, lq, ip, iq in zip(triggers, liq_pre, liq_post, ice_pre, ice_post):
                        add_relative_metric(
                            distributions,
                            denominator_floor_stats,
                            f"singleproc_liq_change_trigger_{code}",
                            lq - lp,
                            lp,
                            denominator_floor,
                        )
                        add_relative_metric(
                            distributions,
                            denominator_floor_stats,
                            f"singleproc_ice_change_trigger_{code}",
                            iq - ip,
                            ip,
                            denominator_floor,
                        )
        if stream in {"all", "selected"}:
            for path in paths:
                masks = ncdump_int_var(path, "sd_event_mask")
                for bit in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
                    event_bit_counts[bit] += sum(1 for value in masks if value & bit)
                sig_masks = ncdump_int_var(path, "sd_event_sig_mask")
                for bit in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
                    event_sig_bit_counts[bit] += sum(1 for value in sig_masks if value & bit)
                diag_masks = ncdump_int_var(path, "sd_diag_mask")
                for bit in (1, 2, 4, 8, 16, 32, 64):
                    diag_bit_counts[bit] += sum(1 for value in diag_masks if value & bit)
        if record_counts_by_rank:
            values = list(record_counts_by_rank.values())
            output["streams"][stream]["records"] = sum(values)
            output["streams"][stream]["records_by_rank"] = dict(sorted(record_counts_by_rank.items()))
            output["streams"][stream]["max_records_per_rank"] = max(values)
            output["streams"][stream]["min_records_per_rank"] = min(values)

    output["trigger_code_counts"] = {str(k): v for k, v in sorted(trigger_counts.items())}
    output["trigger_level_counts"] = {str(k): v for k, v in sorted(trigger_level_counts.items())}
    output["trigger_pair_counts"] = {f"{k[0]}:{k[1]}": v for k, v in sorted(trigger_pair_counts.items())}
    output["sd_event_mask_bit_counts"] = {str(k): v for k, v in sorted(event_bit_counts.items()) if v}
    output["sd_event_sig_mask_bit_counts"] = {str(k): v for k, v in sorted(event_sig_bit_counts.items()) if v}
    output["sd_diag_mask_bit_counts"] = {str(k): v for k, v in sorted(diag_bit_counts.items()) if v}
    output["denominator_floor_stats"] = {
        key: {
            **stats,
            "denominator_below_floor_fraction": (
                stats["denominator_below_floor_count"] / stats["count"] if stats["count"] else 0.0
            ),
        }
        for key, stats in sorted(denominator_floor_stats.items())
    }
    output["relative_mass_change_distributions"] = {
        key.removesuffix("_raw_relative"): summarize_distribution(values)
        for key, values in sorted(distributions.items())
        if key.endswith("_raw_relative")
    }
    output["denominator_filtered_relative_mass_change_distributions"] = {
        key.removesuffix("_denominator_filtered_relative"): summarize_distribution(values)
        for key, values in sorted(distributions.items())
        if key.endswith("_denominator_filtered_relative")
    }
    output["absolute_mass_change_distributions"] = {
        key.removesuffix("_absolute"): summarize_distribution(values)
        for key, values in sorted(distributions.items())
        if key.endswith("_absolute")
    }
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dirs", nargs="+", help="Case output directories to summarize")
    parser.add_argument("--json-out", help="Optional JSON output path")
    parser.add_argument(
        "--denominator-floor",
        type=float,
        default=1.0e-24,
        help=(
            "Denominator floor for filtered relative mass-change percentiles. "
            "Raw relative percentiles still use machine tiny to preserve legacy behavior."
        ),
    )
    args = parser.parse_args()

    summaries = [summarize_case(Path(case_dir), args.denominator_floor) for case_dir in args.case_dirs]
    text = json.dumps(summaries, indent=2, sort_keys=True)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
