#!/usr/bin/env python3
"""Validate cold SDM diagnostic event NetCDF files."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")

REQUIRED_FIELDS = (
    "time",
    "trigger_code",
    "target_reason_mask",
    "x",
    "y",
    "z",
    "sd_event_mask",
    "sd_event_sig_mask",
    "sd_diag_mask",
    "sd_phase_change_flag",
    "sd_liq_radius_max_interval",
    "sd_ice_rvol_max_interval",
    "sd_mixed_rvol_max_interval",
    "sd_rime_mass_max_interval",
    "sd_rime_frac_max_interval",
    "sd_nmono_max_interval",
    "sd_aspect_ratio_max_interval",
)

FORBIDDEN_FIELDS = ("event_type", "event_multiplicity", "num_col")


def ncdump(path: Path) -> str:
    try:
        result = subprocess.run(
            ["ncdump", str(path)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise RuntimeError("ncdump is required for this validator") from None
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"ncdump failed for {path}: {exc.stderr.strip()}") from exc
    return result.stdout


def read_int_variable(text: str, name: str) -> list[int]:
    marker = f" {name} ="
    start = text.find(marker)
    if start < 0:
        marker = f"\n{name} ="
        start = text.find(marker)
    if start < 0:
        raise ValueError(f"variable {name!r} is missing")

    equals = text.find("=", start)
    end = text.find(";", equals)
    if equals < 0 or end < 0:
        raise ValueError(f"variable {name!r} has no complete data block")

    return [int(value) for value in INT_RE.findall(text[equals + 1 : end])]


def require_schema(text: str, path: Path) -> None:
    missing = [name for name in REQUIRED_FIELDS if f" {name}(" not in text]
    if missing:
        raise ValueError(f"{path}: missing required diagnostic fields: {missing}")

    forbidden = [name for name in FORBIDDEN_FIELDS if f" {name}(" in text]
    if forbidden:
        raise ValueError(f"{path}: diagnostic file has forbidden fields: {forbidden}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", required=True, help="Glob for SD_event_diag_NetCDF files")
    parser.add_argument("--expect-none", action="store_true", help="Pass only if no diagnostic event files exist")
    parser.add_argument(
        "--require-trigger-code",
        type=int,
        action="append",
        default=[],
        help="Trigger code that must appear at least --min-count-per-type times",
    )
    parser.add_argument("--min-count-per-type", type=int, default=1)
    args = parser.parse_args()

    paths = sorted(Path(p) for p in glob.glob(args.event_glob))
    if args.expect_none:
        if paths:
            print(f"ERROR: expected no diagnostic event files, found {len(paths)}", file=sys.stderr)
            return 1
        print("diag_event_files=0")
        return 0

    if not paths:
        print(f"ERROR: no diagnostic event files matched {args.event_glob}", file=sys.stderr)
        return 1

    counts: Counter[int] = Counter()
    total_events = 0
    for path in paths:
        try:
            text = ncdump(path)
            require_schema(text, path)
            trigger_code = read_int_variable(text, "trigger_code")
        except (RuntimeError, ValueError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        counts.update(trigger_code)
        total_events += len(trigger_code)

    if total_events == 0:
        print(f"ERROR: matched {len(paths)} files but found zero diagnostic events", file=sys.stderr)
        return 1

    missing = [code for code in args.require_trigger_code if counts[code] < args.min_count_per_type]
    if missing:
        print(
            "ERROR: missing required trigger_code counts: "
            + ", ".join(f"{code}={counts[code]}" for code in missing),
            file=sys.stderr,
        )
        print(f"observed_counts={dict(sorted(counts.items()))}", file=sys.stderr)
        return 1

    print(f"event_files={len(paths)}")
    print(f"total_events={total_events}")
    print("trigger_code_counts=" + ",".join(f"{code}:{count}" for code, count in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
