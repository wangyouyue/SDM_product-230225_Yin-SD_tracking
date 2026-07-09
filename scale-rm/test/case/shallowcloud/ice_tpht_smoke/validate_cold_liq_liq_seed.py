#!/usr/bin/env python3
"""Validate the dedicated cold liquid-liquid collision seed."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"-?\d+")
PHASE_LIQUID = 1
EVENT_COALESCENCE = 1
TRIG_PROC_COALESCENCE = 1
TRIG_LEVEL_OCCURRENCE = 1
TRIG_LEVEL_SIGNIFICANT = 2


def ncdump_var(path: Path, variable: str) -> list[int]:
    result = subprocess.run(
        ["ncdump", "-v", variable, str(path)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise RuntimeError(f"ncdump failed for {path}:{variable}: {result.stderr.strip()}")
    data = result.stdout.split("data:", 1)[-1]
    data = data.split("=", 1)[-1]
    return [int(value) for value in INT_RE.findall(data)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", required=True, help="Glob for SD_event_collision_NetCDF_* files")
    parser.add_argument("--ordinary-glob", help="Optional glob for SD_selected/SD_all output with sd_event_mask")
    parser.add_argument(
        "--require-level",
        type=int,
        action="append",
        choices=(TRIG_LEVEL_OCCURRENCE, TRIG_LEVEL_SIGNIFICANT),
        help="Required coalescence trigger_level. Defaults to both levels.",
    )
    args = parser.parse_args()
    required_levels = args.require_level or [TRIG_LEVEL_OCCURRENCE, TRIG_LEVEL_SIGNIFICANT]

    event_paths = sorted(Path(path) for path in glob.glob(args.event_glob) if not path.endswith(".ids"))
    if not event_paths:
        print(f"ERROR: no event files matched {args.event_glob}", file=sys.stderr)
        return 1

    trigger_counts: Counter[int] = Counter()
    trigger_pair_counts: Counter[tuple[int, int]] = Counter()
    multiplicity_by_pair: Counter[tuple[int, int]] = Counter()
    phase_errors: list[str] = []
    try:
        for path in event_paths:
            triggers = ncdump_var(path, "trigger_code")
            levels = ncdump_var(path, "trigger_level")
            multiplicity = ncdump_var(path, "event_multiplicity")
            p1 = ncdump_var(path, "phase_state1_pre")
            p2 = ncdump_var(path, "phase_state2_pre")
            trigger_counts.update(triggers)
            trigger_pair_counts.update(zip(triggers, levels))
            for code, level, mult, phase1, phase2 in zip(triggers, levels, multiplicity, p1, p2):
                if code == TRIG_PROC_COALESCENCE:
                    multiplicity_by_pair[(code, level)] += mult
                    if phase1 != PHASE_LIQUID or phase2 != PHASE_LIQUID:
                        phase_errors.append(f"{path.name}: trigger_code={code} trigger_level={level} phases={phase1},{phase2}")
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    required_pairs = [(TRIG_PROC_COALESCENCE, level) for level in required_levels]
    missing = [pair for pair in required_pairs if trigger_pair_counts[pair] <= 0]
    if missing:
        print(
            f"ERROR: missing coalescence trigger_code:trigger_level pairs {missing}; "
            f"observed={dict(trigger_pair_counts)}",
            file=sys.stderr,
        )
        return 1
    if phase_errors:
        for error in phase_errors[:20]:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    for pair in required_pairs:
        if multiplicity_by_pair[pair] <= 0:
            print(f"ERROR: trigger pair {pair} has non-positive summed event_multiplicity", file=sys.stderr)
            return 1

    mask_counts: Counter[int] = Counter()
    sig_mask_counts: Counter[int] = Counter()
    if args.ordinary_glob:
        ordinary_paths = sorted(Path(path) for path in glob.glob(args.ordinary_glob) if not path.endswith(".ids"))
        if not ordinary_paths:
            print(f"ERROR: no ordinary files matched {args.ordinary_glob}", file=sys.stderr)
            return 1
        for path in ordinary_paths:
            for mask in ncdump_var(path, "sd_event_mask"):
                if mask & EVENT_COALESCENCE:
                    mask_counts[EVENT_COALESCENCE] += 1
            for mask in ncdump_var(path, "sd_event_sig_mask"):
                if mask & EVENT_COALESCENCE:
                    sig_mask_counts[EVENT_COALESCENCE] += 1
        if mask_counts[EVENT_COALESCENCE] <= 0:
            print(f"ERROR: missing sd_event_mask coalescence bit; counts={dict(mask_counts)}", file=sys.stderr)
            return 1
        if TRIG_LEVEL_SIGNIFICANT in required_levels and sig_mask_counts[EVENT_COALESCENCE] <= 0:
            print(f"ERROR: missing sd_event_sig_mask coalescence bit; counts={dict(sig_mask_counts)}", file=sys.stderr)
            return 1

    print(f"event_files={len(event_paths)}")
    print("trigger_code_counts=" + ",".join(f"{k}:{v}" for k, v in sorted(trigger_counts.items())))
    print("trigger_pair_counts=" + ",".join(f"{k[0]}:{k[1]}:{v}" for k, v in sorted(trigger_pair_counts.items())))
    print("event_multiplicity_sums=" + ",".join(f"{k[0]}:{k[1]}:{v}" for k, v in sorted(multiplicity_by_pair.items())))
    if args.ordinary_glob:
        print("mask_bit_counts=" + ",".join(f"{k}:{v}" for k, v in sorted(mask_counts.items())))
        print("sig_mask_bit_counts=" + ",".join(f"{k}:{v}" for k, v in sorted(sig_mask_counts.items())))
    print("cold_liq_liq_seed=ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
