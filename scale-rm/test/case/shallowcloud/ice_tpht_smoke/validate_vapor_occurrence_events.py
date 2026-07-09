#!/usr/bin/env python3
"""Validate v1.3 vapor-growth occurrence process-event records."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")

PROCESS_BITS = {
    7: 64,   # sublimation
    8: 128,  # condensation
    9: 256,  # evaporation
}


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


def read_mask_values(text: str, variable: str) -> list[int]:
    block_re = re.compile(rf"(?:^|\n)\s*({re.escape(variable)}(?:_\d+)?)\s*=\s*(.*?);", re.S)
    values: list[int] = []
    for _name, data in block_re.findall(text):
        values.extend(int(value) for value in INT_RE.findall(data))
    return values


def parse_required(args: argparse.Namespace) -> set[int]:
    required: set[int] = set()
    if args.require_sublimation_occurrence:
        required.add(7)
    if args.require_condensation_occurrence:
        required.add(8)
    if args.require_evaporation_occurrence:
        required.add(9)
    return required


def validate_events(event_glob: str, required: set[int]) -> tuple[Counter[tuple[int, int]], list[str]]:
    paths = sorted(Path(path) for path in glob.glob(event_glob))
    errors: list[str] = []
    counts: Counter[tuple[int, int]] = Counter()

    if not paths:
        return counts, [f"no event files matched {event_glob}"]

    for path in paths:
        try:
            text = ncdump(path)
            trigger_code = read_int_variable(text, "trigger_code")
            trigger_level = read_int_variable(text, "trigger_level")
        except (RuntimeError, ValueError) as exc:
            errors.append(str(exc))
            continue

        if len(trigger_code) != len(trigger_level):
            errors.append(f"{path}: trigger_code/trigger_level length mismatch")
            continue

        counts.update(zip(trigger_code, trigger_level))
        for code in required:
            if counts[(code, 2)]:
                errors.append(f"{path}: unexpected significant trigger_level=2 for trigger_code={code}")

    for code in sorted(required):
        if counts[(code, 1)] <= 0:
            errors.append(f"missing trigger_code={code}, trigger_level=1")

    return counts, errors


def validate_masks(ordinary_glob: str, required: set[int]) -> list[str]:
    paths = sorted(Path(path) for path in glob.glob(ordinary_glob) if not path.endswith(".ids"))
    errors: list[str] = []
    if not paths:
        return [f"no ordinary/selected files matched {ordinary_glob}"]

    event_masks: list[int] = []
    sig_masks: list[int] = []
    for path in paths:
        try:
            text = ncdump(path)
            event_masks.extend(read_mask_values(text, "sd_event_mask"))
            sig_masks.extend(read_mask_values(text, "sd_event_sig_mask"))
        except (RuntimeError, ValueError) as exc:
            errors.append(str(exc))

    for code in sorted(required):
        bit = PROCESS_BITS[code]
        event_hits = sum(1 for value in event_masks if value & bit)
        sig_hits = sum(1 for value in sig_masks if value & bit)
        if event_hits <= 0:
            errors.append(f"sd_event_mask does not contain bit {bit} for trigger_code={code}")
        if sig_hits != 0:
            errors.append(f"sd_event_sig_mask unexpectedly contains bit {bit} for trigger_code={code}: {sig_hits}")

    return errors


def validate_no_lifecycle(lifecycle_glob: str) -> list[str]:
    paths = sorted(Path(path) for path in glob.glob(lifecycle_glob))
    if not paths:
        return []

    errors: list[str] = []
    for path in paths:
        try:
            text = ncdump(path)
            codes = read_int_variable(text, "lifecycle_code")
        except (RuntimeError, ValueError) as exc:
            errors.append(str(exc))
            continue
        if codes:
            errors.append(f"{path}: lifecycle records are not expected for vapor occurrence seed")
    return errors


def main() -> int:
    if len(sys.argv) == 1 or sys.argv[1:] == ["--self-check"]:
        print("validate_vapor_occurrence_events.py: static import/self-check passed")
        return 0

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", required=True)
    parser.add_argument("--ordinary-glob", default="SD_selected_NetCDF_*.pe*")
    parser.add_argument("--lifecycle-glob", default="SD_lifecycle_NetCDF_*")
    parser.add_argument("--require-sublimation-occurrence", action="store_true")
    parser.add_argument("--require-condensation-occurrence", action="store_true")
    parser.add_argument("--require-evaporation-occurrence", action="store_true")
    parser.add_argument("--check-event-mask", action="store_true")
    parser.add_argument("--forbid-lifecycle-records", action="store_true")
    args = parser.parse_args()

    required = parse_required(args)
    if not required:
        print("ERROR: at least one --require-* occurrence flag is required", file=sys.stderr)
        return 1

    counts, errors = validate_events(args.event_glob, required)
    if args.check_event_mask:
        errors.extend(validate_masks(args.ordinary_glob, required))
    if args.forbid_lifecycle_records:
        errors.extend(validate_no_lifecycle(args.lifecycle_glob))

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        observed = ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(counts.items()))
        print(f"observed_pair_counts={observed}", file=sys.stderr)
        return 1

    observed = ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(counts.items()))
    print(f"event_files={len(glob.glob(args.event_glob))}")
    print(f"trigger_pair_counts={observed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
