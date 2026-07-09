#!/usr/bin/env python3
"""Sweep Cold TPHT event files for schema and trigger-code invariants."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")
VAR_RE = re.compile(
    r"^\s*(?:byte|char|short|int|int64|float|double)\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(",
    re.MULTILINE,
)

FORBIDDEN_VARIABLES = {
    "event_type",
    "transition_code",
    "transition_code1",
    "transition_code2",
    "raw_model_phase_state_pre",
    "raw_model_phase_state_post",
    "raw_model_phase_state1_pre",
    "raw_model_phase_state1_post",
    "raw_model_phase_state2_pre",
    "raw_model_phase_state2_post",
    "canonical_phase_state_pre",
    "canonical_phase_state_post",
    "canonical_phase_state1_pre",
    "canonical_phase_state1_post",
    "canonical_phase_state2_pre",
    "canonical_phase_state2_post",
}

ALLOWED_PHASE_STATES = {0, 1, 10, 11, 99}
OBSOLETE_PHASE_STATES = {2, 3}
ACTIVE_PROCESS_TRIGGER_CODES = set(range(1, 12))
ACTIVE_DIAG_TRIGGER_CODES = set(range(101, 108))
RESERVED_INACTIVE_TRIGGER_CODES = {20, 21, 22, 23, 30, 31, 32, 40, 41}


def ncdump(path: Path, *args: str) -> str:
    try:
        result = subprocess.run(
            ["ncdump", *args, str(path)],
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


def read_int_variable(path: Path, name: str) -> list[int]:
    text = ncdump(path, "-v", name)
    match = re.search(r"\n\s*" + re.escape(name) + r"\s*=\s*(.*?);", text, re.S)
    if match is None:
        raise RuntimeError(f"{path}: variable {name!r} has no data block")
    return [int(value) for value in INT_RE.findall(match.group(1))]


def paths_from_globs(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(Path(path) for path in glob.glob(pattern) if not path.endswith(".ids"))
    return sorted(set(path for path in paths if path.is_file()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", action="append", required=True, help="Event NetCDF glob; can be repeated")
    parser.add_argument(
        "--require-trigger-code",
        type=int,
        action="append",
        default=[],
        help="Trigger code that must appear at least once across all matched files",
    )
    parser.add_argument(
        "--require-trigger-pair",
        action="append",
        default=[],
        metavar="CODE:LEVEL",
        help="Process trigger_code:trigger_level pair that must appear at least once",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Permit zero matched files. Intended only for optional event families.",
    )
    args = parser.parse_args()

    paths = paths_from_globs(args.glob)
    if not paths:
        if args.allow_empty:
            print("event_files=0")
            print("schema_sweep=ok")
            return 0
        print(f"ERROR: no files matched: {args.glob}", file=sys.stderr)
        return 1

    trigger_counts: Counter[int] = Counter()
    trigger_pair_counts: Counter[tuple[int, int]] = Counter()
    phase_counts: Counter[int] = Counter()
    checked_variables = 0

    try:
        for path in paths:
            header = ncdump(path, "-h")
            variables = set(VAR_RE.findall(header))
            forbidden = sorted(FORBIDDEN_VARIABLES & variables)
            if forbidden:
                print(f"ERROR: {path}: forbidden v1.2 variables present: {forbidden}", file=sys.stderr)
                return 1

            if "trigger_code" in variables:
                trigger_values = read_int_variable(path, "trigger_code")
                trigger_counts.update(trigger_values)
                is_diag = "sd_diag_mask" in variables
                is_process = "trigger_level" in variables
                if is_process:
                    trigger_level_values = read_int_variable(path, "trigger_level")
                    if len(trigger_values) != len(trigger_level_values):
                        print(
                            f"ERROR: {path}: trigger_code/trigger_level length mismatch "
                            f"{len(trigger_values)} != {len(trigger_level_values)}",
                            file=sys.stderr,
                        )
                        return 1
                    trigger_pair_counts.update(zip(trigger_values, trigger_level_values))
                    bad_level = sorted(set(trigger_level_values) - {1, 2})
                    if bad_level:
                        print(f"ERROR: {path}: invalid process trigger_level values: {bad_level}", file=sys.stderr)
                        return 1
                    active_codes = ACTIVE_PROCESS_TRIGGER_CODES
                elif is_diag:
                    active_codes = ACTIVE_DIAG_TRIGGER_CODES
                else:
                    active_codes = ACTIVE_PROCESS_TRIGGER_CODES | ACTIVE_DIAG_TRIGGER_CODES
                bad_reserved = sorted(set(trigger_values) & RESERVED_INACTIVE_TRIGGER_CODES)
                if bad_reserved:
                    print(
                        f"ERROR: {path}: inactive reserved trigger codes appeared at runtime: {bad_reserved}",
                        file=sys.stderr,
                    )
                    return 1
                bad_unknown = sorted(set(trigger_values) - active_codes)
                if bad_unknown:
                    print(f"ERROR: {path}: unknown trigger codes appeared: {bad_unknown}", file=sys.stderr)
                    return 1

            for variable in sorted(name for name in variables if name.startswith("phase_state")):
                phase_values = read_int_variable(path, variable)
                checked_variables += 1
                phase_counts.update(phase_values)
                bad_phase = sorted(set(phase_values) - ALLOWED_PHASE_STATES)
                obsolete = sorted(set(phase_values) & OBSOLETE_PHASE_STATES)
                if bad_phase or obsolete:
                    print(
                        f"ERROR: {path}:{variable}: bad phase values={bad_phase}, obsolete={obsolete}",
                        file=sys.stderr,
                    )
                    return 1
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    missing = [code for code in args.require_trigger_code if trigger_counts[code] <= 0]
    if missing:
        print(
            "ERROR: missing required trigger codes: "
            + ", ".join(f"{code}=0" for code in missing),
            file=sys.stderr,
        )
        print(f"observed_trigger_counts={dict(sorted(trigger_counts.items()))}", file=sys.stderr)
        return 1

    required_pairs: list[tuple[int, int]] = []
    for item in args.require_trigger_pair:
        try:
            code_text, level_text = item.split(":", 1)
            required_pairs.append((int(code_text), int(level_text)))
        except ValueError:
            print(f"ERROR: invalid --require-trigger-pair value {item!r}; expected CODE:LEVEL", file=sys.stderr)
            return 1
    missing_pairs = [pair for pair in required_pairs if trigger_pair_counts[pair] <= 0]
    if missing_pairs:
        print(
            "ERROR: missing required trigger_code:trigger_level pairs: "
            + ", ".join(f"{code}:{level}=0" for code, level in missing_pairs),
            file=sys.stderr,
        )
        print(
            "observed_trigger_pair_counts="
            + ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(trigger_pair_counts.items())),
            file=sys.stderr,
        )
        return 1

    print(f"event_files={len(paths)}")
    print("trigger_code_counts=" + ",".join(f"{code}:{count}" for code, count in sorted(trigger_counts.items())))
    if trigger_pair_counts:
        print("trigger_pair_counts=" + ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(trigger_pair_counts.items())))
    print("phase_state_counts=" + ",".join(f"{code}:{count}" for code, count in sorted(phase_counts.items())))
    print(f"phase_variables_checked={checked_variables}")
    print("schema_sweep=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
