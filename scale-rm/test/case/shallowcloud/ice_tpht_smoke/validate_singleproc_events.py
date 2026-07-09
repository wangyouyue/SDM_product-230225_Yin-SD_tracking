#!/usr/bin/env python3
"""Validate cold SDM single-process event NetCDF files."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")

PHASE_DRY_AEROSOL = 0
PHASE_LIQUID = 1
PHASE_ICE = 10
PHASE_MIXED = 11
PHASE_NONE = 99
ALLOWED_PHASE_STATES = {
    PHASE_DRY_AEROSOL,
    PHASE_LIQUID,
    PHASE_ICE,
    PHASE_MIXED,
    PHASE_NONE,
}

REQUIRED_FIELDS = (
    "time",
    "trigger_code",
    "trigger_level",
    "target_reason_mask",
    "phase_state_pre",
    "phase_state_post",
    "x",
    "y",
    "z",
    "hydro_radius_pre",
    "hydro_radius_post",
    "liq_mass_pre",
    "liq_mass_post",
    "ice_mass_pre",
    "ice_mass_post",
)

LEGACY_OR_COLLISION_ONLY_FIELDS = (
    "event_type",
    "event_multiplicity",
    "num_col",
    "transition_code",
    "raw_model_phase_state_pre",
    "raw_model_phase_state_post",
    "canonical_phase_state_pre",
    "canonical_phase_state_post",
)

EXTENDED_GEOMETRY_FIELDS = (
    "ice_re_pre",
    "ice_re_post",
    "ice_rp_pre",
    "ice_rp_post",
    "ice_rho_pre",
    "ice_rho_post",
)

KOHLER_CONTEXT_FIELDS = (
    "kohler_rcrit_pre",
    "kohler_rcrit_post",
    "kohler_margin_pre",
    "kohler_margin_post",
    "activated_state_pre",
    "activated_state_post",
)


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


def read_float_variable(text: str, name: str) -> list[float]:
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

    return [float(value.replace("D", "E").replace("d", "e")) for value in FLOAT_RE.findall(text[equals + 1 : end])]


def require_schema(
    text: str,
    path: Path,
    expect_extended_geometry: bool | None,
    expect_kohler_context: bool | None,
) -> None:
    missing = [name for name in REQUIRED_FIELDS if f" {name}(" not in text]
    if expect_extended_geometry is True:
        missing.extend(name for name in EXTENDED_GEOMETRY_FIELDS if f" {name}(" not in text)
    if expect_kohler_context is True:
        missing.extend(name for name in KOHLER_CONTEXT_FIELDS if f" {name}(" not in text)
    if missing:
        raise ValueError(f"{path}: missing required single-process fields: {missing}")

    forbidden = [name for name in LEGACY_OR_COLLISION_ONLY_FIELDS if f" {name}(" in text]
    if expect_extended_geometry is False:
        forbidden.extend(name for name in EXTENDED_GEOMETRY_FIELDS if f" {name}(" in text)
    if expect_kohler_context is False:
        forbidden.extend(name for name in KOHLER_CONTEXT_FIELDS if f" {name}(" in text)
    if forbidden:
        raise ValueError(f"{path}: single-process file has forbidden fields: {forbidden}")


def validate_phase_values(
    path: Path,
    phase_pre: list[int],
    phase_post: list[int],
    hydro_radius_pre: list[float],
    hydro_radius_post: list[float],
) -> None:
    if not (len(phase_pre) == len(phase_post) == len(hydro_radius_pre) == len(hydro_radius_post)):
        raise ValueError(
            f"{path}: inconsistent phase/radius lengths "
            f"pre={len(phase_pre)} post={len(phase_post)} "
            f"rpre={len(hydro_radius_pre)} rpost={len(hydro_radius_post)}"
        )
    for idx, (pre, post, rpre, rpost) in enumerate(zip(phase_pre, phase_post, hydro_radius_pre, hydro_radius_post)):
        for label, value in (("phase_state_pre", pre), ("phase_state_post", post)):
            if value not in ALLOWED_PHASE_STATES:
                raise ValueError(f"{path}: event {idx}: {label}={value} is outside 0/1/10/11/99")
            if value in (2, 3):
                raise ValueError(f"{path}: event {idx}: {label}={value} uses obsolete canonical ice/mixed code")
        if pre in (PHASE_DRY_AEROSOL, PHASE_NONE) and abs(rpre) > 1.0e-30:
            raise ValueError(f"{path}: event {idx}: dry/missing pre phase has nonzero hydro_radius_pre={rpre}")
        if post in (PHASE_DRY_AEROSOL, PHASE_NONE) and abs(rpost) > 1.0e-30:
            raise ValueError(f"{path}: event {idx}: dry/missing post phase has nonzero hydro_radius_post={rpost}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", required=True, help="Glob for SD_event_singleproc_NetCDF files")
    parser.add_argument(
        "--require-trigger-code",
        type=int,
        action="append",
        default=[],
        help="Trigger code that must appear at least --min-count-per-type times",
    )
    parser.add_argument(
        "--require-trigger-pair",
        action="append",
        default=[],
        metavar="CODE:LEVEL",
        help="Process trigger pair required at least --min-count-per-type times, e.g. 6:2 for significant deposition",
    )
    extended_group = parser.add_mutually_exclusive_group()
    extended_group.add_argument(
        "--require-extended-geometry",
        action="store_true",
        help="Require TRACK_COLD_EVENT_EXTENDED_GEOMETRY single-process ice geometry fields",
    )
    extended_group.add_argument(
        "--forbid-extended-geometry",
        action="store_true",
        help="Require default single-process schema without optional ice geometry fields",
    )
    kohler_group = parser.add_mutually_exclusive_group()
    kohler_group.add_argument(
        "--require-kohler-context",
        action="store_true",
        help="Require TRACK_COLD_OUTPUT_KOHLER_CONTEXT fields",
    )
    kohler_group.add_argument(
        "--forbid-kohler-context",
        action="store_true",
        help="Require schema without optional Kohler context fields",
    )
    parser.add_argument("--min-count-per-type", type=int, default=1)
    args = parser.parse_args()

    expect_extended_geometry: bool | None
    if args.require_extended_geometry:
        expect_extended_geometry = True
    elif args.forbid_extended_geometry:
        expect_extended_geometry = False
    else:
        expect_extended_geometry = None
    if args.require_kohler_context:
        expect_kohler_context: bool | None = True
    elif args.forbid_kohler_context:
        expect_kohler_context = False
    else:
        expect_kohler_context = None

    paths = sorted(Path(p) for p in glob.glob(args.event_glob))
    if not paths:
        print(f"ERROR: no single-process event files matched {args.event_glob}", file=sys.stderr)
        return 1

    counts: Counter[int] = Counter()
    pair_counts: Counter[tuple[int, int]] = Counter()
    total_events = 0
    for path in paths:
        try:
            text = ncdump(path)
            require_schema(text, path, expect_extended_geometry, expect_kohler_context)
            trigger_code = read_int_variable(text, "trigger_code")
            trigger_level = read_int_variable(text, "trigger_level")
            phase_pre = read_int_variable(text, "phase_state_pre")
            phase_post = read_int_variable(text, "phase_state_post")
            hydro_radius_pre = read_float_variable(text, "hydro_radius_pre")
            hydro_radius_post = read_float_variable(text, "hydro_radius_post")
            validate_phase_values(path, phase_pre, phase_post, hydro_radius_pre, hydro_radius_post)
            if len(trigger_code) != len(trigger_level):
                raise ValueError(
                    f"{path}: inconsistent trigger_code/trigger_level lengths "
                    f"{len(trigger_code)} != {len(trigger_level)}"
                )
            for idx, (code, level) in enumerate(zip(trigger_code, trigger_level)):
                if code < 4 or code > 11:
                    raise ValueError(f"{path}: event {idx}: single-process trigger_code={code} outside 4..11")
                if level not in (1, 2):
                    raise ValueError(f"{path}: event {idx}: trigger_level={level} outside 1/2")
                if code in (10, 11) and level == 2:
                    raise ValueError(
                        f"{path}: event {idx}: activation/deactivation significant level is inactive in this patch"
                    )
        except (RuntimeError, ValueError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        counts.update(trigger_code)
        pair_counts.update(zip(trigger_code, trigger_level))
        total_events += len(trigger_code)

    if total_events == 0:
        print(f"ERROR: matched {len(paths)} files but found zero single-process events", file=sys.stderr)
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

    required_pairs: list[tuple[int, int]] = []
    for item in args.require_trigger_pair:
        try:
            code_text, level_text = item.split(":", 1)
            required_pairs.append((int(code_text), int(level_text)))
        except ValueError:
            print(f"ERROR: invalid --require-trigger-pair value {item!r}; expected CODE:LEVEL", file=sys.stderr)
            return 1
    missing_pairs = [pair for pair in required_pairs if pair_counts[pair] < args.min_count_per_type]
    if missing_pairs:
        print(
            "ERROR: missing required trigger_code:trigger_level counts: "
            + ", ".join(f"{code}:{level}={pair_counts[(code, level)]}" for code, level in missing_pairs),
            file=sys.stderr,
        )
        print(
            "observed_pair_counts="
            + ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(pair_counts.items())),
            file=sys.stderr,
        )
        return 1

    print(f"event_files={len(paths)}")
    print(f"total_events={total_events}")
    print("trigger_code_counts=" + ",".join(f"{code}:{count}" for code, count in sorted(counts.items())))
    print("trigger_pair_counts=" + ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(pair_counts.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
