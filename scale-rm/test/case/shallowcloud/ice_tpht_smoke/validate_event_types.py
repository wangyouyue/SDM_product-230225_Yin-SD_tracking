#!/usr/bin/env python3
"""Validate cold SDM collision trigger codes from SD_event_collision_NetCDF files."""

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

REQUIRED_SCHEMA_FIELDS = (
    "trigger_code",
    "trigger_level",
    "target_reason_mask",
    "event_multiplicity",
    "phase_state1_pre",
    "phase_state1_post",
    "phase_state2_pre",
    "phase_state2_post",
    "x",
    "y",
    "z",
    "x1",
    "y1",
    "z1",
    "x2",
    "y2",
    "z2",
    "sd_n1_pre",
    "sd_n2_pre",
    "sd_n1_post",
    "sd_n2_post",
    "hydro_radius1_pre",
    "hydro_radius2_pre",
    "hydro_radius1_post",
    "hydro_radius2_post",
    "hydro_mass1_pre",
    "hydro_mass2_pre",
    "hydro_mass1_post",
    "hydro_mass2_post",
)

ICE_GEOMETRY_FIELDS = (
    "ice_re1_pre",
    "ice_rp1_pre",
    "ice_rho1_pre",
    "ice_re1_post",
    "ice_rp1_post",
    "ice_rho1_post",
    "ice_re2_pre",
    "ice_rp2_pre",
    "ice_rho2_pre",
    "ice_re2_post",
    "ice_rp2_post",
    "ice_rho2_post",
)

RIME_MORPHOLOGY_FIELDS = (
    "rime_mass1_pre",
    "rime_mass1_post",
    "rime_mass2_pre",
    "rime_mass2_post",
    "rime_frac1_pre",
    "rime_frac1_post",
    "rime_frac2_pre",
    "rime_frac2_post",
    "aspect_ratio1_pre",
    "aspect_ratio1_post",
    "aspect_ratio2_pre",
    "aspect_ratio2_post",
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


def require_schema_fields(
    text: str,
    path: Path,
    require_ice_geometry: bool,
    require_rime_morphology: bool,
    forbid_optional_groups: bool,
) -> None:
    missing = [name for name in REQUIRED_SCHEMA_FIELDS if f" {name}(" not in text]
    if require_ice_geometry:
        missing.extend(name for name in ICE_GEOMETRY_FIELDS if f" {name}(" not in text)
    if require_rime_morphology:
        missing.extend(name for name in RIME_MORPHOLOGY_FIELDS if f" {name}(" not in text)
    if missing:
        raise ValueError(f"{path}: missing required cold collision schema fields: {missing}")

    legacy_fields = [
        name
        for name in (
            "event_type",
            "sd_liqice1",
            "sd_liqice2",
            "num_col",
            "transition_code1",
            "transition_code2",
            "raw_model_phase_state1_pre",
            "raw_model_phase_state1_post",
            "raw_model_phase_state2_pre",
            "raw_model_phase_state2_post",
            "canonical_phase_state1_pre",
            "canonical_phase_state1_post",
            "canonical_phase_state2_pre",
            "canonical_phase_state2_post",
        )
        if f" {name}(" in text
    ]
    if forbid_optional_groups:
        legacy_fields.extend(
            name for name in (*ICE_GEOMETRY_FIELDS, *RIME_MORPHOLOGY_FIELDS) if f" {name}(" in text
        )
    if legacy_fields:
        raise ValueError(f"{path}: cold collision file still writes legacy fields: {legacy_fields}")


def parse_event_file(
    path: Path, require_ice_geometry: bool, require_rime_morphology: bool, forbid_optional_groups: bool
) -> tuple[
    list[int],
    list[int],
    list[int],
    list[int],
    list[int],
    list[int],
    list[float],
    list[float],
    list[float],
    list[float],
]:
    text = ncdump(path)
    require_schema_fields(text, path, require_ice_geometry, require_rime_morphology, forbid_optional_groups)
    trigger_code = read_int_variable(text, "trigger_code")
    trigger_level = read_int_variable(text, "trigger_level")
    phase1 = read_int_variable(text, "phase_state1_pre")
    phase2 = read_int_variable(text, "phase_state2_pre")
    phase1_post = read_int_variable(text, "phase_state1_post")
    phase2_post = read_int_variable(text, "phase_state2_post")
    radius1 = read_float_variable(text, "hydro_radius1_pre")
    radius2 = read_float_variable(text, "hydro_radius2_pre")
    radius1_post = read_float_variable(text, "hydro_radius1_post")
    radius2_post = read_float_variable(text, "hydro_radius2_post")
    if not (
        len(trigger_code)
        == len(trigger_level)
        == len(phase1)
        == len(phase2)
        == len(phase1_post)
        == len(phase2_post)
        == len(radius1)
        == len(radius2)
        == len(radius1_post)
        == len(radius2_post)
    ):
        raise ValueError(
            f"{path}: inconsistent lengths trigger_code={len(trigger_code)} "
            f"trigger_level={len(trigger_level)} "
            f"phase_state1_pre={len(phase1)} phase_state2_pre={len(phase2)} "
            f"phase_state1_post={len(phase1_post)} phase_state2_post={len(phase2_post)}"
        )
    return (
        trigger_code,
        trigger_level,
        phase1,
        phase2,
        phase1_post,
        phase2_post,
        radius1,
        radius2,
        radius1_post,
        radius2_post,
    )


def validate_phase_consistency(
    trigger_code: list[int],
    phase1: list[int],
    phase2: list[int],
    phase1_post: list[int],
    phase2_post: list[int],
    radius1: list[float],
    radius2: list[float],
    radius1_post: list[float],
    radius2_post: list[float],
) -> list[str]:
    errors: list[str] = []
    for idx, (code, p1, p2, p1_post, p2_post, r1, r2, r1_post, r2_post) in enumerate(
        zip(trigger_code, phase1, phase2, phase1_post, phase2_post, radius1, radius2, radius1_post, radius2_post)
    ):
        for label, value in (
            ("phase_state1_pre", p1),
            ("phase_state2_pre", p2),
            ("phase_state1_post", p1_post),
            ("phase_state2_post", p2_post),
        ):
            if value not in ALLOWED_PHASE_STATES:
                errors.append(f"event {idx}: {label}={value} is outside 0/1/10/11/99")
            if value in (2, 3):
                errors.append(f"event {idx}: {label}={value} uses obsolete canonical ice/mixed code")
        for label, value, radius in (
            ("phase_state1_pre", p1, r1),
            ("phase_state2_pre", p2, r2),
            ("phase_state1_post", p1_post, r1_post),
            ("phase_state2_post", p2_post, r2_post),
        ):
            if value in (PHASE_DRY_AEROSOL, PHASE_NONE) and abs(radius) > 1.0e-30:
                errors.append(f"event {idx}: {label}={value} has nonzero hydro radius {radius}")
        pair = sorted((p1, p2))
        if code == 2 and pair != [PHASE_LIQUID, PHASE_ICE]:
            errors.append(
                f"event {idx}: trigger_code=2 expects liquid+ice phases, got {p1},{p2}"
            )
        if code == 3 and pair != [PHASE_ICE, PHASE_ICE]:
            errors.append(
                f"event {idx}: trigger_code=3 expects ice+ice phases, got {p1},{p2}"
            )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", help="Glob for SD_event_collision_NetCDF files")
    parser.add_argument("--coal-glob", help="Deprecated alias for --event-glob")
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
        help="Process trigger pair required at least --min-count-per-type times, e.g. 2:2 for significant riming",
    )
    parser.add_argument(
        "--require-event-type",
        type=int,
        action="append",
        default=[],
        help="Deprecated alias for --require-trigger-code",
    )
    parser.add_argument("--min-count-per-type", type=int, default=1)
    parser.add_argument(
        "--require-extended-geometry",
        action="store_true",
        help="Require both optional ice geometry and rime morphology fields",
    )
    parser.add_argument(
        "--require-ice-geometry",
        action="store_true",
        help="Require TRACK_COLD_OUTPUT_ICE_GEOMETRY optional fields",
    )
    parser.add_argument(
        "--require-rime-morphology",
        action="store_true",
        help="Require TRACK_COLD_OUTPUT_RIME_MORPHOLOGY optional fields",
    )
    parser.add_argument(
        "--forbid-optional-groups",
        action="store_true",
        help="Require optional ice/rime fields to be absent",
    )
    args = parser.parse_args()

    event_glob = args.event_glob or args.coal_glob
    if not event_glob:
        print("ERROR: --event-glob is required", file=sys.stderr)
        return 1

    paths = sorted(Path(p) for p in glob.glob(event_glob) if not p.endswith(".ids"))
    if not paths:
        print(f"ERROR: no event output files matched {event_glob}", file=sys.stderr)
        return 1

    counts: Counter[int] = Counter()
    pair_counts: Counter[tuple[int, int]] = Counter()
    context_counts: Counter[tuple[int, int, int, int, int]] = Counter()
    total_events = 0
    phase_errors: list[str] = []
    require_ice_geometry = args.require_extended_geometry or args.require_ice_geometry
    require_rime_morphology = args.require_extended_geometry or args.require_rime_morphology
    for path in paths:
        try:
            (
                trigger_code,
                trigger_level,
                phase1,
                phase2,
                phase1_post,
                phase2_post,
                radius1,
                radius2,
                radius1_post,
                radius2_post,
            ) = parse_event_file(
                path, require_ice_geometry, require_rime_morphology, args.forbid_optional_groups
            )
        except (RuntimeError, ValueError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        counts.update(trigger_code)
        pair_counts.update(zip(trigger_code, trigger_level))
        for idx, (code, level) in enumerate(zip(trigger_code, trigger_level)):
            if code < 1 or code > 3:
                print(f"ERROR: {path.name}: cold collision trigger_code={code} outside 1..3", file=sys.stderr)
                return 1
            if level not in (1, 2):
                print(f"ERROR: {path.name}: trigger_level={level} outside 1/2", file=sys.stderr)
                return 1
            context_counts[(code, idx, level, phase1[idx], phase2[idx])] += 1
        total_events += len(trigger_code)
        phase_errors.extend(
            f"{path.name}: {error}"
            for error in validate_phase_consistency(
                trigger_code, phase1, phase2, phase1_post, phase2_post, radius1, radius2, radius1_post, radius2_post
            )
        )

    if total_events == 0:
        print(f"ERROR: matched {len(paths)} event files but found zero events", file=sys.stderr)
        return 1

    if phase_errors:
        for error in phase_errors[:20]:
            print(f"ERROR: {error}", file=sys.stderr)
        if len(phase_errors) > 20:
            print(f"ERROR: {len(phase_errors) - 20} more phase consistency errors", file=sys.stderr)
        return 1

    required_codes = args.require_trigger_code + args.require_event_type
    missing = [code for code in required_codes if counts[code] < args.min_count_per_type]
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
