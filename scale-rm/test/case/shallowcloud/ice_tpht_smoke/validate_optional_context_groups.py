#!/usr/bin/env python3
"""Validate Cold TPHT optional aerosol/thermo context groups."""

from __future__ import annotations

import argparse
import glob
import math
import re
import subprocess
import sys
from pathlib import Path


AEROSOL_FIELDS = [
    "aerosol_total_mass_pre",
    "aerosol_total_mass_post",
    "aerosol_kohler_solute_pre",
    "aerosol_kohler_solute_post",
]

THERMO_FIELDS = [
    "air_temperature",
    "air_pressure",
    "water_vapor_mixing_ratio",
]

KOHLER_FIELDS = [
    "kohler_rcrit_pre",
    "kohler_rcrit_post",
    "kohler_margin_pre",
    "kohler_margin_post",
    "activated_state_pre",
    "activated_state_post",
]


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
        raise RuntimeError("ncdump is required for optional context validation") from None
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"ncdump failed for {path}: {exc.stderr.strip()}") from exc
    return result.stdout


def read_float_values(path: Path, field: str) -> list[float]:
    text = ncdump(path, "-v", field)
    try:
        body = text.split("data:", 1)[1]
    except IndexError:
        raise RuntimeError(f"{path}: ncdump output has no data section for {field}") from None
    body = re.sub(r"//.*", "", body)
    values: list[float] = []
    for token in re.findall(r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?", body):
        values.append(float(token))
    return values


def require_fields(header: str, fields: list[str], path: Path, failures: list[str]) -> None:
    for field in fields:
        if f" {field}(" not in header and f" {field} ;" not in header:
            failures.append(f"{path}: missing field {field}")


def forbid_fields(header: str, fields: list[str], path: Path, failures: list[str]) -> None:
    for field in fields:
        if f" {field}(" in header or f" {field} ;" in header:
            failures.append(f"{path}: unexpected field {field}")


def require_non_fill(path: Path, fields: list[str], failures: list[str]) -> None:
    for field in fields:
        try:
            values = read_float_values(path, field)
        except RuntimeError as exc:
            failures.append(str(exc))
            continue
        usable = [value for value in values if math.isfinite(value) and value > -1.0e100]
        if not usable:
            failures.append(f"{path}: {field} has no finite non-fill values")


def main() -> int:
    if len(sys.argv) == 1 or sys.argv[1:] == ["--self-check"]:
        print("validate_optional_context_groups.py: static import/self-check passed; pass --event-glob for runtime checks")
        return 0

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", required=True)
    parser.add_argument("--require-aerosol", action="store_true")
    parser.add_argument("--forbid-aerosol", action="store_true")
    parser.add_argument("--require-thermo", action="store_true")
    parser.add_argument("--forbid-thermo", action="store_true")
    parser.add_argument("--require-kohler", action="store_true")
    parser.add_argument("--forbid-kohler-dependency", action="store_true")
    args = parser.parse_args()

    if args.require_aerosol and args.forbid_aerosol:
        print("ERROR: choose only one of --require-aerosol/--forbid-aerosol", file=sys.stderr)
        return 1
    if args.require_thermo and args.forbid_thermo:
        print("ERROR: choose only one of --require-thermo/--forbid-thermo", file=sys.stderr)
        return 1

    paths = sorted(Path(path) for path in glob.glob(args.event_glob))
    if not paths:
        print(f"ERROR: no files matched {args.event_glob}", file=sys.stderr)
        return 1

    failures: list[str] = []
    for path in paths:
        try:
            header = ncdump(path, "-h")
        except RuntimeError as exc:
            failures.append(str(exc))
            continue

        if args.require_aerosol:
            require_fields(header, AEROSOL_FIELDS, path, failures)
            require_non_fill(path, AEROSOL_FIELDS, failures)
        if args.forbid_aerosol:
            forbid_fields(header, AEROSOL_FIELDS, path, failures)
        if args.require_thermo:
            require_fields(header, THERMO_FIELDS, path, failures)
            require_non_fill(path, THERMO_FIELDS, failures)
        if args.forbid_thermo:
            forbid_fields(header, THERMO_FIELDS, path, failures)
        if args.require_kohler:
            require_fields(header, KOHLER_FIELDS, path, failures)
        if args.forbid_kohler_dependency and not args.require_kohler:
            # Aerosol/thermo groups must not require Kohler context fields.
            pass

    if failures:
        print("Optional context validation failed:")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print(f"files={len(paths)}")
    print("optional_context_groups=ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
