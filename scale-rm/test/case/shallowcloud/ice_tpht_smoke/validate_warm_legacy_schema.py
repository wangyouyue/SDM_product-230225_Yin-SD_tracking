#!/usr/bin/env python3
"""Validate warm-mode legacy SDM tracking output names."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")

REQUIRED_WARM_COAL_FIELDS = (
    "event_multiplicity",
    "x",
    "y",
    "z",
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


def ncdump(path: Path, header_only: bool = False) -> str:
    cmd = ["ncdump", "-h", str(path)] if header_only else ["ncdump", str(path)]
    try:
        result = subprocess.run(
            cmd,
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


def require_field(header: str, path: Path, name: str) -> None:
    if f" {name}(" not in header:
        raise ValueError(f"{path}: missing required field {name!r}")


def forbid_field(header: str, path: Path, name: str) -> None:
    if f" {name}(" in header:
        raise ValueError(f"{path}: forbidden cold-only/replacement field {name!r} is present")


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coal-glob", required=True, help="Glob for SD_coal_output_NetCDF files")
    parser.add_argument("--ordinary-glob", required=True, help="Glob for SD_all_NetCDF files")
    args = parser.parse_args()

    coal_paths = sorted(Path(p) for p in glob.glob(args.coal_glob))
    ordinary_paths = sorted(Path(p) for p in glob.glob(args.ordinary_glob))
    if not coal_paths:
        print(f"ERROR: no warm coal files matched {args.coal_glob}", file=sys.stderr)
        return 1
    if not ordinary_paths:
        print(f"ERROR: no warm ordinary files matched {args.ordinary_glob}", file=sys.stderr)
        return 1

    try:
        for path in coal_paths:
            text = ncdump(path)
            for field in REQUIRED_WARM_COAL_FIELDS:
                require_field(text, path, field)
            for forbidden in (
                "num_col",
                "trigger_code",
                "sd_event_mask",
                "sd_diag_mask",
                "sd_phase_change_flag",
            ):
                forbid_field(text, path, forbidden)
            read_int_variable(text, "event_multiplicity")

        for path in ordinary_paths:
            header = ncdump(path, header_only=True)
            require_field(header, path, "if_coal")
            for forbidden in ("sd_event_mask", "sd_diag_mask", "sd_phase_change_flag"):
                forbid_field(header, path, forbidden)
    except (RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"coal_files={len(coal_paths)}")
    print(f"ordinary_files={len(ordinary_paths)}")
    print("warm_legacy_schema=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
