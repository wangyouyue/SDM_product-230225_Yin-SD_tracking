#!/usr/bin/env python3
"""Validate cold-mode ordinary or selected SDM tracking output fields."""

from __future__ import annotations

import argparse
import glob
import subprocess
import sys
from pathlib import Path


REQUIRED_COLD_FIELDS = (
    "sd_event_mask",
    "sd_event_sig_mask",
    "sd_diag_mask",
    "sd_phase_change_flag",
    "sd_spatial_visit_flag",
    "sd_liq_radius_max_interval",
    "sd_ice_rvol_max_interval",
    "sd_mixed_rvol_max_interval",
    "sd_rime_mass_max_interval",
    "sd_rime_frac_max_interval",
    "sd_nmono_max_interval",
    "sd_aspect_ratio_max_interval",
)

FORBIDDEN_LEGACY_FIELDS = ("if_coal",)


def ncdump_header(path: Path) -> str:
    try:
        result = subprocess.run(
            ["ncdump", "-h", str(path)],
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


def validate_file(path: Path) -> None:
    header = ncdump_header(path)
    missing = [name for name in REQUIRED_COLD_FIELDS if f" {name}(" not in header]
    if missing:
        raise ValueError(f"{path}: missing required cold tracking fields: {missing}")

    forbidden = [name for name in FORBIDDEN_LEGACY_FIELDS if f" {name}(" in header]
    if forbidden:
        raise ValueError(f"{path}: forbidden warm legacy fields are present: {forbidden}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Glob for SD_all_NetCDF or SD_selected_NetCDF files")
    parser.add_argument("--label", default="cold_tracking", help="Label printed in the success summary")
    args = parser.parse_args()

    paths = sorted(Path(p) for p in glob.glob(args.glob) if not p.endswith(".ids"))
    if not paths:
        print(f"ERROR: no files matched {args.glob}", file=sys.stderr)
        return 1

    try:
        for path in paths:
            validate_file(path)
    except (RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"{args.label}_files={len(paths)}")
    print(f"{args.label}_schema=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
