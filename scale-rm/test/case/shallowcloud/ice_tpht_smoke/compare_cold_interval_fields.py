#!/usr/bin/env python3
"""Compare cold interval tracking fields across two SD NetCDF output sets."""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path


FIELDS = (
    "sd_event_mask",
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


def netcdf_paths(pattern: str) -> dict[str, Path]:
    paths = [
        Path(path)
        for path in sorted(glob.glob(pattern))
        if not path.endswith(".ids") and os.path.isfile(path)
    ]
    if not paths:
        raise SystemExit(f"no NetCDF files matched: {pattern}")
    return {path.name: path for path in paths}


def dump_field(path: Path, field: str) -> str:
    try:
        result = subprocess.run(
            ["ncdump", "-v", field, str(path)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise SystemExit("ncdump is required for this validator") from None
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"ncdump failed for {path}:{field}: {exc.stderr.strip()}") from exc

    try:
        return result.stdout.split("data:", 1)[1].strip()
    except IndexError:
        raise SystemExit(f"{path}:{field}: ncdump output has no data section") from None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left-glob", required=True)
    parser.add_argument("--right-glob", required=True)
    args = parser.parse_args()

    left = netcdf_paths(args.left_glob)
    right = netcdf_paths(args.right_glob)
    if set(left) != set(right):
        print("file basenames differ", file=sys.stderr)
        print(f"left={sorted(left)}", file=sys.stderr)
        print(f"right={sorted(right)}", file=sys.stderr)
        return 1

    compared = 0
    for name in sorted(left):
        for field in FIELDS:
            left_dump = dump_field(left[name], field)
            right_dump = dump_field(right[name], field)
            if left_dump != right_dump:
                print(f"{name}:{field}: data mismatch", file=sys.stderr)
                return 1
            compared += 1

    print(f"files={len(left)}")
    print(f"fields_compared={compared}")
    print("cold_interval_fields_exact=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
