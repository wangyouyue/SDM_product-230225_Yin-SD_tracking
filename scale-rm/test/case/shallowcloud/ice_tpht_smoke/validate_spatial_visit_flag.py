#!/usr/bin/env python3
"""Validate cold Level-2 spatial-visit interval flags."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from pathlib import Path


def ncdump_variable(path: Path, field: str) -> str:
    try:
        result = subprocess.run(
            ["ncdump", "-v", field, str(path)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise RuntimeError("ncdump is required for this validator") from None
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"ncdump failed for {path}: {exc.stderr.strip()}") from exc

    try:
        return result.stdout.split("data:", 1)[1]
    except IndexError:
        raise RuntimeError(f"{path}: ncdump output has no data section") from None


def read_int_values(path: Path, field: str) -> list[int]:
    text = ncdump_variable(path, field)
    body = re.sub(r"//.*", "", text)
    return [int(value) for value in re.findall(r"-?\d+", body)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Glob for SD_all/SD_selected NetCDF files")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--expect-any", action="store_true", help="Require at least one visited SD")
    group.add_argument("--expect-zero", action="store_true", help="Require all flags to be zero")
    args = parser.parse_args()

    paths = sorted(Path(path) for path in glob.glob(args.glob) if not path.endswith(".ids"))
    if not paths:
        print(f"ERROR: no files matched {args.glob}", file=sys.stderr)
        return 1

    total = 0
    visited = 0
    try:
        for path in paths:
            values = read_int_values(path, "sd_spatial_visit_flag")
            total += len(values)
            visited += sum(1 for value in values if value != 0)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.expect_any and visited <= 0:
        print("ERROR: expected at least one sd_spatial_visit_flag hit", file=sys.stderr)
        return 1
    if args.expect_zero and visited != 0:
        print(f"ERROR: expected zero hits, found {visited}", file=sys.stderr)
        return 1

    print(f"files={len(paths)}")
    print(f"values={total}")
    print(f"visited={visited}")
    print("spatial_visit_flag=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
