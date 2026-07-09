#!/usr/bin/env python3
"""Validate that cold ordinary/selected/history output has expected sd_event_mask bits."""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")
MASK_BLOCK_TEMPLATE = r"(?:^|\n)\s*({name}(?:_\d+)?)\s*=\s*(.*?);"


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


def read_masks(text: str, variable: str) -> list[int]:
    mask_block_re = re.compile(MASK_BLOCK_TEMPLATE.format(name=re.escape(variable)), re.S)
    values: list[int] = []
    for _name, data in mask_block_re.findall(text):
        values.extend(int(value) for value in INT_RE.findall(data))
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glob", required=True, help="Glob for SD_all/SD_selected/SD_*_history files")
    parser.add_argument(
        "--mask-variable",
        default="sd_event_mask",
        choices=("sd_event_mask", "sd_event_sig_mask"),
        help="Interval process mask variable to inspect",
    )
    parser.add_argument(
        "--require-bit",
        type=int,
        action="append",
        default=[],
        help="sd_event_mask bit value that must appear at least --min-count-per-bit times",
    )
    parser.add_argument("--min-count-per-bit", type=int, default=1)
    args = parser.parse_args()

    paths = sorted(Path(p) for p in glob.glob(args.glob) if not p.endswith(".ids"))
    if not paths:
        print(f"ERROR: no files matched {args.glob}", file=sys.stderr)
        return 1

    counts: Counter[int] = Counter()
    mask_values = 0
    for path in paths:
        try:
            masks = read_masks(ncdump(path), args.mask_variable)
        except RuntimeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        mask_values += len(masks)
        for bit in args.require_bit:
            counts[bit] += sum(1 for value in masks if value & bit)

    if mask_values == 0:
        print(f"ERROR: matched {len(paths)} files but found no {args.mask_variable} data", file=sys.stderr)
        return 1

    missing = [bit for bit in args.require_bit if counts[bit] < args.min_count_per_bit]
    if missing:
        print(
            f"ERROR: missing required {args.mask_variable} bits: "
            + ", ".join(f"{bit}={counts[bit]}" for bit in missing),
            file=sys.stderr,
        )
        return 1

    print(f"files={len(paths)}")
    print(f"mask_values={mask_values}")
    print("bit_counts=" + ",".join(f"{bit}:{counts[bit]}" for bit in args.require_bit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
