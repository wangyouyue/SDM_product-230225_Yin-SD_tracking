#!/usr/bin/env python3
"""Create a legacy-format cold SD restart by removing the v1 tracking block."""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path


MARKER_PREFIX = b"SDM_COLD_TRACKING_RESTART_V"


def read_unformatted_records(path: Path) -> list[bytes]:
    records: list[bytes] = []
    with path.open("rb") as handle:
        while True:
            head = handle.read(4)
            if not head:
                break
            if len(head) != 4:
                raise SystemExit(f"{path}: truncated leading record marker")
            (length,) = struct.unpack(">i", head)
            payload = handle.read(length)
            tail = handle.read(4)
            if len(payload) != length or len(tail) != 4:
                raise SystemExit(f"{path}: truncated record payload")
            if tail != head:
                raise SystemExit(f"{path}: record marker mismatch")
            records.append(head + payload + tail)
    return records


def payload(record: bytes) -> bytes:
    (length,) = struct.unpack(">i", record[:4])
    return record[4 : 4 + length]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    records = read_unformatted_records(input_path)

    marker_index = None
    for index in range(len(records) - 1, -1, -1):
        if payload(records[index]).rstrip().startswith(MARKER_PREFIX):
            marker_index = index
            break
    if marker_index is None:
        print(f"{input_path}: cold tracking marker not found", file=sys.stderr)
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        for record in records[:marker_index]:
            handle.write(record)

    print(f"input_records={len(records)}")
    print(f"output_records={marker_index}")
    print(f"legacy_restart={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
