"""Parse TPHT target-set handoff .ids files."""

from __future__ import annotations

import re
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
from typing import Any

from .table_utils import safe_ratio

KEY_VALUE_RE = re.compile(r"([A-Za-z0-9_]+)\s*=\s*([^\s,]+)")
RANK_RE = re.compile(r"pe(\d{6})")
TPHT_META_KEYS = ("PRC_NUM_X", "PRC_NUM_Y", "PRC_nprocs", "TPHT_ID_EPOCH_SEC")


def _parse_meta_value(value: str) -> Any:
    """Parse a metadata value from a TPHT header."""
    try:
        number = float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        return value.strip().strip("\"'")
    return int(number) if number.is_integer() else number


def _parse_header_metadata(text: str) -> dict[str, Any]:
    """Parse key-value or positional TPHT metadata headers."""
    metadata: dict[str, Any] = {}
    for key, value in KEY_VALUE_RE.findall(text):
        metadata[key] = _parse_meta_value(value)
    if metadata:
        return metadata

    parts = text.lstrip("#").split()
    if not parts or parts[0] != "TPHT_META":
        return metadata
    for key, value in zip(TPHT_META_KEYS, parts[1:]):
        metadata[key] = _parse_meta_value(value)
    return metadata


def _parse_id_file_summary(path: Path) -> dict[str, Any]:
    """Read one TPHT .ids file into compact counts and unique pairs."""
    metadata: dict[str, Any] = {}
    warnings: list[str] = []
    unique_pairs: set[tuple[int, int]] = set()
    raw_count = 0
    rank_match = RANK_RE.search(path.name)
    rank = int(rank_match.group(1)) if rank_match else None
    if rank is not None:
        metadata["rank"] = rank

    if not path.exists():
        return {
            "path": path,
            "metadata": metadata,
            "unique_pairs": unique_pairs,
            "raw_count": raw_count,
            "warnings": [f"ID file missing: {path}"],
        }

    with path.open("rb") as handle:
        for raw_line in handle:
            if not raw_line or raw_line[:1] == b"\n":
                continue
            if raw_line[:1] == b"#":
                metadata.update(_parse_header_metadata(raw_line.decode(errors="ignore").strip()))
                continue
            parts = raw_line.split()
            if len(parts) < 2:
                continue
            try:
                pair = (int(parts[0]), int(parts[1]))
            except ValueError:
                warnings.append(f"unparseable ID record in {path.name}: {raw_line[:80]!r}")
                continue
            unique_pairs.add(pair)
            raw_count += 1
    return {
        "path": path,
        "metadata": metadata,
        "unique_pairs": unique_pairs,
        "raw_count": raw_count,
        "warnings": warnings,
    }


def summarize_id_files(paths: list[Path], workers: int = 1) -> tuple[dict[str, Any], list[tuple[int, int]], list[str]]:
    """Summarize a collection of TPHT .ids files."""
    warnings: list[str] = []
    unique_pairs: set[tuple[int, int]] = set()
    raw_id_records = 0
    per_rank: Counter[int] = Counter()
    metadata: dict[str, Any] = {}

    if workers > 1 and len(paths) > 1:
        with Pool(processes=workers) as pool:
            parsed_rows = pool.map(_parse_id_file_summary, paths)
    else:
        parsed_rows = [_parse_id_file_summary(path) for path in paths]

    for parsed in parsed_rows:
        warnings.extend(parsed["warnings"])
        unique_pairs.update(parsed["unique_pairs"])
        raw_id_records += int(parsed.get("raw_count") or 0)
        rank = parsed["metadata"].get("rank")
        if rank is not None:
            per_rank[int(rank)] += int(parsed.get("raw_count") or 0)
        for key, value in parsed["metadata"].items():
            metadata.setdefault(key, value)

    mean_rank_count = sum(per_rank.values()) / len(per_rank) if per_rank else None
    rank_imbalance_ratio = safe_ratio(max(per_rank.values()) if per_rank else None, mean_rank_count)
    summary = {
        "raw_id_records": raw_id_records,
        "unique_pairs": len(unique_pairs),
        "deduplicated_pairs": len(unique_pairs),
        "ids_per_rank": dict(sorted(per_rank.items())),
        "dedup_reduction_ratio": safe_ratio(raw_id_records, len(unique_pairs)),
        "rank_imbalance_ratio": rank_imbalance_ratio,
        "PRC_NUM_X": metadata.get("PRC_NUM_X"),
        "PRC_NUM_Y": metadata.get("PRC_NUM_Y"),
        "PRC_nprocs": metadata.get("PRC_nprocs"),
        "TPHT_ID_EPOCH_SEC": metadata.get("TPHT_ID_EPOCH_SEC"),
    }
    return summary, sorted(unique_pairs), warnings

