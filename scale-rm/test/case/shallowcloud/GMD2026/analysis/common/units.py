"""Unit conversions used by GMD2026 analysis and plotting scripts."""

from __future__ import annotations

BYTES_PER_MIB = 1024.0**2
BYTES_PER_GIB = 1024.0**3
MICROMETER_IN_M = 1.0e-6


def bytes_to_gib(value: float | int | None) -> float | None:
    """Convert bytes to GiB while preserving missing values."""
    return None if value is None else float(value) / BYTES_PER_GIB


def kb_to_mib(value: float | int | None) -> float | None:
    """Convert KiB-like resident-set values to MiB."""
    return None if value is None else float(value) / 1024.0


def meters_to_micrometers(value: float | int | None) -> float | None:
    """Convert meters to micrometers."""
    return None if value is None else float(value) / MICROMETER_IN_M

