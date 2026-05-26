"""Optional NetCDF readers for GMD2026 analysis.

The analysis suite must also run on login nodes where NetCDF Python bindings
are not installed. Import failures therefore become warnings rather than hard
failures unless the caller uses strict mode.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


DIMENSION_NAMES = ("time", "sd", "nsd", "particle", "record", "event", "coal_event")
RADIUS_NAMES = ("sd_r", "r", "radius", "rad")
HEIGHT_NAMES = ("z", "height", "sd_z", "z_sd")
TIME_NAMES = ("time", "t", "time_s")
IF_COAL_NAMES = ("if_coal", "coal_flag", "coalesced")
FW_ID_NAMES = (("dm_id", "sd_id"),)
BW_ID_NAMES = (("pre_dmid", "pre_sdid"), ("pre_dm_id", "pre_sd_id"))


def _dataset_class():
    """Return netCDF4.Dataset or None."""
    try:
        from netCDF4 import Dataset  # type: ignore
    except Exception:
        return None
    return Dataset


def netcdf_available() -> bool:
    """Return True when netCDF4 is importable."""
    return _dataset_class() is not None


def _flatten(values: Any) -> list[float]:
    """Flatten array-like values into finite floats."""
    try:
        import numpy as np

        array = np.asarray(values, dtype=float).ravel()
        return [float(value) for value in array if math.isfinite(float(value))]
    except Exception:
        output: list[float] = []

        def walk(item: Any) -> None:
            try:
                iterator = iter(item)
            except TypeError:
                try:
                    value = float(item)
                except (TypeError, ValueError):
                    return
                if math.isfinite(value):
                    output.append(value)
                return
            for child in iterator:
                walk(child)

        walk(values)
        return output


def _flatten_int(values: Any) -> list[int]:
    """Flatten array-like values into ints."""
    return [int(value) for value in _flatten(values)]


def _limited_paths(paths: list[Path], max_files: int | None = None) -> tuple[list[Path], list[str]]:
    """Apply a conservative file-count limit for foreground debugging."""
    if max_files is None or max_files < 0 or len(paths) <= max_files:
        return paths, []
    return paths[:max_files], [f"NetCDF scan limited to {max_files} files out of {len(paths)}"]


def _netcdf_candidate_paths(paths: Iterable[Path]) -> list[Path]:
    """Return likely NetCDF files while excluding SCALE-SDM sidecar ID files."""
    return sorted(path for path in paths if path.is_file() and path.suffix != ".ids")


def _filename_model_time_s(path: Path) -> float | None:
    """Infer model time from SCALE-SDM per-rank selected-output filenames.

    GMD2026 selected-output files on SQUID may not have a ``.nc`` suffix and
    may not contain an internal time coordinate. In that case the reliable time
    token is the clock field in names such as
    ``SD_selected_NetCDF_00000101-011000.000.pe000034``. The trailing ``pe``
    rank must not be interpreted as model time.
    """
    match = re.search(r"_(\d{8})-(\d{6}(?:\.\d+)?)\.pe\d{6}(?:\.nc)?$", path.name)
    if match is None:
        return None
    clock = match.group(2)
    main, dot, fraction = clock.partition(".")
    main = main.zfill(6)
    try:
        hours = int(main[0:2])
        minutes = int(main[2:4])
        seconds = int(main[4:6])
        fractional_seconds = float(f"0.{fraction}") if dot and fraction else 0.0
    except ValueError:
        return None
    return float(hours * 3600 + minutes * 60 + seconds) + fractional_seconds


def _output_group_key(path: Path) -> str:
    """Return the per-output-time key shared by all ranks."""
    match = re.match(r"^(.*)\.pe\d{6}(?:\.nc)?$", str(path))
    return match.group(1) if match else str(path)


def _variable_length(variable: Any) -> int:
    """Return the first-dimension length for a NetCDF variable."""
    try:
        return int(variable.shape[0])
    except Exception:
        try:
            return len(variable)
        except Exception:
            return 1


def _read_flat_limited(variable: Any, max_records: int | None, chunk_size: int) -> list[float]:
    """Read a variable in first-dimension chunks instead of all at once."""
    length = _variable_length(variable)
    limit = length if max_records is None or max_records < 0 else min(length, max_records)
    output: list[float] = []
    if length <= 1:
        return _flatten(variable[:])[:limit]
    for start in range(0, limit, max(1, chunk_size)):
        stop = min(start + max(1, chunk_size), limit)
        output.extend(_flatten(variable[start:stop]))
        if len(output) >= limit:
            return output[:limit]
    return output


def _iter_first_dim_chunks(length: int, max_records: int | None, chunk_size: int):
    """Yield first-dimension slices up to a record limit."""
    limit = length if max_records is None or max_records < 0 else min(length, max_records)
    for start in range(0, limit, max(1, chunk_size)):
        yield slice(start, min(start + max(1, chunk_size), limit))


def _find_variable(variable_names: Iterable[str], candidates: Iterable[str]) -> str | None:
    """Find the first matching variable name from a candidate list."""
    names = set(variable_names)
    for candidate in candidates:
        if candidate in names:
            return candidate
    return None


def _find_pair_variables(variable_names: Iterable[str], preferred: Iterable[tuple[str, str]]) -> tuple[str, str] | None:
    """Find matching DM/SD ID variables, including numbered variants."""
    names = set(variable_names)
    for dm_name, sd_name in preferred:
        if dm_name in names and sd_name in names:
            return dm_name, sd_name
    numbered_prefixes: list[tuple[str, str]] = []
    if any(dm_name.startswith("pre_") or sd_name.startswith("pre_") for dm_name, sd_name in preferred):
        numbered_prefixes.append(("pre_dmid_", "pre_sdid_"))
    if any(dm_name.startswith("dm_id") or sd_name.startswith("sd_id") for dm_name, sd_name in preferred):
        numbered_prefixes.append(("dm_id_", "sd_id_"))
    for prefix_dm, prefix_sd in numbered_prefixes:
        dm_matches = sorted(name for name in names if name.startswith(prefix_dm))
        sd_matches = sorted(name for name in names if name.startswith(prefix_sd))
        if dm_matches and sd_matches:
            return dm_matches[0], sd_matches[0]
    return None


def choose_first_output_group(paths: list[Path]) -> list[Path]:
    """Select the earliest output-time group across per-rank NetCDF files."""
    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in paths:
        grouped[_output_group_key(path)].append(path)
    if not grouped:
        return []
    return sorted(grouped[sorted(grouped)[0]])


def inspect_netcdf_files(
    paths: list[Path],
    max_files: int | None = None,
    metadata_only: bool = False,
) -> tuple[dict[str, Any], list[str]]:
    """Inspect NetCDF dimensions and variables without requiring local outputs."""
    Dataset = _dataset_class()
    if Dataset is None:
        return {}, ["NetCDF library unavailable; NetCDF-specific metrics skipped"]
    if not paths:
        return {}, ["No NetCDF files were found for inspection"]

    paths, warnings = _limited_paths(paths, max_files)
    variable_names: set[str] = set()
    dimension_max: dict[str, int] = {}
    selected_record_count = 0
    coalescence_event_count = 0
    time_count = 0

    for path in paths:
        try:
            with Dataset(path, "r") as handle:
                variable_names.update(handle.variables.keys())
                for name, dim in handle.dimensions.items():
                    if name in DIMENSION_NAMES:
                        dimension_max[name] = max(dimension_max.get(name, 0), len(dim))
                for name in ("sd", "nsd", "particle", "record"):
                    if name in handle.dimensions:
                        selected_record_count += len(handle.dimensions[name])
                        break
                for name in ("event", "coal_event", "record"):
                    if name in handle.dimensions and "coal" in path.name.lower():
                        coalescence_event_count += len(handle.dimensions[name])
                        break
                time_var = _find_variable(handle.variables.keys(), TIME_NAMES)
                if time_var is not None:
                    if metadata_only:
                        time_count += _variable_length(handle.variables[time_var])
                    else:
                        time_count += len(_read_flat_limited(handle.variables[time_var], None, 100000))
        except Exception as exc:
            warnings.append(f"failed to inspect {path}: {exc}")

    return {
        "netcdf_variable_names": sorted(variable_names),
        "netcdf_dimensions": dimension_max,
        "coalescence_event_count": coalescence_event_count if coalescence_event_count else None,
        "selected_sd_record_count": selected_record_count if selected_record_count else None,
        "netcdf_time_count": time_count if time_count else None,
    }, warnings


def validate_selected_output(paths: list[Path], mode: str) -> tuple[dict[str, Any], list[str]]:
    """Validate FW or BW selected-output tracking ID variables."""
    Dataset = _dataset_class()
    if Dataset is None:
        return {"selected_output_has_required_ids": None}, ["NetCDF library unavailable; selected-output validation skipped"]
    if not paths:
        return {"selected_output_has_required_ids": None}, ["No selected-output NetCDF files found"]
    preferred = FW_ID_NAMES if mode.upper() == "FW" else BW_ID_NAMES
    warnings: list[str] = []
    first_group = choose_first_output_group(paths)
    found_pair: tuple[str, str] | None = None
    for path in first_group:
        try:
            with Dataset(path, "r") as handle:
                found_pair = _find_pair_variables(handle.variables.keys(), preferred)
                if found_pair is not None:
                    break
        except Exception as exc:
            warnings.append(f"failed to validate {path}: {exc}")
    if found_pair is None:
        warnings.append(f"{mode} selected output is missing required ID variables")
    return {
        "selected_output_has_required_ids": found_pair is not None,
        "selected_output_id_variables": ",".join(found_pair) if found_pair else None,
    }, warnings


def read_selected_pairs(
    paths: list[Path],
    mode: str,
    max_files: int | None = None,
    max_records: int | None = None,
    chunk_size: int = 100000,
    metadata_only: bool = False,
    first_group_only: bool = False,
) -> tuple[set[tuple[int, int]], int | None, list[str]]:
    """Read valid selected-output ID pairs from selected-output NetCDF files."""
    Dataset = _dataset_class()
    if Dataset is None:
        return set(), None, ["NetCDF library unavailable; selected-output pairs skipped"]
    if metadata_only:
        return set(), None, ["metadata-only mode skipped selected-output pair reads"]
    candidate_paths = _netcdf_candidate_paths(paths)
    paths_to_read = choose_first_output_group(candidate_paths) if first_group_only else candidate_paths
    if not paths_to_read:
        return set(), None, ["No selected-output NetCDF files found"]
    paths_to_read, limit_warnings = _limited_paths(paths_to_read, max_files)
    preferred = BW_ID_NAMES if mode.upper() == "BW" else FW_ID_NAMES
    warnings: list[str] = list(limit_warnings)
    pairs: set[tuple[int, int]] = set()
    valid_records = 0
    for path in paths_to_read:
        try:
            with Dataset(path, "r") as handle:
                pair_vars = _find_pair_variables(handle.variables.keys(), preferred)
                if pair_vars is None:
                    warnings.append(f"ID variables missing in {path}")
                    continue
                dm_var = handle.variables[pair_vars[0]]
                sd_var = handle.variables[pair_vars[1]]
                length = min(_variable_length(dm_var), _variable_length(sd_var))
                for selection in _iter_first_dim_chunks(length, max_records, chunk_size):
                    dm_values = _flatten_int(dm_var[selection])
                    sd_values = _flatten_int(sd_var[selection])
                    for dm_id, sd_id in zip(dm_values, sd_values):
                        if dm_id < 0 or sd_id < 0:
                            continue
                        valid_records += 1
                        pairs.add((dm_id, sd_id))
        except Exception as exc:
            warnings.append(f"failed to read selected pairs from {path}: {exc}")
    return pairs, valid_records, warnings


def _quantile(values: list[float], fraction: float) -> float | None:
    """Return a linear quantile for finite values."""
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def distribution_summary(values: list[float]) -> dict[str, Any]:
    """Return compact distribution metrics for finite values."""
    if not values:
        return {
            "count": None,
            "mean": None,
            "median": None,
            "q25": None,
            "q75": None,
            "q95": None,
            "maximum": None,
        }
    return {
        "count": len(values),
        "mean": sum(values) / len(values),
        "median": _quantile(values, 0.50),
        "q25": _quantile(values, 0.25),
        "q75": _quantile(values, 0.75),
        "q95": _quantile(values, 0.95),
        "maximum": max(values),
    }


def read_radius_records(
    case_dir: Path,
    max_files: int | None = None,
    max_records: int | None = None,
    chunk_size: int = 100000,
    metadata_only: bool = False,
    skip_heavy_netcdf: bool = False,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Read radius, height, if_coal, and time records from selected-output files."""
    Dataset = _dataset_class()
    if Dataset is None:
        return [], ["NetCDF library unavailable; radius records skipped"]
    if metadata_only:
        return [], ["metadata-only mode skipped radius/height variable reads"]
    if skip_heavy_netcdf:
        return [], ["skip-heavy-netcdf mode skipped radius/height variable reads"]
    paths = sorted(
        _netcdf_candidate_paths(case_dir.rglob("SD_selected_NetCDF_*")),
        key=lambda path: (_filename_model_time_s(path) is None, _filename_model_time_s(path) or 0.0, path.name),
    )
    if not paths:
        return [], ["No selected-output NetCDF files found"]
    paths, warnings = _limited_paths(paths, max_files)
    records: list[dict[str, Any]] = []
    record_order = 0
    for path in paths:
        file_time = _filename_model_time_s(path)
        try:
            with Dataset(path, "r") as handle:
                names = handle.variables.keys()
                radius_name = _find_variable(names, RADIUS_NAMES)
                if radius_name is None:
                    continue
                radius_var = handle.variables[radius_name]
                height_name = _find_variable(names, HEIGHT_NAMES)
                height_var = handle.variables[height_name] if height_name else None
                coal_name = _find_variable(names, IF_COAL_NAMES)
                coal_var = handle.variables[coal_name] if coal_name else None
                time_name = _find_variable(names, TIME_NAMES)
                time_var = handle.variables[time_name] if time_name else None
                pair_vars = _find_pair_variables(names, BW_ID_NAMES) or _find_pair_variables(names, FW_ID_NAMES)
                dm_var = handle.variables[pair_vars[0]] if pair_vars else None
                sd_var = handle.variables[pair_vars[1]] if pair_vars else None
                length = _variable_length(radius_var)
                all_times = None
                if time_var is not None and _variable_length(time_var) <= 10000:
                    all_times = _read_flat_limited(time_var, None, chunk_size)
                for selection in _iter_first_dim_chunks(length, max_records, chunk_size):
                    radii = _flatten(radius_var[selection])
                    heights = _flatten(height_var[selection]) if height_var is not None else []
                    coal = _flatten(coal_var[selection]) if coal_var is not None else []
                    dm_values = _flatten_int(dm_var[selection]) if dm_var is not None else []
                    sd_values = _flatten_int(sd_var[selection]) if sd_var is not None else []
                    time_values = []
                    if time_var is not None:
                        if _variable_length(time_var) == length:
                            time_values = _flatten(time_var[selection])
                        elif all_times is not None:
                            time_values = all_times
                    for index, radius in enumerate(radii):
                        records.append(
                            {
                                "time_s": time_values[index] if index < len(time_values) else (time_values[0] if len(time_values) == 1 else file_time),
                                "radius_m": radius,
                                "height_m": heights[index] if index < len(heights) else None,
                                "if_coal": coal[index] if index < len(coal) else None,
                                "dm_id": dm_values[index] if index < len(dm_values) else None,
                                "sd_id": sd_values[index] if index < len(sd_values) else None,
                                "record_order": record_order,
                                "source_file": path.name,
                            }
                        )
                        record_order += 1
        except Exception as exc:
            warnings.append(f"failed to read radius records from {path}: {exc}")
    if not records:
        warnings.append("No radius variable found in selected-output NetCDF files")
    return records, warnings
