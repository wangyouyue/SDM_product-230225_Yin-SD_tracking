"""Optional NetCDF readers for SCALE-SDM tracking post-processing.

The helpers prefer ``netCDF4`` when it is available.  Local smoke checks often
run on lightweight Python environments, so selected-output ID reads fall back
to ``ncdump`` when Python NetCDF bindings are absent.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Iterable

FW_ID_NAMES = (("dm_id", "sd_id"),)
BW_ID_NAMES = (("pre_dmid", "pre_sdid"), ("pre_dm_id", "pre_sd_id"))
NETCDF_VARIABLE_RE = re.compile(
    r"^\s*(?:byte|char|short|int|int64|float|double)\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(",
    re.MULTILINE,
)
NETCDF_DIMENSION_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\d+)\s*;", re.MULTILINE)


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


def ncdump_available() -> bool:
    """Return True when an ``ncdump`` executable is available."""
    return shutil.which("ncdump") is not None


def _flatten_int(values: Any) -> list[int]:
    """Flatten array-like values into ints."""
    try:
        import numpy as np

        return [int(value) for value in np.asarray(values).ravel()]
    except Exception:
        output: list[int] = []

        def walk(item: Any) -> None:
            try:
                iterator = iter(item)
            except TypeError:
                try:
                    output.append(int(item))
                except (TypeError, ValueError):
                    return
                return
            for child in iterator:
                walk(child)

        walk(values)
        return output


def _variable_length(variable: Any) -> int:
    """Return the first-dimension length for a NetCDF variable."""
    try:
        return int(variable.shape[0])
    except Exception:
        try:
            return len(variable)
        except Exception:
            return 1


def _iter_first_dim_chunks(length: int, max_records: int | None, chunk_size: int):
    """Yield first-dimension slices up to a record limit."""
    limit = length if max_records is None or max_records < 0 else min(length, max_records)
    for start in range(0, limit, max(1, chunk_size)):
        yield slice(start, min(start + max(1, chunk_size), limit))


def _netcdf_candidate_paths(paths: Iterable[Path]) -> list[Path]:
    """Return likely NetCDF files while excluding SCALE-SDM sidecar ID files."""
    return sorted(path for path in paths if path.is_file() and path.suffix != ".ids")


def _filename_model_time_s(path: Path) -> float | None:
    """Infer model time from SCALE-SDM per-rank selected-output filenames."""
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


def choose_output_group(paths: list[Path], latest: bool = False) -> list[Path]:
    """Select the earliest or latest output-time group across per-rank files."""
    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in paths:
        grouped[_output_group_key(path)].append(path)
    if not grouped:
        return []
    keys = sorted(
        grouped,
        key=lambda key: (
            _filename_model_time_s(Path(key + ".pe000000")) is None,
            _filename_model_time_s(Path(key + ".pe000000")) or 0.0,
            key,
        ),
    )
    return sorted(grouped[keys[-1 if latest else 0]])


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


def _ncdump_header(path: Path) -> tuple[set[str], dict[str, int]]:
    """Read variable and dimension names with ncdump."""
    output = subprocess.check_output(["ncdump", "-h", str(path)], text=True)
    variables = set(NETCDF_VARIABLE_RE.findall(output))
    dimensions = {name: int(size) for name, size in NETCDF_DIMENSION_RE.findall(output)}
    return variables, dimensions


def _ncdump_variable(path: Path, variable_name: str) -> list[int]:
    """Read one integer variable with ncdump."""
    output = subprocess.check_output(["ncdump", "-v", variable_name, str(path)], text=True)
    match = re.search(r"\n\s*" + re.escape(variable_name) + r"\s*=\s*(.*?);", output, re.S)
    if match is None:
        raise RuntimeError(f"{variable_name} not found in {path}")
    return [int(value) for value in re.findall(r"-?\d+", match.group(1))]


def inspect_netcdf_files(paths: list[Path], max_files: int | None = None) -> tuple[dict[str, Any], list[str]]:
    """Inspect NetCDF dimensions and variables."""
    paths = _netcdf_candidate_paths(paths)
    if max_files is not None and max_files >= 0:
        paths = paths[:max_files]
    if not paths:
        return {}, ["No NetCDF files were found for inspection"]

    Dataset = _dataset_class()
    variable_names: set[str] = set()
    dimension_max: dict[str, int] = {}
    warnings: list[str] = []

    if Dataset is not None:
        for path in paths:
            try:
                with Dataset(path, "r") as handle:
                    variable_names.update(handle.variables.keys())
                    for name, dim in handle.dimensions.items():
                        dimension_max[name] = max(dimension_max.get(name, 0), len(dim))
            except Exception as exc:
                warnings.append(f"failed to inspect {path}: {exc}")
    elif ncdump_available():
        for path in paths:
            try:
                file_variables, file_dimensions = _ncdump_header(path)
            except Exception as exc:
                warnings.append(f"failed to inspect {path} with ncdump: {exc}")
                continue
            variable_names.update(file_variables)
            for name, size in file_dimensions.items():
                dimension_max[name] = max(dimension_max.get(name, 0), size)
    else:
        warnings.append("NetCDF library and ncdump unavailable; NetCDF-specific metrics skipped")

    return {
        "netcdf_variable_names": sorted(variable_names),
        "netcdf_dimensions": dimension_max,
    }, warnings


def validate_selected_output(paths: list[Path], mode: str) -> tuple[dict[str, Any], list[str]]:
    """Validate FW or BW selected-output tracking ID variables."""
    candidate_paths = _netcdf_candidate_paths(paths)
    if not candidate_paths:
        return {"selected_output_has_required_ids": None}, ["No selected-output NetCDF files found"]
    preferred = FW_ID_NAMES if mode.upper() == "FW" else BW_ID_NAMES
    warnings: list[str] = []
    found_pair: tuple[str, str] | None = None
    first_group = choose_output_group(candidate_paths)
    Dataset = _dataset_class()

    for path in first_group:
        try:
            if Dataset is not None:
                with Dataset(path, "r") as handle:
                    found_pair = _find_pair_variables(handle.variables.keys(), preferred)
            elif ncdump_available():
                variable_names, _ = _ncdump_header(path)
                found_pair = _find_pair_variables(variable_names, preferred)
            else:
                warnings.append("NetCDF library and ncdump unavailable; selected-output validation skipped")
                break
        except Exception as exc:
            warnings.append(f"failed to validate {path}: {exc}")
        if found_pair is not None:
            break

    if found_pair is None and not warnings:
        warnings.append(f"{mode} selected output is missing required ID variables")
    return {
        "selected_output_has_required_ids": found_pair is not None,
        "selected_output_id_variables": ",".join(found_pair) if found_pair else None,
    }, warnings


def _read_selected_pairs_one(args: tuple[Path, tuple[tuple[str, str], ...], int | None, int]) -> tuple[set[tuple[int, int]], int, list[str]]:
    """Read valid selected-output ID pairs from one NetCDF file."""
    path, preferred, max_records, chunk_size = args
    warnings: list[str] = []
    pairs: set[tuple[int, int]] = set()
    valid_records = 0
    Dataset = _dataset_class()

    try:
        if Dataset is not None:
            with Dataset(path, "r") as handle:
                pair_vars = _find_pair_variables(handle.variables.keys(), preferred)
                if pair_vars is None:
                    warnings.append(f"ID variables missing in {path}")
                    return pairs, valid_records, warnings
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
            return pairs, valid_records, warnings

        if not ncdump_available():
            return pairs, valid_records, ["NetCDF library and ncdump unavailable; selected-output pairs skipped"]

        variable_names, _ = _ncdump_header(path)
        pair_vars = _find_pair_variables(variable_names, preferred)
        if pair_vars is None:
            warnings.append(f"ID variables missing in {path}")
            return pairs, valid_records, warnings
        dm_values = _ncdump_variable(path, pair_vars[0])
        sd_values = _ncdump_variable(path, pair_vars[1])
        if max_records is not None and max_records >= 0:
            dm_values = dm_values[:max_records]
            sd_values = sd_values[:max_records]
        for dm_id, sd_id in zip(dm_values, sd_values):
            if dm_id < 0 or sd_id < 0:
                continue
            valid_records += 1
            pairs.add((dm_id, sd_id))
    except Exception as exc:
        warnings.append(f"failed to read selected pairs from {path}: {exc}")
    return pairs, valid_records, warnings


def read_selected_pairs(
    paths: list[Path],
    mode: str,
    max_files: int | None = None,
    max_records: int | None = None,
    chunk_size: int = 100000,
    first_group_only: bool = False,
    latest_group_only: bool = False,
    workers: int = 1,
) -> tuple[set[tuple[int, int]], int | None, list[str]]:
    """Read valid selected-output ID pairs from selected-output NetCDF files."""
    candidate_paths = _netcdf_candidate_paths(paths)
    if latest_group_only:
        paths_to_read = choose_output_group(candidate_paths, latest=True)
    elif first_group_only:
        paths_to_read = choose_output_group(candidate_paths)
    else:
        paths_to_read = candidate_paths
    if not paths_to_read:
        return set(), None, ["No selected-output NetCDF files found"]
    if max_files is not None and max_files >= 0:
        paths_to_read = paths_to_read[:max_files]

    preferred = BW_ID_NAMES if mode.upper() == "BW" else FW_ID_NAMES
    warnings: list[str] = []
    pairs: set[tuple[int, int]] = set()
    valid_records = 0
    worker_args = [(path, preferred, max_records, chunk_size) for path in paths_to_read]
    if workers > 1 and len(worker_args) > 1:
        with Pool(processes=workers) as pool:
            results = pool.map(_read_selected_pairs_one, worker_args)
    else:
        results = [_read_selected_pairs_one(args) for args in worker_args]

    for file_pairs, file_valid_records, file_warnings in results:
        pairs.update(file_pairs)
        valid_records += file_valid_records
        warnings.extend(file_warnings)
    return pairs, valid_records, warnings

