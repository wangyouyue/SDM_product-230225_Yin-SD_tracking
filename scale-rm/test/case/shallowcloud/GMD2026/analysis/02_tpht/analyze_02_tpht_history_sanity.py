#!/usr/bin/env python3
"""Analyze BW ``history.pe*`` fields for TPHT cloud-physics sanity checks.

This GMD2026 script checks whether the meteorological fields in the
``02_tpht_3d_interest_70min/bw_reconstruction`` run look physically plausible.
It is independent of TPHT predecessor-chain reconstruction: the inputs are the
SCALE-SDM history files, and the outputs are compact time series, vertical
profiles, and diagnostic figures.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import analysis_options, build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir, warn_or_raise, warning_text  # noqa: E402
from common.plot_style import PALETTE, configure_matplotlib, figure_size, has_any_number, no_data_panel, read_table, save_figure  # noqa: E402
from common.table_utils import safe_float, write_table_bundle  # noqa: E402


DEFAULT_GROUP = "02_tpht_3d_interest_70min"
DEFAULT_CASE = "bw_reconstruction"
DEFAULT_PROFILE_TIMES_MIN = "1,5,10,30,50,70"
DEFAULT_TIMESERIES_VARIABLES = (
    "QHYD_sd",
    "QC_sd",
    "QR_sd",
    "W",
    "RH",
    "T",
    "QV",
    "LWPT",
    "RAIN_ACC_sd",
    "RAIN",
    "TKE_RS",
    "TKE_SMG",
)
DEFAULT_PROFILE_VARIABLES = (
    "QHYD_sd",
    "QC_sd",
    "QR_sd",
    "W",
    "RH",
    "T",
    "QV",
    "TKE_RS",
    "TKE_SMG",
)
CLOUD_VARIABLE_CANDIDATES = ("QHYD_sd", "QHYD", "QC_sd", "QC")
TIME_NAMES = ("time", "TIME", "t", "time_s")
VERTICAL_NAMES = ("z", "Z", "zh", "ZH", "height", "HEIGHT", "altitude", "lev", "level")
HORIZONTAL_MARKERS = ("x", "y", "lon", "lat", "xi", "eta")
MIXING_RATIO_PREFIXES = ("QHYD", "QC", "QR", "QV", "QTOT")

TIME_SERIES_COLUMNS = [
    "time_s",
    "time_min",
    "history_file_count",
    "cloud_variable",
    "cloud_threshold_kgkg",
    "cloud_fraction",
    "cloud_base_m",
    "cloud_top_m",
    "cloud_depth_m",
    "qhyd_mean_gkg",
    "qhyd_max_gkg",
    "qc_mean_gkg",
    "qc_max_gkg",
    "qr_mean_gkg",
    "qr_max_gkg",
    "w_mean_ms",
    "w_max_ms",
    "w_min_ms",
    "rh_mean",
    "rh_max",
    "t_mean_k",
    "qv_mean_gkg",
    "lwpt_mean",
    "lwpt_max",
    "rain_acc_mean",
    "rain_acc_max",
    "rain_mean",
    "rain_max",
    "tke_rs_max",
    "tke_smg_max",
    "warnings",
]

PROFILE_COLUMNS = [
    "time_s",
    "time_min",
    "profile_time_request_min",
    "variable",
    "z_index",
    "z_m",
    "domain_mean",
    "domain_max",
    "domain_min",
    "domain_count",
    "display_value",
    "display_unit",
    "warnings",
]

TIME_HEIGHT_COLUMNS = [
    "time_s",
    "time_min",
    "z_index",
    "z_m",
    "qhyd_mean_gkg",
    "qhyd_max_gkg",
    "grid_count",
    "warnings",
]

METADATA_COLUMNS = [
    "case_dir",
    "history_file_count",
    "processed_history_file_count",
    "time_count",
    "selected_profile_times_min",
    "available_variables",
    "timeseries_variables_used",
    "profile_variables_used",
    "cloud_variable",
    "cloud_threshold_kgkg",
    "netcdf_backend",
    "warnings",
]


@dataclass
class ScalarAccumulator:
    """Incrementally accumulate scalar statistics without storing full fields."""

    total: float = 0.0
    count: int = 0
    minimum: float | None = None
    maximum: float | None = None

    def update(self, values: Any) -> None:
        """Update the accumulator from an array-like object."""
        np = _numpy()
        array = _finite_array(values, np)
        if array.size == 0:
            return
        self.total += float(array.sum())
        self.count += int(array.size)
        current_min = float(array.min())
        current_max = float(array.max())
        self.minimum = current_min if self.minimum is None else min(self.minimum, current_min)
        self.maximum = current_max if self.maximum is None else max(self.maximum, current_max)

    @property
    def mean(self) -> float | None:
        """Return the accumulated mean."""
        return self.total / self.count if self.count > 0 else None


@dataclass
class ProfileAccumulator:
    """Incrementally accumulate vertical profiles over horizontal/rank chunks."""

    sums: list[float] = field(default_factory=list)
    counts: list[int] = field(default_factory=list)
    maxima: list[float | None] = field(default_factory=list)
    minima: list[float | None] = field(default_factory=list)

    def ensure(self, length: int) -> None:
        """Ensure profile arrays contain at least ``length`` levels."""
        while len(self.sums) < length:
            self.sums.append(0.0)
            self.counts.append(0)
            self.maxima.append(None)
            self.minima.append(None)

    def update(self, sums: list[float], counts: list[int], maxima: list[float | None], minima: list[float | None]) -> None:
        """Update a profile from per-level partial statistics."""
        self.ensure(len(sums))
        for index, (value_sum, value_count, value_max, value_min) in enumerate(zip(sums, counts, maxima, minima)):
            if value_count <= 0:
                continue
            self.sums[index] += value_sum
            self.counts[index] += value_count
            if value_max is not None:
                self.maxima[index] = value_max if self.maxima[index] is None else max(self.maxima[index], value_max)
            if value_min is not None:
                self.minima[index] = value_min if self.minima[index] is None else min(self.minima[index], value_min)

    def means(self) -> list[float | None]:
        """Return the per-level mean profile."""
        return [value_sum / value_count if value_count > 0 else None for value_sum, value_count in zip(self.sums, self.counts)]


def _dataset_class():
    """Return netCDF4.Dataset or None."""
    try:
        from netCDF4 import Dataset  # type: ignore
    except Exception:
        return None
    return Dataset


def _numpy():
    """Import numpy lazily so dry-run and py_compile do not require it."""
    import numpy as np  # type: ignore

    return np


def _finite_array(values: Any, np: Any) -> Any:
    """Convert values to a flat finite numpy array."""
    array = np.asarray(values, dtype=float)
    if hasattr(array, "filled"):
        array = array.filled(np.nan)
    array = np.asarray(array, dtype=float).ravel()
    return array[np.isfinite(array)]


def _as_array(values: Any, np: Any) -> Any:
    """Convert NetCDF values to a float array with missing values as NaN."""
    array = np.asarray(values, dtype=float)
    if hasattr(array, "filled"):
        array = array.filled(np.nan)
    return np.asarray(array, dtype=float)


def _parse_list(text: str | None, default: Iterable[str]) -> list[str]:
    """Parse a comma-separated CLI list."""
    if text is None or not text.strip():
        return list(default)
    return [item.strip() for item in text.split(",") if item.strip()]


def _parse_float_list(text: str) -> list[float]:
    """Parse comma-separated profile times in minutes."""
    values: list[float] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            values.append(float(item))
        except ValueError:
            continue
    return values


def _history_files(case_dir: Path, max_files: int | None) -> tuple[list[Path], list[str]]:
    """Return sorted SCALE history files with optional foreground limit."""
    warnings: list[str] = []
    paths = sorted(path for path in case_dir.glob("history.pe*") if path.is_file())
    if max_files is not None and max_files >= 0 and len(paths) > max_files:
        warnings.append(f"History scan limited to {max_files} files out of {len(paths)}")
        paths = paths[:max_files]
    return paths, warnings


def _dimension_names(handle: Any) -> list[str]:
    """Return dimension names from a NetCDF handle."""
    return list(getattr(handle, "dimensions", {}).keys())


def _variable_names(handle: Any) -> list[str]:
    """Return variable names from a NetCDF handle."""
    return list(getattr(handle, "variables", {}).keys())


def _find_time_variable(handle: Any) -> str | None:
    """Find the history time-coordinate variable."""
    names = set(_variable_names(handle))
    for name in TIME_NAMES:
        if name in names:
            return name
    for name in names:
        lower = name.lower()
        if lower == "time" or lower.endswith("_time"):
            return name
    return None


def _infer_time_axis(variable: Any, time_length: int | None) -> int | None:
    """Infer the time axis for a history variable."""
    dims = tuple(getattr(variable, "dimensions", ()))
    for index, name in enumerate(dims):
        if "time" in name.lower():
            return index
    if time_length is None:
        return None
    shape = tuple(getattr(variable, "shape", ()))
    for index, length in enumerate(shape):
        if int(length) == int(time_length):
            return index
    return None


def _read_times(handle: Any) -> list[float]:
    """Read model/output times in seconds."""
    time_name = _find_time_variable(handle)
    if time_name is not None:
        try:
            np = _numpy()
            values = _finite_array(handle.variables[time_name][:], np)
            return [float(value) for value in values]
        except Exception:
            pass
    for dim_name, dim in getattr(handle, "dimensions", {}).items():
        if "time" in dim_name.lower():
            return [float(index) for index in range(len(dim))]
    return [0.0]


def _axis_after_removing_time(axis: int | None, time_axis: int | None) -> int | None:
    """Return the axis index after a time dimension has been sliced out."""
    if axis is None:
        return None
    if time_axis is None:
        return axis
    if axis == time_axis:
        return None
    return axis - 1 if axis > time_axis else axis


def _infer_vertical_axis(variable: Any, handle: Any, time_axis: int | None) -> int | None:
    """Infer the vertical axis for a history variable."""
    dims = tuple(getattr(variable, "dimensions", ()))
    for index, name in enumerate(dims):
        lower = name.lower()
        if lower in {candidate.lower() for candidate in VERTICAL_NAMES} or lower.startswith(("z", "lev")):
            if index != time_axis:
                return _axis_after_removing_time(index, time_axis)
    shape = tuple(getattr(variable, "shape", ()))
    for index, name in enumerate(dims):
        if index == time_axis:
            continue
        lower = name.lower()
        if any(marker in lower for marker in HORIZONTAL_MARKERS):
            continue
        try:
            if int(shape[index]) > 1:
                return _axis_after_removing_time(index, time_axis)
        except Exception:
            continue
    return None


def _read_vertical_coordinate(handle: Any, variable: Any, vertical_axis: int | None, time_axis: int | None) -> list[float] | None:
    """Read or synthesize a vertical coordinate in meters."""
    if vertical_axis is None:
        return None
    dims = tuple(getattr(variable, "dimensions", ()))
    original_axis = vertical_axis if time_axis is None or vertical_axis < time_axis else vertical_axis + 1
    dim_name = dims[original_axis] if original_axis < len(dims) else None
    names = _variable_names(handle)
    coord_candidates = []
    if dim_name:
        coord_candidates.append(dim_name)
    coord_candidates.extend(VERTICAL_NAMES)
    for name in coord_candidates:
        if name not in names:
            continue
        try:
            np = _numpy()
            values = _finite_array(handle.variables[name][:], np)
            if values.size > 1:
                return [float(value) for value in values]
        except Exception:
            continue
    try:
        length = int(getattr(variable, "shape", ())[original_axis])
    except Exception:
        return None
    return [float(index) for index in range(length)]


def _read_time_slice(variable: Any, time_axis: int | None, time_index: int) -> Any:
    """Read one time slice from a NetCDF variable."""
    if time_axis is None:
        return variable[:]
    slices = [slice(None)] * len(getattr(variable, "shape", ()))
    slices[time_axis] = time_index
    return variable[tuple(slices)]


def _profile_partial_stats(values: Any, vertical_axis: int | None) -> tuple[list[float], list[int], list[float | None], list[float | None]]:
    """Return per-level sums, counts, maxima, and minima for one field chunk."""
    np = _numpy()
    array = _as_array(values, np)
    if vertical_axis is None or array.ndim == 0:
        finite = _finite_array(array, np)
        if finite.size == 0:
            return [], [], [], []
        return [float(finite.sum())], [int(finite.size)], [float(finite.max())], [float(finite.min())]
    moved = np.moveaxis(array, vertical_axis, 0)
    reshaped = moved.reshape((moved.shape[0], -1))
    sums: list[float] = []
    counts: list[int] = []
    maxima: list[float | None] = []
    minima: list[float | None] = []
    for level_values in reshaped:
        finite = level_values[np.isfinite(level_values)]
        if finite.size == 0:
            sums.append(0.0)
            counts.append(0)
            maxima.append(None)
            minima.append(None)
        else:
            sums.append(float(finite.sum()))
            counts.append(int(finite.size))
            maxima.append(float(finite.max()))
            minima.append(float(finite.min()))
    return sums, counts, maxima, minima


def _display_value(variable: str, value: float | None) -> float | None:
    """Convert a raw variable value to plotting/reporting units."""
    if value is None:
        return None
    if variable.startswith(MIXING_RATIO_PREFIXES):
        return value * 1000.0
    return value


def _display_unit(variable: str, raw_unit: str | None = None) -> str:
    """Return a compact display unit for a history variable."""
    if variable.startswith(MIXING_RATIO_PREFIXES):
        return "g kg-1"
    if variable == "W":
        return "m s-1"
    if variable == "T":
        return "K"
    if variable.startswith("TKE"):
        return "m2 s-2"
    return raw_unit or "-"


def _format_time_label(time_min: float) -> str:
    """Format profile time labels without rounding early output to 0 min."""
    if abs(time_min - round(time_min)) < 1.0e-6:
        return f"{time_min:.0f} min"
    if time_min < 1.0:
        return f"{time_min:.2f} min"
    return f"{time_min:.1f} min"


def _column_name(variable: str, statistic: str) -> str:
    """Return a stable time-series column name for a variable and statistic."""
    mapping = {
        ("QHYD_sd", "mean"): "qhyd_mean_gkg",
        ("QHYD_sd", "max"): "qhyd_max_gkg",
        ("QHYD", "mean"): "qhyd_mean_gkg",
        ("QHYD", "max"): "qhyd_max_gkg",
        ("QC_sd", "mean"): "qc_mean_gkg",
        ("QC_sd", "max"): "qc_max_gkg",
        ("QC", "mean"): "qc_mean_gkg",
        ("QC", "max"): "qc_max_gkg",
        ("QR_sd", "mean"): "qr_mean_gkg",
        ("QR_sd", "max"): "qr_max_gkg",
        ("QR", "mean"): "qr_mean_gkg",
        ("QR", "max"): "qr_max_gkg",
        ("W", "mean"): "w_mean_ms",
        ("W", "max"): "w_max_ms",
        ("W", "min"): "w_min_ms",
        ("RH", "mean"): "rh_mean",
        ("RH", "max"): "rh_max",
        ("T", "mean"): "t_mean_k",
        ("QV", "mean"): "qv_mean_gkg",
        ("LWPT", "mean"): "lwpt_mean",
        ("LWPT", "max"): "lwpt_max",
        ("RAIN_ACC_sd", "mean"): "rain_acc_mean",
        ("RAIN_ACC_sd", "max"): "rain_acc_max",
        ("RAIN", "mean"): "rain_mean",
        ("RAIN", "max"): "rain_max",
        ("TKE_RS", "max"): "tke_rs_max",
        ("TKE_SMG", "max"): "tke_smg_max",
    }
    return mapping.get((variable, statistic), f"{variable.lower()}_{statistic}")


def _select_nearest_times(available_times_s: list[float], request_times_min: list[float]) -> dict[float, float]:
    """Map requested profile times in minutes to nearest available model times."""
    if not available_times_s:
        return {}
    output: dict[float, float] = {}
    for request_min in request_times_min:
        request_s = request_min * 60.0
        nearest = min(available_times_s, key=lambda value: abs(value - request_s))
        output[request_min] = nearest
    return output


def _inspect_first_file(paths: list[Path]) -> tuple[list[float], list[str], list[str]]:
    """Inspect times and variable names from the first history file."""
    Dataset = _dataset_class()
    if Dataset is None:
        return [], [], ["NetCDF library unavailable; history fields were not read"]
    if not paths:
        return [], [], ["No history.pe* files found"]
    warnings: list[str] = []
    try:
        with Dataset(paths[0], "r") as handle:
            return _read_times(handle), _variable_names(handle), warnings
    except Exception as exc:
        return [], [], [f"Failed to inspect {paths[0]}: {exc}"]


def _merge_file_times(global_times: list[float], file_times: list[float]) -> list[float]:
    """Merge discovered times while preserving sorted numeric ordering."""
    merged = {round(value, 9): float(value) for value in global_times}
    for value in file_times:
        merged.setdefault(round(value, 9), float(value))
    return [merged[key] for key in sorted(merged)]


def analyze_history(
    case_dir: Path,
    outdir: Path,
    strict: bool,
    options: dict[str, Any],
    timeseries_variables: list[str],
    profile_variables: list[str],
    profile_time_requests_min: list[float],
    cloud_threshold_kgkg: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Read history files and return time-series, profile, and time-height tables."""
    warnings: list[str] = []
    Dataset = _dataset_class()
    max_files = options.get("max_files")
    history_paths, file_warnings = _history_files(case_dir, max_files)
    warnings.extend(file_warnings)
    if Dataset is None:
        warn_or_raise(warnings, "NetCDF library unavailable; history diagnostics skipped", strict=strict)
        metadata = {
            "case_dir": case_dir,
            "history_file_count": 0,
            "processed_history_file_count": 0,
            "time_count": 0,
            "selected_profile_times_min": "NA",
            "available_variables": "NA",
            "timeseries_variables_used": "NA",
            "profile_variables_used": "NA",
            "cloud_variable": "NA",
            "cloud_threshold_kgkg": cloud_threshold_kgkg,
            "netcdf_backend": "unavailable",
            "warnings": warning_text(warnings),
        }
        return [], [], [], metadata
    if not history_paths:
        warn_or_raise(warnings, f"No history.pe* files found in {case_dir}", strict=strict)
        metadata = {
            "case_dir": case_dir,
            "history_file_count": 0,
            "processed_history_file_count": 0,
            "time_count": 0,
            "selected_profile_times_min": "NA",
            "available_variables": "NA",
            "timeseries_variables_used": "NA",
            "profile_variables_used": "NA",
            "cloud_variable": "NA",
            "cloud_threshold_kgkg": cloud_threshold_kgkg,
            "netcdf_backend": "netCDF4",
            "warnings": warning_text(warnings),
        }
        return [], [], [], metadata

    first_times, first_variables, inspect_warnings = _inspect_first_file(history_paths)
    warnings.extend(inspect_warnings)
    available_set = set(first_variables)
    timeseries_used = [name for name in timeseries_variables if name in available_set]
    profile_used = [name for name in profile_variables if name in available_set]
    cloud_variable = next((name for name in CLOUD_VARIABLE_CANDIDATES if name in available_set), None)
    if cloud_variable is None:
        warnings.append("No QHYD/QC cloud variable found; cloud fraction/base/top are NA")
    if not timeseries_used:
        warnings.append("No requested time-series variables were found in history files")
    if not profile_used:
        warnings.append("No requested profile variables were found in history files")

    selected_profile_times = _select_nearest_times(first_times, profile_time_requests_min)
    selected_time_to_request = {round(selected, 9): request for request, selected in selected_profile_times.items()}

    scalar_acc: dict[tuple[float, str], ScalarAccumulator] = {}
    cloud_counts: dict[float, int] = {}
    cloud_totals: dict[float, int] = {}
    cloud_profile: dict[float, ProfileAccumulator] = {}
    profile_acc: dict[tuple[float, str], ProfileAccumulator] = {}
    profile_z: dict[tuple[float, str], list[float] | None] = {}
    qhyd_time_height: dict[float, ProfileAccumulator] = {}
    qhyd_z: list[float] | None = None
    time_file_counts: dict[float, int] = {}
    all_times: list[float] = []

    for path in history_paths:
        try:
            with Dataset(path, "r") as handle:
                file_times = _read_times(handle)
                all_times = _merge_file_times(all_times, file_times)
                variable_names = set(_variable_names(handle))
                file_timeseries = [name for name in timeseries_used if name in variable_names]
                file_profiles = [name for name in profile_used if name in variable_names]
                for time_index, time_s in enumerate(file_times):
                    time_file_counts[time_s] = time_file_counts.get(time_s, 0) + 1
                    for name in file_timeseries:
                        variable = handle.variables[name]
                        time_axis = _infer_time_axis(variable, len(file_times))
                        try:
                            values = _read_time_slice(variable, time_axis, time_index)
                        except Exception as exc:
                            warnings.append(f"Failed reading {name} at t={time_s} in {path.name}: {exc}")
                            continue
                        scalar_acc.setdefault((time_s, name), ScalarAccumulator()).update(values)
                        if name == cloud_variable:
                            np = _numpy()
                            finite = _finite_array(values, np)
                            if finite.size > 0:
                                cloud_counts[time_s] = cloud_counts.get(time_s, 0) + int((finite >= cloud_threshold_kgkg).sum())
                                cloud_totals[time_s] = cloud_totals.get(time_s, 0) + int(finite.size)
                            vertical_axis = _infer_vertical_axis(variable, handle, time_axis)
                            z_values = _read_vertical_coordinate(handle, variable, vertical_axis, time_axis)
                            sums, counts, maxima, minima = _profile_partial_stats(values, vertical_axis)
                            cloud_profile.setdefault(time_s, ProfileAccumulator()).update(sums, counts, maxima, minima)
                            if qhyd_z is None and z_values is not None:
                                qhyd_z = z_values
                            qhyd_time_height.setdefault(time_s, ProfileAccumulator()).update(sums, counts, maxima, minima)
                    if round(time_s, 9) in selected_time_to_request:
                        for name in file_profiles:
                            variable = handle.variables[name]
                            time_axis = _infer_time_axis(variable, len(file_times))
                            try:
                                values = _read_time_slice(variable, time_axis, time_index)
                            except Exception as exc:
                                warnings.append(f"Failed reading profile {name} at t={time_s} in {path.name}: {exc}")
                                continue
                            vertical_axis = _infer_vertical_axis(variable, handle, time_axis)
                            sums, counts, maxima, minima = _profile_partial_stats(values, vertical_axis)
                            key = (time_s, name)
                            profile_acc.setdefault(key, ProfileAccumulator()).update(sums, counts, maxima, minima)
                            if key not in profile_z:
                                profile_z[key] = _read_vertical_coordinate(handle, variable, vertical_axis, time_axis)
        except Exception as exc:
            warn_or_raise(warnings, f"Failed reading history file {path}: {exc}", strict=strict)

    timeseries_rows: list[dict[str, Any]] = []
    for time_s in sorted(all_times):
        row: dict[str, Any] = {
            "time_s": time_s,
            "time_min": time_s / 60.0,
            "history_file_count": time_file_counts.get(time_s),
            "cloud_variable": cloud_variable,
            "cloud_threshold_kgkg": cloud_threshold_kgkg,
            "warnings": None,
        }
        if cloud_totals.get(time_s, 0) > 0:
            row["cloud_fraction"] = cloud_counts.get(time_s, 0) / cloud_totals[time_s]
        cloud_prof = cloud_profile.get(time_s)
        if cloud_prof is not None and qhyd_z is not None:
            cloudy_levels = [
                index
                for index, value in enumerate(cloud_prof.maxima)
                if value is not None and value >= cloud_threshold_kgkg and index < len(qhyd_z)
            ]
            if cloudy_levels:
                base = qhyd_z[min(cloudy_levels)]
                top = qhyd_z[max(cloudy_levels)]
                row["cloud_base_m"] = base
                row["cloud_top_m"] = top
                row["cloud_depth_m"] = top - base
        for name in timeseries_used:
            acc = scalar_acc.get((time_s, name))
            if acc is None:
                continue
            if name in {"W", "RH", "LWPT", "RAIN_ACC_sd", "RAIN"} or name.startswith(("Q", "TKE")) or name == "T":
                mean_column = _column_name(name, "mean")
                max_column = _column_name(name, "max")
                min_column = _column_name(name, "min")
                if mean_column in TIME_SERIES_COLUMNS:
                    row[mean_column] = _display_value(name, acc.mean)
                if max_column in TIME_SERIES_COLUMNS:
                    row[max_column] = _display_value(name, acc.maximum)
                if min_column in TIME_SERIES_COLUMNS:
                    row[min_column] = _display_value(name, acc.minimum)
        timeseries_rows.append(row)

    profile_rows: list[dict[str, Any]] = []
    for (time_s, name), acc in sorted(profile_acc.items(), key=lambda item: (item[0][0], item[0][1])):
        z_values = profile_z.get((time_s, name))
        means = acc.means()
        request_min = selected_time_to_request.get(round(time_s, 9))
        for index, mean_value in enumerate(means):
            z_value = z_values[index] if z_values is not None and index < len(z_values) else index
            row = {
                "time_s": time_s,
                "time_min": time_s / 60.0,
                "profile_time_request_min": request_min,
                "variable": name,
                "z_index": index,
                "z_m": z_value,
                "domain_mean": mean_value,
                "domain_max": acc.maxima[index] if index < len(acc.maxima) else None,
                "domain_min": acc.minima[index] if index < len(acc.minima) else None,
                "domain_count": acc.counts[index] if index < len(acc.counts) else None,
                "display_value": _display_value(name, mean_value),
                "display_unit": _display_unit(name),
                "warnings": None,
            }
            profile_rows.append(row)

    time_height_rows: list[dict[str, Any]] = []
    for time_s, acc in sorted(qhyd_time_height.items()):
        means = acc.means()
        for index, mean_value in enumerate(means):
            z_value = qhyd_z[index] if qhyd_z is not None and index < len(qhyd_z) else index
            qhyd_max = acc.maxima[index] if index < len(acc.maxima) else None
            time_height_rows.append(
                {
                    "time_s": time_s,
                    "time_min": time_s / 60.0,
                    "z_index": index,
                    "z_m": z_value,
                    "qhyd_mean_gkg": _display_value(cloud_variable or "QHYD", mean_value),
                    "qhyd_max_gkg": _display_value(cloud_variable or "QHYD", qhyd_max),
                    "grid_count": acc.counts[index] if index < len(acc.counts) else None,
                    "warnings": None,
                }
            )

    metadata = {
        "case_dir": case_dir,
        "history_file_count": len(sorted(case_dir.glob("history.pe*"))),
        "processed_history_file_count": len(history_paths),
        "time_count": len(all_times),
        "selected_profile_times_min": ";".join(f"{value / 60.0:g}" for value in sorted(set(selected_profile_times.values()))),
        "available_variables": ";".join(first_variables),
        "timeseries_variables_used": ";".join(timeseries_used) if timeseries_used else None,
        "profile_variables_used": ";".join(profile_used) if profile_used else None,
        "cloud_variable": cloud_variable,
        "cloud_threshold_kgkg": cloud_threshold_kgkg,
        "netcdf_backend": "netCDF4",
        "warnings": warning_text(warnings),
    }
    return timeseries_rows, profile_rows, time_height_rows, metadata


def _series(rows: list[dict[str, str]], key: str) -> list[float | None]:
    """Read a numeric series from CSV-like rows."""
    return [safe_float(row.get(key)) for row in rows]


def _times(rows: list[dict[str, str]]) -> list[float | None]:
    """Return time in minutes from CSV-like rows."""
    output = []
    for row in rows:
        value = safe_float(row.get("time_min"))
        if value is None:
            seconds = safe_float(row.get("time_s"))
            value = seconds / 60.0 if seconds is not None else None
        output.append(value)
    return output


def _plot_timeseries(outdir: Path) -> None:
    """Plot BW history time-series sanity checks."""
    rows = read_table(outdir / "tables" / "02_tpht_history_sanity_timeseries.csv")
    plt = configure_matplotlib()
    fig, axes = plt.subplots(3, 2, figsize=figure_size("double", 0.92), sharex=True, constrained_layout=True)
    xs = _times(rows)
    panels = [
        (("cloud_base_m", "cloud_top_m"), "Height (m)", "Cloud layer", (PALETTE[0], PALETTE[1])),
        (("qhyd_max_gkg", "qc_max_gkg", "qr_max_gkg"), "Mixing ratio (g kg-1)", "Hydrometeor maxima", (PALETTE[0], PALETTE[2], PALETTE[3])),
        (("cloud_fraction",), "Cloud fraction (%)", "Cloud fraction", (PALETTE[4],)),
        (("w_max_ms", "w_min_ms"), "Vertical velocity (m s-1)", "Vertical velocity extrema", (PALETTE[2], PALETTE[1])),
        (("rain_acc_max", "rain_max"), "Rain diagnostic", "Rain output", (PALETTE[5], PALETTE[6])),
        (("lwpt_mean", "tke_rs_max", "tke_smg_max"), "History diagnostic", "LWPT/TKE", (PALETTE[0], PALETTE[3], PALETTE[4])),
    ]
    for ax, (keys, ylabel, title, colors) in zip(axes.flat, panels):
        plotted = False
        for key, color in zip(keys, colors):
            ys = _series(rows, key)
            if key == "cloud_fraction":
                ys = [value * 100.0 if value is not None else None for value in ys]
            points = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
            if not points:
                continue
            px, py = zip(*points)
            ax.plot(px, py, label=key, color=color, linewidth=1.0)
            plotted = True
        if plotted:
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.legend(loc="best", frameon=False)
        else:
            no_data_panel(ax, title)
    for ax in axes[-1]:
        ax.set_xlabel("Model time (min)")
    save_figure(fig, outdir / "figures", "02_tpht_history_sanity_timeseries")
    plt.close(fig)


def _plot_profiles(outdir: Path) -> None:
    """Plot selected-time vertical profiles from BW history files."""
    rows = read_table(outdir / "tables" / "02_tpht_history_sanity_profiles.csv")
    preferred = ["QHYD_sd", "QC_sd", "QR_sd", "W", "RH", "TKE_RS"]
    variables = [name for name in preferred if any(row.get("variable") == name for row in rows)]
    variables = variables[:6]
    plt = configure_matplotlib()
    fig, axes = plt.subplots(2, 3, figsize=figure_size("double", 0.78), constrained_layout=True)
    if not variables:
        for ax in axes.flat:
            no_data_panel(ax, "BW history profile")
    for ax, variable in zip(axes.flat, variables):
        var_rows = [row for row in rows if row.get("variable") == variable]
        times = sorted({safe_float(row.get("time_min")) for row in var_rows if safe_float(row.get("time_min")) is not None})
        for index, time_min in enumerate(times):
            time_rows = [row for row in var_rows if safe_float(row.get("time_min")) == time_min]
            values = _series(time_rows, "display_value")
            heights = _series(time_rows, "z_m")
            points = [(value, z) for value, z in zip(values, heights) if value is not None and z is not None]
            if not points:
                continue
            px, py = zip(*points)
            ax.plot(px, py, color=PALETTE[index % len(PALETTE)], linewidth=1.0, label=_format_time_label(time_min))
        unit = next((row.get("display_unit") for row in var_rows if row.get("display_unit")), "-")
        ax.set_xlabel(f"{variable} ({unit})")
        ax.set_ylabel("Height (m)")
        ax.legend(loc="best", frameon=False)
    for ax in list(axes.flat)[len(variables) :]:
        no_data_panel(ax, "BW history profile")
    save_figure(fig, outdir / "figures", "02_tpht_history_sanity_profiles")
    plt.close(fig)


def _plot_time_height(outdir: Path) -> None:
    """Plot height-time hydrometeor envelope from BW history files."""
    rows = read_table(outdir / "tables" / "02_tpht_history_sanity_time_height.csv")
    plt = configure_matplotlib()
    fig, ax = plt.subplots(figsize=figure_size("double", 0.55), constrained_layout=True)
    if not rows:
        no_data_panel(ax, "QHYD height-time envelope")
    else:
        np = _numpy()
        times = sorted({safe_float(row.get("time_min")) for row in rows if safe_float(row.get("time_min")) is not None})
        heights = sorted({safe_float(row.get("z_m")) for row in rows if safe_float(row.get("z_m")) is not None})
        time_index = {value: index for index, value in enumerate(times)}
        height_index = {value: index for index, value in enumerate(heights)}
        grid = np.full((len(heights), len(times)), np.nan)
        for row in rows:
            time_min = safe_float(row.get("time_min"))
            height = safe_float(row.get("z_m"))
            value = safe_float(row.get("qhyd_max_gkg"))
            if time_min is None or height is None or value is None:
                continue
            grid[height_index[height], time_index[time_min]] = value
        if times and heights and np.isfinite(grid).any():
            positive = grid[np.isfinite(grid) & (grid > 0.0)]
            if positive.size > 0:
                from matplotlib.colors import LogNorm

                floor = max(float(positive.min()), 1.0e-9)
                mesh = ax.pcolormesh(times, heights, np.maximum(grid, floor), shading="auto", norm=LogNorm(vmin=floor, vmax=float(positive.max())))
                cbar = fig.colorbar(mesh, ax=ax)
                cbar.set_label("Max QHYD/QC by height (g kg-1)")
            else:
                mesh = ax.pcolormesh(times, heights, grid, shading="auto")
                cbar = fig.colorbar(mesh, ax=ax)
                cbar.set_label("Max QHYD/QC by height (g kg-1)")
            ax.set_xlabel("Model time (min)")
            ax.set_ylabel("Height (m)")
        else:
            no_data_panel(ax, "QHYD height-time envelope")
    save_figure(fig, outdir / "figures", "02_tpht_history_sanity_time_height")
    plt.close(fig)


def write_outputs(
    outdir: Path,
    timeseries_rows: list[dict[str, Any]],
    profile_rows: list[dict[str, Any]],
    time_height_rows: list[dict[str, Any]],
    metadata: dict[str, Any],
) -> None:
    """Write tables and figures for BW history sanity diagnostics."""
    write_table_bundle(timeseries_rows, outdir / "tables" / "02_tpht_history_sanity_timeseries", TIME_SERIES_COLUMNS)
    write_table_bundle(profile_rows, outdir / "tables" / "02_tpht_history_sanity_profiles", PROFILE_COLUMNS)
    write_table_bundle(time_height_rows, outdir / "tables" / "02_tpht_history_sanity_time_height", TIME_HEIGHT_COLUMNS)
    write_table_bundle(metadata, outdir / "tables" / "02_tpht_history_sanity_metadata", METADATA_COLUMNS)
    try:
        _plot_timeseries(outdir)
        _plot_profiles(outdir)
        _plot_time_height(outdir)
    except ModuleNotFoundError as exc:
        # Some login-node Python environments can read tables but lack
        # Matplotlib.  Keep the diagnostic tables rather than failing the scan.
        (outdir / "logs").mkdir(parents=True, exist_ok=True)
        (outdir / "logs" / "02_tpht_history_sanity_plot_warning.txt").write_text(
            f"Skipped history sanity figures because a plotting dependency is unavailable: {exc}\n"
        )


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    parser.add_argument("--case-dir", type=Path, default=None, help="Direct path to the BW case directory containing history.pe* files.")
    parser.add_argument("--group", default=DEFAULT_GROUP, help="GMD2026 group name when --case-dir is not provided.")
    parser.add_argument("--case-name", default=DEFAULT_CASE, help="Case name when --case-dir is not provided.")
    parser.add_argument("--variables", default=None, help="Comma-separated history variables for time-series diagnostics.")
    parser.add_argument("--profile-variables", default=None, help="Comma-separated history variables for vertical profiles.")
    parser.add_argument("--profile-times-min", default=DEFAULT_PROFILE_TIMES_MIN, help="Comma-separated requested profile times in minutes.")
    parser.add_argument("--cloud-threshold-kgkg", type=float, default=1.0e-7, help="Hydrometeor threshold for cloud fraction/base/top.")
    args = parser.parse_args()

    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    case_dir = args.case_dir.resolve() if args.case_dir is not None else root / args.group / args.case_name
    options = analysis_options(args)
    timeseries_variables = _parse_list(args.variables, DEFAULT_TIMESERIES_VARIABLES)
    profile_variables = _parse_list(args.profile_variables, DEFAULT_PROFILE_VARIABLES)
    profile_time_requests_min = _parse_float_list(args.profile_times_min)

    if args.dry_run:
        dry_run_message(
            Path(__file__).name,
            root,
            outdir,
            [
                "tables/02_tpht_history_sanity_timeseries.{csv,md,tex,json}",
                "tables/02_tpht_history_sanity_profiles.{csv,md,tex,json}",
                "tables/02_tpht_history_sanity_time_height.{csv,md,tex,json}",
                "tables/02_tpht_history_sanity_metadata.{csv,md,tex,json}",
                "figures/02_tpht_history_sanity_timeseries.{pdf,svg,png}",
                "figures/02_tpht_history_sanity_profiles.{pdf,svg,png}",
                "figures/02_tpht_history_sanity_time_height.{pdf,svg,png}",
            ],
        )
        print(f"[dry-run] case_dir={case_dir}")
        return

    ensure_output_dirs(outdir)
    timeseries_rows, profile_rows, time_height_rows, metadata = analyze_history(
        case_dir=case_dir,
        outdir=outdir,
        strict=bool(args.strict),
        options=options,
        timeseries_variables=timeseries_variables,
        profile_variables=profile_variables,
        profile_time_requests_min=profile_time_requests_min,
        cloud_threshold_kgkg=float(args.cloud_threshold_kgkg),
    )
    write_outputs(outdir, timeseries_rows, profile_rows, time_height_rows, metadata)


if __name__ == "__main__":
    main()
