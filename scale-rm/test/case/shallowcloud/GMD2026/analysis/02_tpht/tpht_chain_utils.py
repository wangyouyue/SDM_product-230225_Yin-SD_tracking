"""Stepwise BW selected-output chain reconstruction helpers for GMD2026.

The shared SCALE-SDM implementation writes the runtime ``sd_id/dm_id`` arrays
as BW ``pre_dmid/pre_sdid`` variables.  Normal selected output keeps those IDs
stable; history-style output may expose adjacent-output predecessor semantics.
These helpers therefore advance one output level at a time, which correctly
handles stable IDs and still detects broken adjacent links when IDs change.
"""

from __future__ import annotations

import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_netcdf import (  # noqa: E402
    BW_ID_NAMES,
    HEIGHT_NAMES,
    IF_COAL_NAMES,
    RADIUS_NAMES,
    TIME_NAMES,
    _dataset_class,
    _filename_model_time_s,
    _find_pair_variables,
    _find_variable,
    _flatten,
    _flatten_int,
    _iter_first_dim_chunks,
    _limited_paths,
    _variable_length,
)
from common.table_utils import safe_ratio  # noqa: E402


X_NAMES = ("x", "sd_x", "x_sd", "pos_x", "xp")
Y_NAMES = ("y", "sd_y", "y_sd", "pos_y", "yp")
N_NAMES = ("sd_n", "n", "multiplicity")
LARGE_RADIUS_M = 15.0e-6

COAL_EVENT_TIME_NAMES = ("event_time", "time", "time_s")
COAL_PRE1_NAMES = (("pre_dmid1", "pre_sdid1"), ("pre_dm_id1", "pre_sd_id1"))
COAL_PRE2_NAMES = (("pre_dmid2", "pre_sdid2"), ("pre_dm_id2", "pre_sd_id2"))
COAL_RADIUS1_NAMES = ("sd_r1", "r1", "radius1")
COAL_RADIUS2_NAMES = ("sd_r2", "r2", "radius2")
COAL_N1_NAMES = ("sd_n1", "n1", "multiplicity1")
COAL_N2_NAMES = ("sd_n2", "n2", "multiplicity2")
COAL_NUM_COL_NAMES = ("num_col",)

LINK_COLUMNS = [
    "time_s",
    "record_count",
    "matched_parent_records",
    "unmatched_parent_records",
    "valid_link_fraction",
    "unique_output_ids",
    "duplicate_output_ids",
    "active_chain_count",
    "coalescence_event_count_linked",
    "warnings",
]

TREE_COLUMNS = [
    "example_id",
    "target_id",
    "node_id",
    "predecessor_id",
    "time_s",
    "height_m",
    "radius_um",
    "if_coal",
    "node_role",
    "event_time_s",
    "num_col",
    "source_file",
    "notes",
    "warnings",
]

TRAJECTORY_COLUMNS = [
    "target_id",
    "time_s",
    "x_m",
    "y_m",
    "z_m",
    "radius_um",
    "if_coal",
    "category",
    "max_radius_um",
    "first_large_time_s",
    "first_ifcoal_time_s",
    "warnings",
]

IFCOAL_TIMELINE_COLUMNS = [
    "target_id",
    "time_s",
    "if_coal",
    "category",
    "max_radius_um",
    "first_large_time_s",
    "first_ifcoal_time_s",
    "warnings",
]

INTERVAL_COLUMNS = [
    "target_id",
    "time_s",
    "dt_s",
    "r0_um",
    "r1_um",
    "height_m",
    "category",
    "if_coal_interval",
    "sampled_interval",
    "warnings",
]

EVENT_LINK_COLUMNS = [
    "target_id",
    "event_time_s",
    "event_height_m",
    "pre_radius_um",
    "post_radius_um",
    "delta_r3_um3",
    "participant",
    "num_col",
    "category",
    "source_file",
    "warnings",
]


def _finite_or_none(value: float | None) -> float | None:
    """Return a finite value or ``None``."""
    if value is None:
        return None
    return value if math.isfinite(value) else None


def _value_at(values: list[float], index: int) -> float | None:
    """Return an indexed finite value or ``None``."""
    if index >= len(values):
        return None
    return _finite_or_none(values[index])


def _get_optional_values(handle: Any, name: str | None, selection: slice) -> list[float]:
    """Read an optional variable selection as finite floats."""
    if name is None:
        return []
    return _flatten(handle.variables[name][selection])


def _find_pair(variable_names: Iterable[str], preferred: Iterable[tuple[str, str]]) -> tuple[str, str] | None:
    """Find a pair of tracking-ID variables."""
    names = set(variable_names)
    for dm_name, sd_name in preferred:
        if dm_name in names and sd_name in names:
            return dm_name, sd_name
    return None


def _group_paths_by_time(paths: Iterable[Path], max_files: int | None) -> tuple[list[tuple[float | None, list[Path]]], list[str]]:
    """Group per-rank files by output model time."""
    candidate_paths = sorted(
        [path for path in paths if path.is_file() and path.suffix != ".ids"],
        key=lambda path: (_filename_model_time_s(path) is None, _filename_model_time_s(path) or 0.0, path.name),
    )
    candidate_paths, warnings = _limited_paths(candidate_paths, max_files)
    grouped: dict[float | None, list[Path]] = defaultdict(list)
    for path in candidate_paths:
        grouped[_filename_model_time_s(path)].append(path)
    groups = [(time_s, sorted(group)) for time_s, group in grouped.items()]
    groups.sort(key=lambda item: (item[0] is None, item[0] if item[0] is not None else math.inf))
    return groups, warnings


def _new_chain_state(chain_id: str, dm_id: int, sd_id: int) -> dict[str, Any]:
    """Create a lifecycle accumulator for one propagated BW chain."""
    return {
        "chain_id": chain_id,
        "dm_id": dm_id,
        "sd_id": sd_id,
        "record_count": 0,
        "first_time_s": None,
        "last_time_s": None,
        "last_order": -1,
        "final_radius_m": None,
        "max_radius_m": None,
        "first_large_time_s": None,
        "first_large_order": None,
        "first_large_radius_m": None,
        "first_large_height_m": None,
        "first_large_x_m": None,
        "first_large_y_m": None,
        "first_selected_time_s": None,
        "first_selected_order": None,
        "first_selected_category": None,
        "first_coal_time_s": None,
        "first_coal_order": None,
        "first_coal_event_time_s": None,
        "first_coal_event_order": None,
        "duration_ge_15_s": 0.0,
        "duration_ge_20_s": 0.0,
        "previous_time_s": None,
        "previous_radius_m": None,
        "coal_record_count": 0,
        "coal_episode_count_proxy": 0,
        "coal_event_count": 0,
        "coal_num_col_sum": 0.0,
        "previous_coal_active": False,
        "has_large": False,
        "has_coal": False,
        "has_coal_event": False,
    }


def _chain_number(chain_id: str) -> int | None:
    """Return the numeric part of a generated ``chain_########`` identifier."""
    match = re.search(r"(\d+)$", chain_id)
    if not match:
        return None
    return int(match.group(1))


def _sampled_chain(chain_id: str, stride: int) -> bool:
    """Select a deterministic subset of chains for lightweight source tables."""
    number = _chain_number(chain_id)
    return number is not None and number % max(1, stride) == 0


def _state_category(state: dict[str, Any] | None) -> str:
    """Return the final first-interest category for one target chain."""
    if not state:
        return "unknown"
    return state.get("first_selected_category") or "unknown"


def _state_um(state: dict[str, Any] | None, key: str) -> float | None:
    """Return a state radius in micrometers."""
    if not state or state.get(key) is None:
        return None
    return float(state[key]) * 1.0e6


def _first_selected_category(radius_m: float | None, coal_value: float | None) -> str | None:
    """Classify the first record that satisfies a TPHT interest condition."""
    large = radius_m is not None and radius_m >= LARGE_RADIUS_M
    coalesced = coal_value is not None and coal_value > 0.0
    if large and coalesced:
        return "both"
    if large:
        return "radius_only"
    if coalesced:
        return "coal_only"
    return None


def _is_earlier(time_s: float | None, order: int, ref_time: float | None, ref_order: int | None) -> bool:
    """Compare records by model time when available, otherwise by stream order."""
    if ref_order is None:
        return True
    if time_s is not None and ref_time is not None:
        return time_s < ref_time
    return order < ref_order


def _height_bin_center(height_m: float | None) -> float | None:
    """Return a coarse height-bin center used for height-time occurrence tables."""
    if height_m is None:
        return None
    bins = (0.0, 250.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, math.inf)
    for lower, upper in zip(bins[:-1], bins[1:]):
        if lower <= height_m < upper:
            return (lower + upper) * 0.5 if math.isfinite(upper) else lower + 250.0
    return None


def _update_time_accumulator(accumulator: dict[float, dict[str, float]], time_s: float | None, radius_m: float | None, height_m: float | None, coal_value: float | None) -> None:
    """Accumulate selected-output statistics by output time."""
    if time_s is None:
        return
    row = accumulator.setdefault(time_s, {"count": 0.0, "radius_sum": 0.0, "height_sum": 0.0, "height_count": 0.0, "large_count": 0.0, "coal_count": 0.0})
    row["count"] += 1.0
    if radius_m is not None:
        row["radius_sum"] += radius_m
        if radius_m >= LARGE_RADIUS_M:
            row["large_count"] += 1.0
    if height_m is not None:
        row["height_sum"] += height_m
        row["height_count"] += 1.0
    if coal_value is not None and coal_value > 0.0:
        row["coal_count"] += 1.0


def _update_occurrence_accumulator(
    accumulator: dict[tuple[float, float], dict[str, float]],
    time_s: float | None,
    height_m: float | None,
    is_large_target: bool,
    is_ifcoal_target: bool,
) -> None:
    """Accumulate full-history target counts in height-time bins.

    The category flags are chain-level properties computed after following the
    BW predecessor links.  This avoids plotting only the short intervals where a
    target is instantaneously above the radius threshold or instantaneously
    flagged by ``if_coal``.
    """
    height_center = _height_bin_center(height_m)
    if time_s is None or height_center is None:
        return
    row = accumulator.setdefault((time_s, height_center), {"target_count": 0.0, "large_target_count": 0.0, "ifcoal_target_count": 0.0})
    row["target_count"] += 1.0
    if is_large_target:
        row["large_target_count"] += 1.0
    if is_ifcoal_target:
        row["ifcoal_target_count"] += 1.0


def _build_category_occurrence_accumulator(
    selected_groups: list[tuple[float | None, list[Path]]],
    targets: dict[str, dict[str, Any]],
    options: dict[str, Any],
) -> tuple[dict[tuple[float, float], dict[str, float]], list[str]]:
    """Re-scan selected outputs and bin full target histories by final category."""
    by_time_height: dict[tuple[float, float], dict[str, float]] = {}
    warnings: list[str] = []
    previous_id_to_chain: dict[tuple[int, int], str] = {}
    next_chain_index = 0
    for level_index, (time_s, paths) in enumerate(selected_groups):
        records, level_warnings = _read_selected_level(paths, time_s, options)
        warnings.extend(level_warnings)
        current_id_to_chain: dict[tuple[int, int], str] = {}
        for record in records:
            output_id = record["output_id"]
            chain_id = previous_id_to_chain.get(output_id)
            if chain_id is None:
                chain_id = f"chain_{next_chain_index:08d}"
                next_chain_index += 1
            current_id_to_chain.setdefault(output_id, chain_id)
            state = targets.get(chain_id, {})
            _update_occurrence_accumulator(
                by_time_height,
                record.get("time_s"),
                record.get("height_m"),
                bool(state.get("has_large")),
                bool(state.get("has_coal") or state.get("has_coal_event")),
            )
        previous_id_to_chain = current_id_to_chain
    return by_time_height, warnings


def _read_selected_level(paths: list[Path], file_time_s: float | None, options: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    """Read one BW selected-output time level across ranks."""
    Dataset = _dataset_class()
    if Dataset is None:
        return [], ["NetCDF library unavailable; BW selected-output level skipped"]
    records: list[dict[str, Any]] = []
    warnings: list[str] = []
    chunk_size = int(options.get("chunk_size", 100000) or 100000)
    max_records = options.get("max_records")
    record_order = 0
    for path in paths:
        try:
            with Dataset(path, "r") as handle:
                names = handle.variables.keys()
                pair_vars = _find_pair_variables(names, BW_ID_NAMES)
                radius_name = _find_variable(names, RADIUS_NAMES)
                if pair_vars is None:
                    warnings.append(f"BW predecessor ID variables missing in {path.name}")
                    continue
                dm_var = handle.variables[pair_vars[0]]
                sd_var = handle.variables[pair_vars[1]]
                radius_var = handle.variables[radius_name] if radius_name else None
                height_name = _find_variable(names, HEIGHT_NAMES)
                coal_name = _find_variable(names, IF_COAL_NAMES)
                x_name = _find_variable(names, X_NAMES)
                y_name = _find_variable(names, Y_NAMES)
                n_name = _find_variable(names, N_NAMES)
                time_name = _find_variable(names, TIME_NAMES)
                time_var = handle.variables[time_name] if time_name else None
                length = min(_variable_length(dm_var), _variable_length(sd_var))
                if radius_var is not None:
                    length = min(length, _variable_length(radius_var))
                time_matches_records = time_var is not None and _variable_length(time_var) == length
                time_scalar = file_time_s
                if time_var is not None and not time_matches_records and _variable_length(time_var) <= 10000:
                    time_values_all = _flatten(time_var[:])
                    if len(time_values_all) == 1:
                        time_scalar = time_values_all[0]
                for selection in _iter_first_dim_chunks(length, max_records, chunk_size):
                    dm_values = _flatten_int(dm_var[selection])
                    sd_values = _flatten_int(sd_var[selection])
                    radius_values = _flatten(radius_var[selection]) if radius_var is not None else []
                    height_values = _get_optional_values(handle, height_name, selection)
                    coal_values = _get_optional_values(handle, coal_name, selection)
                    x_values = _get_optional_values(handle, x_name, selection)
                    y_values = _get_optional_values(handle, y_name, selection)
                    n_values = _get_optional_values(handle, n_name, selection)
                    time_values = _flatten(time_var[selection]) if time_var is not None and time_matches_records else []
                    n_records = min(len(dm_values), len(sd_values))
                    for index in range(n_records):
                        dm_id = dm_values[index]
                        sd_id = sd_values[index]
                        if dm_id < 0 or sd_id < 0:
                            continue
                        records.append(
                            {
                                "output_id": (dm_id, sd_id),
                                "dm_id": dm_id,
                                "sd_id": sd_id,
                                "time_s": _value_at(time_values, index) if time_matches_records else time_scalar,
                                "radius_m": _value_at(radius_values, index),
                                "height_m": _value_at(height_values, index),
                                "if_coal": _value_at(coal_values, index),
                                "x_m": _value_at(x_values, index),
                                "y_m": _value_at(y_values, index),
                                "multiplicity": _value_at(n_values, index),
                                "source_file": path.name,
                                "record_order": record_order,
                            }
                        )
                        record_order += 1
        except Exception as exc:
            warnings.append(f"failed to read BW selected-output level {path.name}: {exc}")
    return records, warnings


def _read_coal_events(paths: list[Path], file_time_s: float | None, options: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    """Read coalescence events whose BW predecessor IDs match the next output level."""
    Dataset = _dataset_class()
    if Dataset is None:
        return [], ["NetCDF library unavailable; coalescence event reads skipped"]
    events: list[dict[str, Any]] = []
    warnings: list[str] = []
    chunk_size = int(options.get("chunk_size", 100000) or 100000)
    max_records = options.get("max_records")
    for path in paths:
        try:
            with Dataset(path, "r") as handle:
                names = handle.variables.keys()
                pair1 = _find_pair(names, COAL_PRE1_NAMES)
                pair2 = _find_pair(names, COAL_PRE2_NAMES)
                if pair1 is None or pair2 is None:
                    continue
                dm1_var = handle.variables[pair1[0]]
                sd1_var = handle.variables[pair1[1]]
                dm2_var = handle.variables[pair2[0]]
                sd2_var = handle.variables[pair2[1]]
                length = min(_variable_length(dm1_var), _variable_length(sd1_var), _variable_length(dm2_var), _variable_length(sd2_var))
                time_name = _find_variable(names, COAL_EVENT_TIME_NAMES)
                r1_name = _find_variable(names, COAL_RADIUS1_NAMES)
                r2_name = _find_variable(names, COAL_RADIUS2_NAMES)
                n1_name = _find_variable(names, COAL_N1_NAMES)
                n2_name = _find_variable(names, COAL_N2_NAMES)
                num_col_name = _find_variable(names, COAL_NUM_COL_NAMES)
                for selection in _iter_first_dim_chunks(length, max_records, chunk_size):
                    dm1_values = _flatten_int(dm1_var[selection])
                    sd1_values = _flatten_int(sd1_var[selection])
                    dm2_values = _flatten_int(dm2_var[selection])
                    sd2_values = _flatten_int(sd2_var[selection])
                    time_values = _get_optional_values(handle, time_name, selection)
                    r1_values = _get_optional_values(handle, r1_name, selection)
                    r2_values = _get_optional_values(handle, r2_name, selection)
                    n1_values = _get_optional_values(handle, n1_name, selection)
                    n2_values = _get_optional_values(handle, n2_name, selection)
                    num_col_values = _get_optional_values(handle, num_col_name, selection)
                    n_events = min(len(dm1_values), len(sd1_values), len(dm2_values), len(sd2_values))
                    for index in range(n_events):
                        id1 = (dm1_values[index], sd1_values[index])
                        id2 = (dm2_values[index], sd2_values[index])
                        if id1[0] < 0 or id1[1] < 0 or id2[0] < 0 or id2[1] < 0:
                            continue
                        events.append(
                            {
                                "event_time_s": _value_at(time_values, index) if time_values else file_time_s,
                                "id1": id1,
                                "id2": id2,
                                "sd_r1_m": _value_at(r1_values, index),
                                "sd_r2_m": _value_at(r2_values, index),
                                "sd_n1": _value_at(n1_values, index),
                                "sd_n2": _value_at(n2_values, index),
                                "num_col": _value_at(num_col_values, index),
                                "source_file": path.name,
                            }
                        )
        except Exception as exc:
            warnings.append(f"failed to read BW coalescence events from {path.name}: {exc}")
    return events, warnings


def _source_metadata(chain_id: str, targets: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Return commonly used target metadata for optional source rows."""
    state = targets.get(chain_id)
    return {
        "category": _state_category(state),
        "max_radius_um": _state_um(state, "max_radius_m"),
        "first_large_time_s": state.get("first_large_time_s") if state else None,
        "first_ifcoal_time_s": state.get("first_coal_time_s") if state else None,
    }


def build_optional_source_tables(
    selected_files: list[Path],
    coal_files: list[Path] | None,
    targets: dict[str, dict[str, Any]],
    options: dict[str, Any],
) -> tuple[dict[str, list[dict[str, Any]]], list[str]]:
    """Build lightweight optional TPHT source tables for GMD candidate figures.

    The full TPHT target record matrix has about 36 million rows.  To keep the
    manuscript workflow practical, this function writes compact sources:
    deterministic full-run trajectory/interval samples, all positive
    ``if_coal`` records, and target-linked coalescence-event rows.
    """
    if options.get("metadata_only") or options.get("skip_heavy_netcdf"):
        return {}, ["metadata-only or skip-heavy-netcdf mode skipped optional TPHT source tables"]
    selected_groups, selected_warnings = _group_paths_by_time(selected_files, options.get("max_files"))
    coal_groups, coal_warnings = _group_paths_by_time(coal_files or [], options.get("max_files"))
    warnings = selected_warnings + coal_warnings
    if not selected_groups:
        return {}, warnings + ["No BW selected-output files found for optional TPHT source tables"]

    sample_stride = int(options.get("optional_source_stride", 200) or 200)

    trajectory_rows: list[dict[str, Any]] = []
    ifcoal_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []

    previous_id_to_chain: dict[tuple[int, int], str] = {}
    previous_record_by_chain: dict[str, dict[str, Any]] = {}
    previous_selected_time: float | None = None
    coal_group_index = 0
    next_chain_index = 0

    for level_index, (time_s, paths) in enumerate(selected_groups):
        records, level_warnings = _read_selected_level(paths, time_s, options)
        warnings.extend(level_warnings)
        current_id_to_chain: dict[tuple[int, int], str] = {}
        current_record_by_chain: dict[str, dict[str, Any]] = {}

        for record in records:
            output_id = record["output_id"]
            chain_id = previous_id_to_chain.get(output_id)
            if chain_id is None:
                chain_id = f"chain_{next_chain_index:08d}"
                next_chain_index += 1
            current_id_to_chain.setdefault(output_id, chain_id)
            current_record_by_chain[chain_id] = record
            metadata = _source_metadata(chain_id, targets)
            record_time = record.get("time_s")
            radius_um = (record.get("radius_m") * 1.0e6) if record.get("radius_m") is not None else None
            sampled = _sampled_chain(chain_id, sample_stride)
            if sampled and record_time is not None:
                trajectory_rows.append(
                    {
                        "target_id": chain_id,
                        "time_s": record_time,
                        "x_m": record.get("x_m"),
                        "y_m": record.get("y_m"),
                        "z_m": record.get("height_m"),
                        "radius_um": radius_um,
                        "if_coal": record.get("if_coal"),
                        "category": metadata["category"],
                        "max_radius_um": metadata["max_radius_um"],
                        "first_large_time_s": metadata["first_large_time_s"],
                        "first_ifcoal_time_s": metadata["first_ifcoal_time_s"],
                        "warnings": None,
                    }
                )
            coal_active = record.get("if_coal") is not None and record.get("if_coal") > 0.0
            if coal_active:
                ifcoal_rows.append(
                    {
                        "target_id": chain_id,
                        "time_s": record_time,
                        "if_coal": record.get("if_coal"),
                        "category": metadata["category"],
                        "max_radius_um": metadata["max_radius_um"],
                        "first_large_time_s": metadata["first_large_time_s"],
                        "first_ifcoal_time_s": metadata["first_ifcoal_time_s"],
                        "warnings": None,
                    }
                )
            previous_record = previous_record_by_chain.get(chain_id)
            if previous_record is not None:
                previous_time = previous_record.get("time_s")
                previous_radius = previous_record.get("radius_m")
                current_radius = record.get("radius_m")
                if previous_time is not None and record_time is not None and current_radius is not None and previous_radius is not None:
                    dt_s = record_time - previous_time
                    interval_coal = coal_active or (previous_record.get("if_coal") is not None and previous_record.get("if_coal") > 0.0)
                    if dt_s > 0.0 and (sampled or interval_coal):
                        interval_rows.append(
                            {
                                "target_id": chain_id,
                                "time_s": previous_time,
                                "dt_s": dt_s,
                                "r0_um": previous_radius * 1.0e6,
                                "r1_um": current_radius * 1.0e6,
                                "height_m": previous_record.get("height_m"),
                                "category": metadata["category"],
                                "if_coal_interval": 1 if interval_coal else 0,
                                "sampled_interval": 1 if sampled else 0,
                                "warnings": None,
                            }
                        )

        coal_paths: list[Path] = []
        if level_index > 0 and previous_selected_time is not None and time_s is not None:
            while coal_group_index < len(coal_groups):
                coal_time, group_paths = coal_groups[coal_group_index]
                if coal_time is None:
                    break
                if coal_time <= previous_selected_time + 1.0e-9:
                    coal_group_index += 1
                    continue
                if previous_selected_time - 1.0e-9 < coal_time <= time_s + 1.0e-9:
                    coal_paths.extend(group_paths)
                    coal_group_index += 1
                    continue
                break
        if coal_paths:
            events, event_warnings = _read_coal_events(coal_paths, time_s, options)
            warnings.extend(event_warnings)
            for event in events:
                for participant, key, radius_key in (
                    (1, event.get("id1"), "sd_r1_m"),
                    (2, event.get("id2"), "sd_r2_m"),
                ):
                    chain_id = current_id_to_chain.get(key) if key is not None else None
                    if chain_id is None:
                        continue
                    target_record = current_record_by_chain.get(chain_id, {})
                    post_radius_m = target_record.get("radius_m")
                    pre_radius_m = event.get(radius_key)
                    metadata = _source_metadata(chain_id, targets)
                    event_rows.append(
                        {
                            "target_id": chain_id,
                            "event_time_s": event.get("event_time_s"),
                            "event_height_m": target_record.get("height_m"),
                            "pre_radius_um": pre_radius_m * 1.0e6 if pre_radius_m is not None else None,
                            "post_radius_um": post_radius_m * 1.0e6 if post_radius_m is not None else None,
                            "delta_r3_um3": (post_radius_m * 1.0e6) ** 3 - (pre_radius_m * 1.0e6) ** 3 if post_radius_m is not None and pre_radius_m is not None else None,
                            "participant": participant,
                            "num_col": event.get("num_col"),
                            "category": metadata["category"],
                            "source_file": event.get("source_file"),
                            "warnings": None,
                        }
                    )

        previous_id_to_chain = current_id_to_chain
        previous_record_by_chain = current_record_by_chain
        previous_selected_time = time_s

    return (
        {
            "trajectory": trajectory_rows,
            "ifcoal_timeline": ifcoal_rows,
            "interval": interval_rows,
            "event_links": event_rows,
        },
        warnings,
    )


def _update_chain_from_record(state: dict[str, Any], record: dict[str, Any], order: int) -> None:
    """Update one propagated chain with a selected-output record."""
    time_s = record.get("time_s")
    radius_m = record.get("radius_m")
    height_m = record.get("height_m")
    coal_value = record.get("if_coal")
    state["record_count"] += 1
    if state["first_time_s"] is None or (time_s is not None and time_s < state["first_time_s"]):
        state["first_time_s"] = time_s
    if _is_earlier(state.get("last_time_s"), state.get("last_order", -1), time_s, order):
        state["last_time_s"] = time_s
        state["last_order"] = order
        state["final_radius_m"] = radius_m
    if radius_m is not None:
        state["max_radius_m"] = max(radius_m, state["max_radius_m"] or radius_m)
        if radius_m >= LARGE_RADIUS_M:
            state["has_large"] = True
            if _is_earlier(time_s, order, state.get("first_large_time_s"), state.get("first_large_order")):
                state["first_large_time_s"] = time_s
                state["first_large_order"] = order
                state["first_large_radius_m"] = radius_m
                state["first_large_height_m"] = height_m
                state["first_large_x_m"] = record.get("x_m")
                state["first_large_y_m"] = record.get("y_m")
    selected_category = _first_selected_category(radius_m, coal_value)
    if selected_category is not None and _is_earlier(time_s, order, state.get("first_selected_time_s"), state.get("first_selected_order")):
        state["first_selected_time_s"] = time_s
        state["first_selected_order"] = order
        state["first_selected_category"] = selected_category
    coal_active = coal_value is not None and coal_value > 0.0
    if coal_active:
        state["has_coal"] = True
        state["coal_record_count"] += 1
        if _is_earlier(time_s, order, state.get("first_coal_time_s"), state.get("first_coal_order")):
            state["first_coal_time_s"] = time_s
            state["first_coal_order"] = order
        if not state["previous_coal_active"]:
            state["coal_episode_count_proxy"] += 1
    state["previous_coal_active"] = coal_active
    previous_time = state.get("previous_time_s")
    previous_radius = state.get("previous_radius_m")
    if previous_time is not None and previous_radius is not None and time_s is not None and time_s > previous_time:
        delta_t = time_s - previous_time
        if previous_radius >= 15.0e-6:
            state["duration_ge_15_s"] += delta_t
        if previous_radius >= 20.0e-6:
            state["duration_ge_20_s"] += delta_t
    state["previous_time_s"] = time_s
    state["previous_radius_m"] = radius_m


def _update_chain_from_event(state: dict[str, Any], event: dict[str, Any], order: int) -> None:
    """Update one propagated chain with a linked BW coalescence event."""
    event_time = event.get("event_time_s")
    state["has_coal"] = True
    state["has_coal_event"] = True
    state["coal_event_count"] += 1
    state["coal_num_col_sum"] += float(event.get("num_col") or 0.0)
    if event_time is not None and _is_earlier(event_time, order, state.get("first_coal_time_s"), state.get("first_coal_order")):
        state["first_coal_time_s"] = event_time
        state["first_coal_order"] = order
    if event_time is not None and _is_earlier(event_time, order, state.get("first_coal_event_time_s"), state.get("first_coal_event_order")):
        state["first_coal_event_time_s"] = event_time
        state["first_coal_event_order"] = order


def build_stepwise_bw_diagnostics(
    selected_files: list[Path],
    coal_files: list[Path] | None,
    options: dict[str, Any],
) -> tuple[dict[str, dict[str, Any]], dict[float, dict[str, float]], dict[tuple[float, float], dict[str, float]], list[dict[str, Any]], list[str]]:
    """Build BW TPHT diagnostics by following adjacent output-level predecessor links.

    Returns per-chain lifecycle summaries, time aggregates, height-time aggregates,
    per-output link diagnostics, and warnings.  The method streams one output time
    group at a time and keeps only the previous level's ID-to-chain map in memory.
    """
    if options.get("metadata_only") or options.get("skip_heavy_netcdf"):
        return {}, {}, {}, [], ["metadata-only or skip-heavy-netcdf mode skipped stepwise BW chain reconstruction"]
    Dataset = _dataset_class()
    if Dataset is None:
        return {}, {}, {}, [], ["NetCDF library unavailable; stepwise BW chain reconstruction skipped"]

    selected_groups, selected_warnings = _group_paths_by_time(selected_files, options.get("max_files"))
    coal_groups, coal_warnings = _group_paths_by_time(coal_files or [], options.get("max_files"))
    warnings = selected_warnings + coal_warnings
    if not selected_groups:
        return {}, {}, {}, [], warnings + ["No BW selected-output files found for stepwise chain reconstruction"]

    targets: dict[str, dict[str, Any]] = {}
    by_time: dict[float, dict[str, float]] = {}
    by_time_height: dict[tuple[float, float], dict[str, float]] = {}
    link_rows: list[dict[str, Any]] = []
    previous_id_to_chain: dict[tuple[int, int], str] = {}
    previous_selected_time: float | None = None
    coal_group_index = 0
    next_chain_index = 0
    global_order = 0

    for level_index, (time_s, paths) in enumerate(selected_groups):
        records, level_warnings = _read_selected_level(paths, time_s, options)
        warnings.extend(level_warnings)
        current_id_to_chain: dict[tuple[int, int], str] = {}
        duplicate_count = 0
        matched_count = 0
        unmatched_count = 0
        for record in records:
            output_id = record["output_id"]
            chain_id = previous_id_to_chain.get(output_id)
            if chain_id is None:
                if level_index > 0:
                    unmatched_count += 1
                chain_id = f"chain_{next_chain_index:08d}"
                next_chain_index += 1
                targets[chain_id] = _new_chain_state(chain_id, output_id[0], output_id[1])
            else:
                matched_count += 1
            if output_id in current_id_to_chain:
                duplicate_count += 1
            else:
                current_id_to_chain[output_id] = chain_id
            _update_chain_from_record(targets[chain_id], record, global_order)
            _update_time_accumulator(by_time, record.get("time_s"), record.get("radius_m"), record.get("height_m"), record.get("if_coal"))
            global_order += 1

        linked_events = 0
        coal_paths: list[Path] = []
        if level_index > 0 and previous_selected_time is not None and time_s is not None:
            while coal_group_index < len(coal_groups):
                coal_time, group_paths = coal_groups[coal_group_index]
                if coal_time is None:
                    break
                if coal_time <= previous_selected_time + 1.0e-9:
                    coal_group_index += 1
                    continue
                if previous_selected_time - 1.0e-9 < coal_time <= time_s + 1.0e-9:
                    coal_paths.extend(group_paths)
                    coal_group_index += 1
                    continue
                break
        if coal_paths:
            events, event_warnings = _read_coal_events(coal_paths, time_s, options)
            warnings.extend(event_warnings)
            for event in events:
                seen_chains: set[str] = set()
                for key in (event.get("id1"), event.get("id2")):
                    chain_id = current_id_to_chain.get(key)
                    if chain_id is None or chain_id in seen_chains:
                        continue
                    _update_chain_from_event(targets[chain_id], event, global_order)
                    seen_chains.add(chain_id)
                    linked_events += 1
                global_order += 1

        valid_fraction = safe_ratio(matched_count, len(records)) if level_index > 0 else None
        link_rows.append(
            {
                "time_s": time_s,
                "record_count": len(records) if records else None,
                "matched_parent_records": matched_count if level_index > 0 else None,
                "unmatched_parent_records": unmatched_count if level_index > 0 else None,
                "valid_link_fraction": valid_fraction,
                "unique_output_ids": len(current_id_to_chain) if current_id_to_chain else None,
                "duplicate_output_ids": duplicate_count if duplicate_count else 0,
                "active_chain_count": len(set(current_id_to_chain.values())) if current_id_to_chain else None,
                "coalescence_event_count_linked": linked_events if linked_events else 0,
                "warnings": None,
            }
        )
        previous_id_to_chain = current_id_to_chain
        previous_selected_time = time_s

    remaining_coal_groups = len(coal_groups) - coal_group_index
    if remaining_coal_groups > 0:
        warnings.append(
            f"{remaining_coal_groups} coalescence-output time groups were not linked because they are at or after the last selected-output level"
        )

    link_fractions = [row["valid_link_fraction"] for row in link_rows if row.get("valid_link_fraction") is not None]
    if link_fractions and sum(link_fractions) / len(link_fractions) < 0.5:
        warnings.append(
            "Average adjacent-output BW link fraction is below 0.5; verify selected-output cadence and BW time-direction convention before interpreting histories"
        )
    by_time_height, occurrence_warnings = _build_category_occurrence_accumulator(selected_groups, targets, options)
    warnings.extend(occurrence_warnings)
    return targets, by_time, by_time_height, link_rows, warnings


def build_predecessor_tree_examples(
    selected_files: list[Path],
    coal_files: list[Path] | None,
    target_ids: set[str],
    options: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Build compact predecessor-tree rows for selected TPHT target chains.

    The output is intentionally small: it stores the tracked target's selected
    output nodes and short coalescence-flagged partner branches for a few
    representative chains.  Coalescence logs do not carry event height, so event
    branches use the matched target height at the next selected-output level.
    """
    if options.get("metadata_only") or options.get("skip_heavy_netcdf"):
        return [], ["metadata-only or skip-heavy-netcdf mode skipped predecessor-tree example generation"]
    if not target_ids:
        return [], ["No target IDs were provided for predecessor-tree example generation"]
    Dataset = _dataset_class()
    if Dataset is None:
        return [], ["NetCDF library unavailable; predecessor-tree example generation skipped"]

    selected_groups, selected_warnings = _group_paths_by_time(selected_files, options.get("max_files"))
    coal_groups, coal_warnings = _group_paths_by_time(coal_files or [], options.get("max_files"))
    warnings = selected_warnings + coal_warnings
    if not selected_groups:
        return [], warnings + ["No BW selected-output files found for predecessor-tree example generation"]

    rows: list[dict[str, Any]] = []
    previous_id_to_chain: dict[tuple[int, int], str] = {}
    previous_node_by_chain: dict[str, str] = {}
    previous_selected_time: float | None = None
    coal_group_index = 0
    next_chain_index = 0
    event_counter = 0

    for level_index, (time_s, paths) in enumerate(selected_groups):
        records, level_warnings = _read_selected_level(paths, time_s, options)
        warnings.extend(level_warnings)
        current_id_to_chain: dict[tuple[int, int], str] = {}
        current_node_by_chain: dict[str, str] = {}
        current_record_by_chain: dict[str, dict[str, Any]] = {}

        for record_index, record in enumerate(records):
            output_id = record["output_id"]
            chain_id = previous_id_to_chain.get(output_id)
            if chain_id is None:
                chain_id = f"chain_{next_chain_index:08d}"
                next_chain_index += 1
            current_id_to_chain.setdefault(output_id, chain_id)
            if chain_id not in target_ids:
                continue
            node_id = f"{chain_id}:out{level_index:04d}:{record_index:08d}"
            current_node_by_chain[chain_id] = node_id
            current_record_by_chain[chain_id] = record
            rows.append(
                {
                    "example_id": chain_id,
                    "target_id": chain_id,
                    "node_id": node_id,
                    "predecessor_id": previous_node_by_chain.get(chain_id),
                    "time_s": record.get("time_s"),
                    "height_m": record.get("height_m"),
                    "radius_um": (record.get("radius_m") * 1.0e6) if record.get("radius_m") is not None else None,
                    "if_coal": record.get("if_coal"),
                    "node_role": "tracked_target_output",
                    "event_time_s": None,
                    "num_col": None,
                    "source_file": record.get("source_file"),
                    "notes": "tracked target node from BW selected output",
                    "warnings": None,
                }
            )

        coal_paths: list[Path] = []
        if level_index > 0 and previous_selected_time is not None and time_s is not None:
            while coal_group_index < len(coal_groups):
                coal_time, group_paths = coal_groups[coal_group_index]
                if coal_time is None:
                    break
                if coal_time <= previous_selected_time + 1.0e-9:
                    coal_group_index += 1
                    continue
                if previous_selected_time - 1.0e-9 < coal_time <= time_s + 1.0e-9:
                    coal_paths.extend(group_paths)
                    coal_group_index += 1
                    continue
                break

        if coal_paths:
            events, event_warnings = _read_coal_events(coal_paths, time_s, options)
            warnings.extend(event_warnings)
            for event in events:
                for participant_number, key, radius_key in (
                    (1, event.get("id1"), "sd_r1_m"),
                    (2, event.get("id2"), "sd_r2_m"),
                ):
                    if key is None:
                        continue
                    chain_id = current_id_to_chain.get(key)
                    if chain_id not in target_ids:
                        continue
                    target_node = current_node_by_chain.get(chain_id)
                    target_record = current_record_by_chain.get(chain_id, {})
                    event_counter += 1
                    rows.append(
                        {
                            "example_id": chain_id,
                            "target_id": chain_id,
                            "node_id": f"{chain_id}:event{event_counter:06d}:p{participant_number}",
                            "predecessor_id": target_node,
                            "time_s": event.get("event_time_s"),
                            "height_m": target_record.get("height_m"),
                            "radius_um": (event.get(radius_key) * 1.0e6) if event.get(radius_key) is not None else None,
                            "if_coal": 1,
                            "node_role": "coalescence_partner_proxy",
                            "event_time_s": event.get("event_time_s"),
                            "num_col": event.get("num_col"),
                            "source_file": event.get("source_file"),
                            "notes": "partner branch uses target height at next selected-output level; if_coal is a binary flag, not collision count",
                            "warnings": None,
                        }
                    )

        previous_id_to_chain = current_id_to_chain
        previous_node_by_chain.update(current_node_by_chain)
        previous_selected_time = time_s

    if not rows:
        warnings.append("No predecessor-tree rows were generated for the requested target IDs")
    return rows, warnings


def chain_lengths(targets: dict[str, dict[str, Any]]) -> dict[str, int]:
    """Return a compact chain-length distribution."""
    distribution: dict[str, int] = {}
    for state in targets.values():
        count = int(state.get("record_count") or 0)
        if count < 5:
            label = str(count)
        elif count < 20:
            label = "5-19"
        elif count < 100:
            label = "20-99"
        else:
            label = ">=100"
        distribution[label] = distribution.get(label, 0) + 1
    return distribution


def parse_chain_distribution_text(distribution: dict[str, int]) -> str | None:
    """Serialize a compact chain-length distribution for table output."""
    if not distribution:
        return None
    order = ["0", "1", "2", "3", "4", "5-19", "20-99", ">=100"]
    return "; ".join(f"{label}:{distribution[label]}" for label in order if label in distribution)
