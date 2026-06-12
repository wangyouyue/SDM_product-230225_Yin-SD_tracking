#!/usr/bin/env python3
"""Export GMD2026 TPHT QC/QR fields for 3-D Matplotlib rendering.

This SQUID-side post-processing script reads distributed SCALE-SDM
``history.pe*`` files from the TPHT backward-reconstruction case and writes one
compact NetCDF file containing only the Eulerian QC/QR fields needed by the
3-D diagnostic video.  Missing values remain NaN; they are never converted to
zero.  The output is intentionally separated from ``analysis_outputs`` so the
manuscript candidate products can be regenerated without overwriting analysis
tables.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None


DEFAULT_GROUP = "02_tpht_3d_interest_70min"
DEFAULT_CASE = "bw_reconstruction"
TIME_NAMES = ("time", "TIME", "t", "time_s")
QC_CANDIDATES = ("QC_sd", "QC", "QC_hyd", "QHYD_sd", "QHYD")
QR_CANDIDATES = ("QR_sd", "QR", "QR_hyd")
X_DIM_NAMES = ("x", "cx", "i", "ix")
Y_DIM_NAMES = ("y", "cy", "j", "jy")
Z_DIM_NAMES = ("z", "cz", "k", "lev", "level")


@dataclass(frozen=True)
class GridSpec:
    """GMD2026 local/global grid metadata parsed from SCALE namelists."""

    prc_num_x: int
    prc_num_y: int
    imax: int
    jmax: int
    kmax: int
    dx: float
    dy: float
    dz: float

    @property
    def nx(self) -> int:
        return self.prc_num_x * self.imax

    @property
    def ny(self) -> int:
        return self.prc_num_y * self.jmax


@dataclass(frozen=True)
class AxisMap:
    """Axis mapping for one history variable after the time dimension is sliced."""

    z_axis: int
    y_axis: int
    x_axis: int
    time_axis: int | None


def parse_args() -> argparse.Namespace:
    """Return command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Path to the GMD2026 directory.")
    parser.add_argument("--case-dir", type=Path, default=None, help="Direct path to the BW case directory.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output root for exported 3-D fields.")
    parser.add_argument("--history-pattern", default="history.pe*", help="History file glob inside --case-dir.")
    parser.add_argument("--qc-variable", default=None, help="Override source QC variable name.")
    parser.add_argument("--qr-variable", default=None, help="Override source QR variable name.")
    parser.add_argument("--stride-x", type=int, default=1, help="Global x stride for exported fields.")
    parser.add_argument("--stride-y", type=int, default=1, help="Global y stride for exported fields.")
    parser.add_argument("--stride-z", type=int, default=1, help="Global z stride for exported fields.")
    parser.add_argument("--time-stride", type=int, default=1, help="Output-time stride.")
    parser.add_argument("--compression-level", type=int, default=4, help="NetCDF zlib compression level.")
    parser.add_argument("--max-files", type=int, default=None, help="Debugging limit on history.pe* file count.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work without reading field arrays.")
    return parser.parse_args()


def positive_stride(value: int, name: str) -> int:
    """Validate a positive stride value."""
    stride = int(value)
    if stride < 1:
        raise ValueError(f"{name} must be >= 1")
    return stride


def strip_namelist_comments(text: str) -> str:
    """Remove simple Fortran namelist comments before regex parsing."""
    lines = []
    for line in text.splitlines():
        lines.append(line.split("!", 1)[0])
    return "\n".join(lines)


def parse_namelist_number(text: str, name: str, default: float | None = None) -> float:
    """Parse a scalar numeric namelist value."""
    clean = strip_namelist_comments(text)
    pattern = re.compile(rf"\b{re.escape(name)}\s*=\s*([^,\n/]+)", re.IGNORECASE)
    match = pattern.search(clean)
    if not match:
        if default is None:
            raise ValueError(f"Required namelist value {name} was not found")
        return float(default)
    value = match.group(1).strip().strip("'\"").replace("D", "E").replace("d", "e")
    return float(value)


def case_config_text(case_dir: Path) -> str:
    """Return concatenated run/init config text for a GMD2026 case."""
    parts = []
    for name in ("run.conf", "init.conf"):
        path = case_dir / name
        if path.exists():
            parts.append(path.read_text(encoding="utf-8", errors="ignore"))
    return "\n".join(parts)


def read_grid_spec(case_dir: Path) -> GridSpec:
    """Read the TPHT case grid from SCALE namelist files."""
    text = case_config_text(case_dir)
    if not text:
        raise FileNotFoundError(f"No run.conf or init.conf found in {case_dir}")
    return GridSpec(
        prc_num_x=int(parse_namelist_number(text, "PRC_NUM_X")),
        prc_num_y=int(parse_namelist_number(text, "PRC_NUM_Y")),
        imax=int(parse_namelist_number(text, "IMAX")),
        jmax=int(parse_namelist_number(text, "JMAX")),
        kmax=int(parse_namelist_number(text, "KMAX")),
        dx=float(parse_namelist_number(text, "DX")),
        dy=float(parse_namelist_number(text, "DY")),
        dz=float(parse_namelist_number(text, "DZ")),
    )


def rank_id_from_history(path: Path) -> int:
    """Parse the rank id from a SCALE ``history.peXXXXXX`` filename."""
    match = re.search(r"\.pe(\d+)", path.name)
    if not match:
        raise ValueError(f"Cannot parse rank id from {path.name}")
    return int(match.group(1))


def history_files(case_dir: Path, pattern: str, max_files: int | None) -> list[Path]:
    """Return sorted history files, optionally truncated for debugging."""
    files = sorted(path for path in case_dir.glob(pattern) if path.is_file())
    if max_files is not None and max_files >= 0:
        files = files[:max_files]
    return files


def find_time_variable(dataset: Dataset) -> str | None:
    """Find a likely time coordinate variable."""
    names = set(dataset.variables)
    for name in TIME_NAMES:
        if name in names:
            return name
    for name in names:
        lower = name.lower()
        if lower == "time" or lower.endswith("_time"):
            return name
    return None


def read_times(dataset: Dataset) -> np.ndarray:
    """Read output times in seconds from one history file."""
    time_name = find_time_variable(dataset)
    if time_name is not None:
        values = np.asarray(dataset.variables[time_name][:], dtype=float)
        if hasattr(values, "filled"):
            values = values.filled(np.nan)
        finite = values[np.isfinite(values)]
        if finite.size:
            return np.asarray(finite, dtype=float)
    for dim_name, dim in dataset.dimensions.items():
        if "time" in dim_name.lower():
            return np.arange(len(dim), dtype=float)
    return np.asarray([0.0], dtype=float)


def choose_variable(dataset: Dataset, override: str | None, candidates: tuple[str, ...], label: str) -> str:
    """Choose a source variable from a candidate list."""
    if override:
        if override not in dataset.variables:
            raise KeyError(f"Requested {label} variable {override!r} is absent")
        return override
    for name in candidates:
        if name in dataset.variables:
            return name
    available = ", ".join(sorted(dataset.variables)[:80])
    raise KeyError(f"No {label} variable found from {candidates}. Available variables start with: {available}")


def infer_time_axis(variable: Any, time_length: int) -> int | None:
    """Infer the time axis for a NetCDF variable."""
    dims = tuple(getattr(variable, "dimensions", ()))
    for index, name in enumerate(dims):
        if "time" in name.lower():
            return index
    shape = tuple(getattr(variable, "shape", ()))
    for index, length in enumerate(shape):
        if int(length) == int(time_length):
            return index
    return None


def _match_dim(name: str, candidates: tuple[str, ...]) -> bool:
    lower = name.lower()
    return lower in candidates or any(lower.endswith(candidate) for candidate in candidates)


def _nearest_unused_axis(shape: tuple[int, ...], expected: int, used: set[int]) -> int | None:
    candidates = [(abs(int(length) - expected), index) for index, length in enumerate(shape) if index not in used and int(length) > 1]
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][1]


def infer_axis_map(variable: Any, time_axis: int | None, grid: GridSpec) -> AxisMap:
    """Infer z/y/x axes after slicing out the time dimension."""
    original_dims = tuple(getattr(variable, "dimensions", ()))
    original_shape = tuple(int(value) for value in getattr(variable, "shape", ()))
    dims = [name for index, name in enumerate(original_dims) if index != time_axis]
    shape = tuple(length for index, length in enumerate(original_shape) if index != time_axis)
    if len(shape) < 3:
        raise ValueError(f"Expected a 3-D field after time slicing, got shape {shape} for {getattr(variable, 'name', 'variable')}")

    axis_x = axis_y = axis_z = None
    for index, name in enumerate(dims):
        if axis_x is None and _match_dim(name, X_DIM_NAMES):
            axis_x = index
        elif axis_y is None and _match_dim(name, Y_DIM_NAMES):
            axis_y = index
        elif axis_z is None and _match_dim(name, Z_DIM_NAMES):
            axis_z = index

    used = {axis for axis in (axis_x, axis_y, axis_z) if axis is not None}
    if axis_z is None:
        axis_z = _nearest_unused_axis(shape, grid.kmax, used)
        if axis_z is not None:
            used.add(axis_z)
    if axis_y is None and len(shape) >= 3:
        axis_y = 1 if 1 not in used else _nearest_unused_axis(shape, grid.jmax, used)
        if axis_y is not None:
            used.add(axis_y)
    if axis_x is None and len(shape) >= 3:
        axis_x = 2 if 2 not in used else _nearest_unused_axis(shape, grid.imax, used)
    if axis_x is None or axis_y is None or axis_z is None:
        raise ValueError(f"Failed to infer z/y/x axes for {getattr(variable, 'name', 'variable')} dims={original_dims} shape={original_shape}")
    return AxisMap(z_axis=int(axis_z), y_axis=int(axis_y), x_axis=int(axis_x), time_axis=time_axis)


def physical_indices(axis_size: int, expected_count: int) -> np.ndarray:
    """Return centered physical-domain indices without treating missing halos as zero."""
    axis_size = int(axis_size)
    expected_count = int(expected_count)
    if expected_count > 0 and axis_size >= expected_count:
        start = max(0, (axis_size - expected_count) // 2)
        return np.arange(start, start + expected_count, dtype=int)
    if axis_size > 2:
        return np.arange(1, axis_size - 1, dtype=int)
    return np.arange(axis_size, dtype=int)


def read_field_block(variable: Any, axis_map: AxisMap, time_index: int, grid: GridSpec) -> np.ndarray:
    """Read one local rank/time block as ``(z, y, x)``."""
    slices = [slice(None)] * len(getattr(variable, "shape", ()))
    if axis_map.time_axis is not None:
        slices[axis_map.time_axis] = time_index
    raw = variable[tuple(slices)]
    array = np.asarray(raw, dtype=np.float32)
    if hasattr(raw, "filled"):
        array = raw.filled(np.nan).astype(np.float32)
    array = np.squeeze(array)
    moved = np.moveaxis(array, (axis_map.z_axis, axis_map.y_axis, axis_map.x_axis), (0, 1, 2))
    z_indices = physical_indices(moved.shape[0], grid.kmax)
    y_indices = physical_indices(moved.shape[1], grid.jmax)
    x_indices = physical_indices(moved.shape[2], grid.imax)
    block = moved[np.ix_(z_indices, y_indices, x_indices)]
    if block.shape != (grid.kmax, grid.jmax, grid.imax):
        raise ValueError(f"Unexpected cropped block shape {block.shape}; expected {(grid.kmax, grid.jmax, grid.imax)}")
    return block


def create_output(path: Path, case_dir: Path, grid: GridSpec, times: np.ndarray, strides: tuple[int, int, int], time_stride: int, compression_level: int) -> tuple[Dataset, dict[str, Any]]:
    """Create the compact rendering NetCDF output file."""
    sx, sy, sz = strides
    selected_times = times[::time_stride]
    x = (np.arange(grid.nx, dtype=float)[::sx] + 0.5) * grid.dx
    y = (np.arange(grid.ny, dtype=float)[::sy] + 0.5) * grid.dy
    z = (np.arange(grid.kmax, dtype=float)[::sz] + 0.5) * grid.dz

    path.parent.mkdir(parents=True, exist_ok=True)
    dataset = Dataset(path, "w", format="NETCDF4")
    dataset.createDimension("time", selected_times.size)
    dataset.createDimension("z", z.size)
    dataset.createDimension("y", y.size)
    dataset.createDimension("x", x.size)
    dataset.case_name = DEFAULT_CASE
    dataset.source_case_path = str(case_dir)
    dataset.dx = float(grid.dx * sx)
    dataset.dy = float(grid.dy * sy)
    dataset.dz = float(grid.dz * sz)
    dataset.stride_x = sx
    dataset.stride_y = sy
    dataset.stride_z = sz
    dataset.time_stride = time_stride
    dataset.note = "GMD2026 SCALE-SDM TPHT QC/QR export for 3-D diagnostic video"

    for name, values, units, long_name in (
        ("time", selected_times, "s", "model time"),
        ("x", x, "m", "x coordinate"),
        ("y", y, "m", "y coordinate"),
        ("z", z, "m", "height AGL"),
    ):
        var = dataset.createVariable(name, "f8", (name,))
        var[:] = values
        var.units = units
        var.long_name = long_name

    chunksizes = (1, max(1, min(z.size, 64)), max(1, min(y.size, 32)), max(1, min(x.size, 32)))
    variables = {}
    for name, long_name in (("QC_hyd", "cloud water or hydrometeor mixing ratio"), ("QR_hyd", "rain water mixing ratio")):
        var = dataset.createVariable(
            name,
            "f4",
            ("time", "z", "y", "x"),
            zlib=True,
            complevel=compression_level,
            chunksizes=chunksizes,
            fill_value=np.float32(np.nan),
        )
        var.units = "kg kg-1"
        var.long_name = long_name
        variables[name] = var
    return dataset, variables


def output_indices_for_rank(rank_id: int, grid: GridSpec, strides: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return local and output indices for one rank block."""
    sx, sy, sz = strides
    rank_x = rank_id % grid.prc_num_x
    rank_y = rank_id // grid.prc_num_x
    global_x = rank_x * grid.imax + np.arange(grid.imax, dtype=int)
    global_y = rank_y * grid.jmax + np.arange(grid.jmax, dtype=int)
    global_z = np.arange(grid.kmax, dtype=int)
    local_x = np.flatnonzero(global_x % sx == 0)
    local_y = np.flatnonzero(global_y % sy == 0)
    local_z = np.flatnonzero(global_z % sz == 0)
    out_x = global_x[local_x] // sx
    out_y = global_y[local_y] // sy
    out_z = global_z[local_z] // sz
    return local_z, local_y, local_x, out_z, out_y, out_x


def write_rank_block(variable: Any, output_time_index: int, block: np.ndarray, local_indices: tuple[np.ndarray, np.ndarray, np.ndarray], output_indices: tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    """Write one possibly strided rank block into the global output variable."""
    local_z, local_y, local_x = local_indices
    out_z, out_y, out_x = output_indices
    selected = block[np.ix_(local_z, local_y, local_x)]
    for z_index, global_z in enumerate(out_z):
        for y_index, global_y in enumerate(out_y):
            variable[output_time_index, int(global_z), int(global_y), out_x] = selected[z_index, y_index, :]


def validate_export(path: Path) -> None:
    """Fail fast if the exported fields contain no finite information."""
    with Dataset(path, "r") as dataset:
        problems: list[str] = []
        for name in ("QC_hyd", "QR_hyd"):
            variable = dataset.variables[name]
            finite_count = 0
            for time_index in sorted({0, len(dataset.dimensions["time"]) // 2, len(dataset.dimensions["time"]) - 1}):
                values = np.asarray(variable[time_index], dtype=float)
                finite_count += int(np.count_nonzero(np.isfinite(values)))
            if finite_count == 0:
                problems.append(f"{name} has no finite values in sampled frames")
        if problems:
            raise RuntimeError("; ".join(problems))


def export_case(args: argparse.Namespace) -> Path:
    """Export the BW-reconstruction QC/QR fields."""
    root = args.root.resolve()
    case_dir = (args.case_dir or root / DEFAULT_GROUP / DEFAULT_CASE).resolve()
    outdir = (args.outdir or root / "analysis_outputs_for_GMD" / "video_sources" / "tpht_3d").resolve()
    strides = (
        positive_stride(args.stride_x, "stride-x"),
        positive_stride(args.stride_y, "stride-y"),
        positive_stride(args.stride_z, "stride-z"),
    )
    time_stride = positive_stride(args.time_stride, "time-stride")
    files = history_files(case_dir, args.history_pattern, args.max_files)
    grid = read_grid_spec(case_dir)
    output_path = outdir / "tpht_bw_qc_qr_3d_fullgrid.nc"
    if args.dry_run:
        print("[dry-run] GMD2026 TPHT 3-D QC/QR export")
        print(f"[dry-run] case_dir={case_dir}")
        print(f"[dry-run] history_files={len(files)}")
        print(f"[dry-run] grid={grid}")
        print("[dry-run] qc_variable=auto qr_variable=auto")
        print(f"[dry-run] output={output_path}")
        return output_path

    if Dataset is None:
        raise RuntimeError("netCDF4 is required. Activate the GMD2026 sdm_env environment on SQUID.")
    if not files:
        raise FileNotFoundError(f"No history files matched {case_dir / args.history_pattern}")

    with Dataset(files[0], "r") as first:
        times = read_times(first)
        qc_name = choose_variable(first, args.qc_variable, QC_CANDIDATES, "QC")
        qr_name = choose_variable(first, args.qr_variable, QR_CANDIDATES, "QR")
        qc_axis = infer_axis_map(first.variables[qc_name], infer_time_axis(first.variables[qc_name], len(times)), grid)
        qr_axis = infer_axis_map(first.variables[qr_name], infer_time_axis(first.variables[qr_name], len(times)), grid)

    dataset, out_vars = create_output(output_path, case_dir, grid, times, strides, time_stride, args.compression_level)
    dataset.variables["QC_hyd"].source_variable = qc_name
    dataset.variables["QR_hyd"].source_variable = qr_name
    selected_time_indices = np.arange(times.size, dtype=int)[::time_stride]
    try:
        for file_index, path in enumerate(files, start=1):
            rank_id = rank_id_from_history(path)
            if rank_id >= grid.prc_num_x * grid.prc_num_y:
                print(f"[warning] skipping rank outside configured decomposition: {path}", flush=True)
                continue
            local_z, local_y, local_x, out_z, out_y, out_x = output_indices_for_rank(rank_id, grid, strides)
            with Dataset(path, "r") as source:
                qc_var = source.variables[qc_name]
                qr_var = source.variables[qr_name]
                for output_time_index, source_time_index in enumerate(selected_time_indices):
                    qc_block = read_field_block(qc_var, qc_axis, int(source_time_index), grid)
                    qr_block = read_field_block(qr_var, qr_axis, int(source_time_index), grid)
                    write_rank_block(out_vars["QC_hyd"], output_time_index, qc_block, (local_z, local_y, local_x), (out_z, out_y, out_x))
                    write_rank_block(out_vars["QR_hyd"], output_time_index, qr_block, (local_z, local_y, local_x), (out_z, out_y, out_x))
            if file_index == 1 or file_index % 8 == 0 or file_index == len(files):
                print(f"[export_3d] processed {file_index}/{len(files)} history files", flush=True)
    finally:
        dataset.close()

    validate_export(output_path)
    metadata = {
        "output": str(output_path),
        "case_dir": str(case_dir),
        "history_file_count": len(files),
        "grid": grid.__dict__,
        "qc_variable": qc_name,
        "qr_variable": qr_name,
        "stride_x": strides[0],
        "stride_y": strides[1],
        "stride_z": strides[2],
        "time_stride": time_stride,
        "time_first_s": float(times[0]) if times.size else None,
        "time_last_s": float(times[-1]) if times.size else None,
        "time_count": int(times[::time_stride].size),
    }
    (outdir / "tpht_bw_qc_qr_3d_fullgrid_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"[export_3d] wrote {output_path}", flush=True)
    return output_path


def main() -> None:
    """CLI entry point."""
    export_case(parse_args())


if __name__ == "__main__":
    main()
