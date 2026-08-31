#!/usr/bin/env python3
"""Extract stable FW/BW super-droplet trajectories from SCALE-SDM snapshots.

The legacy per-case script aligned records by their array column.  That is not
safe after MPI migration, particle loss, or a missing output level.  This
implementation discovers the actual output files and joins records by the
model identity pair ``(dm_id, sd_id)`` (FW) or ``(pre_dmid, pre_sdid)`` (BW).
"""

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

try:
    import numpy as np
    from netCDF4 import Dataset
except ImportError as exc:  # pragma: no cover - exercised by deployment checks
    raise SystemExit(
        "sd_output.py requires numpy and netCDF4. Activate the project Python "
        "environment before running it."
    ) from exc


OUTPUT_RE = re.compile(
    r"^(SD_(selected|all)_NetCDF_)(\d{8})-(\d{6}(?:\.\d+)?)"
    r"\.pe(\d{6})(?:\.nc)?$"
)
FLOAT_ALIASES = {
    "sd_x": ("sd_x", "x", "x_sd", "pos_x"),
    "sd_y": ("sd_y", "y", "y_sd", "pos_y"),
    "sd_z": ("sd_z", "z", "height", "z_sd"),
    "sd_r": ("sd_r", "r", "radius", "rad"),
}
INT_ALIASES = {
    "sd_n": ("sd_n", "n", "multiplicity"),
    "if_coal": ("if_coal", "coal_flag", "coalesced"),
}
FW_IDS = (("dm_id", "sd_id"),)
BW_IDS = (("pre_dmid", "pre_sdid"), ("pre_dm_id", "pre_sd_id"))


def parse_args(argv=None):
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="..", help="directory containing SD NetCDF snapshots")
    parser.add_argument("--output", default="tracking_trajectories.nc", help="output trajectory NetCDF file")
    parser.add_argument("--stream", choices=("auto", "selected", "all"), default="auto")
    parser.add_argument("--mode", choices=("auto", "forward", "backward"), default="auto")
    parser.add_argument(
        "--anchor",
        choices=("auto", "first", "last"),
        default="auto",
        help="time level used to choose targets; auto uses first for FW and last for BW",
    )
    parser.add_argument(
        "--target",
        action="append",
        default=[],
        metavar="DM_ID:SD_ID",
        help="extract an exact identity pair; may be repeated",
    )
    parser.add_argument("--target-file", help="text/CSV file containing dm_id and sd_id pairs")
    parser.add_argument(
        "--max-trajectories",
        type=int,
        default=100,
        help="maximum automatically selected trajectories; 0 keeps every anchor identity",
    )
    parser.add_argument(
        "--selection",
        choices=("largest", "identity"),
        default="largest",
        help="deterministic automatic target ordering at the anchor time",
    )
    parser.add_argument("--strict", action="store_true", help="fail on duplicate identities or weak adjacent overlap")
    return parser.parse_args(argv)


def _clock_seconds(clock):
    """Convert SCALE's HHMMSS.sss token to seconds."""
    main, dot, fraction = clock.partition(".")
    main = main.zfill(6)
    seconds = int(main[0:2]) * 3600 + int(main[2:4]) * 60 + int(main[4:6])
    return float(seconds) + (float("0." + fraction) if dot and fraction else 0.0)


def _date_ordinal(token):
    """Return an ordinal for SCALE's YYYYMMDD token, accepting year 0000."""
    year = max(1, int(token[0:4]))
    return date(year, int(token[4:6]), int(token[6:8])).toordinal()


def parse_output_path(path):
    """Return stream, group key, rank, and clock metadata for one snapshot."""
    match = OUTPUT_RE.match(path.name)
    if match is None:
        return None
    stream = match.group(2)
    date_token = match.group(3)
    clock_token = match.group(4)
    return {
        "stream": stream,
        "key": (date_token, clock_token),
        "rank": int(match.group(5)),
        "clock_s": _clock_seconds(clock_token),
        "ordinal": _date_ordinal(date_token),
    }


def discover_groups(input_dir, requested_stream):
    """Discover real output levels rather than inventing a time sequence."""
    files_by_stream = {"selected": [], "all": []}
    for path in input_dir.iterdir():
        if not path.is_file():
            continue
        metadata = parse_output_path(path)
        if metadata is not None:
            files_by_stream[metadata["stream"]].append((path, metadata))

    if requested_stream == "auto":
        stream = "selected" if files_by_stream["selected"] else "all"
    else:
        stream = requested_stream
    if not files_by_stream[stream]:
        raise RuntimeError("no SD_{0}_NetCDF_* files found in {1}".format(stream, input_dir))

    grouped = defaultdict(list)
    metadata_by_key = {}
    for path, metadata in files_by_stream[stream]:
        grouped[metadata["key"]].append(path)
        metadata_by_key[metadata["key"]] = metadata
    keys = sorted(grouped)
    first_ordinal = metadata_by_key[keys[0]]["ordinal"]
    groups = []
    for time_index, key in enumerate(keys):
        metadata = metadata_by_key[key]
        time_s = metadata["clock_s"] + 86400.0 * (metadata["ordinal"] - first_ordinal)
        groups.append(
            {
                "time_index": time_index,
                "time_s": time_s,
                "label": "{0}-{1}".format(*key),
                "paths": sorted(grouped[key]),
            }
        )
    return stream, groups


def _find_variable(names, aliases, prefixes=()):
    """Find an exact variable alias, then a numbered compatibility alias."""
    name_set = set(names)
    for alias in aliases:
        if alias in name_set:
            return alias
    for prefix in prefixes:
        matches = sorted(name for name in name_set if name.startswith(prefix))
        if matches:
            return matches[0]
    return None


def _find_id_pair(names, mode):
    """Find the correct model identity variables for FW or BW output."""
    preferred = BW_IDS if mode == "backward" else FW_IDS
    for dm_name, sd_name in preferred:
        if dm_name in names and sd_name in names:
            return dm_name, sd_name
    prefixes = (("pre_dmid_", "pre_sdid_"),) if mode == "backward" else (("dm_id_", "sd_id_"),)
    for dm_prefix, sd_prefix in prefixes:
        dm_matches = sorted(name for name in names if name.startswith(dm_prefix))
        sd_matches = sorted(name for name in names if name.startswith(sd_prefix))
        if dm_matches and sd_matches:
            return dm_matches[0], sd_matches[0]
    return None


def detect_mode(path, requested_mode):
    """Detect FW/BW semantics from variable names when requested."""
    if requested_mode != "auto":
        return requested_mode
    with Dataset(path, "r") as handle:
        names = set(handle.variables)
    if _find_id_pair(names, "backward") is not None:
        return "backward"
    if _find_id_pair(names, "forward") is not None:
        return "forward"
    raise RuntimeError("could not find FW dm_id/sd_id or BW pre_dmid/pre_sdid in {0}".format(path))


def _flat(variable, integer=False):
    """Read one variable as a plain one-dimensional NumPy array."""
    values = np.ma.asarray(variable[:]).reshape(-1)
    fill = -999 if integer else np.nan
    values = np.ma.filled(values, fill)
    return np.asarray(values, dtype=np.int64 if integer else np.float64)


def read_rank_records(path, mode, target_pairs=None):
    """Read valid identity records from one rank snapshot."""
    metadata = parse_output_path(path)
    source_rank = metadata["rank"] if metadata is not None else -1
    with Dataset(path, "r") as handle:
        names = set(handle.variables)
        pair_names = _find_id_pair(names, mode)
        if pair_names is None:
            raise RuntimeError("required {0} identity variables are missing in {1}".format(mode, path.name))
        dm_values = _flat(handle.variables[pair_names[0]], integer=True)
        sd_values = _flat(handle.variables[pair_names[1]], integer=True)
        length = min(len(dm_values), len(sd_values))
        dm_values = dm_values[:length]
        sd_values = sd_values[:length]
        valid = (dm_values >= 0) & (sd_values > 0)
        if target_pairs is not None:
            target_mask = np.fromiter(
                ((int(dm_values[i]), int(sd_values[i])) in target_pairs for i in range(length)),
                dtype=bool,
                count=length,
            )
            valid &= target_mask
        indices = np.nonzero(valid)[0]

        arrays = {}
        for output_name, aliases in FLOAT_ALIASES.items():
            variable_name = _find_variable(names, aliases)
            arrays[output_name] = _flat(handle.variables[variable_name])[:length] if variable_name else None
        for output_name, aliases in INT_ALIASES.items():
            prefixes = ("if_coal_",) if output_name == "if_coal" else ()
            variable_name = _find_variable(names, aliases, prefixes)
            arrays[output_name] = _flat(handle.variables[variable_name], integer=True)[:length] if variable_name else None

    records = []
    for index in indices:
        z_value = arrays["sd_z"][index] if arrays["sd_z"] is not None else np.nan
        if math.isfinite(float(z_value)) and float(z_value) < 0.0:
            continue
        record = {
            "pair": (int(dm_values[index]), int(sd_values[index])),
            "source_rank": source_rank,
            "source_index": int(index),
        }
        for name in FLOAT_ALIASES:
            record[name] = float(arrays[name][index]) if arrays[name] is not None else math.nan
        record["sd_n"] = int(arrays["sd_n"][index]) if arrays["sd_n"] is not None else -1
        record["if_coal"] = int(arrays["if_coal"][index]) if arrays["if_coal"] is not None else -1
        records.append(record)
    return records


def _parse_pair(text):
    """Parse one DM:SD or two-column identity token."""
    fields = re.split(r"[:,\s]+", text.strip())
    if len(fields) < 2:
        raise ValueError("invalid identity pair: {0}".format(text))
    dm_id, sd_id = int(fields[0]), int(fields[1])
    if dm_id < 0 or sd_id <= 0:
        raise ValueError("identity must satisfy dm_id >= 0 and sd_id > 0: {0}".format(text))
    return dm_id, sd_id


def read_target_file(path):
    """Read identity pairs from a TPHT .ids file or a CSV table."""
    target_path = Path(path)
    if not target_path.is_file():
        raise RuntimeError("target file does not exist: {0}".format(target_path))
    pairs = set()
    if target_path.suffix.lower() == ".csv":
        with target_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "dm_id" not in reader.fieldnames or "sd_id" not in reader.fieldnames:
                raise RuntimeError("target CSV requires dm_id and sd_id columns")
            for row in reader:
                pairs.add(_parse_pair("{0}:{1}".format(row["dm_id"], row["sd_id"])))
    else:
        for raw_line in target_path.read_text(errors="replace").splitlines():
            line = raw_line.strip()
            if line and not line.startswith("#"):
                pairs.add(_parse_pair(line))
    return pairs


def requested_targets(args):
    """Return exact targets requested on the command line or in a file."""
    explicit = set(_parse_pair(value) for value in args.target)
    if args.target_file:
        explicit.update(read_target_file(args.target_file))
    return explicit


def ranked_anchor_records(anchor_records, selection):
    """Deduplicate and deterministically rank anchor-time candidates."""
    anchor_by_pair = {}
    for record in anchor_records:
        anchor_by_pair.setdefault(record["pair"], record)
    records = list(anchor_by_pair.values())
    if selection == "largest":
        records.sort(
            key=lambda item: (
                -item["sd_r"] if math.isfinite(item["sd_r"]) else math.inf,
                item["pair"],
            )
        )
    else:
        records.sort(key=lambda item: item["pair"])
    return records


def choose_targets(anchor_records, args, explicit):
    """Choose a deterministic target set from an anchor output level."""
    anchor_by_pair = {record["pair"]: record for record in anchor_records}
    if explicit:
        missing = sorted(explicit.difference(anchor_by_pair))
        if missing:
            raise RuntimeError("explicit targets absent from anchor output: {0}".format(missing[:10]))
        return explicit

    records = ranked_anchor_records(anchor_records, args.selection)
    if args.max_trajectories > 0:
        records = records[: args.max_trajectories]
    return set(record["pair"] for record in records)


def nearest_run_conf(input_dir):
    """Find the closest run.conf for provenance and x-period inference."""
    for directory in (input_dir,) + tuple(input_dir.parents):
        candidate = directory / "run.conf"
        if candidate.is_file():
            return candidate
    return None


def read_namelist_number(path, key):
    """Read one simple numeric namelist assignment without evaluating Fortran."""
    if path is None:
        return None
    pattern = re.compile(r"^\s*{0}\s*=\s*([^,!]+)".format(re.escape(key)), re.IGNORECASE)
    for raw_line in path.read_text(errors="replace").splitlines():
        match = pattern.match(raw_line)
        if match:
            token = re.sub(r"[dD]", "e", match.group(1).strip())
            try:
                return float(token)
            except ValueError:
                return None
    return None


def infer_x_period(run_conf):
    """Infer periodic x-domain length from DX and IMAXG when available."""
    dx = read_namelist_number(run_conf, "DX")
    imaxg = read_namelist_number(run_conf, "IMAXG")
    if dx is None or imaxg is None or dx <= 0.0 or imaxg <= 0.0:
        return math.nan
    return dx * imaxg


def write_output(output_path, mode, stream, groups, target_pairs, rows, warnings, run_conf):
    """Write a tidy trajectory NetCDF plus compact CSV/JSON summaries."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs = sorted(target_pairs)
    pair_to_index = {pair: index for index, pair in enumerate(pairs)}
    rows.sort(key=lambda row: (row["time_index"], pair_to_index[row["pair"]]))
    rows_by_pair = defaultdict(list)
    for row in rows:
        rows_by_pair[row["pair"]].append(row)

    with Dataset(output_path, "w", format="NETCDF4") as handle:
        handle.createDimension("record", len(rows))
        handle.createDimension("trajectory", len(pairs))
        handle.createDimension("time", len(groups))
        handle.setncattr("tracking_mode", mode)
        handle.setncattr("source_stream", stream)
        handle.setncattr("identity_semantics", "dm_id:sd_id" if mode == "forward" else "pre_dmid:pre_sdid")
        handle.setncattr("join_method", "stable identity pair across actual output levels")
        handle.setncattr("source_run_conf", str(run_conf) if run_conf else "")
        handle.setncattr("x_period_m", infer_x_period(run_conf))
        handle.setncattr("warning_count", len(warnings))

        time_var = handle.createVariable("time", "f8", ("time",))
        time_var.units = "seconds since SCALE model clock origin"
        time_var[:] = np.asarray([group["time_s"] for group in groups], dtype=np.float64)

        def int_var(name, values, dimensions=("record",), fill=-999):
            variable = handle.createVariable(name, "i8", dimensions, fill_value=fill)
            variable[:] = np.asarray(values, dtype=np.int64)
            return variable

        def float_var(name, values, dimensions=("record",)):
            variable = handle.createVariable(name, "f8", dimensions, fill_value=np.nan)
            variable[:] = np.asarray(values, dtype=np.float64)
            return variable

        int_var("trajectory_index", [pair_to_index[row["pair"]] for row in rows])
        int_var("particle_index", [pair_to_index[row["pair"]] for row in rows])
        int_var("time_index", [row["time_index"] for row in rows])
        float_var("time_s", [row["time_s"] for row in rows])
        int_var("dm_id", [row["pair"][0] for row in rows])
        int_var("sd_id", [row["pair"][1] for row in rows])
        int_var("source_rank", [row["source_rank"] for row in rows])
        int_var("source_index", [row["source_index"] for row in rows])
        for name in FLOAT_ALIASES:
            variable = float_var(name, [row[name] for row in rows])
            if name == "sd_r":
                variable.units = "m"
            elif name in ("sd_x", "sd_y", "sd_z"):
                variable.units = "m"
        int_var("sd_n", [row["sd_n"] for row in rows])
        int_var("if_coal", [row["if_coal"] for row in rows], fill=-999)

        int_var("trajectory_dm_id", [pair[0] for pair in pairs], ("trajectory",))
        int_var("trajectory_sd_id", [pair[1] for pair in pairs], ("trajectory",))
        int_var("trajectory_record_count", [len(rows_by_pair[pair]) for pair in pairs], ("trajectory",), fill=0)
        float_var(
            "trajectory_max_radius_m",
            [
                max((row["sd_r"] for row in rows_by_pair[pair] if math.isfinite(row["sd_r"])), default=math.nan)
                for pair in pairs
            ],
            ("trajectory",),
        )

    summary_rows = []
    for index, pair in enumerate(pairs):
        pair_rows = rows_by_pair[pair]
        finite_radii = [row["sd_r"] for row in pair_rows if math.isfinite(row["sd_r"])]
        summary_rows.append(
            {
                "trajectory_index": index,
                "dm_id": pair[0],
                "sd_id": pair[1],
                "record_count": len(pair_rows),
                "first_time_s": min((row["time_s"] for row in pair_rows), default=None),
                "last_time_s": max((row["time_s"] for row in pair_rows), default=None),
                "max_radius_um": max(finite_radii) * 1.0e6 if finite_radii else None,
            }
        )
    csv_path = output_path.with_name(output_path.stem + "_summary.csv")
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]) if summary_rows else ["trajectory_index", "dm_id", "sd_id"])
        writer.writeheader()
        writer.writerows(summary_rows)
    json_path = output_path.with_name(output_path.stem + "_summary.json")
    json_path.write_text(
        json.dumps(
            {
                "tracking_mode": mode,
                "source_stream": stream,
                "time_level_count": len(groups),
                "target_count": len(pairs),
                "record_count": len(rows),
                "warnings": warnings,
                "trajectories": summary_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return csv_path, json_path


def main(argv=None):
    """Run trajectory extraction."""
    args = parse_args(argv)
    if args.max_trajectories < 0:
        raise SystemExit("--max-trajectories must be non-negative")
    input_dir = Path(args.input_dir).resolve()
    output_path = Path(args.output).resolve()
    try:
        stream, groups = discover_groups(input_dir, args.stream)
        mode = detect_mode(groups[0]["paths"][0], args.mode)
        anchor_name = args.anchor
        if anchor_name == "auto":
            anchor_name = "last" if mode == "backward" else "first"
        anchor_group = groups[-1] if anchor_name == "last" else groups[0]
        explicit_targets = requested_targets(args)
        anchor_records = []
        for path in anchor_group["paths"]:
            anchor_records.extend(
                read_rank_records(path, mode, explicit_targets if explicit_targets else None)
            )
            if not explicit_targets and args.max_trajectories > 0:
                anchor_records = ranked_anchor_records(anchor_records, args.selection)[
                    : args.max_trajectories
                ]
        target_pairs = choose_targets(anchor_records, args, explicit_targets)
        if not target_pairs:
            raise RuntimeError("no valid trajectory targets were selected")

        rows = []
        warnings = []
        previous_pairs = None
        duplicate_total = 0
        for group in groups:
            level_records = []
            for path in group["paths"]:
                level_records.extend(read_rank_records(path, mode, target_pairs))
            unique = {}
            for record in level_records:
                if record["pair"] in unique:
                    duplicate_total += 1
                    continue
                unique[record["pair"]] = record
            current_pairs = set(unique)
            if previous_pairs is not None and current_pairs:
                overlap = len(previous_pairs.intersection(current_pairs)) / float(len(current_pairs))
                if overlap < 0.5:
                    warnings.append(
                        "adjacent identity overlap is {0:.3f} at {1}; verify output cadence and ID semantics".format(
                            overlap, group["label"]
                        )
                    )
            for record in unique.values():
                record["time_index"] = group["time_index"]
                record["time_s"] = group["time_s"]
                rows.append(record)
            previous_pairs = current_pairs
        if duplicate_total:
            warnings.append("ignored {0} duplicate identity records within output levels".format(duplicate_total))
        if args.strict and warnings:
            raise RuntimeError("strict trajectory checks failed: " + "; ".join(warnings))

        run_conf = nearest_run_conf(input_dir)
        csv_path, json_path = write_output(
            output_path, mode, stream, groups, target_pairs, rows, warnings, run_conf
        )
    except (OSError, RuntimeError, ValueError) as exc:
        raise SystemExit("ERROR: {0}".format(exc))

    print("tracking mode: {0}".format(mode))
    print("source stream: SD_{0}_NetCDF_*".format(stream))
    print("output levels: {0}".format(len(groups)))
    print("target trajectories: {0}".format(len(target_pairs)))
    print("trajectory records: {0}".format(len(rows)))
    print("wrote: {0}".format(output_path))
    print("wrote: {0}".format(csv_path))
    print("wrote: {0}".format(json_path))
    for warning in warnings:
        print("WARNING: {0}".format(warning), file=sys.stderr)


if __name__ == "__main__":
    main()
