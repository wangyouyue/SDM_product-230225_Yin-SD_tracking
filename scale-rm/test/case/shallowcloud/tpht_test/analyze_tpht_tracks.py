#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
from collections import defaultdict
from netCDF4 import Dataset


FILE_RE = re.compile(r"^(?P<prefix>.+?)_(?P<stamp>\d{8}-\d{6}\.\d{3})\.pe(?P<rank>\d{6})(?:\..*)?$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Post-process TPHT BW outputs using (pre_dmid, pre_sdid) as the stable "
            "trajectory key. The script summarizes trajectory evolution and, when "
            "available, coalescence-event history from SD_coal_output_NetCDF_*."
        )
    )
    parser.add_argument(
        "--bw-glob",
        required=True,
        help="Glob for BW SD_selected_NetCDF outputs.",
    )
    parser.add_argument(
        "--coal-glob",
        default="",
        help="Optional glob for SD_coal_output_NetCDF outputs.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where CSV/JSON analysis outputs will be written.",
    )
    parser.add_argument(
        "--sample-count",
        type=int,
        default=20,
        help="Number of representative TPHT particle pairs to export in trajectory_samples.csv.",
    )
    return parser.parse_args()


def parse_stamp_rank(path: str):
    base = os.path.basename(path)
    match = FILE_RE.match(base)
    if not match:
        raise SystemExit(f"Could not parse timestamp/rank from path: {path}")
    return match.group("stamp"), int(match.group("rank"))


def group_paths_by_stamp(paths: list[str]):
    grouped = defaultdict(dict)
    for path in sorted(paths):
        stamp, rank = parse_stamp_rank(path)
        if rank in grouped[stamp]:
            raise SystemExit(f"Duplicate rank {rank} for timestamp {stamp}: {path}")
        grouped[stamp][rank] = path
    ordered_stamps = sorted(grouped)
    return ordered_stamps, grouped


def read_tracking_var(nc_handle, names, prefixes=None):
    variable_names = set(nc_handle.variables.keys())
    for name in names:
        if name in variable_names:
            return nc_handle.variables[name][:]
    if prefixes is not None:
        all_names = sorted(variable_names)
        for prefix in prefixes:
            matched = [name for name in all_names if name.startswith(prefix)]
            if matched:
                return nc_handle.variables[matched[-1]][:]
    raise SystemExit(
        "Failed to find any of the expected variables "
        f"{names} / prefixes {prefixes} in file variables {sorted(variable_names)[:20]}"
    )


def read_bw_group(stamp: str, rank_to_path: dict[int, str]):
    records = {}
    ranks = sorted(rank_to_path)

    tracked_count = 0
    nonlocal_count = 0
    if_coal_count = 0
    sd_r_values = []
    sd_z_values = []
    sd_x_values = []
    sd_y_values = []

    for rank in ranks:
        path = rank_to_path[rank]
        with Dataset(path, "r") as nc_handle:
            sd_x = nc_handle.variables["sd_x"][:]
            sd_y = nc_handle.variables["sd_y"][:]
            sd_z = nc_handle.variables["sd_z"][:]
            sd_r = nc_handle.variables["sd_r"][:]
            sd_n = nc_handle.variables["sd_n"][:]
            if_coal = read_tracking_var(nc_handle, ["if_coal"], ["if_coal_"])
            pre_sdid = read_tracking_var(nc_handle, ["pre_sdid", "sd_id"], ["pre_sdid_", "sd_id_"])
            pre_dmid = read_tracking_var(nc_handle, ["pre_dmid", "dm_id"], ["pre_dmid_", "dm_id_"])

            for idx in range(len(pre_sdid)):
                pair = (int(pre_dmid[idx]), int(pre_sdid[idx]))
                if pair[0] < 0 or pair[1] < 0:
                    continue
                if pair in records:
                    raise SystemExit(
                        f"Duplicate TPHT pair {pair} at timestamp {stamp}; "
                        f"previous source={records[pair]['source_path']}, duplicate source={path}"
                    )

                record = {
                    "pair_dm": pair[0],
                    "pair_sd": pair[1],
                    "file_rank": rank,
                    "sd_x": float(sd_x[idx]),
                    "sd_y": float(sd_y[idx]),
                    "sd_z": float(sd_z[idx]),
                    "sd_r": float(sd_r[idx]),
                    "sd_n": int(sd_n[idx]),
                    "if_coal": int(if_coal[idx]),
                    "source_path": path,
                }
                records[pair] = record

                tracked_count += 1
                if pair[0] != rank:
                    nonlocal_count += 1
                if record["if_coal"] > 0:
                    if_coal_count += 1

                sd_r_values.append(record["sd_r"])
                sd_z_values.append(record["sd_z"])
                sd_x_values.append(record["sd_x"])
                sd_y_values.append(record["sd_y"])

    if tracked_count == 0:
        raise SystemExit(f"No valid TPHT records found at timestamp {stamp}")

    summary = {
        "time_label": stamp,
        "tracked_count": tracked_count,
        "unique_pairs": len(records),
        "nonlocal_pairs": nonlocal_count,
        "if_coal_count": if_coal_count,
        "sd_r_mean": sum(sd_r_values) / len(sd_r_values),
        "sd_r_min": min(sd_r_values),
        "sd_r_max": max(sd_r_values),
        "sd_z_min": min(sd_z_values),
        "sd_z_max": max(sd_z_values),
        "sd_x_min": min(sd_x_values),
        "sd_x_max": max(sd_x_values),
        "sd_y_min": min(sd_y_values),
        "sd_y_max": max(sd_y_values),
    }
    return records, summary


def read_coalescence_group(stamp: str, rank_to_path: dict[int, str]):
    pair_stats = defaultdict(
        lambda: {
            "event_count": 0,
            "num_col_sum": 0,
            "partners": set(),
            "first_collision_time": None,
            "last_collision_time": None,
        }
    )
    event_rows = []

    for path in (rank_to_path[rank] for rank in sorted(rank_to_path)):
        _, file_rank = parse_stamp_rank(path)
        with Dataset(path, "r") as nc_handle:
            pre_sdid1 = read_tracking_var(nc_handle, ["pre_sdid1", "sd_id1"], ["pre_sdid1_", "sd_id1_"])
            pre_dmid1 = read_tracking_var(nc_handle, ["pre_dmid1", "dm_id1"], ["pre_dmid1_", "dm_id1_"])
            pre_sdid2 = read_tracking_var(nc_handle, ["pre_sdid2", "sd_id2"], ["pre_sdid2_", "sd_id2_"])
            pre_dmid2 = read_tracking_var(nc_handle, ["pre_dmid2", "dm_id2"], ["pre_dmid2_", "dm_id2_"])
            num_col = nc_handle.variables["num_col"][:]

            for idx in range(len(num_col)):
                pair1 = (int(pre_dmid1[idx]), int(pre_sdid1[idx]))
                pair2 = (int(pre_dmid2[idx]), int(pre_sdid2[idx]))
                ncol = int(num_col[idx])

                event_rows.append(
                    {
                        "time_label": stamp,
                        "file_rank": file_rank,
                        "pre_dmid1": pair1[0],
                        "pre_sdid1": pair1[1],
                        "pre_dmid2": pair2[0],
                        "pre_sdid2": pair2[1],
                        "num_col": ncol,
                    }
                )

                if pair1[0] >= 0 and pair1[1] >= 0:
                    pair_stats[pair1]["event_count"] += 1
                    pair_stats[pair1]["num_col_sum"] += ncol
                    if pair_stats[pair1]["first_collision_time"] is None:
                        pair_stats[pair1]["first_collision_time"] = stamp
                    pair_stats[pair1]["last_collision_time"] = stamp
                    if pair2[0] >= 0 and pair2[1] >= 0:
                        pair_stats[pair1]["partners"].add(pair2)

                if pair2[0] >= 0 and pair2[1] >= 0:
                    pair_stats[pair2]["event_count"] += 1
                    pair_stats[pair2]["num_col_sum"] += ncol
                    if pair_stats[pair2]["first_collision_time"] is None:
                        pair_stats[pair2]["first_collision_time"] = stamp
                    pair_stats[pair2]["last_collision_time"] = stamp
                    if pair1[0] >= 0 and pair1[1] >= 0:
                        pair_stats[pair2]["partners"].add(pair1)

    return pair_stats, event_rows


def stamp_to_seconds(stamp: str, start_stamp: str):
    def parse_time_part(label: str):
        try:
            time_part = label.split("-", 1)[1]
            hour = int(time_part[0:2])
            minute = int(time_part[2:4])
            second = int(time_part[4:6])
            millisecond = int(time_part[7:10])
        except Exception as exc:
            raise SystemExit(f"Failed to parse TPHT timestamp '{label}': {exc}")
        return hour * 3600.0 + minute * 60.0 + second + millisecond / 1000.0

    return parse_time_part(stamp) - parse_time_part(start_stamp)


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def main():
    args = parse_args()
    ensure_dir(args.output_dir)

    bw_paths = sorted(glob.glob(args.bw_glob))
    if not bw_paths:
        raise SystemExit(f"No BW outputs matched: {args.bw_glob}")

    bw_stamps, bw_groups = group_paths_by_stamp(bw_paths)

    reference_pairs = None
    per_time_records = {}
    time_summaries = []
    pair_overview = defaultdict(
        lambda: {
            "n_times_present": 0,
            "first_time": None,
            "last_time": None,
            "max_radius": -1.0,
            "if_coal_time_count": 0,
        }
    )

    for stamp in bw_stamps:
        records, summary = read_bw_group(stamp, bw_groups[stamp])
        current_pairs = set(records)
        if reference_pairs is None:
            reference_pairs = current_pairs
        missing_vs_ref = len(reference_pairs - current_pairs)
        extra_vs_ref = len(current_pairs - reference_pairs)
        summary["missing_vs_first_group"] = missing_vs_ref
        summary["extra_vs_first_group"] = extra_vs_ref
        summary["time_seconds_from_start"] = stamp_to_seconds(stamp, bw_stamps[0])

        per_time_records[stamp] = records
        time_summaries.append(summary)

        for pair, record in records.items():
            overview = pair_overview[pair]
            overview["n_times_present"] += 1
            overview["first_time"] = stamp if overview["first_time"] is None else overview["first_time"]
            overview["last_time"] = stamp
            overview["max_radius"] = max(overview["max_radius"], record["sd_r"])
            if record["if_coal"] > 0:
                overview["if_coal_time_count"] += 1

    all_pairs = sorted(reference_pairs if reference_pairs is not None else [])
    all_pair_set = set(all_pairs)

    coal_paths = sorted(glob.glob(args.coal_glob)) if args.coal_glob else []
    coal_pair_stats = defaultdict(lambda: {"event_count": 0, "num_col_sum": 0, "partners": set()})
    collision_events = []
    coal_stamps = []
    if coal_paths:
        coal_stamps, coal_groups = group_paths_by_stamp(coal_paths)
        for stamp in coal_stamps:
            current_stats, current_events = read_coalescence_group(stamp, coal_groups[stamp])
            collision_events.extend(current_events)
            for pair, stats in current_stats.items():
                coal_pair_stats[pair]["event_count"] += stats["event_count"]
                coal_pair_stats[pair]["num_col_sum"] += stats["num_col_sum"]
                coal_pair_stats[pair]["partners"].update(stats["partners"])
                if coal_pair_stats[pair].get("first_collision_time") is None:
                    coal_pair_stats[pair]["first_collision_time"] = stats["first_collision_time"]
                coal_pair_stats[pair]["last_collision_time"] = stats["last_collision_time"]

    collision_time_summary = []
    if collision_events:
        grouped_events = defaultdict(list)
        for row in collision_events:
            pair1 = (row["pre_dmid1"], row["pre_sdid1"])
            pair2 = (row["pre_dmid2"], row["pre_sdid2"])
            row["pair1_tracked"] = pair1 in all_pair_set
            row["pair2_tracked"] = pair2 in all_pair_set
            row["tracked_any"] = row["pair1_tracked"] or row["pair2_tracked"]
            row["tracked_both"] = row["pair1_tracked"] and row["pair2_tracked"]
            grouped_events[row["time_label"]].append(row)

        for stamp in sorted(grouped_events):
            rows = grouped_events[stamp]
            tracked_any_pairs = set()
            tracked_both_pairs = set()
            num_col_sum_all = 0
            num_col_sum_tracked_any = 0
            num_col_sum_tracked_both = 0

            for row in rows:
                num_col_sum_all += row["num_col"]
                pair1 = (row["pre_dmid1"], row["pre_sdid1"])
                pair2 = (row["pre_dmid2"], row["pre_sdid2"])
                if row["tracked_any"]:
                    num_col_sum_tracked_any += row["num_col"]
                    if row["pair1_tracked"]:
                        tracked_any_pairs.add(pair1)
                    if row["pair2_tracked"]:
                        tracked_any_pairs.add(pair2)
                if row["tracked_both"]:
                    num_col_sum_tracked_both += row["num_col"]
                    tracked_both_pairs.add(pair1)
                    tracked_both_pairs.add(pair2)

            collision_time_summary.append(
                {
                    "time_label": stamp,
                    "time_seconds_from_start": stamp_to_seconds(stamp, bw_stamps[0]),
                    "event_count_all": len(rows),
                    "event_count_tracked_any": sum(1 for row in rows if row["tracked_any"]),
                    "event_count_tracked_both": sum(1 for row in rows if row["tracked_both"]),
                    "num_col_sum_all": num_col_sum_all,
                    "num_col_sum_tracked_any": num_col_sum_tracked_any,
                    "num_col_sum_tracked_both": num_col_sum_tracked_both,
                    "tracked_pair_count_any": len(tracked_any_pairs),
                    "tracked_pair_count_both": len(tracked_both_pairs),
                }
            )

    full_timeseries_csv = os.path.join(args.output_dir, "tpht_pair_timeseries.csv")
    with open(full_timeseries_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_label",
                "time_seconds_from_start",
                "pre_dmid",
                "pre_sdid",
                "file_rank",
                "sd_x",
                "sd_y",
                "sd_z",
                "sd_r",
                "sd_n",
                "if_coal",
            ]
        )
        for summary in time_summaries:
            stamp = summary["time_label"]
            records = per_time_records[stamp]
            for pair in all_pairs:
                if pair not in records:
                    continue
                record = records[pair]
                writer.writerow(
                    [
                        stamp,
                        summary["time_seconds_from_start"],
                        pair[0],
                        pair[1],
                        record["file_rank"],
                        record["sd_x"],
                        record["sd_y"],
                        record["sd_z"],
                        record["sd_r"],
                        record["sd_n"],
                        record["if_coal"],
                    ]
                )

    pair_overview_rows = []
    for pair in all_pairs:
        overview = pair_overview[pair]
        coal_stats = coal_pair_stats[pair]
        pair_overview_rows.append(
            {
                "pre_dmid": pair[0],
                "pre_sdid": pair[1],
                "n_times_present": overview["n_times_present"],
                "first_time": overview["first_time"],
                "last_time": overview["last_time"],
                "first_sd_r": per_time_records[overview["first_time"]][pair]["sd_r"],
                "last_sd_r": per_time_records[overview["last_time"]][pair]["sd_r"],
                "delta_sd_r": per_time_records[overview["last_time"]][pair]["sd_r"]
                - per_time_records[overview["first_time"]][pair]["sd_r"],
                "first_sd_z": per_time_records[overview["first_time"]][pair]["sd_z"],
                "last_sd_z": per_time_records[overview["last_time"]][pair]["sd_z"],
                "delta_sd_z": per_time_records[overview["last_time"]][pair]["sd_z"]
                - per_time_records[overview["first_time"]][pair]["sd_z"],
                "first_file_rank": per_time_records[overview["first_time"]][pair]["file_rank"],
                "last_file_rank": per_time_records[overview["last_time"]][pair]["file_rank"],
                "max_radius": overview["max_radius"],
                "if_coal_time_count": overview["if_coal_time_count"],
                "event_count": coal_stats["event_count"],
                "num_col_sum": coal_stats["num_col_sum"],
                "partner_pair_count": len(coal_stats["partners"]),
                "first_collision_time": coal_stats.get("first_collision_time"),
                "last_collision_time": coal_stats.get("last_collision_time"),
            }
        )

    overview_csv = os.path.join(args.output_dir, "tpht_pair_overview.csv")
    with open(overview_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pre_dmid",
                "pre_sdid",
                "n_times_present",
                "first_time",
                "last_time",
                "first_sd_r",
                "last_sd_r",
                "delta_sd_r",
                "first_sd_z",
                "last_sd_z",
                "delta_sd_z",
                "first_file_rank",
                "last_file_rank",
                "max_radius",
                "if_coal_time_count",
                "event_count",
                "num_col_sum",
                "partner_pair_count",
                "first_collision_time",
                "last_collision_time",
            ],
        )
        writer.writeheader()
        writer.writerows(pair_overview_rows)

    collision_rows = [
        row
        for row in pair_overview_rows
        if row["if_coal_time_count"] > 0 or row["event_count"] > 0 or row["num_col_sum"] > 0
    ]
    collision_csv = os.path.join(args.output_dir, "tpht_collision_summary.csv")
    with open(collision_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pre_dmid",
                "pre_sdid",
                "if_coal_time_count",
                "event_count",
                "num_col_sum",
                "partner_pair_count",
            ],
        )
        writer.writeheader()
        for row in collision_rows:
            writer.writerow(
                {
                    "pre_dmid": row["pre_dmid"],
                    "pre_sdid": row["pre_sdid"],
                    "if_coal_time_count": row["if_coal_time_count"],
                    "event_count": row["event_count"],
                    "num_col_sum": row["num_col_sum"],
                    "partner_pair_count": row["partner_pair_count"],
                }
            )

    collision_events_csv = os.path.join(args.output_dir, "tpht_collision_events.csv")
    with open(collision_events_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "time_label",
                "file_rank",
                "pre_dmid1",
                "pre_sdid1",
                "pre_dmid2",
                "pre_sdid2",
                "num_col",
                "pair1_tracked",
                "pair2_tracked",
                "tracked_any",
                "tracked_both",
            ],
        )
        writer.writeheader()
        writer.writerows(collision_events)

    collision_time_csv = os.path.join(args.output_dir, "tpht_collision_time_summary.csv")
    with open(collision_time_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "time_label",
                "time_seconds_from_start",
                "event_count_all",
                "event_count_tracked_any",
                "event_count_tracked_both",
                "num_col_sum_all",
                "num_col_sum_tracked_any",
                "num_col_sum_tracked_both",
                "tracked_pair_count_any",
                "tracked_pair_count_both",
            ],
        )
        writer.writeheader()
        writer.writerows(collision_time_summary)

    sample_rows = sorted(
        pair_overview_rows,
        key=lambda row: (-row["max_radius"], row["pre_dmid"], row["pre_sdid"]),
    )[: args.sample_count]
    sample_pairs = {(row["pre_dmid"], row["pre_sdid"]) for row in sample_rows}

    sample_csv = os.path.join(args.output_dir, "tpht_trajectory_samples.csv")
    with open(sample_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_label",
                "time_seconds_from_start",
                "pre_dmid",
                "pre_sdid",
                "sd_x",
                "sd_y",
                "sd_z",
                "sd_r",
                "sd_n",
                "if_coal",
            ]
        )
        for summary in time_summaries:
            stamp = summary["time_label"]
            records = per_time_records[stamp]
            for pair in sorted(sample_pairs):
                if pair not in records:
                    continue
                record = records[pair]
                writer.writerow(
                    [
                        stamp,
                        summary["time_seconds_from_start"],
                        pair[0],
                        pair[1],
                        record["sd_x"],
                        record["sd_y"],
                        record["sd_z"],
                        record["sd_r"],
                        record["sd_n"],
                        record["if_coal"],
                    ]
                )

    summary_json = os.path.join(args.output_dir, "tpht_analysis_summary.json")
    summary_payload = {
        "bw_glob": args.bw_glob,
        "coal_glob": args.coal_glob,
        "output_dir": args.output_dir,
        "time_labels": bw_stamps,
        "time_count": len(bw_stamps),
        "reference_pair_count": len(all_pairs),
        "pair_set_stable_across_times": all(
            summary["missing_vs_first_group"] == 0 and summary["extra_vs_first_group"] == 0
            for summary in time_summaries
        ),
        "coal_file_count": len(coal_paths),
        "collision_event_count": len(collision_events),
        "collision_pair_count": len(collision_rows),
        "time_summaries": time_summaries,
        "collision_time_summaries": collision_time_summary,
        "sample_pairs": sample_rows,
    }
    with open(summary_json, "w") as handle:
        json.dump(summary_payload, handle, indent=2)

    print(f"bw_time_count={len(bw_stamps)}")
    print(f"reference_pair_count={len(all_pairs)}")
    print(
        "pair_set_stable_across_times="
        f"{summary_payload['pair_set_stable_across_times']}"
    )
    print(f"coal_file_count={len(coal_paths)}")
    print(f"collision_event_count={len(collision_events)}")
    print(f"collision_pair_count={len(collision_rows)}")
    print(f"time_summary_csv={os.path.join(args.output_dir, 'tpht_time_summary.csv')}")
    print(f"pair_timeseries_csv={full_timeseries_csv}")
    print(f"pair_overview_csv={overview_csv}")
    print(f"collision_summary_csv={collision_csv}")
    print(f"collision_events_csv={collision_events_csv}")
    print(f"collision_time_summary_csv={collision_time_csv}")
    print(f"trajectory_samples_csv={sample_csv}")
    print(f"summary_json={summary_json}")

    time_summary_csv = os.path.join(args.output_dir, "tpht_time_summary.csv")
    with open(time_summary_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "time_label",
                "time_seconds_from_start",
                "tracked_count",
                "unique_pairs",
                "nonlocal_pairs",
                "if_coal_count",
                "missing_vs_first_group",
                "extra_vs_first_group",
                "sd_r_mean",
                "sd_r_min",
                "sd_r_max",
                "sd_z_min",
                "sd_z_max",
                "sd_x_min",
                "sd_x_max",
                "sd_y_min",
                "sd_y_max",
            ],
        )
        writer.writeheader()
        writer.writerows(time_summaries)


if __name__ == "__main__":
    main()
