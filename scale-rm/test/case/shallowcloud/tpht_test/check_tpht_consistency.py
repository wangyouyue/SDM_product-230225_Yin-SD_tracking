#!/usr/bin/env python3
import argparse
import glob
import os
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fw-ids",
        default="./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids",
    )
    parser.add_argument(
        "--bw-glob",
        default="./bt_interest_id_baseline/bw_output/SD_selected_NetCDF_*.pe*.nc",
    )
    parser.add_argument(
        "--max-report",
        type=int,
        default=20,
    )
    return parser.parse_args()


def read_fw_ids(file_path):
    pairs = set()
    with open(file_path, "r") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.split()
            if len(parts) < 2:
                continue
            pairs.add((int(parts[0]), int(parts[1])))
    return pairs


def choose_first_bw_group(file_paths):
    grouped = {}
    for file_path in sorted(file_paths):
        match = re.match(r"^(.*)\.pe\d{6}\.nc$", file_path)
        group_key = match.group(1) if match else file_path
        grouped.setdefault(group_key, []).append(file_path)
    first_key = sorted(grouped.keys())[0]
    return first_key, sorted(grouped[first_key])


def find_tracking_variables(nc_handle):
    variable_names = set(nc_handle.variables.keys())
    if "pre_dmid" in variable_names and "pre_sdid" in variable_names:
        return "pre_dmid", "pre_sdid"

    pre_dm = sorted(
        name for name in variable_names if re.match(r"^pre_dmid_\d+$", name)
    )
    pre_sd = sorted(
        name for name in variable_names if re.match(r"^pre_sdid_\d+$", name)
    )
    if pre_dm and pre_sd:
        return pre_dm[0], pre_sd[0]

    if "dm_id" in variable_names and "sd_id" in variable_names:
        return "dm_id", "sd_id"

    dm_time = sorted(
        name for name in variable_names if re.match(r"^dm_id_\d+$", name)
    )
    sd_time = sorted(
        name for name in variable_names if re.match(r"^sd_id_\d+$", name)
    )
    if dm_time and sd_time:
        return dm_time[0], sd_time[0]

    raise RuntimeError("No BW tracking ID variables found in NetCDF file")


def read_bw_pairs(file_paths):
    from netCDF4 import Dataset

    pairs = set()
    total_valid_records = 0
    selected_dm_var = None
    selected_sd_var = None

    for file_path in file_paths:
        with Dataset(file_path, "r") as nc_handle:
            dm_var, sd_var = find_tracking_variables(nc_handle)
            selected_dm_var = dm_var
            selected_sd_var = sd_var
            dm_values = nc_handle.variables[dm_var][:]
            sd_values = nc_handle.variables[sd_var][:]

            for dm_id, sd_id in zip(dm_values, sd_values):
                dm_id = int(dm_id)
                sd_id = int(sd_id)
                if dm_id < 0 or sd_id < 0:
                    continue
                total_valid_records += 1
                pairs.add((dm_id, sd_id))

    return pairs, total_valid_records, selected_dm_var, selected_sd_var


def print_sample(title, items, max_report):
    print(title)
    for pair in list(sorted(items))[:max_report]:
        print("  {0} {1}".format(pair[0], pair[1]))


def main():
    args = parse_args()

    if not os.path.exists(args.fw_ids):
        raise SystemExit("FW ID file not found: {0}".format(args.fw_ids))

    bw_candidates = sorted(glob.glob(args.bw_glob))
    if not bw_candidates:
        raise SystemExit("No BW NetCDF files matched: {0}".format(args.bw_glob))

    group_key, bw_files = choose_first_bw_group(bw_candidates)
    fw_pairs = read_fw_ids(args.fw_ids)
    bw_pairs, bw_record_count, dm_var, sd_var = read_bw_pairs(bw_files)

    missing_in_bw = fw_pairs - bw_pairs
    extra_in_bw = bw_pairs - fw_pairs

    print("fw_ids={0}".format(args.fw_ids))
    print("bw_group={0}".format(group_key))
    print("bw_file_count={0}".format(len(bw_files)))
    print("bw_tracking_variables={0},{1}".format(dm_var, sd_var))
    print("fw_unique_pairs={0}".format(len(fw_pairs)))
    print("bw_valid_records={0}".format(bw_record_count))
    print("bw_unique_pairs={0}".format(len(bw_pairs)))
    print("missing_in_bw={0}".format(len(missing_in_bw)))
    print("extra_in_bw={0}".format(len(extra_in_bw)))

    if missing_in_bw:
        print_sample("missing_pairs:", missing_in_bw, args.max_report)
    if extra_in_bw:
        print_sample("extra_pairs:", extra_in_bw, args.max_report)

    if missing_in_bw or extra_in_bw:
        return 1

    print("TPHT consistency check passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
