#!/usr/bin/env python
from __future__ import print_function

import argparse
from glob import glob
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-glob",
        default="./ft_interest_id_baseline/fw_tracking/tracking_interest_ids.pe*.nc",
    )
    parser.add_argument(
        "--output",
        default="./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.nc",
    )
    return parser.parse_args()


def collect_records(file_paths):
    from netCDF4 import Dataset

    merged = {}
    tracking_id_output_basename = ""
    radius_datatype = "f8"

    for file_path in file_paths:
        with Dataset(file_path, "r") as nc:
            if not tracking_id_output_basename and "tracking_id_output_basename" in nc.ncattrs():
                tracking_id_output_basename = nc.getncattr("tracking_id_output_basename")

            radius_datatype = nc.variables["sd_r"].datatype
            dm_id = nc.variables["dm_id"][:]
            sd_id = nc.variables["sd_id"][:]
            first_time = nc.variables["first_time"][:]
            sd_r = nc.variables["sd_r"][:]
            flag_radius = nc.variables["flag_radius"][:]
            flag_coalescence = nc.variables["flag_coalescence"][:]

            for idx in range(len(dm_id)):
                key = (int(dm_id[idx]), int(sd_id[idx]))
                record = merged.get(key)
                if record is None:
                    merged[key] = {
                        "first_time": float(first_time[idx]),
                        "sd_r": float(sd_r[idx]),
                        "flag_radius": int(flag_radius[idx]),
                        "flag_coalescence": int(flag_coalescence[idx]),
                    }
                    continue

                record["flag_radius"] = max(record["flag_radius"], int(flag_radius[idx]))
                record["flag_coalescence"] = max(record["flag_coalescence"], int(flag_coalescence[idx]))

                if float(first_time[idx]) < record["first_time"]:
                    record["first_time"] = float(first_time[idx])
                    record["sd_r"] = float(sd_r[idx])

    return merged, tracking_id_output_basename, radius_datatype


def write_output(output_path, merged, tracking_id_output_basename, radius_datatype, source_count):
    from netCDF4 import Dataset

    sorted_items = sorted(
        merged.items(),
        key=lambda item: (item[1]["first_time"], item[0][0], item[0][1]),
    )

    dm_id = [item[0][0] for item in sorted_items]
    sd_id = [item[0][1] for item in sorted_items]
    first_time = [item[1]["first_time"] for item in sorted_items]
    sd_r = [item[1]["sd_r"] for item in sorted_items]
    flag_radius = [item[1]["flag_radius"] for item in sorted_items]
    flag_coalescence = [item[1]["flag_coalescence"] for item in sorted_items]

    output_dir = output_path.rsplit("/", 1)[0] if "/" in output_path else "."
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    with Dataset(output_path, "w", format="NETCDF4") as nc:
        nc.createDimension("tracked_id", len(dm_id))
        nc.createVariable("dm_id", "i4", ("tracked_id",))[:] = dm_id
        nc.createVariable("sd_id", "i4", ("tracked_id",))[:] = sd_id
        nc.createVariable("first_time", "f8", ("tracked_id",))[:] = first_time
        nc.createVariable("sd_r", radius_datatype, ("tracked_id",))[:] = sd_r
        nc.createVariable("flag_radius", "i4", ("tracked_id",))[:] = flag_radius
        nc.createVariable("flag_coalescence", "i4", ("tracked_id",))[:] = flag_coalescence
        if tracking_id_output_basename:
            nc.setncattr("tracking_id_output_basename", tracking_id_output_basename)
        nc.setncattr("merge_source_file_count", source_count)


def main():
    args = parse_args()
    input_files = sorted(glob(args.input_glob))
    if not input_files:
        raise SystemExit("No files matched: {0}".format(args.input_glob))

    merged, tracking_id_output_basename, radius_datatype = collect_records(input_files)
    output_path = args.output
    write_output(output_path, merged, tracking_id_output_basename, radius_datatype, len(input_files))

    print("input_files={0}".format(len(input_files)))
    print("unique_pairs={0}".format(len(merged)))
    print("output={0}".format(output_path))


if __name__ == "__main__":
    main()
