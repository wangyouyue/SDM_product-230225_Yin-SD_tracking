#!/usr/bin/env python
from __future__ import print_function

import argparse
from glob import glob
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-glob",
        default="./ft_interest_id_baseline/fw_tracking/tracking_interest_ids.pe*.ids",
    )
    parser.add_argument(
        "--output",
        default="./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_merged.ids",
    )
    parser.add_argument(
        "--rank-bucket-basename",
        default="",
        help=(
            "Optional basename for deduplicated per-rank ID outputs. "
            "Example: ./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_dedup "
            "writes *.pe000000.ids, *.pe000001.ids, ..."
        ),
    )
    return parser.parse_args()


def collect_records_nc(file_paths):
    from netCDF4 import Dataset

    merged = {}
    tracking_id_output_basename = ""
    radius_datatype = "f8"
    tpht_meta = None

    for file_path in file_paths:
        with Dataset(file_path, "r") as nc:
            if not tracking_id_output_basename and "tracking_id_output_basename" in nc.ncattrs():
                tracking_id_output_basename = nc.getncattr("tracking_id_output_basename")
            if tpht_meta is None:
                attrs = nc.ncattrs()
                if (
                    "tpht_prc_num_x" in attrs
                    and "tpht_prc_num_y" in attrs
                    and "tpht_prc_nprocs" in attrs
                ):
                    tpht_meta = (
                        int(nc.getncattr("tpht_prc_num_x")),
                        int(nc.getncattr("tpht_prc_num_y")),
                        int(nc.getncattr("tpht_prc_nprocs")),
                    )

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

    return merged, tracking_id_output_basename, radius_datatype, tpht_meta


def collect_records_ids(file_paths):
    merged = set()
    tpht_meta = None

    for file_path in file_paths:
        with open(file_path, "r") as ids_in:
            for raw_line in ids_in:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 5 and parts[0] == "#" and parts[1] == "TPHT_META":
                        current_meta = (int(parts[2]), int(parts[3]), int(parts[4]))
                        if tpht_meta is None:
                            tpht_meta = current_meta
                        elif tpht_meta != current_meta:
                            raise SystemExit(
                                "Inconsistent TPHT_META across input files: {0} vs {1}".format(
                                    tpht_meta, current_meta
                                )
                            )
                    continue

                parts = line.split()
                if len(parts) < 2:
                    continue
                merged.add((int(parts[0]), int(parts[1])))

    if tpht_meta is None:
        raise SystemExit("Missing TPHT_META header in input ID files.")

    return merged, tpht_meta


def write_ids_output(output_ids_path, sorted_pairs, tpht_meta):
    output_dir = output_ids_path.rsplit("/", 1)[0] if "/" in output_ids_path else "."
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    with open(output_ids_path, "w") as ids_out:
        ids_out.write("# TPHT_META {0} {1} {2}\n".format(tpht_meta[0], tpht_meta[1], tpht_meta[2]))
        for dm_id, sd_id in sorted_pairs:
            ids_out.write("{0} {1}\n".format(dm_id, sd_id))

    return output_ids_path


def write_rank_bucket_ids_output(rank_bucket_basename, sorted_pairs, tpht_meta):
    if not rank_bucket_basename:
        return []

    if rank_bucket_basename.endswith(".ids"):
        rank_bucket_basename = rank_bucket_basename[:-4]

    output_dir = rank_bucket_basename.rsplit("/", 1)[0] if "/" in rank_bucket_basename else "."
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    nprocs = int(tpht_meta[2])
    buckets = {rank: [] for rank in range(nprocs)}
    for dm_id, sd_id in sorted_pairs:
        if dm_id < 0 or dm_id >= nprocs:
            raise SystemExit(
                "dm_id {0} is outside TPHT_META nprocs={1}".format(dm_id, nprocs)
            )
        buckets[dm_id].append((dm_id, sd_id))

    output_paths = []
    for rank in range(nprocs):
        output_path = "{0}.pe{1:06d}.ids".format(rank_bucket_basename, rank)
        with open(output_path, "w") as ids_out:
            ids_out.write("# TPHT_META {0} {1} {2}\n".format(tpht_meta[0], tpht_meta[1], tpht_meta[2]))
            for dm_id, sd_id in buckets[rank]:
                ids_out.write("{0} {1}\n".format(dm_id, sd_id))
        output_paths.append(output_path)

    return output_paths


def write_output_nc(output_path, merged, tracking_id_output_basename, radius_datatype, source_count, tpht_meta):
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
        if tpht_meta is not None:
            nc.setncattr("tpht_prc_num_x", tpht_meta[0])
            nc.setncattr("tpht_prc_num_y", tpht_meta[1])
            nc.setncattr("tpht_prc_nprocs", tpht_meta[2])

    if output_path.endswith(".nc"):
        output_ids_path = output_path[:-3] + ".ids"
    else:
        output_ids_path = output_path + ".ids"

    sorted_pairs = list(zip(dm_id, sd_id))
    write_ids_output(output_ids_path, sorted_pairs, tpht_meta)

    return output_ids_path


def main():
    args = parse_args()
    input_files = sorted(glob(args.input_glob))
    if not input_files:
        raise SystemExit("No files matched: {0}".format(args.input_glob))

    first_ext = os.path.splitext(input_files[0])[1].lower()
    output_path = args.output

    if first_ext == ".ids":
        merged, tpht_meta = collect_records_ids(input_files)
        sorted_pairs = sorted(merged)
        if output_path.endswith(".nc"):
            output_ids_path = output_path[:-3] + ".ids"
        else:
            output_ids_path = output_path
            if not output_ids_path.endswith(".ids"):
                output_ids_path = output_ids_path + ".ids"
        write_ids_output(output_ids_path, sorted_pairs, tpht_meta)
        rank_bucket_paths = write_rank_bucket_ids_output(
            args.rank_bucket_basename, sorted_pairs, tpht_meta
        )
        print("input_files={0}".format(len(input_files)))
        print("unique_pairs={0}".format(len(merged)))
        print("output_ids={0}".format(output_ids_path))
        if rank_bucket_paths:
            print("rank_bucket_file_count={0}".format(len(rank_bucket_paths)))
            print("rank_bucket_first={0}".format(rank_bucket_paths[0]))
        return

    merged, tracking_id_output_basename, radius_datatype, tpht_meta = collect_records_nc(input_files)
    output_ids_path = write_output_nc(output_path, merged, tracking_id_output_basename, radius_datatype, len(input_files), tpht_meta)
    sorted_pairs = sorted(merged)
    rank_bucket_paths = write_rank_bucket_ids_output(
        args.rank_bucket_basename, sorted_pairs, tpht_meta
    )

    print("input_files={0}".format(len(input_files)))
    print("unique_pairs={0}".format(len(merged)))
    print("output_nc={0}".format(output_path))
    print("output_ids={0}".format(output_ids_path))
    if rank_bucket_paths:
        print("rank_bucket_file_count={0}".format(len(rank_bucket_paths)))
        print("rank_bucket_first={0}".format(rank_bucket_paths[0]))


if __name__ == "__main__":
    main()
