#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import os
import re
import struct
import sys
from array import array


RANK_RE = re.compile(r"\.pe(\d{6})(?:\..*)?$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Verify that TPHT FW handoff IDs and BW restart-selected IDs are consistent. "
            "The script reads deduplicated per-rank .ids files and SD restart files, "
            "then compares (dm_id, sd_id) sets rank by rank."
        )
    )
    parser.add_argument(
        "--ids-basename",
        required=True,
        help=(
            "Basename of deduplicated per-rank TPHT IDs, for example "
            "./ft_interest_id_baseline/fw_tracking/tracking_interest_ids_dedup"
        ),
    )
    parser.add_argument(
        "--bw-restart-glob",
        required=True,
        help=(
            "Glob for BW superdroplet restart files, for example "
            "./bt_interest_id_baseline/bw_output/superdroplet_restart_*"
        ),
    )
    parser.add_argument(
        "--fw-restart-glob",
        default="",
        help=(
            "Optional glob for FW superdroplet restart files. When provided, the script "
            "also checks whether the TPHT target IDs exist in the FW restart archive."
        ),
    )
    parser.add_argument(
        "--record-marker-bytes",
        type=int,
        default=4,
        choices=(4, 8),
        help="Fortran sequential unformatted record-marker size in bytes. Default: 4.",
    )
    parser.add_argument(
        "--endian",
        default="<",
        choices=("<", ">"),
        help="Binary endianness for restart files: '<' little, '>' big. Default: '<'.",
    )
    parser.add_argument(
        "--show-samples",
        type=int,
        default=5,
        help="How many missing/extra/sample IDs to print per mismatching rank. Default: 5.",
    )
    return parser.parse_args()


def rank_from_path(path: str) -> int:
    match = RANK_RE.search(path)
    if not match:
        raise SystemExit(f"Could not infer rank from path: {path}")
    return int(match.group(1))


def ids_path_from_basename(basename: str, rank: int) -> str:
    if basename.endswith(".ids"):
        basename = basename[:-4]
    return f"{basename}.pe{rank:06d}.ids"


def read_ids_file(path: str):
    tpht_meta = None
    pairs = set()
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("# TPHT_META"):
                    parts = line.split()
                    if len(parts) >= 5:
                        tpht_meta = tuple(map(int, parts[2:5]))
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            dm_id = int(parts[0])
            sd_id = int(parts[1])
            pairs.add((dm_id, sd_id))
    return tpht_meta, pairs


def read_marker(fid, marker_bytes: int, endian: str) -> int:
    raw = fid.read(marker_bytes)
    if not raw:
        raise EOFError
    if len(raw) != marker_bytes:
        raise SystemExit("Unexpected EOF while reading Fortran record marker")
    fmt = endian + ("I" if marker_bytes == 4 else "Q")
    return struct.unpack(fmt, raw)[0]


def read_fortran_record(fid, marker_bytes: int, endian: str) -> bytes:
    nbytes = read_marker(fid, marker_bytes, endian)
    payload = fid.read(nbytes)
    if len(payload) != nbytes:
        raise SystemExit("Unexpected EOF while reading Fortran record payload")
    nbytes_end = read_marker(fid, marker_bytes, endian)
    if nbytes != nbytes_end:
        raise SystemExit(
            f"Fortran record marker mismatch: begin={nbytes}, end={nbytes_end}"
        )
    return payload


def parse_restart_header(payload: bytes, endian: str):
    if len(payload) == 28:
        otime, rp_kind, dp_kind, sd_num, sd_numasl, sdfmnum = struct.unpack(
            endian + "diiiii", payload
        )
    elif len(payload) == 24:
        otime, rp_kind, dp_kind, sd_num, sd_numasl, sdfmnum = struct.unpack(
            endian + "fiiiii", payload
        )
    else:
        raise SystemExit(
            f"Unexpected restart header size {len(payload)} bytes; expected 24 or 28"
        )
    return {
        "otime": float(otime),
        "rp_kind": int(rp_kind),
        "dp_kind": int(dp_kind),
        "sd_num": int(sd_num),
        "sd_numasl": int(sd_numasl),
        "sdfmnum": int(sdfmnum),
    }


def parse_int32_array(payload: bytes, count: int, endian: str):
    if len(payload) != count * 4:
        raise SystemExit(
            f"Unexpected integer payload size {len(payload)} bytes for count={count}"
        )
    arr = array("i")
    arr.frombytes(payload)
    if sys.byteorder == "big" and endian == "<":
        arr.byteswap()
    if sys.byteorder == "little" and endian == ">":
        arr.byteswap()
    return arr


def read_restart_valid_pairs(path: str, marker_bytes: int, endian: str):
    with open(path, "rb") as fid:
        header = parse_restart_header(read_fortran_record(fid, marker_bytes, endian), endian)

        # Skip records written before sdid/dmid:
        # sdn, sdrk, sdx, sdy, sdz, sdr, sdu, sdv, sdvz, sdasl, sdliqice
        for _ in range(11):
            read_fortran_record(fid, marker_bytes, endian)

        sdid = parse_int32_array(
            read_fortran_record(fid, marker_bytes, endian), header["sd_num"], endian
        )
        dmid = parse_int32_array(
            read_fortran_record(fid, marker_bytes, endian), header["sd_num"], endian
        )

    valid_pairs = set()
    for sd_id, dm_id in zip(sdid, dmid):
        if (sd_id >= 0 or sd_id <= -1000) and dm_id >= 0:
            valid_pairs.add((int(dm_id), int(sd_id)))

    return header, valid_pairs


def collect_restart_sets(glob_pattern: str, marker_bytes: int, endian: str):
    paths = sorted(glob.glob(glob_pattern))
    if not paths:
        raise SystemExit(f"No restart files matched: {glob_pattern}")
    by_rank = {}
    headers = {}
    for path in paths:
        rank = rank_from_path(path)
        header, pairs = read_restart_valid_pairs(path, marker_bytes, endian)
        by_rank[rank] = pairs
        headers[rank] = header
    return by_rank, headers


def collect_ids_sets(ids_basename: str, ranks):
    by_rank = {}
    metas = {}
    for rank in ranks:
        path = ids_path_from_basename(ids_basename, rank)
        if not os.path.exists(path):
            raise SystemExit(f"Missing TPHT ID file for rank {rank}: {path}")
        meta, pairs = read_ids_file(path)
        by_rank[rank] = pairs
        metas[rank] = meta
    return by_rank, metas


def summarize_rank_diff(label: str, rank: int, expected: set, actual: set, sample_n: int):
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    print(
        f"{label} rank={rank:06d} expected={len(expected)} actual={len(actual)} "
        f"missing={len(missing)} extra={len(extra)}"
    )
    if missing:
        print(f"  missing_sample={missing[:sample_n]}")
    if extra:
        print(f"  extra_sample={extra[:sample_n]}")
    return len(missing) == 0 and len(extra) == 0


def main():
    args = parse_args()

    bw_sets, bw_headers = collect_restart_sets(
        args.bw_restart_glob, args.record_marker_bytes, args.endian
    )
    ranks = sorted(bw_sets)
    ids_sets, ids_metas = collect_ids_sets(args.ids_basename, ranks)

    print("=== TPHT BW vs deduplicated handoff ===")
    bw_ok = True
    for rank in ranks:
        meta = ids_metas[rank]
        if meta is not None:
            prc_num_x, prc_num_y, nprocs = meta
            if rank == 0:
                print(
                    f"TPHT_META rank0: PRC_NUM_X={prc_num_x} PRC_NUM_Y={prc_num_y} PRC_nprocs={nprocs}"
                )
        ok = summarize_rank_diff(
            "BW_MATCH", rank, ids_sets[rank], bw_sets[rank], args.show_samples
        )
        bw_ok = bw_ok and ok

    if args.fw_restart_glob:
        fw_sets, _ = collect_restart_sets(
            args.fw_restart_glob, args.record_marker_bytes, args.endian
        )
        print("=== TPHT FW restart coverage ===")
        fw_ok = True
        for rank in ranks:
            if rank not in fw_sets:
                print(f"FW_COVER rank={rank:06d} missing restart file")
                fw_ok = False
                continue
            expected = ids_sets[rank]
            actual = fw_sets[rank]
            missing = sorted(expected - actual)
            print(
                f"FW_COVER rank={rank:06d} expected={len(expected)} "
                f"present_in_restart={len(expected) - len(missing)} missing={len(missing)}"
            )
            if missing:
                print(f"  missing_sample={missing[:args.show_samples]}")
                fw_ok = False
    else:
        fw_ok = True

    print("=== Restart header sample ===")
    rank0 = ranks[0]
    print(f"BW header rank0: {bw_headers[rank0]}")

    if bw_ok and fw_ok:
        print("RESULT: PASS")
        return 0

    print("RESULT: FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
