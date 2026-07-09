#!/usr/bin/env python3
"""Test-only harness for Cold TPHT lifecycle helper semantics.

This script intentionally does not run the production SD-number adjustment path.
It mirrors the lifecycle helper contracts with minimal in-memory SD slots and
creates a harness-only SD_lifecycle_NetCDF_* file via ncgen when available.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict


REPO_ROOT = pathlib.Path(__file__).resolve().parents[5]
DEFAULT_OUTDIR = pathlib.Path("/private/tmp/cold_tpht_lifecycle_harness")

TRACK_ID_INVALID = -999
TRACK_SD_ID_DYNAMIC_START = -1000

REQUIRED_FIELDS = [
    "time",
    "lifecycle_code",
    "n_parent",
    "n_child",
    "parent_sd_id1",
    "parent_dm_id1",
    "parent_sd_id2",
    "parent_dm_id2",
    "child_sd_id1",
    "child_dm_id1",
    "child_sd_id2",
    "child_dm_id2",
    "child_sd_id3",
    "child_dm_id3",
    "child_sd_id4",
    "child_dm_id4",
    "x",
    "y",
    "z",
]


@dataclass
class SlotState:
    sd_id: int
    dm_id: int
    sd_event_mask: int
    sd_event_sig_mask: int
    sd_diag_mask: int
    sd_phase_change_flag: int
    sd_spatial_visit_flag: int
    liq_radius_max: float
    ice_rvol_max: float
    mixed_rvol_max: float
    rime_mass_max: float
    rime_frac_max: float
    nmono_max: float
    aspect_ratio_max: float


class HarnessFailure(Exception):
    pass


def valid_sd_id(value: int) -> bool:
    return value >= 0 or value <= TRACK_SD_ID_DYNAMIC_START


def invalid_sd_id(value: int) -> bool:
    return value == TRACK_ID_INVALID


def valid_dm_id(value: int) -> bool:
    return value >= 0


def reset_invalid_slot(slot: SlotState) -> None:
    slot.sd_id = TRACK_ID_INVALID
    slot.dm_id = TRACK_ID_INVALID
    slot.sd_event_mask = 0
    slot.sd_event_sig_mask = 0
    slot.sd_diag_mask = 0
    slot.sd_phase_change_flag = 0
    slot.sd_spatial_visit_flag = 0
    slot.liq_radius_max = 0.0
    slot.ice_rvol_max = 0.0
    slot.mixed_rvol_max = 0.0
    slot.rime_mass_max = 0.0
    slot.rime_frac_max = 0.0
    slot.nmono_max = 0.0
    slot.aspect_ratio_max = 0.0


def assign_dynamic_id(slot: SlotState, next_dynamic_sd_id: int, mype: int) -> int:
    slot.sd_id = next_dynamic_sd_id
    slot.dm_id = mype
    return next_dynamic_sd_id - 1


def copy_interval_state(src: SlotState, dst: SlotState) -> None:
    dst.sd_event_mask = src.sd_event_mask
    dst.sd_event_sig_mask = src.sd_event_sig_mask
    dst.sd_diag_mask = src.sd_diag_mask
    dst.sd_phase_change_flag = src.sd_phase_change_flag
    dst.sd_spatial_visit_flag = src.sd_spatial_visit_flag
    dst.liq_radius_max = src.liq_radius_max
    dst.ice_rvol_max = src.ice_rvol_max
    dst.mixed_rvol_max = src.mixed_rvol_max
    dst.rime_mass_max = src.rime_mass_max
    dst.rime_frac_max = src.rime_frac_max
    dst.nmono_max = src.nmono_max
    dst.aspect_ratio_max = src.aspect_ratio_max


def lifecycle_record(
    *,
    time: float,
    lifecycle_code: int,
    n_parent: int,
    n_child: int,
    parent_ids: list[tuple[int, int]] | None = None,
    child_ids: list[tuple[int, int]] | None = None,
    x: float,
    y: float,
    z: float,
) -> dict[str, int | float]:
    parent_ids = list(parent_ids or [])
    child_ids = list(child_ids or [])
    parent = parent_ids[:2] + [(TRACK_ID_INVALID, TRACK_ID_INVALID)] * (2 - len(parent_ids[:2]))
    child = child_ids[:4] + [(TRACK_ID_INVALID, TRACK_ID_INVALID)] * (4 - len(child_ids[:4]))
    record: dict[str, int | float] = {
        "time": time,
        "lifecycle_code": lifecycle_code,
        "n_parent": n_parent,
        "n_child": n_child,
        "parent_sd_id1": parent[0][0],
        "parent_dm_id1": parent[0][1],
        "parent_sd_id2": parent[1][0],
        "parent_dm_id2": parent[1][1],
        "child_sd_id1": child[0][0],
        "child_dm_id1": child[0][1],
        "child_sd_id2": child[1][0],
        "child_dm_id2": child[1][1],
        "child_sd_id3": child[2][0],
        "child_dm_id3": child[2][1],
        "child_sd_id4": child[3][0],
        "child_dm_id4": child[3][1],
        "x": x,
        "y": y,
        "z": z,
    }
    return record


def parse_source_constants() -> dict[str, int]:
    tracking_text = (REPO_ROOT / "contrib/SDM/sdm_tracking_cold.f90").read_text()
    common_text = (REPO_ROOT / "contrib/SDM/sdm_common.f90").read_text()
    constants = {}
    for name in (
        "LIFE_SDADD_SPLIT",
        "LIFE_SDREMOVE_INVALIDATION",
        "LIFE_DOMAIN_ENTRY",
        "LIFE_GLOBAL_HALO_ENTRY",
        "LIFE_SEEDING_ENTRY",
        "LIFE_ASLFORM_NEW_SD",
    ):
        match = re.search(rf"integer,\s*parameter\s*::\s*{name}\s*=\s*(-?\d+)\b", tracking_text)
        if not match:
            raise HarnessFailure(f"missing lifecycle constant {name}")
        constants[name] = int(match.group(1))
    match = re.search(r"integer,\s*parameter\s*::\s*TRACK_SD_ID_DYNAMIC_START\s*=\s*(-?\d+)\b", common_text)
    if not match:
        raise HarnessFailure("missing TRACK_SD_ID_DYNAMIC_START")
    constants["TRACK_SD_ID_DYNAMIC_START"] = int(match.group(1))
    match = re.search(r"integer,\s*parameter\s*::\s*TRACK_ID_INVALID\s*=\s*INVALID_i4", common_text)
    if not match:
        raise HarnessFailure("TRACK_ID_INVALID is not tied to INVALID_i4")
    return constants


def assert_slot_reset(slot: SlotState) -> None:
    if slot.sd_id != TRACK_ID_INVALID or slot.dm_id != TRACK_ID_INVALID:
        raise HarnessFailure(f"slot IDs not reset: {slot}")
    for key, value in asdict(slot).items():
        if key in {"sd_id", "dm_id"}:
            continue
        if value != 0:
            raise HarnessFailure(f"slot field {key} not reset: {value}")


def run_harness() -> tuple[list[dict[str, int | float]], dict[str, object]]:
    constants = parse_source_constants()
    if constants["TRACK_SD_ID_DYNAMIC_START"] != TRACK_SD_ID_DYNAMIC_START:
        raise HarnessFailure("dynamic ID start changed without updating harness")

    nonzero = SlotState(
        sd_id=10,
        dm_id=0,
        sd_event_mask=7,
        sd_event_sig_mask=2,
        sd_diag_mask=3,
        sd_phase_change_flag=1,
        sd_spatial_visit_flag=1,
        liq_radius_max=1.0,
        ice_rvol_max=2.0,
        mixed_rvol_max=3.0,
        rime_mass_max=4.0,
        rime_frac_max=5.0,
        nmono_max=6.0,
        aspect_ratio_max=7.0,
    )

    reset_only = SlotState(**asdict(nonzero))
    reset_invalid_slot(reset_only)
    assert_slot_reset(reset_only)

    remove_parent = SlotState(**asdict(nonzero))
    remove_record = lifecycle_record(
        time=1.0,
        lifecycle_code=constants["LIFE_SDREMOVE_INVALIDATION"],
        n_parent=1,
        n_child=0,
        parent_ids=[(remove_parent.sd_id, remove_parent.dm_id)],
        child_ids=[],
        x=100.0,
        y=200.0,
        z=300.0,
    )
    reset_invalid_slot(remove_parent)
    assert_slot_reset(remove_parent)

    parent = SlotState(**asdict(nonzero))
    parent.sd_id = 5
    parent.dm_id = 0
    child = SlotState(
        sd_id=TRACK_ID_INVALID,
        dm_id=TRACK_ID_INVALID,
        sd_event_mask=0,
        sd_event_sig_mask=0,
        sd_diag_mask=0,
        sd_phase_change_flag=0,
        sd_spatial_visit_flag=0,
        liq_radius_max=0.0,
        ice_rvol_max=0.0,
        mixed_rvol_max=0.0,
        rime_mass_max=0.0,
        rime_frac_max=0.0,
        nmono_max=0.0,
        aspect_ratio_max=0.0,
    )
    next_dynamic_id = TRACK_SD_ID_DYNAMIC_START
    next_dynamic_id = assign_dynamic_id(child, next_dynamic_id, mype=0)
    copy_interval_state(parent, child)
    if parent.sd_id != 5 or parent.dm_id != 0:
        raise HarnessFailure("split parent ID changed")
    if not valid_sd_id(child.sd_id) or child.sd_id > TRACK_SD_ID_DYNAMIC_START:
        raise HarnessFailure(f"split child dynamic sd_id invalid: {child.sd_id}")
    if child.sd_id == parent.sd_id or not valid_dm_id(child.dm_id):
        raise HarnessFailure("split child ID pair invalid")
    for key in (
        "sd_event_mask",
        "sd_event_sig_mask",
        "sd_diag_mask",
        "sd_phase_change_flag",
        "sd_spatial_visit_flag",
        "liq_radius_max",
        "ice_rvol_max",
        "mixed_rvol_max",
        "rime_mass_max",
        "rime_frac_max",
        "nmono_max",
        "aspect_ratio_max",
    ):
        if getattr(child, key) != getattr(parent, key):
            raise HarnessFailure(f"split child did not copy {key}")
    split_record = lifecycle_record(
        time=2.0,
        lifecycle_code=constants["LIFE_SDADD_SPLIT"],
        n_parent=1,
        n_child=1,
        parent_ids=[(parent.sd_id, parent.dm_id)],
        child_ids=[(child.sd_id, child.dm_id)],
        x=110.0,
        y=210.0,
        z=310.0,
    )

    parentless_records = []
    for offset, code_name in enumerate(
        ("LIFE_DOMAIN_ENTRY", "LIFE_GLOBAL_HALO_ENTRY", "LIFE_SEEDING_ENTRY"), start=1
    ):
        test_sd_id = TRACK_SD_ID_DYNAMIC_START - offset
        if not valid_sd_id(test_sd_id) or not invalid_sd_id(TRACK_ID_INVALID):
            raise HarnessFailure("ID predicate mismatch in parentless harness")
        parentless_records.append(
            lifecycle_record(
                time=2.0 + offset,
                lifecycle_code=constants[code_name],
                n_parent=0,
                n_child=1,
                parent_ids=[],
                child_ids=[(test_sd_id, 0)],
                x=10.0 * offset,
                y=20.0 * offset,
                z=30.0 * offset,
            )
        )

    records = [remove_record, split_record] + parentless_records
    for record in records:
        for field in REQUIRED_FIELDS:
            if field not in record:
                raise HarnessFailure(f"record missing field {field}")
        if record["n_parent"] == 0:
            for field in ("parent_sd_id1", "parent_dm_id1", "parent_sd_id2", "parent_dm_id2"):
                if record[field] != TRACK_ID_INVALID:
                    raise HarnessFailure(f"parentless record has non-fill {field}: {record[field]}")
        if record["n_child"] == 0:
            for field in (
                "child_sd_id1",
                "child_dm_id1",
                "child_sd_id2",
                "child_dm_id2",
                "child_sd_id3",
                "child_dm_id3",
                "child_sd_id4",
                "child_dm_id4",
            ):
                if record[field] != TRACK_ID_INVALID:
                    raise HarnessFailure(f"childless record has non-fill {field}: {record[field]}")

    summary = {
        "reset_only": "pass",
        "remove_record": "pass",
        "split_record": "pass",
        "parentless_entry_records": len(parentless_records),
        "next_dynamic_id_after_split": next_dynamic_id,
        "records": len(records),
    }
    return records, summary


def cdl_value_list(records: list[dict[str, int | float]], field: str) -> str:
    values = [record[field] for record in records]
    if field in {"time", "x", "y", "z"}:
        return ", ".join(f"{float(value):.10g}" for value in values)
    return ", ".join(str(int(value)) for value in values)


def write_cdl(records: list[dict[str, int | float]], cdl_path: pathlib.Path) -> None:
    lines = [
        "netcdf SD_lifecycle_NetCDF_harness {",
        "dimensions:",
        f"    event = {len(records)} ;",
        "variables:",
        "    double time(event) ;",
        "    int lifecycle_code(event) ;",
        "    int n_parent(event) ;",
        "    int n_child(event) ;",
        "    int parent_sd_id1(event) ;",
        "    int parent_dm_id1(event) ;",
        "    int parent_sd_id2(event) ;",
        "    int parent_dm_id2(event) ;",
        "    int child_sd_id1(event) ;",
        "    int child_dm_id1(event) ;",
        "    int child_sd_id2(event) ;",
        "    int child_dm_id2(event) ;",
        "    int child_sd_id3(event) ;",
        "    int child_dm_id3(event) ;",
        "    int child_sd_id4(event) ;",
        "    int child_dm_id4(event) ;",
        "    double x(event) ;",
        "    double y(event) ;",
        "    double z(event) ;",
        "data:",
    ]
    for field in REQUIRED_FIELDS:
        lines.append(f"    {field} = {cdl_value_list(records, field)} ;")
    lines.append("}")
    cdl_path.write_text("\n".join(lines) + "\n")


def write_optional_netcdf(records: list[dict[str, int | float]], outdir: pathlib.Path) -> pathlib.Path | None:
    ncgen = shutil.which("ncgen")
    ncdump = shutil.which("ncdump")
    if not ncgen or not ncdump:
        return None
    cdl_path = outdir / "SD_lifecycle_NetCDF_harness.cdl"
    nc_path = outdir / "SD_lifecycle_NetCDF_harness.nc"
    write_cdl(records, cdl_path)
    subprocess.run([ncgen, "-o", str(nc_path), str(cdl_path)], check=True)
    dump = subprocess.run([ncdump, str(nc_path)], check=True, text=True, capture_output=True).stdout
    for field in REQUIRED_FIELDS:
        if field not in dump:
            raise HarnessFailure(f"ncdump output missing field {field}")
    for needle in ("lifecycle_code = 2, 1, 3, 4, 5", "n_parent = 1, 1, 0, 0, 0", "n_child = 0, 1, 1, 1, 1"):
        if needle not in dump:
            raise HarnessFailure(f"ncdump output missing expected values: {needle}")
    (outdir / "SD_lifecycle_NetCDF_harness.ncdump.txt").write_text(dump)
    return nc_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    try:
        args.outdir.mkdir(parents=True, exist_ok=True)
        records, summary = run_harness()
        records_path = args.outdir / "lifecycle_harness_records.json"
        summary_path = args.outdir / "lifecycle_harness_summary.json"
        records_path.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
        nc_path = write_optional_netcdf(records, args.outdir)
        summary["records_json"] = str(records_path)
        summary["netcdf_file"] = str(nc_path) if nc_path else "not written; ncgen/ncdump unavailable"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    except (HarnessFailure, subprocess.CalledProcessError) as exc:
        print(f"Lifecycle harness validation failed: {exc}")
        return 1

    print("Lifecycle harness validation passed")
    print(f"  records: {summary['records']}")
    print(f"  records_json: {summary['records_json']}")
    print(f"  netcdf_file: {summary['netcdf_file']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
