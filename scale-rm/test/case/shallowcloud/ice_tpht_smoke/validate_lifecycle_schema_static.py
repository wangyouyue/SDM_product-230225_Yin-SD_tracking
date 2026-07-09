#!/usr/bin/env python3
"""Static checks for the Cold TPHT lifecycle sidecar skeleton."""

import pathlib
import re
import sys


REPO_ROOT = pathlib.Path(__file__).resolve().parents[5]

SCHEMA_FIELDS = [
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

LIFECYCLE_CODES = {
    "LIFE_SDADD_SPLIT": 1,
    "LIFE_SDREMOVE_INVALIDATION": 2,
    "LIFE_DOMAIN_ENTRY": 3,
    "LIFE_GLOBAL_HALO_ENTRY": 4,
    "LIFE_SEEDING_ENTRY": 5,
    "LIFE_ASLFORM_NEW_SD": 6,
}


def read(rel_path):
    return (REPO_ROOT / rel_path).read_text()


def fail(failures, label, detail=""):
    failures.append(f"{label}: {detail}" if detail else label)


def extract_subroutine(text, name):
    pattern = re.compile(
        rf"subroutine\s+{re.escape(name)}\b(?P<body>.*?)end\s+subroutine\s+{re.escape(name)}",
        re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(text)
    return match.group(0) if match else ""


def check_schema(failures):
    io_text = read("contrib/SDM/sdm_io.f90")
    if "public :: sdm_lifecycle_outnetcdf" not in io_text:
        fail(failures, "sdm_lifecycle_outnetcdf is not public")

    body = extract_subroutine(io_text, "sdm_lifecycle_outnetcdf")
    if not body:
        fail(failures, "missing sdm_lifecycle_outnetcdf subroutine")
        return

    if "SD_lifecycle" not in body or "_NetCDF_" not in body:
        fail(failures, "lifecycle writer does not construct SD_lifecycle_NetCDF_* name")
    if 'nf90_def_dim(ncid, "event", NF90_UNLIMITED' not in body:
        fail(failures, "lifecycle writer missing unlimited event dimension")
    if "integer :: parent_sd(2), parent_dm(2), child_sd(4), child_dm(4)" not in body:
        fail(failures, "lifecycle writer does not use fixed-width 2-parent/4-child buffers")
    if "parent_sd(:) = TRACK_ID_INVALID" not in body or "child_sd(:) = TRACK_ID_INVALID" not in body:
        fail(failures, "lifecycle writer does not fill unused ID slots with TRACK_ID_INVALID")
    if "ncopy = min(2, max(0, n_parent)" not in body:
        fail(failures, "lifecycle writer missing parent clamp to 2")
    if "ncopy = min(4, max(0, n_child)" not in body:
        fail(failures, "lifecycle writer missing child clamp to 4")

    for field in SCHEMA_FIELDS:
        if field == "time":
            expected = f'def_time("{field}"'
        elif field in {"x", "y", "z"}:
            expected = f'def_real("{field}"'
        else:
            expected = f'def_int("{field}"'
        if expected not in body:
            fail(failures, "missing lifecycle schema definition", field)
        if field != "time" and field not in {"x", "y", "z"}:
            varid = f"{field}_id"
            if varid not in body:
                fail(failures, "missing lifecycle variable id", varid)

    for helper, code_name in (
        ("sdm_tracking_lifecycle_write_domain_entry", "LIFE_DOMAIN_ENTRY"),
        ("sdm_tracking_lifecycle_write_global_halo_entry", "LIFE_GLOBAL_HALO_ENTRY"),
        ("sdm_tracking_lifecycle_write_seeding_entry", "LIFE_SEEDING_ENTRY"),
    ):
        if f"public :: {helper}" not in io_text:
            fail(failures, "lifecycle entry helper is not public", helper)
        helper_body = extract_subroutine(io_text, helper)
        if not helper_body:
            fail(failures, "missing lifecycle entry helper", helper)
            continue
        if code_name not in helper_body:
            fail(failures, "lifecycle entry helper missing code", f"{helper}/{code_name}")
        if f"call sdm_lifecycle_outnetcdf(otime, {code_name}, 0, 1" not in helper_body:
            fail(failures, "lifecycle entry helper must write n_parent=0 n_child=1", helper)
        if "TRACK_ID_INVALID" not in helper_body:
            fail(failures, "lifecycle entry helper must fill parent IDs with invalid sentinel", helper)


def check_codes(failures):
    tracking_text = read("contrib/SDM/sdm_tracking_cold.f90")
    for name, value in LIFECYCLE_CODES.items():
        pattern = re.compile(rf"integer,\s*parameter\s*::\s*{name}\s*=\s*{value}\b")
        if not pattern.search(tracking_text):
            fail(failures, "missing or changed lifecycle code", f"{name}={value}")
        public_lines = [line for line in tracking_text.splitlines() if "public ::" in line and name in line]
        if not public_lines:
            fail(failures, "lifecycle code is not public", name)

    driver_text = read("contrib/SDM/scale_atmos_phy_mp_sdm.F90")
    if "LIFE_SDREMOVE_INVALIDATION" not in driver_text:
        fail(failures, "driver missing LIFE_SDREMOVE_INVALIDATION reference")
    if "LIFE_SDADD_SPLIT" not in driver_text:
        fail(failures, "driver missing LIFE_SDADD_SPLIT reference")
    if "sdm_lifecycle_outnetcdf" not in driver_text:
        fail(failures, "driver missing lifecycle writer call")


def check_restart_counter(failures):
    driver_text = read("contrib/SDM/scale_atmos_phy_mp_sdm.F90")
    required = [
        "SDM_COLD_TRACKING_RESTART_V4",
        "cold_tracking_restart_version_v4 = 4",
        "read(fid_sd_i, iostat=ierr) tracking_next_dynamic_sd_id",
        "write(fid_sd_o) tracking_next_dynamic_sd_id",
        "tracking_next_dynamic_sd_id = TRACK_SD_ID_DYNAMIC_START",
    ]
    for needle in required:
        if needle not in driver_text:
            fail(failures, "missing dynamic-ID restart counter handling", needle)


def check_no_active_adjsdnum(failures):
    driver_text = read("contrib/SDM/scale_atmos_phy_mp_sdm.F90")
    for lineno, line in enumerate(driver_text.splitlines(), 1):
        stripped = line.strip().lower()
        if stripped.startswith("!") or not stripped:
            continue
        if "call sdm_adjsdnum" in stripped:
            fail(failures, "sdm_adjsdnum production path appears active", f"line {lineno}: {line.strip()}")


def main():
    failures = []
    check_schema(failures)
    check_codes(failures)
    check_restart_counter(failures)
    check_no_active_adjsdnum(failures)

    if failures:
        print("Lifecycle schema static validation failed:")
        for item in failures:
            print(f"  {item}")
        return 1

    print("Lifecycle schema static validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
