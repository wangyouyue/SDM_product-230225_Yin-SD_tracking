#!/usr/bin/env python3
"""Static checks for terminal precipitation/outflow tracking cleanup."""

from __future__ import annotations

import argparse
import pathlib
import re
import sys


REPO_ROOT = pathlib.Path(__file__).resolve().parents[5]
DRIVER = REPO_ROOT / "contrib/SDM/scale_atmos_phy_mp_sdm.F90"
TRACKING = REPO_ROOT / "contrib/SDM/sdm_tracking_cold.f90"
IO = REPO_ROOT / "contrib/SDM/sdm_io.f90"


FORBIDDEN_FILES = [
    "contrib/SDM/sdm_boundary.f90",
    "contrib/SDM/sdm_coalescence_cold.f90",
    "contrib/SDM/sdm_sd2fluid.f90",
]

CLEANUP_RE = re.compile(
    r"COLD_TPHT_CLEANUP\s+"
    r"(?P<label>\S+)\s+"
    r"rank=\s*(?P<rank>-?\d+)\s+"
    r"candidate=\s*(?P<candidate>\d+)\s+"
    r"invalid_after=\s*(?P<invalid_after>\d+)\s+"
    r"newly_invalid=\s*(?P<newly_invalid>\d+)\s+"
    r"reset=\s*(?P<reset>\d+)\s+"
    r"lifecycle=\s*(?P<lifecycle>\d+)"
)


def fail(failures: list[str], message: str) -> None:
    failures.append(message)


def parse_cleanup_markers(path: pathlib.Path) -> list[dict[str, int | str]]:
    records: list[dict[str, int | str]] = []
    text = path.read_text(errors="replace")
    for match in CLEANUP_RE.finditer(text):
        item: dict[str, int | str] = {"label": match.group("label")}
        for key in ("rank", "candidate", "invalid_after", "newly_invalid", "reset", "lifecycle"):
            item[key] = int(match.group(key))
        records.append(item)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=pathlib.Path, help="optional runtime log containing COLD_TPHT_CLEANUP markers")
    parser.add_argument(
        "--require-label",
        action="append",
        default=[],
        help="cleanup label that must have newly_invalid>0 and reset==newly_invalid",
    )
    args = parser.parse_args()

    failures: list[str] = []
    driver = DRIVER.read_text()
    tracking = TRACKING.read_text()
    io_text = IO.read_text()

    if "subroutine sdm_tracking_cleanup_invalidated_slots" not in driver:
        fail(failures, "driver missing centralized invalidated-slot cleanup helper")
    if driver.count("call sdm_tracking_cleanup_invalidated_slots") < 2:
        fail(failures, "driver should call cleanup helper for vertical outflow and sdm_sd2prec paths")
    if "call sdm_jdginvdv" not in driver or "call sdm_sd2prec" not in driver:
        fail(failures, "driver missing expected vertical/sd2prec call sites")
    if "sd_itmp3(:) = 0" not in driver:
        fail(failures, "driver cleanup candidate marker is missing")
    if "if( sd_rk(n) >= VALID2INVALID ) sd_itmp3(n) = 1" not in driver:
        fail(failures, "vertical cleanup candidate marker is missing")
    if "if( sd_rk(n) < VALID2INVALID .and. sd_rk(n) > PREC2INVALID ) sd_itmp3(n) = 1" not in driver:
        fail(failures, "sdm_sd2prec precipitation candidate marker is missing")
    if "LIFE_SDREMOVE_INVALIDATION" not in driver:
        fail(failures, "cleanup helper does not write/remain aware of lifecycle invalidation code")
    if "sdm_tracking_reset_invalid_slot" not in driver:
        fail(failures, "cleanup helper does not call reset helper")
    if "sd_event_sig_mask" not in tracking:
        fail(failures, "reset helper/source does not include sd_event_sig_mask")
    if "sd_event_sig_mask" not in driver:
        fail(failures, "driver cleanup/copy path does not include sd_event_sig_mask")

    for rel_path in FORBIDDEN_FILES:
        text = (REPO_ROOT / rel_path).read_text()
        if "sdm_tracking_cleanup_invalidated_slots" in text or "sdm_tracking_reset_invalid_slot" in text:
            fail(failures, f"forbidden cleanup helper call in {rel_path}")

    for name, value in {
        "LIFE_DOMAIN_ENTRY": 3,
        "LIFE_GLOBAL_HALO_ENTRY": 4,
        "LIFE_SEEDING_ENTRY": 5,
    }.items():
        if not re.search(rf"integer,\s*parameter\s*::\s*{name}\s*=\s*{value}\b", tracking):
            fail(failures, f"missing lifecycle entry constant {name}={value}")

    for helper in (
        "sdm_tracking_lifecycle_write_domain_entry",
        "sdm_tracking_lifecycle_write_global_halo_entry",
        "sdm_tracking_lifecycle_write_seeding_entry",
    ):
        if f"public :: {helper}" not in io_text:
            fail(failures, f"{helper} is not public")
        if f"subroutine {helper}" not in io_text:
            fail(failures, f"{helper} subroutine is missing")
        driver_calls = [line for line in driver.splitlines() if helper in line and not line.strip().startswith("!")]
        if driver_calls:
            fail(failures, f"{helper} should remain inactive in driver path")

    if args.log is not None:
        if not args.log.exists():
            fail(failures, f"cleanup runtime log not found: {args.log}")
        else:
            records = parse_cleanup_markers(args.log)
            if not records:
                fail(failures, f"no COLD_TPHT_CLEANUP markers found in {args.log}")
            for label in args.require_label:
                matching = [item for item in records if item["label"] == label]
                if not matching:
                    fail(failures, f"missing cleanup marker label {label}")
                    continue
                newly_invalid = sum(int(item["newly_invalid"]) for item in matching)
                reset_count = sum(int(item["reset"]) for item in matching)
                lifecycle_count = sum(int(item["lifecycle"]) for item in matching)
                if newly_invalid <= 0:
                    fail(failures, f"{label} did not report newly_invalid>0")
                if reset_count != newly_invalid:
                    fail(failures, f"{label} reset count {reset_count} != newly_invalid {newly_invalid}")
                if lifecycle_count > reset_count:
                    fail(failures, f"{label} lifecycle count {lifecycle_count} exceeds reset count {reset_count}")

    if failures:
        print("Invalid cleanup static validation failed:")
        for item in failures:
            print(f"  {item}")
        return 1

    if args.log is not None:
        print(f"Invalid cleanup runtime markers validated: {args.log}")
    else:
        print("Invalid cleanup static validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
