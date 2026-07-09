#!/usr/bin/env python3
"""Static/unit checks for Cold TPHT Level-3 spatial-visit tracking."""

from __future__ import annotations

import pathlib
import re
import sys


REPO_ROOT = pathlib.Path(__file__).resolve().parents[5]
TRACKING = REPO_ROOT / "contrib/SDM/sdm_tracking_cold.f90"
DRIVER = REPO_ROOT / "contrib/SDM/scale_atmos_phy_mp_sdm.F90"


def segment_intersects_box(p0, p1, box) -> bool:
    tmin = 0.0
    tmax = 1.0
    for axis, (lo, hi) in enumerate(((box[0], box[1]), (box[2], box[3]), (box[4], box[5]))):
        d = p1[axis] - p0[axis]
        if abs(d) <= 1.0e-300:
            if not lo <= p0[axis] <= hi:
                return False
            continue
        t1 = (lo - p0[axis]) / d
        t2 = (hi - p0[axis]) / d
        tlo = min(t1, t2)
        thi = max(t1, t2)
        tmin = max(tmin, tlo)
        tmax = min(tmax, thi)
        if tmax < tmin:
            return False
    return True


def run_geometry_tests(failures: list[str]) -> None:
    box = (0.49, 0.51, -1.0, 1.0, -1.0, 1.0)
    cases = [
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), True, "thin x crossing"),
        ((0.0, 2.0, 0.0), (1.0, 2.0, 0.0), False, "parallel outside y"),
        ((0.5, 0.0, 0.0), (0.5, 0.5, 0.0), True, "endpoint starts inside"),
        ((0.0, 0.0, 2.0), (1.0, 0.0, 2.0), False, "outside z slab"),
    ]
    for p0, p1, expected, label in cases:
        observed = segment_intersects_box(p0, p1, box)
        if observed != expected:
            failures.append(f"geometry unit failed for {label}: expected {expected}, got {observed}")


def check_static(failures: list[str]) -> None:
    tracking_text = TRACKING.read_text()
    driver_text = DRIVER.read_text()

    required_tracking = [
        "public :: sdm_cold_tracking_update_spatial_visit_segment",
        "public :: sdm_cold_segment_intersects_box",
        "subroutine sdm_cold_tracking_update_spatial_visit_segment",
        "logical function sdm_cold_segment_intersects_box",
        "call update_axis(x0, x1",
        "call update_axis(y0, y1",
        "call update_axis(z0, z1",
    ]
    for needle in required_tracking:
        if needle not in tracking_text:
            failures.append(f"missing Level-3 tracking source: {needle}")

    required_driver = [
        "track_spatial_segments = sdm_cold .and. tracking_spatial_visit_enable",
        "allocate(spatial_prev_x(1:sd_num))",
        "spatial_prev_x(:) = sd_x(:)",
        "call sdm_cold_tracking_update_spatial_visit_segment",
    ]
    for needle in required_driver:
        if needle not in driver_text:
            failures.append(f"missing Level-3 driver source: {needle}")

    match = re.search(
        r"call\s+sdm_cold_tracking_update_spatial_visit_segment\b.*?call\s+sdm_boundary\b",
        driver_text,
        re.IGNORECASE | re.DOTALL,
    )
    if not match:
        failures.append("Level-3 segment check is not before the horizontal boundary call")

    if "tracking_spatial_region_mask" in driver_text or "sd_spatial_region_mask" in driver_text:
        failures.append("multiple-region spatial mask appeared in v1.2+ Level-3 implementation")


def main() -> int:
    failures: list[str] = []
    run_geometry_tests(failures)
    check_static(failures)
    if failures:
        print("Spatial Level-3 validation failed:")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("Spatial Level-3 static/unit validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
