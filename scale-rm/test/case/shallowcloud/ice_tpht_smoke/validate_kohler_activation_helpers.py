#!/usr/bin/env python3
"""Static and helper-level checks for derived Kohler activation events."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[5]


def read_rel(path: str) -> str:
    return (REPO_ROOT / path).read_text()


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def kohler_rcrit(curve_f: float, asl_ff: float, temperature: float, terms: list[tuple[float, float, float]]) -> float:
    coef_a = curve_f / temperature
    coef_b = asl_ff * sum(mass * ion / mol_weight for mass, ion, mol_weight in terms)
    if coef_a <= 0.0 or coef_b <= 0.0:
        return math.inf
    return math.sqrt(3.0 * coef_b / coef_a)


def activated_state(valid: bool, phase: str, radius: float, rcrit: float, eps: float = 0.0) -> bool:
    return valid and phase == "liquid" and radius > rcrit * (1.0 + eps)


def crossing(pre: bool, post: bool) -> str:
    if (not pre) and post:
        return "activation"
    if pre and (not post):
        return "deactivation"
    return "none"


def check_helper_math(failures: list[str]) -> None:
    rcrit = kohler_rcrit(
        curve_f=3.3e-7,
        asl_ff=4.3e-6,
        temperature=280.0,
        terms=[(2.0e-18, 3.0, 0.13214)],
    )
    require(math.isfinite(rcrit) and rcrit > 0.0, "Python reference rcrit should be finite positive", failures)

    require(
        crossing(
            activated_state(True, "liquid", 0.95 * rcrit, rcrit),
            activated_state(True, "liquid", 1.05 * rcrit, rcrit),
        )
        == "activation",
        "helper semantics should identify inactive->active crossing",
        failures,
    )
    require(
        crossing(
            activated_state(True, "liquid", 1.05 * rcrit, rcrit),
            activated_state(True, "liquid", 0.95 * rcrit, rcrit),
        )
        == "deactivation",
        "helper semantics should identify active->inactive crossing",
        failures,
    )
    require(
        crossing(
            activated_state(True, "liquid", 0.90 * rcrit, rcrit),
            activated_state(True, "liquid", 0.95 * rcrit, rcrit),
        )
        == "none",
        "helper semantics should not report false->false as crossing",
        failures,
    )
    require(
        crossing(
            activated_state(True, "liquid", 1.05 * rcrit, rcrit),
            activated_state(True, "liquid", 1.10 * rcrit, rcrit),
        )
        == "none",
        "helper semantics should not report true->true as crossing",
        failures,
    )
    require(
        not activated_state(False, "liquid", 2.0 * rcrit, rcrit),
        "invalid SD must not be considered activated",
        failures,
    )
    require(
        not activated_state(True, "dry", 2.0 * rcrit, rcrit),
        "dry aerosol / phase 0 must not be considered activated",
        failures,
    )
    require(
        not activated_state(True, "ice", 2.0 * rcrit, rcrit),
        "ice must not be considered activated",
        failures,
    )
    require(
        not activated_state(True, "mixed", 2.0 * rcrit, rcrit),
        "mixed phase must not be considered activated in first implementation",
        failures,
    )
    require(
        not activated_state(True, "liquid", 1.005 * rcrit, rcrit, eps=0.01),
        "hysteresis should suppress within-band activation when eps is nonzero",
        failures,
    )


def check_source_contracts(failures: list[str]) -> None:
    tracking = read_rel("contrib/SDM/sdm_tracking_cold.f90")
    idutil = read_rel("contrib/SDM/sdm_idutil.f90")
    driver = read_rel("contrib/SDM/scale_atmos_phy_mp_sdm.F90")
    io_src = read_rel("contrib/SDM/sdm_io.f90")
    common = read_rel("contrib/SDM/sdm_common.f90")

    require("public :: sdm_kohler_critical_radius" in tracking, "sdm_kohler_critical_radius must be public", failures)
    require("public :: sdm_kohler_activated_state" in tracking, "sdm_kohler_activated_state must be public", failures)
    require("public :: sdm_kohler_solute_params" in tracking, "sdm_kohler_solute_params must be public", failures)
    require("function sdm_kohler_critical_radius" in tracking, "sdm_kohler_critical_radius function is missing", failures)
    require("function sdm_kohler_activated_state" in tracking, "sdm_kohler_activated_state function is missing", failures)
    require("subroutine sdm_kohler_solute_params" in tracking, "sdm_kohler_solute_params subroutine is missing", failures)
    require(
        "sqrt(3.0_RP*coef_b/coef_a)" not in idutil,
        "selected-output activated criterion should call shared helper, not inline sqrt(3*coef_b/coef_a)",
        failures,
    )
    require(
        "sdm_kohler_solute_params" in idutil
        and "sdm_kohler_critical_radius" in idutil
        and "sdm_kohler_activated_state" in idutil,
        "sdm_copy_selected_sd should use shared Kohler solute and activation helper functions",
        failures,
    )
    require(
        "sdm_kohler_solute_params" in driver and "sd_kohler_rcrit_pre" in driver and "sd_kohler_active_pre" in driver,
        "condevp wrapper should use shared solute setup and capture pre/post Kohler activation state",
        failures,
    )
    require(
        "TRIG_PROC_ACTIVATION" in io_src and "TRIG_PROC_DEACTIVATION" in io_src,
        "singleproc writer should import activation/deactivation trigger codes",
        failures,
    )
    require(
        "TRACK_COLD_OUTPUT_KOHLER_CONTEXT" in common and "TRACK_COLD_OUTPUT_KOHLER_CONTEXT" in io_src,
        "optional Kohler context switch should be declared and used by singleproc writer",
        failures,
    )
    for name in (
        "kohler_rcrit_pre",
        "kohler_rcrit_post",
        "kohler_margin_pre",
        "kohler_margin_post",
        "activated_state_pre",
        "activated_state_post",
    ):
        require(name in io_src, f"optional singleproc Kohler context variable is missing: {name}", failures)

    require(
        "TRIG_PROC_ACTIVATION" not in read_rel("contrib/SDM/sdm_io.f90").split("subroutine sdm_lifecycle_outnetcdf", 1)[-1],
        "activation/deactivation must not be written through lifecycle sidecar",
        failures,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    failures: list[str] = []
    check_helper_math(failures)
    check_source_contracts(failures)

    if failures:
        print("Kohler activation helper validation failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("Kohler activation helper validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
