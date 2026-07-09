#!/usr/bin/env python3
"""Validate Cold TPHT spatial rank-subdomain pruning.

Without runtime arguments this performs static checks. Optional runtime
arguments compare sd_spatial_visit_flag counts between conservative and
rank-pruned outputs and inspect LOG.pe* active-rank diagnostics.
"""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[5]
COMMON = REPO_ROOT / "contrib/SDM/sdm_common.f90"
DRIVER = REPO_ROOT / "contrib/SDM/scale_atmos_phy_mp_sdm.F90"
SCALIB_DRIVER = REPO_ROOT / "scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_sdm.F90"
CONFIG_GENERATOR = (
    REPO_ROOT
    / "scale-rm/test/case/shallowcloud/ice_tpht_smoke/squid_validation/make_v12_validation_configs.py"
)


def ncdump_variable(path: Path, field: str) -> str:
    try:
        result = subprocess.run(
            ["ncdump", "-v", field, str(path)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise RuntimeError("ncdump is required for runtime spatial pruning validation") from None
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"ncdump failed for {path}: {exc.stderr.strip()}") from exc
    try:
        return result.stdout.split("data:", 1)[1]
    except IndexError:
        raise RuntimeError(f"{path}: ncdump output has no data section") from None


def read_int_values(path: Path, field: str) -> list[int]:
    text = ncdump_variable(path, field)
    body = re.sub(r"//.*", "", text)
    return [int(value) for value in re.findall(r"-?\d+", body)]


def count_spatial_flags(pattern: str) -> tuple[int, int, int]:
    paths = sorted(Path(path) for path in glob.glob(pattern) if not path.endswith(".ids"))
    if not paths:
        raise RuntimeError(f"no NetCDF files matched {pattern}")
    total = 0
    visited = 0
    for path in paths:
        values = read_int_values(path, "sd_spatial_visit_flag")
        total += len(values)
        visited += sum(1 for value in values if value != 0)
    return len(paths), total, visited


def parse_active_rank_logs(pattern: str) -> tuple[int, int]:
    paths = sorted(Path(path) for path in glob.glob(pattern))
    if not paths:
        raise RuntimeError(f"no LOG files matched {pattern}")
    active = 0
    inactive = 0
    for path in paths:
        text = path.read_text(errors="replace")
        for match in re.finditer(
            r"tracking_spatial_active_rank\s*=\s*([TF])\s+pruning\s*=\s*([TF])",
            text,
        ):
            if match.group(2) != "T":
                continue
            if match.group(1) == "T":
                active += 1
            else:
                inactive += 1
        for match in re.finditer(
            r"tracking_spatial_rank_pruning\s+mype\s+active\s+margin\s*=\s*"
            r"\d+\s+([TF])\s+[-+0-9.Ee]+",
            text,
        ):
            if match.group(1) == "T":
                active += 1
            else:
                inactive += 1
    return active, inactive


def check_static() -> list[str]:
    failures: list[str] = []
    common_text = COMMON.read_text()
    driver_text = DRIVER.read_text()
    scalib_text = SCALIB_DRIVER.read_text()
    generator_text = CONFIG_GENERATOR.read_text()

    required_common = [
        "tracking_spatial_rank_pruning_enable = .false.",
        "tracking_spatial_rank_pruning_enable, &",
        "tracking_spatial_rank_margin, &",
    ]
    for needle in required_common:
        if needle not in common_text:
            failures.append(f"missing common namelist/default source: {needle}")

    required_driver = [
        "if( tracking_spatial_rank_pruning_enable ) then",
        "spatial_margin = max(0.0_RP, tracking_spatial_rank_margin)",
        "rank_x_min = min(GRID_FX(IS-1), GRID_FX(IE))",
        "rank_y_min = min(GRID_FY(JS-1), GRID_FY(JE))",
        "rank_z_min = minval(REAL_FZ(KS-1:KE,IS:IE,JS:JE))",
        "tracking_spatial_active_rank = &",
        "tracking_spatial_active_rank = .true.",
    ]
    for needle in required_driver:
        if needle not in driver_text:
            failures.append(f"missing driver pruning source: {needle}")
        if needle not in scalib_text:
            failures.append(f"missing scalelib driver pruning source: {needle}")

    required_configs = [
        "run_v13_spatial_rank_pruning_off_rank2.conf",
        "run_v13_spatial_rank_pruning_on_rank2.conf",
        "run_v13_spatial_rank_pruning_all_overlap_rank2.conf",
        "run_v13_spatial_rank_pruning_margin_rank2.conf",
        "tracking_spatial_rank_pruning_enable = .true.,",
    ]
    for needle in required_configs:
        if needle not in generator_text:
            failures.append(f"missing generated pruning config source: {needle}")

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-glob", help="Conservative SD_all/SD_selected NetCDF glob")
    parser.add_argument("--candidate-glob", help="Pruned SD_all/SD_selected NetCDF glob")
    parser.add_argument("--expect-any", action="store_true", help="Require candidate visited count > 0")
    parser.add_argument("--log-glob", help="LOG.pe* glob for pruning active-rank diagnostics")
    parser.add_argument("--expect-active-ranks", type=int)
    parser.add_argument("--expect-inactive-ranks", type=int)
    args = parser.parse_args()

    failures = check_static()
    if args.baseline_glob or args.candidate_glob:
        if not args.baseline_glob or not args.candidate_glob:
            failures.append("--baseline-glob and --candidate-glob must be used together")
        else:
            try:
                b_files, b_total, b_visited = count_spatial_flags(args.baseline_glob)
                c_files, c_total, c_visited = count_spatial_flags(args.candidate_glob)
                if b_total != c_total or b_visited != c_visited:
                    failures.append(
                        "spatial pruning changed visit flags: "
                        f"baseline files={b_files} total={b_total} visited={b_visited}; "
                        f"candidate files={c_files} total={c_total} visited={c_visited}"
                    )
                if args.expect_any and c_visited <= 0:
                    failures.append("candidate output has zero spatial visits")
            except RuntimeError as exc:
                failures.append(str(exc))

    if args.log_glob:
        try:
            active, inactive = parse_active_rank_logs(args.log_glob)
            if args.expect_active_ranks is not None and active != args.expect_active_ranks:
                failures.append(f"expected {args.expect_active_ranks} active ranks, observed {active}")
            if args.expect_inactive_ranks is not None and inactive != args.expect_inactive_ranks:
                failures.append(f"expected {args.expect_inactive_ranks} inactive ranks, observed {inactive}")
        except RuntimeError as exc:
            failures.append(str(exc))

    if failures:
        print("Spatial rank-pruning validation failed:")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print("Spatial rank-pruning validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
