#!/usr/bin/env python3
"""Static checks for lifecycle-safe TPHT tracking IDs."""

import pathlib
import re
import sys


REPO_ROOT = pathlib.Path(__file__).resolve().parents[5]

CHECK_FILES = [
    "contrib/SDM/sdm_idutil.f90",
    "contrib/SDM/sdm_io.f90",
    "contrib/SDM/scale_atmos_phy_mp_sdm.F90",
    "scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_sdm.F90",
    "scale-rm/test/case/shallowcloud/ice_tpht_smoke/check_tpht_consistency.py",
    "scale-rm/test/case/shallowcloud/ice_tpht_smoke/analyze_tpht_tracks.py",
    "scale-rm/test/case/shallowcloud/tpht_test/check_tpht_consistency.py",
    "scale-rm/test/case/shallowcloud/tpht_test/analyze_tpht_tracks.py",
    "scale-rm/test/case/shallowcloud/tpht_test/verify_tpht_restart_consistency.py",
]

DANGEROUS_PATTERNS = [
    (re.compile(r"\bsd_id\b\s*<\s*0"), "generic sd_id < 0 invalid check"),
    (re.compile(r"\bsd_id\b\s*<=\s*INVALID_i4"), "generic sd_id <= INVALID_i4 invalid check"),
    (re.compile(r"\bsd_id\b\s*>\s*INVALID_i4"), "generic sd_id > INVALID_i4 valid check"),
    (re.compile(r"\bsd_id\(n\)\s*<=\s*INVALID_i4"), "per-slot sd_id <= INVALID_i4 invalid check"),
    (re.compile(r"\bsd_id\(n\)\s*>\s*INVALID_i4"), "per-slot sd_id > INVALID_i4 valid check"),
    (re.compile(r"\bpair\[1\]\s*<\s*0"), "generic pair[1] < 0 sd_id invalid check"),
]


def iter_hits(path, text):
    for lineno, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("!") or stripped.startswith("#"):
            continue
        for pattern, label in DANGEROUS_PATTERNS:
            if pattern.search(line):
                yield path, lineno, label, stripped


def main():
    failures = []
    for rel_path in CHECK_FILES:
        path = REPO_ROOT / rel_path
        if not path.exists():
            failures.append((rel_path, 0, "missing checked file", ""))
            continue
        text = path.read_text()
        failures.extend(iter_hits(rel_path, text))

    idutil_text = (REPO_ROOT / "contrib/SDM/sdm_idutil.f90").read_text()
    if "sdm_tracking_find_slot_by_id" not in idutil_text:
        failures.append(("contrib/SDM/sdm_idutil.f90", 0, "missing BW ID lookup helper", ""))
    if "target_sd_id(m) >= 1" in idutil_text:
        failures.append(
            (
                "contrib/SDM/sdm_idutil.f90",
                0,
                "BW lookup still has slot-index-only target_sd_id path",
                "target_sd_id(m) >= 1",
            )
        )

    reset_hits = []
    for rel_path in CHECK_FILES:
        path = REPO_ROOT / rel_path
        if path.exists():
            for lineno, line in enumerate(path.read_text().splitlines(), 1):
                if "call sdm_tracking_reset_invalid_slot" in line:
                    reset_hits.append((rel_path, lineno, line.strip()))
    approved = {
        "contrib/SDM/scale_atmos_phy_mp_sdm.F90",
        "scalelib/src/atmos-physics/microphysics/scale_atmos_phy_mp_sdm.F90",
    }
    for rel_path, lineno, line in reset_hits:
        if rel_path not in approved:
            failures.append((rel_path, lineno, "reset helper call outside approved lifecycle sites", line))

    if failures:
        print("Dynamic ID predicate validation failed:")
        for rel_path, lineno, label, line in failures:
            location = f"{rel_path}:{lineno}" if lineno else rel_path
            print(f"  {location}: {label}: {line}")
        return 1

    print("Dynamic ID predicate validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
