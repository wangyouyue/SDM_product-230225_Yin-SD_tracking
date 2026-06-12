#!/usr/bin/env python3
"""Audit GMD2026 physical restart chains and TPHT ID handoff settings."""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
GMD_ROOT = SCRIPT_DIR.parent
REPO_ROOT = GMD_ROOT.parents[4]
SPEC_PATH = SCRIPT_DIR / "expected_gmd2026_cases.json"

COLUMNS = [
    "case_group",
    "case_name",
    "policy",
    "status",
    "failures",
    "warnings",
    "ATMOS_RESTART_IN_BASENAME",
    "ATMOS_RESTART_OUTPUT",
    "ATMOS_RESTART_OUT_BASENAME",
    "SD_IN_BASENAME",
    "SD_OUT_BASENAME",
    "TIME_STARTDATE",
    "TIME_STARTMS",
    "TIME_DURATION",
    "TIME_DT_ATMOS_RESTART",
    "tracking_id_input_basename",
    "squid_shared_init_policy",
    "physical_restart_time_check",
    "supersaturation_cap_check",
]

INIT_RESTART_BASENAME = "./init_dycom_rf02_00000101-000000.000"
BASE_RESTART_OUT = "./restart_output/base_restart"
EXPECTED_RESTART_KEYS = {
    "ATMOS_RESTART_IN_BASENAME",
    "ATMOS_RESTART_OUTPUT",
    "ATMOS_RESTART_OUT_BASENAME",
}


def strip_comment(line: str) -> str:
    """Remove comments outside quotes."""
    quote: str | None = None
    result: list[str] = []
    for char in line:
        if char in {"'", '"'}:
            quote = None if quote == char else char if quote is None else quote
        if char == "!" and quote is None:
            break
        result.append(char)
    return "".join(result)


def parse_conf(path: Path) -> dict[str, str]:
    """Parse SCALE namelist assignments."""
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = strip_comment(raw_line).strip()
        if not line or "=" not in line or line.startswith("&"):
            continue
        key, value = line.split("=", 1)
        value = value.strip()
        if value.endswith(",") and value.count(",") == 1:
            value = value[:-1].strip()
        values[key.strip()] = value
    return values


def parse_squid_exports(path: Path) -> dict[str, str]:
    """Parse simple export KEY=value lines from a SQUID wrapper."""
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = raw_line.strip()
        if not line.startswith("export ") or "=" not in line:
            continue
        key, value = line[len("export ") :].split("=", 1)
        values[key.strip()] = clean_string(value.strip())
    return values


def clean_string(value: Any) -> str:
    """Normalize quoted values."""
    if value is None:
        return "NA"
    text = str(value).strip()
    if (text.startswith('"') and text.endswith('"')) or (text.startswith("'") and text.endswith("'")):
        return text[1:-1]
    return text


def to_float(value: Any) -> float | None:
    """Parse Fortran numbers."""
    try:
        return float(clean_string(value).replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def to_bool(value: Any) -> bool | None:
    """Parse Fortran logicals."""
    text = clean_string(value).lower()
    if text in {".true.", "true", "t"}:
        return True
    if text in {".false.", "false", "f"}:
        return False
    return None


def close_float(actual: float | None, expected: float) -> bool:
    """Compare restart times."""
    return actual is not None and math.isclose(actual, expected, rel_tol=1e-8, abs_tol=1e-10)


def model_cases(root: Path) -> list[Path]:
    """Return run.conf paths for model cases exactly two levels below GMD2026."""
    return sorted(path for path in root.glob("*/*/run.conf") if "postprocess" not in path.parts)


def discover_restart_keys() -> set[str]:
    """Discover restart-related keys from existing shallowcloud template conventions."""
    shallowcloud_root = GMD_ROOT.parent
    keys: set[str] = set()
    for conf_path in shallowcloud_root.glob("*/*.conf"):
        if "GMD2026" in conf_path.parts:
            continue
        keys.update(key for key in parse_conf(conf_path) if "RESTART" in key)
    return keys


def source_has_first_hour_cap(repo_root: Path) -> tuple[bool, list[str]]:
    """Find the DYCOMS-II first-hour supersaturation cap implementation."""
    candidates = [
        repo_root / "contrib" / "SDM" / "sdm_condensation_water.f90",
        repo_root / "scale-rm" / "src" / "contrib" / "SDM" / "sdm_condensation_water.f90",
    ]
    evidence: list[str] = []
    for path in candidates:
        if not path.exists():
            continue
        text = path.read_text(errors="ignore")
        if "TIME_NOWSEC<3600" in text.replace(" ", "") or "TIME_NOWSEC < 3600" in text:
            evidence.append(str(path.relative_to(repo_root)))
    return bool(evidence), evidence


class RestartAudit:
    """Accumulate restart-chain findings for one case."""

    def __init__(self, group: str, case: str, conf: dict[str, str], cap_found: bool, case_dir: Path) -> None:
        self.group = group
        self.case = case
        self.conf = conf
        self.cap_found = cap_found
        self.case_dir = case_dir
        self.squid_exports = parse_squid_exports(case_dir / "squid_run.sh")
        self.policy = ""
        self.failures: list[str] = []
        self.warnings: list[str] = []
        self.squid_shared_init_policy = "NA"
        self.physical_restart_time_check = "NA"
        self.supersaturation_cap_check = "NA"

    def require_string(self, key: str, expected: str) -> None:
        """Require a string setting."""
        actual = clean_string(self.conf.get(key))
        if actual != expected:
            self.failures.append(f'{key}="{actual}" expected "{expected}"')

    def require_bool(self, key: str, expected: bool) -> None:
        """Require a logical setting."""
        actual = to_bool(self.conf.get(key))
        if actual is not expected:
            self.failures.append(f"{key}={clean_string(self.conf.get(key))} expected {expected}")

    def require_float(self, key: str, expected: float) -> None:
        """Require a float setting."""
        actual = to_float(self.conf.get(key))
        if not close_float(actual, expected):
            self.failures.append(f"{key}={clean_string(self.conf.get(key))} expected {expected}")

    def check_no_base_overwrite(self) -> None:
        """Ensure branch cases do not write into the common restart basename."""
        out_base = clean_string(self.conf.get("ATMOS_RESTART_OUT_BASENAME"))
        if out_base == BASE_RESTART_OUT or "00_base_restart_3d" in out_base:
            self.failures.append("branch case must not overwrite the common base restart output")

    def reject_fw_discovery_restart_input(self) -> None:
        """Reject physical restart inputs produced by a FW discovery run."""
        for key in ("ATMOS_RESTART_IN_BASENAME", "SD_IN_BASENAME"):
            basename = clean_string(self.conf.get(key))
            if basename == "NA" or basename == "":
                continue
            if "fw_discovery" in basename or "fw_output" in basename or "superdroplet_restart" in basename:
                self.failures.append(
                    f"{key} must not read FW discovery output; TPHT BW must start from the same initial physical state as FW"
                )

    def require_squid_shared_init(self, expected_use_shared: bool, expected_dir: str = "") -> None:
        """Check that the SQUID wrapper matches the intended shared-init policy."""
        use_shared = self.squid_exports.get("GMD_USE_SHARED_INIT", "NA")
        run_init = self.squid_exports.get("GMD_RUN_INIT", "NA")
        shared_dir = self.squid_exports.get("GMD_SHARED_INIT_DIR", "")
        expected_use = "1" if expected_use_shared else "0"
        expected_run_init = "0" if expected_use_shared else "1"
        self.squid_shared_init_policy = (
            f"use_shared={use_shared}, run_init={run_init}, shared_dir={shared_dir or 'NA'}"
        )
        if use_shared != expected_use:
            self.failures.append(f"GMD_USE_SHARED_INIT={use_shared} expected {expected_use}")
        if run_init != expected_run_init:
            self.failures.append(f"GMD_RUN_INIT={run_init} expected {expected_run_init}")
        if expected_use_shared and shared_dir != expected_dir:
            self.failures.append(f'GMD_SHARED_INIT_DIR="{shared_dir}" expected "{expected_dir}"')

    def check_cold_start_window(self, expected_duration: float, label: str) -> None:
        """Check that a case starts from the initial physical state for the requested window."""
        self.require_string("ATMOS_RESTART_IN_BASENAME", INIT_RESTART_BASENAME)
        self.require_bool("ATMOS_RESTART_OUTPUT", False)
        self.require_float("TIME_DURATION", expected_duration)
        self.physical_restart_time_check = label
        if self.cap_found:
            self.supersaturation_cap_check = (
                "source cap uses TIME_NOWSEC < 3600; cold-start cases intentionally include the first-hour regime"
            )
        else:
            self.supersaturation_cap_check = "first-hour supersaturation cap source pattern not found"
            self.warnings.append("could not statically verify first-hour supersaturation cap behavior")

    def apply_policy(self) -> None:
        """Apply the expected restart policy for this group/case."""
        if self.group == "00_base_restart_3d" and self.case == "base_restart_3d":
            self.policy = "cold start, writes 60 min restart"
            self.require_squid_shared_init(False)
            self.require_string("ATMOS_RESTART_IN_BASENAME", INIT_RESTART_BASENAME)
            self.require_bool("ATMOS_RESTART_OUTPUT", True)
            self.require_string("ATMOS_RESTART_OUT_BASENAME", BASE_RESTART_OUT)
            self.require_float("TIME_DURATION", 3600.0)
            self.physical_restart_time_check = "cold start spin-up 0-60 min"
            self.supersaturation_cap_check = "first-hour cap intentionally applies during base restart spin-up"
            return

        if self.group == "01_bench_3d_samp_30min":
            self.policy = "cold start, 30 min controlled benchmark"
            self.require_squid_shared_init(self.case != "nt_nolog", "../nt_nolog" if self.case != "nt_nolog" else "")
            self.check_cold_start_window(1800.0, "cold-start benchmark window 0-30 min")
            self.check_no_base_overwrite()
            return

        if self.group == "05_outint_io_lite_10min":
            self.policy = "cold start, 10 min output-interval benchmark"
            self.require_squid_shared_init(self.case != "nt_nolog_10min", "../nt_nolog_10min" if self.case != "nt_nolog_10min" else "")
            self.check_cold_start_window(600.0, "cold-start output-interval window 0-10 min")
            self.check_no_base_overwrite()
            return

        if self.group == "02_tpht_3d_interest_70min":
            if self.case == "fw_discovery":
                self.policy = "cold start TPHT FW discovery"
                self.require_squid_shared_init(False)
            elif self.case == "bw_reconstruction":
                self.policy = "cold start TPHT BW reconstruction + TPHT ids"
                self.require_squid_shared_init(True, "../fw_discovery")
            else:
                self.policy = "unknown TPHT case"
                self.failures.append("unknown TPHT case")
                return
            self.check_cold_start_window(4200.0, "cold-start TPHT window 0-70 min")
            self.check_no_base_overwrite()
            if self.case == "bw_reconstruction":
                ids_base = clean_string(self.conf.get("tracking_id_input_basename"))
                if ids_base != "../fw_discovery/fw_tracking/tracking_interest_ids_dedup":
                    self.failures.append("TPHT BW must read deduplicated .ids basename")
                self.reject_fw_discovery_restart_input()
                self.require_string("SD_OUT_BASENAME", "./bw_output/superdroplet_restart")
            return

        if self.group == "03_fw_rep_2d_600s":
            self.policy = "cold start, shared init, no 60 min restart"
            self.require_squid_shared_init(
                self.case != "ref_full_stratified",
                "../ref_full_stratified" if self.case != "ref_full_stratified" else "",
            )
            self.require_string("ATMOS_RESTART_IN_BASENAME", INIT_RESTART_BASENAME)
            self.require_bool("ATMOS_RESTART_OUTPUT", False)
            self.require_float("TIME_DURATION", 600.0)
            self.physical_restart_time_check = "cold-start 2D sampling verification with shared initial state"
            self.supersaturation_cap_check = "short cold-start sampling test, not a 60-70 min branch"
            return

        if self.group == "04_sdnc_scaling_lite_30min":
            self.policy = "cold start, 30 min SDNC implementation-scaling benchmark"
            match = re.match(r"^(sdnc\d+)_(nt|fw005|bw005)$", self.case)
            if match:
                self.require_squid_shared_init(match.group(2) != "nt", f"../{match.group(1)}_nt")
            self.check_cold_start_window(1800.0, "cold-start SDNC scaling window 0-30 min")
            return

        if self.group == "06_fortran_diagnostics_hook_smoke_squid":
            self.policy = "short diagnostics-hook smoke test with shared init"
            self.require_squid_shared_init(self.case != "shared_init_hook_2min", "../shared_init_hook_2min" if self.case != "shared_init_hook_2min" else "")
            self.check_cold_start_window(120.0, "short smoke-test window 0-2 min")
            self.check_no_base_overwrite()
            if self.case == "tpht_bw_hook_2min":
                ids_base = clean_string(self.conf.get("tracking_id_input_basename"))
                if ids_base != "../tpht_ids_hook_2min/fw_tracking/tracking_interest_ids_dedup":
                    self.failures.append("TPHT BW smoke test must read deduplicated .ids basename")
                self.reject_fw_discovery_restart_input()
            if self.case == "tpht_bw_missing_ids_failfast":
                ids_base = clean_string(self.conf.get("tracking_id_input_basename"))
                if ids_base != "./missing_ids/tracking_interest_ids_dedup":
                    self.failures.append("missing-ID failfast case must point to the intentionally absent .ids basename")
                self.reject_fw_discovery_restart_input()
            return

        self.policy = "unknown"
        self.failures.append("unknown restart policy group")

    def finalize(self) -> dict[str, str]:
        """Return one output row."""
        status = "FAIL" if self.failures else "WARN" if self.warnings else "PASS"
        return {
            "case_group": self.group,
            "case_name": self.case,
            "policy": self.policy,
            "status": status,
            "failures": "; ".join(self.failures),
            "warnings": "; ".join(self.warnings),
            "ATMOS_RESTART_IN_BASENAME": clean_string(self.conf.get("ATMOS_RESTART_IN_BASENAME")),
            "ATMOS_RESTART_OUTPUT": clean_string(self.conf.get("ATMOS_RESTART_OUTPUT")),
            "ATMOS_RESTART_OUT_BASENAME": clean_string(self.conf.get("ATMOS_RESTART_OUT_BASENAME")),
            "SD_IN_BASENAME": clean_string(self.conf.get("SD_IN_BASENAME")),
            "SD_OUT_BASENAME": clean_string(self.conf.get("SD_OUT_BASENAME")),
            "TIME_STARTDATE": clean_string(self.conf.get("TIME_STARTDATE")),
            "TIME_STARTMS": clean_string(self.conf.get("TIME_STARTMS")),
            "TIME_DURATION": clean_string(self.conf.get("TIME_DURATION")),
            "TIME_DT_ATMOS_RESTART": clean_string(self.conf.get("TIME_DT_ATMOS_RESTART")),
            "tracking_id_input_basename": clean_string(self.conf.get("tracking_id_input_basename")),
            "squid_shared_init_policy": self.squid_shared_init_policy,
            "physical_restart_time_check": self.physical_restart_time_check,
            "supersaturation_cap_check": self.supersaturation_cap_check,
        }


def write_outputs(rows: list[dict[str, str]], output_base: Path) -> None:
    """Write restart audit outputs."""
    with output_base.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    output_base.with_suffix(".json").write_text(json.dumps(rows, indent=2))
    md_lines = [
        "| " + " | ".join(COLUMNS) + " |",
        "| " + " | ".join(["---"] * len(COLUMNS)) + " |",
    ]
    for row in rows:
        md_lines.append("| " + " | ".join(str(row.get(column, "NA")) for column in COLUMNS) + " |")
    output_base.with_suffix(".md").write_text("\n".join(md_lines) + "\n")


def main() -> int:
    spec = json.loads(SPEC_PATH.read_text())
    discovered_restart_keys = discover_restart_keys()
    missing_convention_keys = sorted(EXPECTED_RESTART_KEYS - discovered_restart_keys)
    cap_found, cap_evidence = source_has_first_hour_cap(REPO_ROOT)
    rows: list[dict[str, str]] = []

    tpht_restarts: dict[str, str] = {}
    audits: list[RestartAudit] = []
    for run_conf in model_cases(GMD_ROOT):
        case_dir = run_conf.parent
        audit = RestartAudit(case_dir.parent.name, case_dir.name, parse_conf(run_conf), cap_found, case_dir)
        audit.apply_policy()
        if missing_convention_keys:
            audit.failures.append("restart keys not found in existing templates: " + ", ".join(missing_convention_keys))
        if audit.group == "02_tpht_3d_interest_70min":
            tpht_restarts[audit.case] = clean_string(audit.conf.get("ATMOS_RESTART_IN_BASENAME"))
        audits.append(audit)

    fw_restart = tpht_restarts.get("fw_discovery")
    bw_restart = tpht_restarts.get("bw_reconstruction")
    if fw_restart and bw_restart and fw_restart != bw_restart:
        for audit in audits:
            if audit.group == "02_tpht_3d_interest_70min":
                audit.failures.append("TPHT FW and BW must read the same initial physical state")

    if cap_found:
        for audit in audits:
            if audit.group == "00_base_restart_3d":
                audit.warnings.append("first-hour cap evidence for spin-up: " + ", ".join(cap_evidence))

    rows = [audit.finalize() for audit in sorted(audits, key=lambda item: (item.group, item.case))]
    write_outputs(rows, GMD_ROOT / "restart_audit")
    pass_count = sum(row["status"] == "PASS" for row in rows)
    fail_count = sum(row["status"] == "FAIL" for row in rows)
    warn_count = sum(row["status"] == "WARN" for row in rows)
    print(f"PASS: {pass_count} cases")
    print(f"FAIL: {fail_count} cases")
    print(f"WARN: {warn_count} cases")
    return 1 if fail_count else 0


if __name__ == "__main__":
    sys.exit(main())
