#!/usr/bin/env python3
"""Audit GMD2026 namelist settings against the expected experiment design."""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
GMD_ROOT = SCRIPT_DIR.parent
SPEC_PATH = SCRIPT_DIR / "expected_gmd2026_cases.json"

AUDIT_KEYS = [
    "dimension",
    "domain_size",
    "resolution",
    "TIME_DT",
    "TIME_DT_ATMOS_DYN",
    "TIME_DT_ATMOS_PHY_MP",
    "TIME_DT_over_TIME_DT_ATMOS_DYN",
    "ATMOS_DYN_TYPE",
    "ATMOS_DYN_TINTEG_SHORT_TYPE",
    "ATMOS_DYN_wdamp_height",
    "ATMOS_DYN_wdamp_tau",
    "tracking_mode",
    "tracking_selection_mode",
    "tracking_fraction",
    "max_tracked_sds",
    "tracking_height_min",
    "tracking_height_max",
    "tracking_radius_min",
    "tracking_radius_max",
    "coalescence_output_enable",
    "gmd_benchmark_diag_enable",
    "sdm_dmpvar",
    "sdm_dmpitvl",
    "HISTORY_DEFAULT_TINTERVAL",
    "HISTORY_OUTPUT_STEP0",
    "tracking_id_output_basename",
    "tracking_id_input_basename",
    "SD_IN_BASENAME",
    "SD_OUT_BASENAME",
    "tracking_interest_radius_enable",
    "tracking_interest_radius_threshold",
    "tracking_interest_coalescence_enable",
    "random_perturbation_enable",
    "random_perturbation_amp",
    "tracking_sampling_seed",
]

COLUMNS = ["case_group", "case_name", "status", "failures", "warnings"] + AUDIT_KEYS


def strip_comment(line: str) -> str:
    """Remove Fortran namelist comments while preserving quoted strings."""
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
    """Parse scalar namelist assignments from a SCALE-style config file."""
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = strip_comment(raw_line).strip()
        if not line or "=" not in line or line.startswith("&"):
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if value.endswith(",") and value.count(",") == 1:
            value = value[:-1].strip()
        values[key] = value
    return values


def clean_string(value: Any) -> str:
    """Normalize quoted namelist strings."""
    if value is None:
        return "NA"
    text = str(value).strip()
    if (text.startswith('"') and text.endswith('"')) or (text.startswith("'") and text.endswith("'")):
        text = text[1:-1]
    return text


def to_float(value: Any) -> float | None:
    """Parse Fortran D-notation numbers."""
    if value is None:
        return None
    text = clean_string(value).replace("D", "E").replace("d", "e")
    try:
        return float(text)
    except ValueError:
        return None


def to_int(value: Any) -> int | None:
    """Parse integer-like namelist values, including leading-zero flags."""
    number = to_float(value)
    if number is None:
        return None
    return int(round(number))


def to_bool(value: Any) -> bool | None:
    """Parse Fortran logical values."""
    text = clean_string(value).lower()
    if text in {".true.", "true", "t"}:
        return True
    if text in {".false.", "false", "f"}:
        return False
    return None


def value_for_row(conf: dict[str, str], key: str) -> str:
    """Return a display value for an audited key."""
    return clean_string(conf[key]) if key in conf else "NA"


def close_float(actual: float | None, expected: float, rel_tol: float = 1e-8, abs_tol: float = 1e-12) -> bool:
    """Compare floating-point values with tight tolerances for namelists."""
    return actual is not None and math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol)


class CaseAudit:
    """Accumulate failures and warnings for one model case."""

    def __init__(self, group: str, case: str, conf: dict[str, str], init_conf: dict[str, str]) -> None:
        self.group = group
        self.case = case
        self.conf = conf
        self.init_conf = init_conf
        self.failures: list[str] = []
        self.warnings: list[str] = []
        self.row: dict[str, Any] = {
            "case_group": group,
            "case_name": case,
            "status": "PASS",
            "failures": "",
            "warnings": "",
        }

    def require_key(self, key: str) -> bool:
        """Fail if a required key is absent."""
        if key in self.conf:
            return True
        self.failures.append(f"missing required key {key}")
        return False

    def require_float(self, key: str, expected: float) -> None:
        """Require a float setting."""
        if not self.require_key(key):
            return
        actual = to_float(self.conf.get(key))
        if not close_float(actual, expected):
            self.failures.append(f"{key}={clean_string(self.conf.get(key))} expected {expected}")

    def require_int(self, key: str, expected: int) -> None:
        """Require an integer setting."""
        if not self.require_key(key):
            return
        actual = to_int(self.conf.get(key))
        if actual != expected:
            self.failures.append(f"{key}={clean_string(self.conf.get(key))} expected {expected}")

    def require_bool(self, key: str, expected: bool) -> None:
        """Require a logical setting."""
        if not self.require_key(key):
            return
        actual = to_bool(self.conf.get(key))
        if actual is not expected:
            self.failures.append(f"{key}={clean_string(self.conf.get(key))} expected {expected}")

    def require_string(self, key: str, expected: str) -> None:
        """Require a string setting."""
        if not self.require_key(key):
            return
        actual = clean_string(self.conf.get(key))
        if actual != expected:
            self.failures.append(f'{key}="{actual}" expected "{expected}"')

    def finalize(self) -> dict[str, Any]:
        """Build the output row."""
        for key in AUDIT_KEYS:
            if key in {"dimension", "domain_size", "resolution", "TIME_DT_over_TIME_DT_ATMOS_DYN"}:
                continue
            self.row[key] = value_for_row(self.conf, key)

        dx = to_float(self.conf.get("DX"))
        dy = to_float(self.conf.get("DY"))
        dz = to_float(self.conf.get("DZ"))
        imax = to_int(self.conf.get("IMAX"))
        jmax = to_int(self.conf.get("JMAX"))
        kmax = to_int(self.conf.get("KMAX"))
        prc_x = to_int(self.conf.get("PRC_NUM_X"))
        prc_y = to_int(self.conf.get("PRC_NUM_Y"))
        if None not in {dx, dy, dz, imax, jmax, kmax, prc_x, prc_y}:
            assert dx is not None and dy is not None and dz is not None
            assert imax is not None and jmax is not None and kmax is not None
            assert prc_x is not None and prc_y is not None
            self.row["domain_size"] = f"{prc_x * imax * dx:g} x {prc_y * jmax * dy:g} x {kmax * dz:g} m"
            self.row["resolution"] = f"{dx:g} x {dy:g} x {dz:g} m"
        else:
            self.row["domain_size"] = "NA"
            self.row["resolution"] = "NA"

        dt = to_float(self.conf.get("TIME_DT"))
        dt_dyn = to_float(self.conf.get("TIME_DT_ATMOS_DYN"))
        self.row["TIME_DT_over_TIME_DT_ATMOS_DYN"] = "NA" if not dt or not dt_dyn else f"{dt / dt_dyn:.12g}"
        self.row["dimension"] = self.row.get("dimension", "NA")
        self.row["status"] = "FAIL" if self.failures else "WARN" if self.warnings else "PASS"
        self.row["failures"] = "; ".join(self.failures)
        self.row["warnings"] = "; ".join(self.warnings)
        return self.row


def model_cases(root: Path) -> list[Path]:
    """Return run.conf paths for model cases exactly two levels below GMD2026."""
    return sorted(
        path
        for path in root.glob("*/*/run.conf")
        if "postprocess" not in path.parts and path.parent.parent.name != "common"
    )


def apply_global_3d(audit: CaseAudit, spec: dict[str, Any]) -> None:
    """Apply required common 3D settings."""
    expected = spec["global_3d"]
    audit.row["dimension"] = "3d"
    audit.require_string("ATMOS_DYN_TYPE", expected["ATMOS_DYN_TYPE"])
    if clean_string(audit.conf.get("ATMOS_DYN_TYPE", "")).upper() == "HEVE":
        audit.failures.append('ATMOS_DYN_TYPE="HEVE" is rejected; use "HEVI"')
    audit.require_string("ATMOS_DYN_TINTEG_SHORT_TYPE", expected["ATMOS_DYN_TINTEG_SHORT_TYPE"])
    audit.require_float("TIME_DT", expected["TIME_DT"])
    audit.require_float("TIME_DT_ATMOS_DYN", expected["TIME_DT_ATMOS_DYN"])
    audit.require_float("TIME_DT_ATMOS_PHY_MP", expected["TIME_DT_ATMOS_PHY_MP"])
    dt = to_float(audit.conf.get("TIME_DT"))
    dt_dyn = to_float(audit.conf.get("TIME_DT_ATMOS_DYN"))
    if dt is None or dt_dyn in {None, 0.0} or not close_float(dt / dt_dyn, expected["TIME_DT_ratio"]):
        audit.failures.append("TIME_DT / TIME_DT_ATMOS_DYN must be 2")
    audit.require_float("ATMOS_DYN_wdamp_height", expected["ATMOS_DYN_wdamp_height"])
    audit.require_float("ATMOS_DYN_wdamp_tau", expected["ATMOS_DYN_wdamp_tau"])
    audit.require_bool("random_perturbation_enable", expected["random_perturbation_enable"])
    audit.require_float("random_perturbation_amp", expected["random_perturbation_amp"])

    dx = to_float(audit.conf.get("DX"))
    dy = to_float(audit.conf.get("DY"))
    dz = to_float(audit.conf.get("DZ"))
    imax = to_int(audit.conf.get("IMAX"))
    jmax = to_int(audit.conf.get("JMAX"))
    kmax = to_int(audit.conf.get("KMAX"))
    prc_x = to_int(audit.conf.get("PRC_NUM_X"))
    prc_y = to_int(audit.conf.get("PRC_NUM_Y"))
    values = [dx, dy, dz, imax, jmax, kmax, prc_x, prc_y]
    if any(value is None for value in values):
        audit.failures.append("missing grid/domain keys for 3D domain check")
        return
    assert dx is not None and dy is not None and dz is not None
    assert imax is not None and jmax is not None and kmax is not None
    assert prc_x is not None and prc_y is not None
    if [prc_x * imax, prc_y * jmax, kmax] != expected["grid_count"]:
        audit.failures.append(f"grid_count={[prc_x * imax, prc_y * jmax, kmax]} expected {expected['grid_count']}")
    domain = [prc_x * imax * dx, prc_y * jmax * dy, kmax * dz]
    if any(not close_float(actual, target) for actual, target in zip(domain, expected["domain_m"])):
        audit.failures.append(f"domain_m={domain} expected {expected['domain_m']}")
    resolution = [dx, dy, dz]
    if any(not close_float(actual, target) for actual, target in zip(resolution, expected["resolution_m"])):
        audit.failures.append(f"resolution_m={resolution} expected {expected['resolution_m']}")


def apply_expected_settings(audit: CaseAudit, expected: dict[str, Any]) -> None:
    """Apply typed expected namelist settings."""
    for key, value in expected.items():
        if key in {
            "dimension",
            "restart_policy",
            "cases",
            "sampled_fractions",
            "sampled_modes",
            "seed_count_per_mode_fraction",
            "seed_key",
            "seed_namelist",
            "sdnc_values",
            "output_intervals",
        }:
            continue
        if isinstance(value, bool):
            audit.require_bool(key, value)
        elif isinstance(value, int):
            audit.require_int(key, value)
        elif isinstance(value, float):
            audit.require_float(key, value)
        elif isinstance(value, str):
            audit.require_string(key, value)


def apply_03_checks(audit: CaseAudit, group_spec: dict[str, Any]) -> tuple[str, float, int] | None:
    """Check one representativeness case and return sampled grouping info."""
    audit.row["dimension"] = "2d"
    apply_expected_settings(audit, group_spec)
    audit.require_string("ATMOS_RESTART_IN_BASENAME", "./init_dycom_rf02_00000101-000000.000")
    audit.require_bool("ATMOS_RESTART_OUTPUT", False)
    if "RANDOM_NUMBER_SEED" in audit.init_conf:
        audit.failures.append("representativeness cases must not change PARAM_SDMRANDOM/RANDOM_NUMBER_SEED")
    if "RANDOM_SEED_SCALE" in audit.conf:
        audit.failures.append("representativeness cases must not use PARAM_RANDOM/RANDOM_SEED_SCALE")
    if "tracking_sampling_seed" not in audit.conf:
        audit.failures.append("representativeness case missing tracking_sampling_seed")
    run_conf_path = GMD_ROOT / audit.group / audit.case / "run.conf"
    if run_conf_path.exists() and "&HISTITEM" in run_conf_path.read_text(errors="ignore"):
        audit.failures.append("representativeness cases must not contain HISTITEM blocks when history output is disabled")
    case = audit.case
    if case.startswith("ref_full_"):
        mode = case.replace("ref_full_", "")
        audit.require_string("tracking_selection_mode", mode)
        audit.require_float("tracking_fraction", 1.0)
        return None
    match = re.match(r"^(stratified|random)_f(\d{3})_s(\d{2})$", case)
    if not match:
        audit.failures.append("representativeness case name does not match expected pattern")
        return None
    mode, fraction_code, seed_code = match.groups()
    fraction = int(fraction_code) / 100.0
    seed_index = int(seed_code)
    audit.require_string("tracking_selection_mode", mode)
    audit.require_float("tracking_fraction", fraction)
    if fraction not in group_spec["sampled_fractions"]:
        audit.failures.append(f"tracking_fraction {fraction} not in expected sampled fractions")
    return (mode, fraction, seed_index)


def apply_04_checks(audit: CaseAudit, group_spec: dict[str, Any]) -> None:
    """Check one optional SDNC scaling fallback case."""
    apply_expected_settings(
        audit,
        {
            "TIME_DURATION": group_spec["TIME_DURATION"],
            "HISTORY_DEFAULT_TINTERVAL": group_spec["HISTORY_DEFAULT_TINTERVAL"],
            "coalescence_output_enable": group_spec["coalescence_output_enable"],
        },
    )
    match = re.match(r"^sdnc(\d+)_(nt|fw005|bw005)$", audit.case)
    if not match:
        audit.failures.append("SDNC scaling case name does not match expected pattern")
        return
    sdnc = float(match.group(1))
    mode_name = match.group(2)
    audit.require_float("sdm_inisdnc", sdnc)
    if sdnc not in group_spec["sdnc_values"]:
        audit.failures.append(f"sdm_inisdnc {sdnc} not in expected SDNC values")
    expected_mode = {"nt": 0, "fw005": 1, "bw005": 2}[mode_name]
    audit.require_int("tracking_mode", expected_mode)
    if mode_name == "nt":
        audit.require_int("sdm_dmpvar", 0)
    else:
        audit.require_int("sdm_dmpvar", 100)
        audit.require_string("tracking_selection_mode", group_spec["tracking_selection_mode"])
        audit.require_float("tracking_fraction", group_spec["tracking_fraction"])
    audit.require_float("tracking_height_max", 800.0)


def apply_05_checks(audit: CaseAudit, group_spec: dict[str, Any]) -> None:
    """Check one optional output-interval scaling case."""
    apply_expected_settings(
        audit,
        {
            "TIME_DURATION": group_spec["TIME_DURATION"],
            "HISTORY_DEFAULT_TINTERVAL": group_spec["HISTORY_DEFAULT_TINTERVAL"],
            "coalescence_output_enable": group_spec["coalescence_output_enable"],
        },
    )
    if audit.case == "nt_nolog_10min":
        apply_expected_settings(audit, group_spec.get("cases", {}).get(audit.case, {}))
        audit.require_float("tracking_fraction", group_spec["tracking_fraction"])
        audit.require_float("tracking_height_max", 800.0)
        return
    match = re.match(r"^out(\d+)_(fw005|bw005)$", audit.case)
    if not match:
        audit.failures.append("output-interval case name does not match expected pattern")
        return
    interval = float(match.group(1))
    mode_name = match.group(2)
    expected_mode = {"fw005": 1, "bw005": 2}[mode_name]
    audit.require_int("tracking_mode", expected_mode)
    audit.require_float("sdm_dmpitvl", interval)
    audit.require_int("sdm_dmpvar", 100)
    audit.require_string("tracking_selection_mode", group_spec["tracking_selection_mode"])
    audit.require_float("tracking_fraction", group_spec["tracking_fraction"])
    audit.require_float("tracking_height_max", 800.0)
    if interval not in group_spec["output_intervals"]:
        audit.failures.append(f"output interval {interval} not in expected output intervals")


def write_table(rows: list[dict[str, Any]], path: Path) -> None:
    """Write CSV, JSON, and Markdown audit outputs."""
    csv_path = path.with_suffix(".csv")
    json_path = path.with_suffix(".json")
    md_path = path.with_suffix(".md")
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(rows, indent=2))
    md_lines = [
        "| " + " | ".join(COLUMNS) + " |",
        "| " + " | ".join(["---"] * len(COLUMNS)) + " |",
    ]
    for row in rows:
        md_lines.append("| " + " | ".join(str(row.get(column, "NA")) for column in COLUMNS) + " |")
    md_path.write_text("\n".join(md_lines) + "\n")


def main() -> int:
    spec = json.loads(SPEC_PATH.read_text())
    rows: list[dict[str, Any]] = []
    audits_by_case: dict[tuple[str, str], CaseAudit] = {}
    sampled_cases: defaultdict[tuple[str, float], list[tuple[int, CaseAudit, str]]] = defaultdict(list)

    for run_conf in model_cases(GMD_ROOT):
        case_dir = run_conf.parent
        group = case_dir.parent.name
        case = case_dir.name
        group_spec = spec["groups"].get(group)
        if group_spec is None:
            continue
        audit = CaseAudit(group, case, parse_conf(run_conf), parse_conf(case_dir / "init.conf"))
        audits_by_case[(group, case)] = audit
        if group_spec.get("dimension") == "3d":
            apply_global_3d(audit, spec)

        if group == "00_base_restart_3d":
            apply_expected_settings(audit, group_spec["cases"][case])
        elif group == "01_bench_3d_samp_30min":
            apply_expected_settings(
                audit,
                {
                    "TIME_DURATION": group_spec["TIME_DURATION"],
                    "HISTORY_DEFAULT_TINTERVAL": group_spec["HISTORY_DEFAULT_TINTERVAL"],
                },
            )
            apply_expected_settings(audit, group_spec["cases"].get(case, {}))
            audit.require_string("ATMOS_RESTART_IN_BASENAME", spec["initial_restart_basename"])
        elif group == "02_tpht_3d_interest_70min":
            apply_expected_settings(audit, {"TIME_DURATION": group_spec["TIME_DURATION"]})
            apply_expected_settings(audit, group_spec["cases"].get(case, {}))
            audit.require_string("ATMOS_RESTART_IN_BASENAME", spec["initial_restart_basename"])
            if case == "bw_reconstruction":
                input_basename = clean_string(audit.conf.get("tracking_id_input_basename", ""))
                if input_basename:
                    selection = clean_string(audit.conf.get("tracking_selection_mode", "NA"))
                    fraction = to_float(audit.conf.get("tracking_fraction"))
                    if selection not in {"", "NA", "none"} or (fraction is not None and not close_float(fraction, 1.0)):
                        audit.warnings.append("tracking_id_input_basename should take priority over sampling settings")
                else:
                    audit.failures.append("BW reconstruction must define tracking_id_input_basename and fail if .ids are absent")
        elif group == "03_fw_rep_2d_600s":
            sampled_info = apply_03_checks(audit, group_spec)
            if sampled_info is not None:
                mode, fraction, seed_index = sampled_info
                seed_value = clean_string(audit.conf.get("tracking_sampling_seed", "NA"))
                sampled_cases[(mode, fraction)].append((seed_index, audit, seed_value))
        elif group == "04_sdnc_scaling_lite_30min":
            apply_04_checks(audit, group_spec)
        elif group == "05_outint_io_lite_10min":
            apply_05_checks(audit, group_spec)
        else:
            audit.failures.append("unknown GMD2026 group")

        if group.startswith(("00_", "01_", "02_", "03_", "04_", "05_")):
            audit.require_bool("gmd_benchmark_diag_enable", True)

        history_interval = to_float(audit.conf.get("HISTORY_DEFAULT_TINTERVAL"))
        history_step0 = to_bool(audit.conf.get("HISTORY_OUTPUT_STEP0"))
        if history_interval is not None and history_interval > 0.0 and history_step0 is not True:
            audit.failures.append("History-output cases must set HISTORY_OUTPUT_STEP0=.true.")
        if history_interval is not None and history_interval <= 0.0 and history_step0 is True:
            audit.failures.append("History-disabled cases should not request step-0 History output")

    rep_spec = spec["groups"]["03_fw_rep_2d_600s"]
    for mode in rep_spec["sampled_modes"]:
        for fraction in rep_spec["sampled_fractions"]:
            entries = sampled_cases.get((mode, fraction), [])
            if len(entries) != rep_spec["seed_count_per_mode_fraction"]:
                message = f"{mode} fraction {fraction} has {len(entries)} seeds; expected {rep_spec['seed_count_per_mode_fraction']}"
                for _seed_index, audit, _seed_value in entries:
                    audit.failures.append(message)
                if not entries:
                    placeholder = CaseAudit("03_fw_rep_2d_600s", f"{mode}_f{int(fraction * 100):03d}_missing", {}, {})
                    placeholder.failures.append(message)
                    audits_by_case[(placeholder.group, placeholder.case)] = placeholder
            seed_indices = [seed_index for seed_index, _audit, _seed_value in entries]
            seed_values = [seed_value for _seed_index, _audit, seed_value in entries]
            if len(seed_indices) != len(set(seed_indices)) or len(seed_values) != len(set(seed_values)):
                for _seed_index, audit, _seed_value in entries:
                    audit.failures.append(f"{mode} fraction {fraction} has duplicate deterministic seeds")

    for key in sorted(audits_by_case):
        rows.append(audits_by_case[key].finalize())

    write_table(rows, GMD_ROOT / "config_audit")
    pass_count = sum(row["status"] == "PASS" for row in rows)
    fail_count = sum(row["status"] == "FAIL" for row in rows)
    warn_count = sum(row["status"] == "WARN" for row in rows)
    print(f"PASS: {pass_count} cases")
    print(f"FAIL: {fail_count} cases")
    print(f"WARN: {warn_count} cases")
    return 1 if fail_count else 0


if __name__ == "__main__":
    sys.exit(main())
