#!/usr/bin/env python3
"""Create SQUID-first smoke-test cases for the GMD diagnostics hook.

The generated cases are lightweight engineering tests, not manuscript
experiments.  They use one shared initialization so hooked/unhooked model runs
start from exactly the same SD population and random-number state.
"""

from __future__ import annotations

import re
import shutil
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GMD_ROOT = ROOT.parent

CASE_TEMPLATES = {
    "shared_init_hook_2min": GMD_ROOT / "01_bench_3d_samp_30min" / "nt_nolog",
    "nt_nolog_hook_2min": GMD_ROOT / "01_bench_3d_samp_30min" / "nt_nolog",
    "fw_all_hook_2min": GMD_ROOT / "01_bench_3d_samp_30min" / "fw005_coallog",
    "fw_all_unhooked_2min": GMD_ROOT / "01_bench_3d_samp_30min" / "fw005_coallog",
    "bw_selected_hook_2min": GMD_ROOT / "01_bench_3d_samp_30min" / "bw005_coallog",
    "tpht_ids_hook_2min": GMD_ROOT / "02_tpht_3d_interest_70min" / "fw_discovery",
    "tpht_bw_hook_2min": GMD_ROOT / "02_tpht_3d_interest_70min" / "bw_reconstruction",
    "tpht_bw_missing_ids_failfast": GMD_ROOT / "02_tpht_3d_interest_70min" / "bw_reconstruction",
}

CASE_ORDER = tuple(CASE_TEMPLATES)

OUTPUT_PATTERNS = (
    "LOG.pe*",
    "*.o*",
    "*.e*",
    "job_metrics.json",
    "SD_selected_NetCDF_*",
    "SD_all_NetCDF_*",
    "SD_coal_output_NetCDF_*",
    "history*.nc",
    "time_scale_rm*.log",
    "fw_tracking",
    "fw_output",
    "bw_output",
    "restart_output",
    "dycom_restart_rf02*",
    "init_dycom_rf02*",
    "random_number_init.pe*",
    "random_number_output.pe*",
)


def ignore_generated(_directory: str, names: list[str]) -> set[str]:
    """Skip generated outputs and forbidden helper files when copying templates."""
    ignored: set[str] = {"particle_tracer_opt.f90", "random_traj.py", "__pycache__"}
    for name in names:
        for pattern in OUTPUT_PATTERNS:
            if Path(name).match(pattern):
                ignored.add(name)
                break
    return ignored


def replace_assignment(text: str, key: str, value: str) -> str:
    """Replace one scalar namelist assignment while preserving comments."""
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^!\n]*?)(\s*(?:!.*)?$)", re.MULTILINE)

    def repl(match: re.Match[str]) -> str:
        old_value = match.group(2).rstrip()
        comma = "," if old_value.endswith(",") else ""
        return f"{match.group(1)}{value}{comma}{match.group(3)}"

    new_text, count = pattern.subn(repl, text, count=1)
    if count == 0:
        raise KeyError(f"Missing required namelist key: {key}")
    return new_text


def replace_inline_comment(text: str, key: str, comment: str) -> str:
    """Replace the inline comment for one namelist assignment."""
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*[^!\n]*)(?:\s*!.*)?$", re.MULTILINE)
    new_text, count = pattern.subn(rf"\1 ! {comment}", text, count=1)
    if count == 0:
        raise KeyError(f"Missing required namelist key: {key}")
    return new_text


def set_assignment(text: str, key: str, value: str) -> str:
    """Replace an assignment if present, otherwise append it to the SDM namelist."""
    try:
        return replace_assignment(text, key, value)
    except KeyError:
        insertion = f"{key} = {value},\n"
        marker = re.search(r"^\s*coalescence_output_enable\s*=", text, flags=re.MULTILINE)
        if marker:
            return text[: marker.start()] + insertion + text[marker.start() :]
        return re.sub(r"^\s*/\s*$", insertion + "/", text, count=1, flags=re.MULTILINE)


def patch_common_run_conf(text: str) -> str:
    """Apply common small-domain, no-history smoke-test settings."""
    replacements = {
        "PRC_NUM_X": "2",
        "PRC_NUM_Y": "2",
        "IMAX": "5",
        "JMAX": "5",
        "KMAX": "300",
        "TIME_DURATION": "120.0D0",
        "TIME_DT_ATMOS_RESTART": "120.0D0",
        "ATMOS_RESTART_IN_BASENAME": '"./init_dycom_rf02_00000101-000000.000"',
        "ATMOS_RESTART_OUTPUT": ".false.",
        "sdm_dmpitva": "10.0d0",
        "sdm_dmpitvb": "10.0d0",
        "sdm_dmpitvl": "10.0d0",
        "tracking_height_min": "0.0d0",
        "tracking_height_max": "800.0d0",
        "tracking_radius_min": "1.0d-9",
        "HISTORY_DEFAULT_TINTERVAL": "0.0D0",
        "HISTORY_OUTPUT_STEP0": ".false.",
    }
    for key, value in replacements.items():
        text = replace_assignment(text, key, value)
    text = set_assignment(text, "gmd_benchmark_diag_enable", ".true.")
    text = text.replace(
        "               ! when sampling is enabled (tracking_fraction < 1 or max_tracked_sds > 0),\n"
        "               ! set sdm_dmpvar = 100 to output sampled trajectories into SD_selected_NetCDF_*\n",
        "               ! sdm_dmpvar=010 writes SD_all_NetCDF_*; sdm_dmpvar=100 writes SD_selected_NetCDF_*.\n",
    )
    text = replace_inline_comment(text, "sdm_dmpitva", "Text-output interval; text output remains disabled in these smoke cases")
    text = replace_inline_comment(text, "sdm_dmpitvb", "SD_all_NetCDF_* interval when sdm_dmpvar=010")
    text = replace_inline_comment(text, "sdm_dmpitvl", "SD_selected_NetCDF_* interval when sdm_dmpvar=100")
    text = replace_inline_comment(text, "tracking_height_min", "Height filter for random/stratified sampling; ignored when selection_mode='none'")
    text = replace_inline_comment(text, "tracking_height_max", "Height filter for random/stratified sampling; ignored when selection_mode='none'")
    text = replace_inline_comment(text, "tracking_radius_min", "Radius filter for random/stratified sampling; ignored when selection_mode='none'")
    text = replace_inline_comment(text, "tracking_fraction", "Sampling fraction; ignored when tracking_selection_mode='none'")
    text = replace_inline_comment(text, "max_tracked_sds", "Sampling cap; ignored when tracking_selection_mode='none'")
    return text


def patch_common_init_conf(text: str) -> str:
    """Match the four-rank smoke-test layout in init.conf."""
    replacements = {
        "PRC_NUM_X": "2",
        "PRC_NUM_Y": "2",
        "IMAX": "5",
        "JMAX": "5",
        "KMAX": "300",
    }
    for key, value in replacements.items():
        text = replace_assignment(text, key, value)
    return text


def patch_makefile(text: str) -> str:
    """Use the generated four-rank decomposition."""
    return re.sub(r"^TPROC\s*=.*$", "TPROC     = 4", text, flags=re.MULTILINE)


def patch_tracking_defaults(text: str) -> str:
    """Normalize optional tracking keys across templates."""
    text = set_assignment(text, "tracking_id_output_basename", '""')
    text = set_assignment(text, "tracking_id_input_basename", '""')
    text = set_assignment(text, "tracking_interest_radius_enable", ".false.")
    text = set_assignment(text, "tracking_interest_coalescence_enable", ".false.")
    return text


def write_squid_run(case_dir: Path, case_name: str, run_init: bool, expect_failure: bool) -> None:
    """Write a compact SQUID job script that delegates to the common runner."""
    run_init_value = "1" if run_init else "0"
    run_main_value = "0" if run_init else "1"
    use_shared_init = "0" if run_init else "1"
    shared_init_dir = "" if run_init else "../shared_init_hook_2min"
    expect_failure_value = "1" if expect_failure else "0"
    script = f"""\
#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=4
#PBS -l elapstim_req=00:30:00
#PBS -T intmpi
#PBS -N {case_name}
set -euo pipefail

cd "${{PBS_O_WORKDIR:-$(pwd)}}"
export GMD_CASE_GROUP="06_fortran_diagnostics_hook_smoke_squid"
export GMD_CASE_NAME="{case_name}"
export GMD_RUN_INIT={run_init_value}
export GMD_RUN_MAIN={run_main_value}
export GMD_USE_SHARED_INIT={use_shared_init}
export GMD_SHARED_INIT_DIR="{shared_init_dir}"
export GMD_EXPECT_FAILURE={expect_failure_value}
source ../../common/squid_model_job.sh
"""
    (case_dir / "squid_run.sh").write_text(textwrap.dedent(script))
    (case_dir / "squid_run.sh").chmod(0o755)


def patch_case(case_dir: Path, case_name: str) -> None:
    """Patch copied template files for one smoke-test case."""
    run_conf = case_dir / "run.conf"
    init_conf = case_dir / "init.conf"
    makefile = case_dir / "Makefile"

    run_text = patch_tracking_defaults(patch_common_run_conf(run_conf.read_text()))
    expect_failure = case_name == "tpht_bw_missing_ids_failfast"

    if case_name == "shared_init_hook_2min":
        run_text = replace_assignment(run_text, "tracking_mode", "0")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"none"')
        run_text = replace_assignment(run_text, "sdm_dmpvar", "000")
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".false.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "No SD dump is needed during the shared initialization run")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Inactive because tracking_mode=0")
    elif case_name == "nt_nolog_hook_2min":
        run_text = replace_assignment(run_text, "tracking_mode", "0")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"none"')
        run_text = replace_assignment(run_text, "sdm_dmpvar", "000")
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".false.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "No SD dump in the no-tracking baseline")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Inactive because tracking_mode=0")
    elif case_name in {"fw_all_hook_2min", "fw_all_unhooked_2min"}:
        run_text = replace_assignment(run_text, "tracking_mode", "1")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"none"')
        run_text = replace_assignment(run_text, "tracking_fraction", "1.0d0")
        run_text = replace_assignment(run_text, "max_tracked_sds", "0")
        run_text = replace_assignment(run_text, "sdm_dmpvar", "010")
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".true.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "All valid SDs are written to SD_all_NetCDF_* for selection_mode='none'")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Track all valid SDs; height/radius/fraction/max filters are ignored")
        if case_name == "fw_all_unhooked_2min":
            run_text = replace_assignment(run_text, "gmd_benchmark_diag_enable", ".false.")
    elif case_name == "bw_selected_hook_2min":
        run_text = replace_assignment(run_text, "tracking_mode", "2")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"stratified"')
        run_text = replace_assignment(run_text, "tracking_fraction", "0.05d0")
        run_text = replace_assignment(run_text, "max_tracked_sds", "0")
        run_text = replace_assignment(run_text, "sdm_dmpvar", "100")
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".true.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "Selected BW targets are written to SD_selected_NetCDF_*")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Stratified sampling defines the BW target set in this smoke case")
    elif case_name == "tpht_ids_hook_2min":
        run_text = replace_assignment(run_text, "tracking_mode", "1")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"none"')
        run_text = replace_assignment(run_text, "tracking_fraction", "1.0d0")
        run_text = replace_assignment(run_text, "max_tracked_sds", "0")
        run_text = replace_assignment(run_text, "sdm_dmpvar", "000")
        run_text = replace_assignment(run_text, "tracking_id_output_basename", '"./fw_tracking/tracking_interest_ids"')
        run_text = replace_assignment(run_text, "tracking_interest_radius_enable", ".true.")
        run_text = replace_assignment(run_text, "tracking_interest_radius_threshold", "1.0d-6")
        run_text = replace_assignment(run_text, "tracking_interest_coalescence_enable", ".true.")
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".true.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "No SD NetCDF dump; TPHT FW writes interest IDs through tracking_id_output_basename")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Track all valid SDs before applying TPHT interest filters")
    elif case_name == "tpht_bw_hook_2min":
        run_text = replace_assignment(run_text, "tracking_mode", "2")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"none"')
        run_text = replace_assignment(run_text, "tracking_fraction", "1.0d0")
        run_text = replace_assignment(run_text, "max_tracked_sds", "0")
        run_text = replace_assignment(run_text, "sdm_dmpvar", "100")
        run_text = replace_assignment(run_text, "tracking_id_input_basename", '"../tpht_ids_hook_2min/fw_tracking/tracking_interest_ids_dedup"')
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".true.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "TPHT BW target set is written to SD_selected_NetCDF_*")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Ignored because tracking_id_input_basename defines BW targets")
    elif case_name == "tpht_bw_missing_ids_failfast":
        run_text = replace_assignment(run_text, "tracking_mode", "2")
        run_text = replace_assignment(run_text, "tracking_selection_mode", '"none"')
        run_text = replace_assignment(run_text, "tracking_fraction", "1.0d0")
        run_text = replace_assignment(run_text, "max_tracked_sds", "0")
        run_text = replace_assignment(run_text, "sdm_dmpvar", "100")
        run_text = replace_assignment(run_text, "tracking_id_input_basename", '"./missing_ids/tracking_interest_ids_dedup"')
        run_text = replace_assignment(run_text, "coalescence_output_enable", ".true.")
        run_text = replace_inline_comment(run_text, "sdm_dmpvar", "Negative-control TPHT BW case should fail before producing selected output")
        run_text = replace_inline_comment(run_text, "tracking_selection_mode", "Ignored because tracking_id_input_basename should define BW targets")
    else:
        raise ValueError(case_name)

    run_conf.write_text(run_text)
    init_conf.write_text(patch_common_init_conf(init_conf.read_text()))
    makefile.write_text(patch_makefile(makefile.read_text()))
    write_squid_run(case_dir, case_name, run_init=(case_name == "shared_init_hook_2min"), expect_failure=expect_failure)


def main() -> None:
    """Generate all smoke-test model directories."""
    for case_name, template in CASE_TEMPLATES.items():
        destination = ROOT / case_name
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(template, destination, ignore=ignore_generated)
        patch_case(destination, case_name)
        (destination / "logs").mkdir(exist_ok=True)
        if case_name == "tpht_ids_hook_2min":
            (destination / "fw_tracking").mkdir(exist_ok=True)
            (destination / "fw_output").mkdir(exist_ok=True)
        if "bw" in case_name:
            (destination / "bw_output").mkdir(exist_ok=True)
    print("Prepared hook smoke-test cases:")
    for case_name in CASE_ORDER:
        print(f"  - {case_name}")


if __name__ == "__main__":
    main()
