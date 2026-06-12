#!/usr/bin/env python3
"""Generate the 600 s 2D forward-sampling representativeness cases.

The checked-in suite already contains the generated cases. This script is kept
as a reproducible case generator for future fractions or seeds. It copies a
local generated reference case, skips forbidden legacy helper files, updates
the tracking namelist settings, and preserves the GMD2026 Makefile depth.
"""

from __future__ import annotations

import argparse
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


FORBIDDEN = {"particle_tracer_opt.f90", "random_traj.py"}
OUTPUT_PATTERNS = {
    "LOG.pe*",
    "*.o*",
    "*.e*",
    "job_metrics.json",
    "SD_selected_NetCDF_*",
    "SD_all_NetCDF_*",
    "SD_coal_output_NetCDF_*",
    "history*.nc",
}


@dataclass(frozen=True)
class CaseSpec:
    """One representativeness case definition."""

    name: str
    selection_mode: str
    tracking_fraction: float
    seed: int


def build_specs() -> list[CaseSpec]:
    """Return the fixed GMD2026 representativeness case list."""
    specs = [
        CaseSpec("ref_full_stratified", "stratified", 1.0, 1001),
        CaseSpec("ref_full_random", "random", 1.0, 2001),
    ]
    for selection_mode, seed_base in (("stratified", 10000), ("random", 20000)):
        for label, fraction in (("f005", 0.05), ("f010", 0.10), ("f020", 0.20)):
            for seed_index in range(10):
                seed = seed_base + int(fraction * 1000) * 10 + seed_index
                specs.append(CaseSpec(f"{selection_mode}_{label}_s{seed_index:02d}", selection_mode, fraction, seed))
    return specs


def should_ignore(path: Path) -> bool:
    """Return True for forbidden or generated files that should not be copied."""
    if path.name in FORBIDDEN:
        return True
    for pattern in OUTPUT_PATTERNS:
        if path.match(pattern):
            return True
    return path.name == "__pycache__"


def copy_template(template_dir: Path, target_dir: Path) -> None:
    """Copy a template case while skipping forbidden and generated files."""
    for source in template_dir.rglob("*"):
        relative = source.relative_to(template_dir)
        if any(should_ignore(part) for part in [source, *source.parents]):
            continue
        target = target_dir / relative
        if source.is_dir():
            target.mkdir(parents=True, exist_ok=True)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def replace_assignment(text: str, key: str, value: str) -> str:
    """Replace a simple Fortran namelist assignment."""
    pattern = re.compile(r"^(\s*" + re.escape(key) + r"\s*=\s*)(.*?)(,?\s*(?:!.*)?)$", re.MULTILINE)
    if not pattern.search(text):
        raise ValueError(f"Missing assignment: {key}")
    return pattern.sub(lambda match: f"{match.group(1)}{value}{match.group(3)}", text, count=1)


def configure_run_conf(run_conf: Path, spec: CaseSpec) -> None:
    """Update the run.conf tracking settings for one case."""
    text = run_conf.read_text()
    replacements = {
        "tracking_mode": "1",
        "TIME_DURATION": "600.0D0",
        "TIME_DT_ATMOS_RESTART": "600.0D0",
        "tracking_selection_mode": f'"{spec.selection_mode}"',
        "tracking_fraction": "1.0d0" if spec.tracking_fraction >= 1.0 else f"{spec.tracking_fraction:.2f}d0",
        "max_tracked_sds": "0",
        "tracking_sampling_seed": str(spec.seed),
        "sdm_dmpvar": "100",
        "sdm_dmpitvl": "60.0d0",
        "HISTORY_DEFAULT_TINTERVAL": "0.0D0",
        "HISTORY_OUTPUT_STEP0": ".false.",
        "coalescence_output_enable": ".true.",
        "gmd_benchmark_diag_enable": ".true.",
        "random_perturbation_enable": ".false.",
        "random_perturbation_amp": "0.0d0",
    }
    for key, value in replacements.items():
        if key == "gmd_benchmark_diag_enable" and key not in text:
            text = re.sub(
                r"(coalescence_output_enable\s*=\s*\.(?:true|false)\.,[^\n]*\n)",
                r"\1gmd_benchmark_diag_enable = .true., ! Emit GMD_BENCH_DIAG/GMD_IO_DIAG for benchmark post-processing\n",
                text,
                count=1,
                flags=re.IGNORECASE,
            )
        else:
            text = replace_assignment(text, key, value)
    text = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("&HISTITEM")) + "\n"
    run_conf.write_text(text)


def write_squid_run(case_dir: Path, case_name: str) -> None:
    """Write a case-local SQUID wrapper using shared init where appropriate."""
    use_shared_init = case_name != "ref_full_stratified"
    run_init = "0" if use_shared_init else "1"
    shared_init = "../ref_full_stratified" if use_shared_init else ""
    use_shared = "1" if use_shared_init else "0"
    script = f"""#!/usr/bin/env bash
#PBS -q SQUID
#PBS --group=hp250136
#PBS -m b
#PBS -b 1
#PBS -l cpunum_job=20
#PBS -l elapstim_req=02:00:00
#PBS -T intmpi
#PBS -N {case_name}
set -euo pipefail

cd "${{PBS_O_WORKDIR:-$(pwd)}}"
export GMD_CASE_GROUP="03_fw_rep_2d_600s"
export GMD_CASE_NAME="{case_name}"
export GMD_RUN_INIT={run_init}
export GMD_RUN_MAIN=1
export GMD_USE_SHARED_INIT={use_shared}
export GMD_SHARED_INIT_DIR="{shared_init}"
export GMD_EXPECT_FAILURE=0
source ../../common/squid_model_job.sh
"""
    target = case_dir / "squid_run.sh"
    target.write_text(script)
    target.chmod(0o755)


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="Replace existing generated case directories.")
    parser.add_argument(
        "--template",
        type=Path,
        default=Path("ref_full_stratified"),
        help="Template case directory relative to this group.",
    )
    return parser.parse_args()


def main() -> None:
    """Generate all 2D representativeness cases."""
    args = parse_args()
    group_dir = Path(__file__).resolve().parent
    template_dir = args.template if args.template.is_absolute() else group_dir / args.template
    if not template_dir.is_dir():
        raise SystemExit(f"Template directory not found: {template_dir}")

    for spec in build_specs():
        target = group_dir / spec.name
        if target.exists():
            if not args.overwrite:
                print(f"skip existing {spec.name}")
                continue
            shutil.rmtree(target)
        copy_template(template_dir, target)
        configure_run_conf(target / "run.conf", spec)
        write_squid_run(target, spec.name)
        print(f"generated {spec.name}")


if __name__ == "__main__":
    main()
