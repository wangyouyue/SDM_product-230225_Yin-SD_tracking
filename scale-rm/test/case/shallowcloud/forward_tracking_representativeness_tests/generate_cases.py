from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil


@dataclass(frozen=True)
class CaseSpec:
    case_name: str
    mode: str
    fraction: float
    seed_scale: int


ROOT = Path(__file__).resolve().parent
SHALLOWCLOUD_ROOT = ROOT.parent
TEMPLATE = SHALLOWCLOUD_ROOT / "forward_tracking_sampling_tests" / "ft_stratified_baseline"
FRACTIONS = (0.05, 0.10, 0.20)
MODES = ("stratified", "random")
SEEDS_PER_GROUP = 10


def _replace_key_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(^\s*{re.escape(key)}\s*=\s*)([^,\n]+)(\s*,.*$)", re.MULTILINE)
    return pattern.sub(rf"\g<1>{value}\g<3>", text)


def _ensure_param_random(text: str, seed_scale: int) -> str:
    """Inject PARAM_RANDOM into run.conf for reproducible random seeds."""
    block = (
        "&PARAM_RANDOM\n"
        " RANDOM_FIX = .true.,\n"
        f" RANDOM_SEED_SCALE = {seed_scale},\n"
        "/\n\n"
    )
    if "&PARAM_RANDOM" in text:
        return re.sub(r"&PARAM_RANDOM[\s\S]*?/\n", block, text, count=1)
    match = re.search(r"&PARAM_TIME[\s\S]*?/\n", text)
    if match:
        return text[: match.end()] + "\n" + block + text[match.end() :]
    return block + text


def _apply_runconf(path: Path, mode: str, fraction: float, seed_scale: int) -> None:
    """Apply tracking controls and random seed to a copied run.conf."""
    text = path.read_text()
    text = re.sub(r"^\s*tracking_sample_initialized\s*=.*$\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\s*coalescence_output_enable\s*=.*$\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\s*random_perturbation_enable\s*=.*$\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\s*random_perturbation_amp\s*=.*$\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\s*tracked_grid_coal_only\s*=.*$\n?", "", text, flags=re.MULTILINE)
    values = {
        "tracking_mode": "1",
        "tracking_selection_mode": f'"{mode}"',
        "tracking_fraction": "1.d0" if fraction == 1.0 else f"{fraction:.2f}d0",
        "max_tracked_sds": "0",
        "tracking_height_min": "400.d0",
        "tracking_height_max": "800.d0",
        "tracking_radius_min": "0.d0",
        "tracking_radius_max": "0.d0",
        "tracking_nz_bin": "8",
        "tracking_nr_bin": "10",
        "tracking_min_per_bin": "1",
        "tracking_fallback_to_random": ".true.",
        "sdm_noise_amp": "0.d0",
    }
    for key, value in values.items():
        text = _replace_key_value(text, key, value)
    text = _ensure_param_random(text, seed_scale)
    text = _remove_tracked_radius_threshold(text)
    path.write_text(text)


def _apply_job_name(path: Path, name: str) -> None:
    if not path.exists():
        return
    text = path.read_text()
    text = re.sub(r"^#PBS -N .*$", f"#PBS -N {name}", text, flags=re.MULTILINE)
    path.write_text(text)


def _rewrite_tracking_comment_block(path: Path) -> None:
    text = path.read_text()
    lines = text.splitlines(keepends=True)
    start = None
    end = None
    for idx, line in enumerate(lines):
        if re.match(r"^\s*(tracking_mode|tracking_fraction)\s*=", line):
            start = idx
            break
    for idx, line in enumerate(lines):
        if re.match(r"^\s*tracking_fallback_to_random\s*=", line):
            end = idx
            break
    if start is None or end is None or end < start:
        return

    def _get(key: str, default: str) -> str:
        match = re.search(rf"^\s*{key}\s*=\s*([^,\n]+)", text, flags=re.MULTILINE)
        return match.group(1).strip() if match else default

    mode = _get("tracking_selection_mode", '"stratified"').strip('"')
    block = (
        "tracking_mode = 1,                ! Tracking mode: 0=no tracking, 1=forward tracking, 2=backward tracking\n"
        f"tracking_fraction = {_get('tracking_fraction', '0.10d0')},          ! Fraction of SDs to track (0-1]\n"
        "max_tracked_sds = 0,                ! Maximum number of tracked SDs (0 = unlimited)\n"
        f"tracking_selection_mode = \"{mode}\", ! \"random\" or \"stratified\"\n"
        f"tracking_height_min = {_get('tracking_height_min', '400.d0')},       ! Minimum height for stratified selection [m]\n"
        f"tracking_height_max = {_get('tracking_height_max', '800.d0')},       ! Maximum height for stratified selection [m]\n"
        f"tracking_radius_min = {_get('tracking_radius_min', '0.d0')},         ! Minimum radius for stratified selection [m]\n"
        f"tracking_radius_max = {_get('tracking_radius_max', '0.d0')},         ! Maximum radius for stratified selection [m], <= tracking_radius_min uses runtime maximum radius\n"
        f"tracking_nz_bin = {_get('tracking_nz_bin', '8')},                ! Number of vertical bins for stratified selection\n"
        f"tracking_nr_bin = {_get('tracking_nr_bin', '10')},               ! Number of radius bins for stratified selection\n"
        f"tracking_min_per_bin = {_get('tracking_min_per_bin', '1')},           ! Minimum selected SDs per non-empty bin\n"
        "tracking_fallback_to_random = .true., ! Fallback to random selection if stratified candidates are insufficient\n"
        "coalescence_output_enable = .true., ! Master switch for SD_coal_output_NetCDF_* (default ON)\n"
        "random_perturbation_enable = .false., ! Master switch for random perturbation in SD motion (default OFF)\n"
        "random_perturbation_amp = 0.d0, ! Random perturbation amplitude [m^1.5 * s^-0.5]\n"
    )
    lines[start : end + 1] = [block]
    path.write_text("".join(lines))


def _collect_binaries(bin_dir: Path) -> list[Path]:
    if not bin_dir.exists() or not bin_dir.is_dir():
        return []
    return sorted([p for p in bin_dir.glob("scale-rm*") if p.is_file()])


def _copy_binaries_to_case(case_dir: Path, binaries: list[Path]) -> None:
    for src in binaries:
        shutil.copy2(src, case_dir / src.name)


def _remove_tracked_radius_threshold(text: str) -> str:
    return re.sub(r"^.*tracked_radius_threshold\s*=.*$\n?", "", text, flags=re.MULTILINE)


def _clean_tracked_radius_threshold_in_tree(root_dir: Path) -> int:
    updated = 0
    for runconf in root_dir.glob("*/run.conf"):
        text = runconf.read_text()
        cleaned = _remove_tracked_radius_threshold(text)
        if cleaned != text:
            runconf.write_text(cleaned)
            updated += 1
    return updated


def _resolve_bin_dir(user_bin_dir: str | None) -> Path | None:
    if user_bin_dir:
        return Path(user_bin_dir).expanduser().resolve()
    env_bin_dir = os.environ.get("TRACKING_BIN_DIR", "").strip()
    if env_bin_dir:
        return Path(env_bin_dir).expanduser().resolve()
    default_in_repo = ROOT.parents[4] / "bin"
    if default_in_repo.exists():
        return default_in_repo
    for sibling in ROOT.parent.parent.iterdir():
        if sibling.is_dir() and sibling.name.startswith("SDM_product"):
            candidate = sibling / "bin"
            if candidate.exists():
                return candidate
    return None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=str, default=None)
    parser.add_argument("--clean-tracked-radius-threshold-only", action="store_true")
    return parser.parse_args()


def build_case_specs() -> list[CaseSpec]:
    specs: list[CaseSpec] = [
        CaseSpec("ref_stratified_f100_seed1001", "stratified", 1.0, 1001),
        CaseSpec("ref_random_f100_seed2001", "random", 1.0, 2001),
    ]
    for mode_idx, mode in enumerate(MODES):
        for fraction in FRACTIONS:
            frac_tag = int(round(fraction * 100))
            for i in range(1, SEEDS_PER_GROUP + 1):
                seed = (mode_idx + 1) * 10000 + frac_tag * 100 + i
                specs.append(
                    CaseSpec(
                        case_name=f"sample_{mode}_f{frac_tag:03d}_seed{seed:05d}",
                        mode=mode,
                        fraction=fraction,
                        seed_scale=seed,
                    )
                )
    return specs


def generate_cases(bin_dir_arg: str | None = None) -> None:
    """Generate all representativeness cases under this directory."""
    specs = build_case_specs()
    bin_dir = _resolve_bin_dir(bin_dir_arg)
    binaries = _collect_binaries(bin_dir) if bin_dir is not None else []
    if bin_dir is None:
        print("No bin directory resolved, skip binary copy")
    elif binaries:
        print(f"Found {len(binaries)} binaries in {bin_dir}")
    else:
        print(f"No scale-rm* binaries found in {bin_dir}, skip binary copy")
    if not TEMPLATE.exists():
        raise FileNotFoundError(f"Template not found: {TEMPLATE}")
    keep_files = {
        "generate_cases.py",
        "evaluate_representativeness.py",
        "submit_all_uoh.sh",
        "submit_all_squid.sh",
        "ensemble_index.csv",
    }
    for child in ROOT.iterdir():
        if child.name in keep_files:
            continue
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()
    for spec in specs:
        case_dir = ROOT / spec.case_name
        shutil.copytree(TEMPLATE, case_dir)
        _apply_runconf(case_dir / "run.conf", spec.mode, spec.fraction, spec.seed_scale)
        _rewrite_tracking_comment_block(case_dir / "run.conf")
        _apply_job_name(case_dir / "UoH_run.pbs", spec.case_name)
        _apply_job_name(case_dir / "squid_run.sh", spec.case_name)
        _copy_binaries_to_case(case_dir, binaries)
    index_path = ROOT / "ensemble_index.csv"
    with index_path.open("w") as f:
        f.write("case_name,mode,tracking_fraction,random_seed_scale\n")
        for spec in specs:
            f.write(f"{spec.case_name},{spec.mode},{spec.fraction:.2f},{spec.seed_scale}\n")


if __name__ == "__main__":
    args = _parse_args()
    if args.clean_tracked_radius_threshold_only:
        changed = _clean_tracked_radius_threshold_in_tree(ROOT)
        print(f"Updated {changed} run.conf files in {ROOT}")
    else:
        generate_cases(args.bin_dir)
