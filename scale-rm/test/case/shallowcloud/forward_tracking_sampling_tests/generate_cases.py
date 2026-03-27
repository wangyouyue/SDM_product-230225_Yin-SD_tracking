from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent


CASE_SPECS = {
    "ft_stratified_baseline": {
        "mode": "stratified",
        "fraction": "0.10d0",
        "height_min": "400.d0",
        "height_max": "800.d0",
        "radius_min": "0.d0",
        "radius_max": "0.d0",
        "nz_bin": "8",
        "nr_bin": "10",
        "min_per_bin": "1",
    },
    "ft_random_mode": {
        "mode": "random",
        "fraction": "0.10d0",
        "height_min": "400.d0",
        "height_max": "800.d0",
        "radius_min": "0.d0",
        "radius_max": "0.d0",
        "nz_bin": "8",
        "nr_bin": "10",
        "min_per_bin": "1",
    },
    "ft_fraction_0p2": {
        "mode": "stratified",
        "fraction": "0.20d0",
        "height_min": "400.d0",
        "height_max": "800.d0",
        "radius_min": "0.d0",
        "radius_max": "0.d0",
        "nz_bin": "8",
        "nr_bin": "10",
        "min_per_bin": "1",
    },
    "ft_nzbin_5": {
        "mode": "stratified",
        "fraction": "0.10d0",
        "height_min": "400.d0",
        "height_max": "800.d0",
        "radius_min": "0.d0",
        "radius_max": "0.d0",
        "nz_bin": "5",
        "nr_bin": "10",
        "min_per_bin": "1",
    },
    "ft_nrbin_6": {
        "mode": "stratified",
        "fraction": "0.10d0",
        "height_min": "400.d0",
        "height_max": "800.d0",
        "radius_min": "0.d0",
        "radius_max": "0.d0",
        "nz_bin": "8",
        "nr_bin": "6",
        "min_per_bin": "1",
    },
    "ft_height_0_1000": {
        "mode": "stratified",
        "fraction": "0.10d0",
        "height_min": "0.d0",
        "height_max": "1000.d0",
        "radius_min": "0.d0",
        "radius_max": "0.d0",
        "nz_bin": "8",
        "nr_bin": "10",
        "min_per_bin": "1",
    },
    "ft_radius_min_1e6": {
        "mode": "stratified",
        "fraction": "0.10d0",
        "height_min": "400.d0",
        "height_max": "800.d0",
        "radius_min": "1.d-6",
        "radius_max": "0.d0",
        "nz_bin": "8",
        "nr_bin": "10",
        "min_per_bin": "1",
    },
}


def _replace_key_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(^\s*{re.escape(key)}\s*=\s*)([^,\n]+)(\s*,.*$)", re.MULTILINE)
    return pattern.sub(rf"\g<1>{value}\g<3>", text)


def _remove_legacy_keys(text: str) -> str:
    text = re.sub(r"^\s*tracking_sample_initialized\s*=.*$\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\s*coalescence_output_enable\s*=.*$\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\s*tracked_grid_coal_only\s*=.*$\n?", "", text, flags=re.MULTILINE)
    return text


def _rewrite_tracking_block(text: str, spec: dict[str, str]) -> str:
    start = re.search(r"^\s*forward_tracking_enable\s*=.*$", text, flags=re.MULTILINE)
    if start is None:
        start = re.search(r"^\s*backward_tracking_enable\s*=.*$", text, flags=re.MULTILINE)
    end = re.search(r"^\s*tracking_fallback_to_random\s*=.*$", text, flags=re.MULTILINE)
    if start is None or end is None or end.start() < start.start():
        return text
    block = (
        "forward_tracking_enable = .true.,  ! Master switch for forward tracking\n"
        "backward_tracking_enable = .false., ! Master switch for backward tracking\n"
        f"tracking_fraction = {spec['fraction']},          ! Fraction of SDs to track (0-1]\n"
        "max_tracked_sds = 0,                ! Maximum number of tracked SDs (0 = unlimited)\n"
        f"tracking_selection_mode = \"{spec['mode']}\", ! \"random\" or \"stratified\"\n"
        f"tracking_height_min = {spec['height_min']},       ! Minimum height for stratified selection [m]\n"
        f"tracking_height_max = {spec['height_max']},       ! Maximum height for stratified selection [m]\n"
        f"tracking_radius_min = {spec['radius_min']},         ! Minimum radius for stratified selection [m]\n"
        f"tracking_radius_max = {spec['radius_max']},         ! Maximum radius for stratified selection [m], <= tracking_radius_min uses runtime maximum radius\n"
        f"tracking_nz_bin = {spec['nz_bin']},                ! Number of vertical bins for stratified selection\n"
        f"tracking_nr_bin = {spec['nr_bin']},               ! Number of radius bins for stratified selection\n"
        f"tracking_min_per_bin = {spec['min_per_bin']},           ! Minimum selected SDs per non-empty bin\n"
        "tracking_fallback_to_random = .true., ! Fallback to random selection if stratified candidates are insufficient\n"
    )
    return text[: start.start()] + block + text[end.end() :]


def _apply_case(case_dir: Path, spec: dict[str, str]) -> None:
    runconf = case_dir / "run.conf"
    text = runconf.read_text()
    text = _remove_legacy_keys(text)
    text = _replace_key_value(text, "forward_tracking_enable", ".true.")
    text = _replace_key_value(text, "backward_tracking_enable", ".false.")
    text = _replace_key_value(text, "max_tracked_sds", "0")
    text = _replace_key_value(text, "sdm_noise_amp", "0.d0")
    text = _replace_key_value(text, "coal_output", "0")
    text = _rewrite_tracking_block(text, spec)
    runconf.write_text(text)

    for job_name in ("UoH_run.pbs", "squid_run.sh"):
        job_script = case_dir / job_name
        if job_script.exists():
            job_text = job_script.read_text()
            job_text = re.sub(r"^#PBS -N .*$", f"#PBS -N {case_dir.name}", job_text, flags=re.MULTILINE)
            job_script.write_text(job_text)


def main() -> None:
    for case_name, spec in CASE_SPECS.items():
        case_dir = ROOT / case_name
        if not case_dir.exists():
            raise FileNotFoundError(f"missing case directory: {case_dir}")
        _apply_case(case_dir, spec)
    print(f"updated {len(CASE_SPECS)} forward sampling cases in {ROOT}")


if __name__ == "__main__":
    main()
