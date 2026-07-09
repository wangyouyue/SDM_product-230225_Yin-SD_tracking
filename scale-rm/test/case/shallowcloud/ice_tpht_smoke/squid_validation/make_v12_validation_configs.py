#!/usr/bin/env python3
"""Generate disposable SQUID validation namelists for Cold TPHT v1.3."""

from __future__ import annotations

from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent
SMOKE_DIR = THIS_DIR.parent
COLD_CASE = SMOKE_DIR / "fw_event_type_smoke"
WARM_CASE = SMOKE_DIR.parent / "tpht_test" / "ft_interest_id_baseline"


def insert_tracking_entries(text: str, entries: list[str]) -> str:
    marker = "&PARAM_ATMOS_PHY_MP_SDM_TRACKING\n"
    start = text.find(marker)
    if start < 0:
        raise SystemExit("PARAM_ATMOS_PHY_MP_SDM_TRACKING block is missing")
    insert_at = start + len(marker)
    body = "".join(f"{entry}\n" for entry in entries)
    return text[:insert_at] + body + text[insert_at:]


def replace_line(text: str, prefix: str, replacement: str) -> str:
    lines = text.splitlines()
    for idx, line in enumerate(lines):
        if line.strip().startswith(prefix):
            lines[idx] = replacement
            return "\n".join(lines) + "\n"
    raise SystemExit(f"could not find line starting with {prefix!r}")


def write_config(source_name: str, dest_name: str, entries: list[str]) -> None:
    source = COLD_CASE / source_name
    dest = COLD_CASE / dest_name
    text = source.read_text()
    dest.write_text(insert_tracking_entries(text, entries))
    print(f"wrote {dest.relative_to(COLD_CASE)}")


def write_config_from_text(source_name: str, dest_name: str, text: str, entries: list[str]) -> None:
    dest = COLD_CASE / dest_name
    dest.write_text(insert_tracking_entries(text, entries))
    print(f"wrote {dest.relative_to(COLD_CASE)}")


def write_coalescence_seed_configs() -> None:
    occurrence_text = (COLD_CASE / "run_cold_liq_liq_seed.conf").read_text()
    occurrence_text = replace_line(
        occurrence_text,
        "tracking_sig_liq_liq_coal_enable",
        "tracking_sig_liq_liq_coal_enable = .false.,",
    )
    occurrence_text = replace_line(
        occurrence_text,
        "tracking_sig_liq_liq_coal_relmass_threshold",
        "tracking_sig_liq_liq_coal_relmass_threshold = 1.0d300,",
    )
    write_config_from_text(
        "run_cold_liq_liq_seed.conf",
        "run_v13_coalescence_seed_occurrence.conf",
        occurrence_text,
        [
            "! v1.3 process schema: coalescence occurrence coverage",
            "tracking_evt_coalescence_enable = .true.,",
            "tracking_sig_coalescence_enable = .false.,",
            "tracking_sig_coalescence_relmass_threshold = 1.0d300,",
        ],
    )

    significant_text = (COLD_CASE / "run_cold_liq_liq_seed.conf").read_text()
    write_config_from_text(
        "run_cold_liq_liq_seed.conf",
        "run_v13_coalescence_seed_significant.conf",
        significant_text,
        [
            "! v1.3 process schema: coalescence significant coverage",
            "tracking_evt_coalescence_enable = .true.,",
            "tracking_sig_coalescence_enable = .true.,",
            "tracking_sig_coalescence_relmass_threshold = 0.0d0,",
        ],
    )


def write_vapor_occurrence_config() -> None:
    text = (COLD_CASE / "run_singleproc.conf").read_text()
    for name in ("deposition", "sublimation"):
        text = replace_line(
            text,
            f"tracking_sig_{name}_enable",
            f"tracking_sig_{name}_enable = .false.,",
        )
        text = replace_line(
            text,
            f"tracking_sig_{name}_relmass_threshold",
            f"tracking_sig_{name}_relmass_threshold = 1.0d300,",
        )
    write_config_from_text(
        "run_singleproc.conf",
        "run_v12_vapor_occurrence.conf",
        text,
        [
            "tracking_evt_deposition_enable = .true.,",
            "tracking_evt_sublimation_enable = .true.,",
            "tracking_evt_condensation_enable = .true.,",
            "tracking_evt_evaporation_enable = .true.,",
        ],
    )
    write_config_from_text(
        "run_singleproc.conf",
        "run_vapor_occurrence_sublimation_seed.conf",
        text,
        [
            "! v1.3+ controlled sublimation occurrence seed",
            "! This lowers QV before the real sdm_subldep call; it does not write events directly.",
            "tracking_evt_deposition_enable = .false.,",
            "tracking_evt_sublimation_enable = .true.,",
            "tracking_evt_condensation_enable = .false.,",
            "tracking_evt_evaporation_enable = .false.,",
            "tracking_vapor_sublimation_smoke_enable = .true.,",
            "tracking_vapor_sublimation_smoke_qv_factor = 1.0d-2,",
        ],
    )


def write_rank4_configs() -> None:
    init_text = (COLD_CASE / "init_rank2.conf").read_text()
    run_text = (COLD_CASE / "run_rank2.conf").read_text()
    init_text = replace_line(init_text, "PRC_NUM_Y", " PRC_NUM_Y       = 2,")
    run_text = replace_line(run_text, "PRC_NUM_Y", " PRC_NUM_Y       = 2,")
    (COLD_CASE / "init_rank4.conf").write_text(init_text)
    (COLD_CASE / "run_rank4.conf").write_text(run_text)
    print("wrote init_rank4.conf")
    print("wrote run_rank4.conf")


def write_warm_rank1_configs() -> None:
    init_text = (WARM_CASE / "init.conf").read_text()
    run_text = (WARM_CASE / "run.conf").read_text()
    init_text = replace_line(init_text, "PRC_NUM_X", " PRC_NUM_X       = 1,")
    run_text = replace_line(run_text, "PRC_NUM_X", " PRC_NUM_X       = 1,")
    run_text = replace_line(
        run_text,
        "sdm_dmpvar",
        "sdm_dmpvar  = 010,    ! Rank-1 warm regression writes ordinary SD_all_NetCDF_* with legacy if_coal",
    )
    (WARM_CASE / "init_rank1_warm.conf").write_text(init_text)
    (WARM_CASE / "run_rank1_warm.conf").write_text(run_text)
    print("wrote warm init_rank1_warm.conf")
    print("wrote warm run_rank1_warm.conf")


def write_rank2_gap_configs(spatial_entries: list[str]) -> None:
    run_rank2_text = (COLD_CASE / "run_rank2.conf").read_text()
    write_config_from_text(
        "run_rank2.conf",
        "run_v12_spatial_rank2.conf",
        run_rank2_text,
        spatial_entries,
    )

    for source, dest in (
        ("run_restart_continuous.conf", "run_v12_restart_rank2_continuous.conf"),
        ("run_restart_part1.conf", "run_v12_restart_rank2_part1.conf"),
        ("run_restart_part2.conf", "run_v12_restart_rank2_part2.conf"),
    ):
        text = (COLD_CASE / source).read_text()
        text = replace_line(text, "PRC_NUM_X", " PRC_NUM_X       = 2,")
        text = replace_line(text, "PRC_NUM_Y", " PRC_NUM_Y       = 1,")
        write_config_from_text(source, dest, text, spatial_entries)


def write_spatial_level3_config() -> None:
    text = (COLD_CASE / "run.conf").read_text()
    text = replace_line(text, "domovement", "domovement       = .true.,")
    text = replace_line(text, "sdm_dmpvar", "sdm_dmpvar  = 010,")
    text = replace_line(text, "TIME_DURATION", " TIME_DURATION              = 1.D0,")
    entries = [
        "tracking_spatial_visit_enable = .true.,",
        "tracking_spatial_x_min = 85.0d0,",
        "tracking_spatial_x_max = 86.0d0,",
        "tracking_spatial_y_min = -1.0d9,",
        "tracking_spatial_y_max = 1.0d9,",
        "tracking_spatial_z_min = -1.0d9,",
        "tracking_spatial_z_max = 1.0d9,",
        "tracking_spatial_level3_smoke_enable = .true.,",
        "tracking_spatial_level3_smoke_dx = 20.0d0,",
        "tracking_spatial_level3_smoke_dy = 0.0d0,",
    ]
    write_config_from_text("run.conf", "run_v12_spatial_level3.conf", text, entries)


def write_spatial_rank_pruning_configs() -> None:
    base_text = (COLD_CASE / "run.conf").read_text()
    base_text = replace_line(base_text, "PRC_NUM_X", " PRC_NUM_X       = 2,")
    base_text = replace_line(base_text, "PRC_NUM_Y", " PRC_NUM_Y       = 1,")
    base_text = replace_line(base_text, "domovement", "domovement       = .true.,")
    base_text = replace_line(base_text, "sdm_dmpvar", "sdm_dmpvar  = 010,")
    base_text = replace_line(base_text, "TIME_DURATION", " TIME_DURATION              = 1.D0,")
    subset_entries = [
        "tracking_spatial_visit_enable = .true.,",
        "tracking_spatial_x_min = 85.0d0,",
        "tracking_spatial_x_max = 86.0d0,",
        "tracking_spatial_y_min = -1.0d9,",
        "tracking_spatial_y_max = 1.0d9,",
        "tracking_spatial_z_min = -1.0d9,",
        "tracking_spatial_z_max = 1.0d9,",
        "tracking_spatial_level3_smoke_enable = .true.,",
        "tracking_spatial_level3_smoke_dx = 20.0d0,",
        "tracking_spatial_level3_smoke_dy = 0.0d0,",
    ]
    write_config_from_text(
        "run.conf",
        "run_v13_spatial_rank_pruning_off_rank2.conf",
        base_text,
        subset_entries + [
            "tracking_spatial_rank_pruning_enable = .false.,",
        ],
    )
    write_config_from_text(
        "run.conf",
        "run_v13_spatial_rank_pruning_on_rank2.conf",
        base_text,
        subset_entries + [
            "tracking_spatial_rank_pruning_enable = .true.,",
            "tracking_spatial_rank_margin = 0.0d0,",
        ],
    )
    write_config_from_text(
        "run.conf",
        "run_v13_spatial_rank_pruning_all_overlap_rank2.conf",
        base_text,
        [
            "tracking_spatial_visit_enable = .true.,",
            "tracking_spatial_x_min = -1.0d9,",
            "tracking_spatial_x_max = 1.0d9,",
            "tracking_spatial_y_min = -1.0d9,",
            "tracking_spatial_y_max = 1.0d9,",
            "tracking_spatial_z_min = -1.0d9,",
            "tracking_spatial_z_max = 1.0d9,",
            "tracking_spatial_level3_smoke_enable = .true.,",
            "tracking_spatial_level3_smoke_dx = 20.0d0,",
            "tracking_spatial_level3_smoke_dy = 0.0d0,",
            "tracking_spatial_rank_pruning_enable = .true.,",
            "tracking_spatial_rank_margin = 0.0d0,",
        ],
    )
    write_config_from_text(
        "run.conf",
        "run_v13_spatial_rank_pruning_margin_rank2.conf",
        base_text,
        [
            "tracking_spatial_visit_enable = .true.,",
            "tracking_spatial_x_min = 151.0d0,",
            "tracking_spatial_x_max = 152.0d0,",
            "tracking_spatial_y_min = -1.0d9,",
            "tracking_spatial_y_max = 1.0d9,",
            "tracking_spatial_z_min = -1.0d9,",
            "tracking_spatial_z_max = 1.0d9,",
            "tracking_spatial_level3_smoke_enable = .true.,",
            "tracking_spatial_level3_smoke_dx = 20.0d0,",
            "tracking_spatial_level3_smoke_dy = 0.0d0,",
            "tracking_spatial_rank_pruning_enable = .true.,",
            "tracking_spatial_rank_margin = 2.0d0,",
        ],
    )


def write_kohler_activation_configs() -> None:
    base_text = (COLD_CASE / "run_singleproc_significant.conf").read_text()
    base_text = replace_line(base_text, "docondensation", "docondensation   = .true.,")
    # Cold event-file output is still gated by the historical
    # coalescence_output_enable/coal_output master switch, and the driver turns
    # that off when doautoconversion is false. Keep autoconversion enabled so
    # the derived Kohler process stream is active; the event validator checks
    # only the single-process activation/deactivation records.
    base_text = replace_line(base_text, "doautoconversion", "doautoconversion = .true.,")
    base_text = replace_line(base_text, "domeltfreeze", "domeltfreeze     = .false., ! disabled for isolated Kohler crossing smoke")
    base_text = replace_line(base_text, "dosublimation", "dosublimation    = .false., ! disabled for isolated Kohler crossing smoke")
    for name in ("condensation", "evaporation"):
        base_text = replace_line(
            base_text,
            f"tracking_sig_{name}_enable",
            f"tracking_sig_{name}_enable = .false.,",
        )
        base_text = replace_line(
            base_text,
            f"tracking_sig_{name}_relmass_threshold",
            f"tracking_sig_{name}_relmass_threshold = 1.0d300,",
        )

    write_config_from_text(
        "run_singleproc_significant.conf",
        "run_kohler_activation_occurrence.conf",
        base_text,
        [
            "! v1.3+ derived Kohler activation occurrence crossing smoke",
            "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            "tracking_evt_condensation_enable = .true.,",
            "tracking_evt_evaporation_enable = .false.,",
            "tracking_evt_activation_enable = .true.,",
            "tracking_evt_deactivation_enable = .false.,",
            "tracking_sig_activation_enable = .false.,",
            "tracking_sig_deactivation_enable = .false.,",
            "tracking_kohler_activation_smoke_enable = .true.,",
            "tracking_kohler_deactivation_smoke_enable = .false.,",
            "tracking_kohler_smoke_asl_mass = 1.0d-16,",
            "tracking_kohler_smoke_margin_fraction = 1.0d-10,",
        ],
    )

    write_config_from_text(
        "run_singleproc_significant.conf",
        "run_kohler_deactivation_occurrence.conf",
        base_text,
        [
            "! v1.3+ derived Kohler deactivation occurrence crossing smoke",
            "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            "tracking_evt_condensation_enable = .false.,",
            "tracking_evt_evaporation_enable = .true.,",
            "tracking_evt_activation_enable = .false.,",
            "tracking_evt_deactivation_enable = .true.,",
            "tracking_sig_activation_enable = .false.,",
            "tracking_sig_deactivation_enable = .false.,",
            "tracking_kohler_activation_smoke_enable = .false.,",
            "tracking_kohler_deactivation_smoke_enable = .true.,",
            "tracking_kohler_smoke_asl_mass = 1.0d-16,",
            "tracking_kohler_smoke_margin_fraction = 1.0d-10,",
        ],
    )

    write_config_from_text(
        "run_singleproc_significant.conf",
        "run_kohler_activation_significant.conf",
        base_text,
        [
            "! v1.3+ derived Kohler activation significant crossing smoke",
            "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            "tracking_evt_condensation_enable = .true.,",
            "tracking_evt_evaporation_enable = .false.,",
            "tracking_evt_activation_enable = .true.,",
            "tracking_evt_deactivation_enable = .false.,",
            "tracking_sig_activation_enable = .true.,",
            "tracking_sig_activation_radius_threshold = 0.0d0,",
            "tracking_sig_deactivation_enable = .false.,",
            "tracking_kohler_activation_smoke_enable = .true.,",
            "tracking_kohler_deactivation_smoke_enable = .false.,",
            "tracking_kohler_smoke_asl_mass = 1.0d-16,",
            "tracking_kohler_smoke_margin_fraction = 1.0d-10,",
        ],
    )

    write_config_from_text(
        "run_singleproc_significant.conf",
        "run_kohler_deactivation_significant.conf",
        base_text,
        [
            "! v1.3+ derived Kohler deactivation significant crossing smoke",
            "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            "tracking_evt_condensation_enable = .false.,",
            "tracking_evt_evaporation_enable = .true.,",
            "tracking_evt_activation_enable = .false.,",
            "tracking_evt_deactivation_enable = .true.,",
            "tracking_sig_activation_enable = .false.,",
            "tracking_sig_deactivation_enable = .true.,",
            "tracking_sig_deactivation_radius_threshold = 0.0d0,",
            "tracking_kohler_activation_smoke_enable = .false.,",
            "tracking_kohler_deactivation_smoke_enable = .true.,",
            "tracking_kohler_smoke_asl_mass = 1.0d-16,",
            "tracking_kohler_smoke_margin_fraction = 1.0d-10,",
        ],
    )


def write_optional_context_configs() -> None:
    activation_text = (COLD_CASE / "run_kohler_activation_occurrence.conf").read_text()
    write_config_from_text(
        "run_kohler_activation_occurrence.conf",
        "run_v13_aerosol_context.conf",
        activation_text,
        [
            "TRACK_COLD_OUTPUT_AEROSOL_CONTEXT = .true.,",
            "TRACK_COLD_OUTPUT_THERMO_CONTEXT = .false.,",
        ],
    )
    write_config_from_text(
        "run_kohler_activation_occurrence.conf",
        "run_v13_thermo_context.conf",
        activation_text,
        [
            "TRACK_COLD_OUTPUT_AEROSOL_CONTEXT = .false.,",
            "TRACK_COLD_OUTPUT_THERMO_CONTEXT = .true.,",
        ],
    )


def main() -> int:
    write_warm_rank1_configs()
    write_config(
        "run_singleproc_significant.conf",
        "run_v12_ice_geometry.conf",
        ["TRACK_COLD_OUTPUT_ICE_GEOMETRY = .true.,"],
    )
    write_config(
        "run.conf",
        "run_v12_rime_morphology.conf",
        ["TRACK_COLD_OUTPUT_RIME_MORPHOLOGY = .true.,"],
    )
    write_vapor_occurrence_config()
    write_coalescence_seed_configs()
    spatial_entries = [
        "tracking_spatial_visit_enable = .true.,",
        "tracking_spatial_x_min = -1.0d9,",
        "tracking_spatial_x_max = 1.0d9,",
        "tracking_spatial_y_min = -1.0d9,",
        "tracking_spatial_y_max = 1.0d9,",
        "tracking_spatial_z_min = -1.0d9,",
        "tracking_spatial_z_max = 1.0d9,",
    ]
    write_config(
        "run_restart_continuous.conf",
        "run_v12_spatial_restart_continuous.conf",
        spatial_entries,
    )
    write_config(
        "run_restart_part1.conf",
        "run_v12_spatial_restart_part1.conf",
        spatial_entries,
    )
    write_config(
        "run_restart_part2.conf",
        "run_v12_spatial_restart_part2.conf",
        spatial_entries,
    )
    write_rank2_gap_configs(spatial_entries)
    write_rank4_configs()
    write_spatial_level3_config()
    write_spatial_rank_pruning_configs()
    write_kohler_activation_configs()
    write_optional_context_configs()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
