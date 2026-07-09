#!/usr/bin/env python3
"""Generate longer v1.3 Cold TPHT production-calibration namelists."""

from __future__ import annotations

import os
from pathlib import Path

import make_v12_production_io_configs as base


LONG_DURATION = float(os.environ.get("COLD_TPHT_LONG_DURATION_SEC", "180.0"))
LOW_DIAG_DURATION = float(os.environ.get("COLD_TPHT_LONG_LOW_DIAG_DURATION_SEC", "30.0"))
RANKS = int(os.environ.get("COLD_TPHT_LONG_RANKS", os.environ.get("COLD_TPHT_PROD_RANKS", "4")))
RANK_TO_PRC = {
    1: (1, 1),
    2: (2, 1),
    4: (2, 2),
}


def with_rank_layout(text: str) -> str:
    try:
        prc_num_x, prc_num_y = RANK_TO_PRC[RANKS]
    except KeyError as exc:
        supported = ", ".join(str(rank) for rank in sorted(RANK_TO_PRC))
        raise SystemExit(f"unsupported COLD_TPHT_LONG_RANKS={RANKS}; supported: {supported}") from exc
    text = base.replace_line(text, "PRC_NUM_X", f" PRC_NUM_X       = {prc_num_x},")
    text = base.replace_line(text, "PRC_NUM_Y", f" PRC_NUM_Y       = {prc_num_y},")
    return text


def with_duration(text: str, duration: float) -> str:
    return base.replace_line(text, "TIME_DURATION", f" TIME_DURATION              = {duration:.6g}D0,")


def replace_existing_tracking_line(text: str, prefix: str, replacement: str) -> str:
    try:
        return base.replace_line(text, prefix, replacement)
    except SystemExit:
        return base.insert_tracking_entries(text, [replacement])


def write_config(source_name: str, dest_name: str, duration: float, entries: list[str]) -> None:
    text = (base.COLD_CASE / source_name).read_text()
    text = with_rank_layout(with_duration(text, duration))
    text = base.insert_tracking_entries(text, entries)
    (base.COLD_CASE / dest_name).write_text(text)
    print(f"wrote {dest_name}")


def write_moderate_diag_config() -> None:
    text = (base.COLD_CASE / "run_diag.conf").read_text()
    text = with_rank_layout(with_duration(text, LONG_DURATION))
    replacements = [
        ("tracking_diag_liq_radius_threshold", "tracking_diag_liq_radius_threshold = 5.0d-6,"),
        ("tracking_diag_ice_rvol_threshold", "tracking_diag_ice_rvol_threshold = 5.0d-6,"),
        ("tracking_diag_mixed_rvol_threshold", "tracking_diag_mixed_rvol_threshold = 5.0d-6,"),
        ("tracking_diag_rime_mass_threshold", "tracking_diag_rime_mass_threshold = 1.0d-15,"),
        ("tracking_diag_rime_frac_threshold", "tracking_diag_rime_frac_threshold = 0.05d0,"),
        ("tracking_diag_nmono_threshold", "tracking_diag_nmono_threshold = 2.0d0,"),
        ("tracking_diag_aspect_ratio_threshold", "tracking_diag_aspect_ratio_threshold = 0.05d0,"),
    ]
    for prefix, replacement in replacements:
        text = replace_existing_tracking_line(text, prefix, replacement)
    text = base.insert_tracking_entries(
        text,
        [
            "! v1.3 longer production calibration: moderate diagnostic thresholds",
        ],
    )
    (base.COLD_CASE / "run_v13_long_moderate_diag.conf").write_text(text)
    print("wrote run_v13_long_moderate_diag.conf")


def main() -> int:
    write_config(
        "run.conf",
        "run_v13_long_conservative.conf",
        LONG_DURATION,
        [
            "! v1.3 longer production calibration: conservative/default-disabled high-volume triggers",
            "tracking_evt_liq_liq_coal_enable = .true.,",
            "tracking_evt_riming_enable = .true.,",
            "tracking_evt_aggregation_enable = .true.,",
            "tracking_evt_freezing_enable = .true.,",
            "tracking_evt_melting_enable = .true.,",
        ],
    )
    write_config(
        "run_singleproc_significant.conf",
        "run_v13_long_low_significant.conf",
        LONG_DURATION,
        [
            "! v1.3 longer production calibration: low significant thresholds",
            "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            "tracking_sig_deposition_enable = .true.,",
            "tracking_sig_deposition_relmass_threshold = 0.0d0,",
            "tracking_sig_sublimation_enable = .true.,",
            "tracking_sig_sublimation_relmass_threshold = 0.0d0,",
            "tracking_sig_condensation_enable = .true.,",
            "tracking_sig_condensation_relmass_threshold = 0.0d0,",
            "tracking_sig_evaporation_enable = .true.,",
            "tracking_sig_evaporation_relmass_threshold = 0.0d0,",
            "tracking_sig_freezing_enable = .true.,",
            "tracking_sig_freezing_relmass_threshold = 0.0d0,",
            "tracking_sig_melting_enable = .true.,",
            "tracking_sig_melting_relmass_threshold = 0.0d0,",
            "tracking_sig_riming_enable = .true.,",
            "tracking_sig_riming_relmass_threshold = 0.0d0,",
            "tracking_sig_aggregation_enable = .true.,",
            "tracking_sig_aggregation_relmass_threshold = 0.0d0,",
            "tracking_sig_coalescence_enable = .true.,",
            "tracking_sig_coalescence_relmass_threshold = 0.0d0,",
            "tracking_sig_activation_enable = .true.,",
            "tracking_sig_activation_radius_threshold = 0.0d0,",
            "tracking_sig_deactivation_enable = .true.,",
            "tracking_sig_deactivation_radius_threshold = 0.0d0,",
        ],
    )
    write_moderate_diag_config()
    write_config(
        "run_diag.conf",
        "run_v13_long_low_diag_cap.conf",
        LOW_DIAG_DURATION,
        [
            "! v1.3 longer production calibration: capped near-zero diagnostic stress",
        ],
    )
    write_config(
        "run_kohler_activation_occurrence.conf",
        "run_v13_long_context_on.conf",
        LONG_DURATION,
        [
            "! v1.3 longer production calibration: optional context overhead",
            "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            "TRACK_COLD_OUTPUT_AEROSOL_CONTEXT = .true.,",
            "TRACK_COLD_OUTPUT_THERMO_CONTEXT = .true.,",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
