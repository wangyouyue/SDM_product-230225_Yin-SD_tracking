#!/usr/bin/env python3
"""Generate Cold TPHT production-style I/O calibration namelists."""

from __future__ import annotations

import os
from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent
SMOKE_DIR = THIS_DIR.parent
COLD_CASE = SMOKE_DIR / "fw_event_type_smoke"
DURATION = float(os.environ.get("COLD_TPHT_PROD_DURATION_SEC", "30.0"))
RANKS = int(os.environ.get("COLD_TPHT_PROD_RANKS", "4"))
RANK_TO_PRC = {
    1: (1, 1),
    2: (2, 1),
    4: (2, 2),
}


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


def with_duration(text: str) -> str:
    return replace_line(text, "TIME_DURATION", f" TIME_DURATION              = {DURATION:.6g}D0,")


def with_rank_layout(text: str) -> str:
    try:
        prc_num_x, prc_num_y = RANK_TO_PRC[RANKS]
    except KeyError as exc:
        supported = ", ".join(str(rank) for rank in sorted(RANK_TO_PRC))
        raise SystemExit(f"unsupported COLD_TPHT_PROD_RANKS={RANKS}; supported: {supported}") from exc
    text = replace_line(text, "PRC_NUM_X", f" PRC_NUM_X       = {prc_num_x},")
    text = replace_line(text, "PRC_NUM_Y", f" PRC_NUM_Y       = {prc_num_y},")
    return text


def write_config(source_name: str, dest_name: str, entries: list[str]) -> None:
    text = with_rank_layout(with_duration((COLD_CASE / source_name).read_text()))
    (COLD_CASE / dest_name).write_text(insert_tracking_entries(text, entries))
    print(f"wrote {dest_name}")


def main() -> int:
    write_config(
        "run.conf",
        "run_v12_prod_conservative.conf",
        [
            "! production I/O calibration: conservative/default-disabled significant and diagnostic triggers",
        ],
    )
    write_config(
        "run_singleproc_significant.conf",
        "run_v12_prod_low_significant.conf",
        [
            "! production I/O calibration: low significant thresholds",
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
        ],
    )
    write_config(
        "run_diag.conf",
        "run_v12_prod_low_diag.conf",
        [
            "! production I/O calibration: low diagnostic thresholds",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
