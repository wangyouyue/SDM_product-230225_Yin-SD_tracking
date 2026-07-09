#!/usr/bin/env python3
"""Generate v1.3 Cold TPHT production-style I/O calibration namelists."""

from __future__ import annotations

from pathlib import Path

import make_v12_production_io_configs as base


def copy_config(old_name: str, new_name: str) -> None:
    source = base.COLD_CASE / old_name
    target = base.COLD_CASE / new_name
    text = source.read_text()
    text = text.replace("production I/O calibration:", "v1.3 production I/O calibration:")
    if new_name == "run_v13_prod_low_significant.conf":
        text = base.insert_tracking_entries(
            text,
            [
                "! v1.3 production I/O calibration: optional Kohler context overhead",
                "TRACK_COLD_OUTPUT_KOHLER_CONTEXT = .true.,",
            ],
        )
    target.write_text(text)
    print(f"wrote {new_name}")


def main() -> int:
    base.main()
    copy_config("run_v12_prod_conservative.conf", "run_v13_prod_conservative.conf")
    copy_config("run_v12_prod_low_significant.conf", "run_v13_prod_low_significant.conf")
    copy_config("run_v12_prod_low_diag.conf", "run_v13_prod_low_diag.conf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
