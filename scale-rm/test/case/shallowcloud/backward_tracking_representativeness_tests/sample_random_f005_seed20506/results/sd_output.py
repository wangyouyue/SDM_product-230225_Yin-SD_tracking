#!/usr/bin/env python3
"""Compatibility entry point for the shared trajectory extractor."""

import runpy
import sys
from pathlib import Path


current = Path(__file__).resolve()
for parent in current.parents:
    maintained = parent / "tracking_postprocess" / "sd_output.py"
    if maintained.is_file() and maintained.resolve() != current:
        sys.argv[0] = str(maintained)
        runpy.run_path(str(maintained), run_name="__main__")
        break
else:
    raise SystemExit("Cannot locate shallowcloud/tracking_postprocess/sd_output.py")
