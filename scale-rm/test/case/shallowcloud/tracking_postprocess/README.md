# Tracking trajectory post-processing

This directory contains the maintained warm-SDM FW/BW snapshot trajectory
workflow. It adapts the identity-pair and actual-output-time handling used by
the GMD TPHT analysis to the legacy per-case snapshot workflow. The Bitbucket
cold-TPHT event-stream scripts solve a different problem and are therefore not
used as a drop-in replacement. The small `results/sd_output.py` and
`results/traj_plot.py` files in individual cases are compatibility entry points
that dispatch here.

The extractor groups the files that actually exist, accepts MPI rank 0 as a
valid `dm_id`, treats `if_coal` as optional, and joins records by the model
identity pair rather than by a mutable NetCDF array index.  Missing outputs and
periodic-boundary crossings therefore create line breaks instead of false
trajectory segments.

The required Python packages are `numpy`, `netCDF4`, and `matplotlib`.

From a case's `results/` directory, run:

```bash
python3 sd_output.py --input-dir .. --output tracking_trajectories.nc
python3 traj_plot.py --input tracking_trajectories.nc \
  --style gmd --figure-width double \
  --output-base particle_trajectories
```

`gmd` is the default style. It follows the GMD2026 figure contract: 85 or
170 mm column width, compact 7--8 pt typography, inward ticks, color-blind-safe
endpoint colors, a perceptually uniform radius scale, embedded vector text,
and PDF/SVG/600 dpi PNG output. Use `--style diagnostic` when a larger title
and screen-oriented layout are more useful than a manuscript panel. Manuscript
titles should normally be supplied in the caption; use `--title` only when an
in-panel title is needed, and `--panel-label '(a)'` when assembling panels.
The default view zooms to the selected trajectories. Add `--full-domain` when
the complete periodic x domain and the distance from the surface are part of
the scientific comparison; `--z-min-m` and `--z-max-m` provide explicit,
reproducible height limits for matched panels.

For short smoke tests or no-wind controls, absolute positions can hide small
but valid motion. Add `--relative-position` to plot each trajectory as
horizontal and vertical displacement from its own first record. This mode is
especially useful for checking that a no-UV experiment does not contain the
kilometre-scale horizontal jumps produced by the obsolete array-index script.

Start/end markers use increasing model time. Add `--reverse-time` to the
plotting command only when a BW figure should be shown in reconstruction
direction.

The extractor automatically uses `SD_selected_NetCDF_*` when present and
otherwise uses `SD_all_NetCDF_*`.  It automatically prefers FW `dm_id/sd_id` or
BW `pre_dmid/pre_sdid` according to the file schema.  Automatic extraction is
limited to 100 deterministic anchor-time targets; use
`--max-trajectories 0` only when the resulting dataset size is acceptable.
Use `--stream selected` or `--stream all` when both output streams coexist and
the automatic preference is not appropriate. `if_coal` is optional; it is
stored as `-1` when event logging or physical coalescence is disabled.

For one exact trajectory, use its model identity pair:

```bash
python3 sd_output.py --input-dir .. --target 0:1473897
python3 traj_plot.py --target 0:1473897
```

Do not use the compatibility `particle_index` as a physical identifier.  The
authoritative identifier is always the pair `(dm_id, sd_id)` for FW or the
corresponding `(pre_dmid, pre_sdid)` fields written by BW.
