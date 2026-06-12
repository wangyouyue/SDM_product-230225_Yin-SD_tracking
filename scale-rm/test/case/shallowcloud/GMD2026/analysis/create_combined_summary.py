#!/usr/bin/env python3
"""Create combined GMD2026 manuscript tables and reviewer-coverage notes."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ANALYSIS_DIR))

from common.paths import build_arg_parser, dry_run_message, ensure_output_dirs, resolve_outdir  # noqa: E402
from common.table_utils import NA, read_csv_rows, write_table_bundle  # noqa: E402


def _first(rows: list[dict], key: str, value: str) -> dict:
    """Return the first row matching a key/value pair."""
    for row in rows:
        if row.get(key) == value:
            return row
    return {}


def _value(row: dict, key: str) -> str:
    """Return a compact string value for markdown."""
    return str(row.get(key) or NA)


def analyze(outdir: Path) -> None:
    """Write combined summary and manuscript/reviewer markdown outputs."""
    tables = outdir / "tables"
    bench = read_csv_rows(tables / "01_benchmark_summary.csv")
    tpht = read_csv_rows(tables / "02_tpht_summary.csv")
    sampling = read_csv_rows(tables / "03_sampling_seed_statistics.csv")
    scaling = read_csv_rows(tables / "04_scaling_slopes.csv")
    outint = read_csv_rows(tables / "05_outint_io_summary.csv")

    combined_rows = []
    for row in bench:
        combined_rows.append(
            {
                "analysis_group": "01_benchmark",
                "case_or_metric": row.get("case_name"),
                "wallclock_s": row.get("wallclock_s"),
                "core_hours": row.get("core_hours"),
                "memory_mib": row.get("peak_memory_rank_max_mib"),
                "output_bytes": row.get("scientific_output_bytes") or row.get("total_output_bytes"),
                "key_result": row.get("wallclock_relative_to_nt_nolog"),
            }
        )
    if tpht:
        row = tpht[0]
        combined_rows.append(
            {
                "analysis_group": "02_tpht",
                "case_or_metric": "target-restricted reconstruction",
                "wallclock_s": f"FW={_value(row, 'fw_wallclock_s')}; BW={_value(row, 'bw_wallclock_s')}",
                "core_hours": row.get("total_tpht_core_hours"),
                "memory_mib": f"FW={_value(row, 'fw_peak_memory_rank_max_mib')}; BW={_value(row, 'bw_peak_memory_rank_max_mib')}",
                "output_bytes": row.get("tpht_total_reconstruction_bytes"),
                "key_result": f"status={_value(row, 'consistency_status')}; target_reduction_ratio={_value(row, 'target_reduction_ratio')}",
            }
        )
    for row in scaling:
        combined_rows.append(
            {
                "analysis_group": "04_sdnc_scaling",
                "case_or_metric": f"{row.get('tracking_label')} {row.get('metric')} slope",
                "wallclock_s": NA,
                "core_hours": NA,
                "memory_mib": NA,
                "output_bytes": NA,
                "key_result": row.get("log_log_slope_b"),
            }
        )
    for row in outint:
        if row.get("case_name") == "nt_nolog_10min":
            continue
        combined_rows.append(
            {
                "analysis_group": "05_outint_io",
                "case_or_metric": row.get("case_name"),
                "wallclock_s": row.get("wallclock_s"),
                "core_hours": row.get("core_hours"),
                "memory_mib": row.get("peak_memory_rank_max_mib"),
                "output_bytes": row.get("sd_selected_output_bytes"),
                "key_result": f"write_time_total_s={_value(row, 'sd_output_write_time_total_s')}",
            }
        )
    if not combined_rows:
        combined_rows.append(
            {
                "analysis_group": NA,
                "case_or_metric": NA,
                "wallclock_s": NA,
                "core_hours": NA,
                "memory_mib": NA,
                "output_bytes": NA,
                "key_result": "Run analysis scripts first.",
            }
        )
    write_table_bundle(combined_rows, tables / "GMD2026_combined_performance_summary")

    bench_fw = _first(bench, "case_name", "fw005_coallog")
    bench_bw = _first(bench, "case_name", "bw005_coallog")
    tpht_row = tpht[0] if tpht else {}
    sample_note = sampling[0] if sampling else {}
    slope_wall = _first(scaling, "metric", "wallclock")
    out30 = _first(outint, "case_name", "out30_fw005")
    key_lines = [
        "# GMD2026 Key Results for Manuscript",
        "",
        "- 01 benchmark overhead: FW005 relative wallclock = "
        f"{_value(bench_fw, 'wallclock_relative_to_nt_nolog')}; BW005 relative wallclock = "
        f"{_value(bench_bw, 'wallclock_relative_to_nt_nolog')}.",
        "- 02 TPHT consistency and target reduction: consistency_status = "
        f"{_value(tpht_row, 'consistency_status')}; target_reduction_ratio = {_value(tpht_row, 'target_reduction_ratio')}; "
        f"estimated_storage_reduction_factor = {_value(tpht_row, 'estimated_storage_reduction_factor')}.",
        "- 03 sampling verification: short 2D sampling-procedure verification only; first seed-statistics row seed_count = "
        f"{_value(sample_note, 'seed_count')}.",
        "- 04 SDNC scaling slopes: first available wallclock slope = "
        f"{_value(slope_wall, 'log_log_slope_b')}; see 04_scaling_slopes for NT/FW005/BW005-specific slopes.",
        "- 05 output interval I/O sensitivity: out30_fw005 selected output bytes = "
        f"{_value(out30, 'sd_selected_output_bytes')}; see 05_outint_io_summary for all intervals.",
        "",
        "Interpretation guardrails: 03_fw_rep_2d_600s is not evidence of long-time 3D representativeness, and "
        "01_bench_3d_samp_30min is a controlled cold-start computational benchmark rather than a mature-cloud physical benchmark.",
    ]
    (tables / "GMD2026_key_results_for_manuscript.md").write_text("\n".join(key_lines) + "\n")

    coverage_lines = [
        "# GMD2026 Reviewer Coverage",
        "",
        "| Reviewer concern | Evidence |",
        "|---|---|",
        "| no-tracking baseline | 01 nt_nolog and 05 nt_nolog_10min |",
        "| coalescence log overhead | 01 nt_nolog vs nt_coallog |",
        "| FW/BW fair comparison | 01 fw005_coallog vs bw005_coallog |",
        "| memory footprint | peak_memory_rank_max_mib, peak_memory_rank_sum_mib, tracking_id_memory_bytes, if_coal_memory_bytes |",
        "| scaling with super-droplet number | 04 SDNC scaling slopes |",
        "| I/O/storage burden | 01 output size and 05 interval sensitivity |",
        "| 3D backward practicality | 02 TPHT target-restricted reconstruction |",
        "| sampling representativeness | 03 short 2D sampling verification only; manuscript wording must remain conservative |",
        "",
        "These analyses do not by themselves solve:",
        "- forward tracking conceptual novelty;",
        "- scientific novelty overclaim;",
        "- literature review gaps;",
        "- long-time 3D representativeness.",
        "",
        "Those issues require manuscript text revision.",
    ]
    (tables / "GMD2026_reviewer_coverage.md").write_text("\n".join(coverage_lines) + "\n")


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(Path(__file__).name, root, outdir, ["tables/GMD2026_combined_performance_summary.*", "tables/GMD2026_key_results_for_manuscript.md", "tables/GMD2026_reviewer_coverage.md"])
        return
    ensure_output_dirs(outdir)
    analyze(outdir)


if __name__ == "__main__":
    main()
