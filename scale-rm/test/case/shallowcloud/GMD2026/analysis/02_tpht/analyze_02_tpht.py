#!/usr/bin/env python3
"""Analyze the GMD2026 3D TPHT interest-restricted reconstruction workflow."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

from common.parse_ids import summarize_id_files  # noqa: E402
from common.parse_logs import collect_case_metrics  # noqa: E402
from common.parse_netcdf import read_selected_pairs  # noqa: E402
from common.paths import analysis_options, build_arg_parser, case_path, dry_run_message, ensure_output_dirs, iter_cases, load_config, resolve_outdir, warning_text  # noqa: E402
from common.table_utils import safe_float, safe_ratio, write_table_bundle  # noqa: E402


COLUMNS = [
    "fw_wallclock_s",
    "bw_wallclock_s",
    "fw_core_hours",
    "bw_core_hours",
    "total_tpht_core_hours",
    "fw_peak_memory_rank_max_mib",
    "bw_peak_memory_rank_max_mib",
    "fw_peak_memory_rank_sum_mib",
    "bw_peak_memory_rank_sum_mib",
    "fw_chain_count",
    "bw_chain_count",
    "raw_id_file_count",
    "raw_id_bytes",
    "merged_id_bytes",
    "dedup_id_file_count",
    "dedup_id_bytes",
    "fw_raw_id_records",
    "fw_unique_pairs",
    "deduplicated_target_pairs",
    "bw_valid_record_rows",
    "bw_valid_records",
    "bw_unique_pairs",
    "missing_in_bw",
    "extra_in_bw",
    "target_reduction_ratio",
    "target_reduction_ratio_is_approximate",
    "tpht_id_write_time_total_s",
    "tpht_id_write_count",
    "tpht_id_records_written_total",
    "bw_sd_selected_output_bytes",
    "bw_total_output_bytes",
    "fw_total_output_bytes",
    "bw_scientific_output_bytes",
    "fw_scientific_output_bytes",
    "bw_scientific_output_file_count",
    "fw_scientific_output_file_count",
    "bw_auxiliary_file_bytes",
    "fw_auxiliary_file_bytes",
    "coalescence_log_bytes_fw",
    "coalescence_log_bytes_bw",
    "tpht_handoff_bytes",
    "tpht_total_reconstruction_bytes",
    "bytes_per_target",
    "core_hours_per_target",
    "estimated_full_bw_output_bytes",
    "estimated_storage_reduction_factor",
    "id_epoch_sec_fw",
    "id_epoch_sec_bw",
    "mpi_decomposition_match",
    "fallback_to_sampling",
    "consistency_status",
    "warnings",
]


RANK_COLUMNS = [
    "rank",
    "raw_records_per_rank",
    "unique_pairs_per_rank",
    "dedup_targets_per_rank",
    "dedup_reduction_ratio",
    "rank_imbalance_raw",
    "rank_imbalance_dedup",
]

FEATURE_STATUS_COLUMNS = [
    "feature_id",
    "feature",
    "status",
    "reason",
    "output_tables",
]


def _sum_existing(paths: list[Path]) -> int:
    """Return total bytes for existing paths."""
    return sum(path.stat().st_size for path in paths if path.exists())


def _combine_fallback_to_sampling(*values: object) -> bool | None:
    """Preserve explicit false fallback diagnostics instead of treating them as missing."""
    saw_false = False
    for value in values:
        if value is True:
            return True
        if value is False or value == 0:
            saw_false = True
    return False if saw_false else None


def build_summary(root: Path, options: dict | None = None) -> tuple[dict, list[str]]:
    """Build the TPHT summary row without writing outputs."""
    config = load_config()
    cases = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}
    fw_case = cases["fw_discovery"]
    bw_case = cases["bw_reconstruction"]
    options = options or {}
    fw_row, fw_files, _fw_records, fw_warnings = collect_case_metrics(root, fw_case, options)
    bw_row, bw_files, _bw_records, bw_warnings = collect_case_metrics(root, bw_case, options)
    warnings = fw_warnings + bw_warnings

    fw_dir = case_path(root, fw_case)
    bw_dir = case_path(root, bw_case)
    fw_tracking = fw_dir / "fw_tracking"
    raw_files = sorted(path for path in fw_tracking.glob("tracking_interest_ids.pe*.ids") if path.is_file())
    merged_file = fw_tracking / "tracking_interest_ids_merged.ids"
    dedup_files = sorted(path for path in fw_tracking.glob("tracking_interest_ids_dedup.pe*.ids") if path.is_file())

    raw_summary, raw_pairs, raw_warnings = summarize_id_files(raw_files)
    warnings.extend(raw_warnings)
    merged_summary, merged_pairs, merged_warnings = summarize_id_files([merged_file] if merged_file.exists() else [])
    warnings.extend(merged_warnings)
    dedup_summary, dedup_pairs, dedup_warnings = summarize_id_files(dedup_files)
    warnings.extend(dedup_warnings)

    if not raw_files:
        warnings.append("fw_discovery must have raw .ids files")
    if safe_float(fw_row.get("sd_selected_output_bytes")) not in (None, 0.0):
        warnings.append("fw_discovery should not have SD_selected output")
    if (safe_float(bw_row.get("sd_selected_output_bytes")) or 0.0) <= 0.0:
        warnings.append("bw_reconstruction must have SD_selected output")
        bw_row["sd_selected_output_bytes"] = None
    if (safe_float(fw_row.get("coalescence_log_bytes")) or 0.0) <= 0.0:
        fw_row["coalescence_log_bytes"] = None
    if (safe_float(bw_row.get("coalescence_log_bytes")) or 0.0) <= 0.0:
        bw_row["coalescence_log_bytes"] = None

    expected_pairs = set(dedup_pairs or merged_pairs or raw_pairs)
    bw_pairs, bw_valid_record_rows, bw_pair_warnings = read_selected_pairs(
        bw_files.get("sd_selected_output_bytes", []),
        "BW",
        max_files=options.get("max_files"),
        max_records=options.get("max_records"),
        chunk_size=options.get("chunk_size", 100000),
        metadata_only=options.get("metadata_only", False),
        first_group_only=True,
    )
    warnings.extend(bw_pair_warnings)
    bw_unique_pair_count = len(bw_pairs) if bw_valid_record_rows is not None else None
    missing_in_bw = len(expected_pairs - bw_pairs) if expected_pairs and bw_unique_pair_count is not None else None
    extra_in_bw = len(bw_pairs - expected_pairs) if expected_pairs and bw_unique_pair_count is not None else None

    full_valid_sd_count = safe_float(fw_row.get("tracking_chain_count"))
    target_count = len(expected_pairs) if expected_pairs else None
    target_reduction_ratio = safe_ratio(target_count, full_valid_sd_count)
    target_reduction_is_approx = full_valid_sd_count is not None

    fw_core = safe_float(fw_row.get("core_hours"))
    bw_core = safe_float(bw_row.get("core_hours"))
    total_core = fw_core + bw_core if fw_core is not None and bw_core is not None else None
    handoff_bytes = _sum_existing(raw_files) + _sum_existing(dedup_files)
    bw_selected_bytes = safe_float(bw_row.get("sd_selected_output_bytes"))
    total_reconstruction_bytes = handoff_bytes + bw_selected_bytes if bw_selected_bytes is not None else None

    estimated_full_bw = None
    estimated_storage_factor = None
    if bw_selected_bytes is not None and target_reduction_ratio is not None and target_reduction_ratio > 0.0:
        estimated_full_bw = bw_selected_bytes / target_reduction_ratio
        if total_reconstruction_bytes:
            estimated_storage_factor = estimated_full_bw / total_reconstruction_bytes

    id_epoch_fw = fw_row.get("id_epoch_sec") or raw_summary.get("TPHT_ID_EPOCH_SEC")
    id_epoch_bw = bw_row.get("id_epoch_sec") or dedup_summary.get("TPHT_ID_EPOCH_SEC")
    mpi_decomp_match = None
    if raw_summary.get("PRC_nprocs") is not None and dedup_summary.get("PRC_nprocs") is not None:
        mpi_decomp_match = raw_summary.get("PRC_nprocs") == dedup_summary.get("PRC_nprocs")
    fallback = _combine_fallback_to_sampling(fw_row.get("fallback_to_sampling"), bw_row.get("fallback_to_sampling"))

    fail_conditions = [
        missing_in_bw not in (None, 0),
        extra_in_bw not in (None, 0),
        bool(fallback) is True,
        id_epoch_fw is not None and id_epoch_bw is not None and id_epoch_fw != id_epoch_bw,
    ]
    unknown_conditions = [missing_in_bw is None, extra_in_bw is None, id_epoch_fw is None or id_epoch_bw is None]
    if any(fail_conditions):
        consistency_status = "FAIL"
    elif any(unknown_conditions):
        consistency_status = None
    else:
        consistency_status = "PASS"

    row = {
        "fw_wallclock_s": fw_row.get("wallclock_s"),
        "bw_wallclock_s": bw_row.get("wallclock_s"),
        "fw_core_hours": fw_row.get("core_hours"),
        "bw_core_hours": bw_row.get("core_hours"),
        "total_tpht_core_hours": total_core,
        "fw_peak_memory_rank_max_mib": fw_row.get("peak_memory_rank_max_mib"),
        "bw_peak_memory_rank_max_mib": bw_row.get("peak_memory_rank_max_mib"),
        "fw_peak_memory_rank_sum_mib": fw_row.get("peak_memory_rank_sum_mib"),
        "bw_peak_memory_rank_sum_mib": bw_row.get("peak_memory_rank_sum_mib"),
        "fw_chain_count": fw_row.get("tracking_chain_count"),
        "bw_chain_count": bw_row.get("tracking_chain_count"),
        "raw_id_file_count": len(raw_files) if fw_tracking.exists() else None,
        "raw_id_bytes": _sum_existing(raw_files) if fw_tracking.exists() else None,
        "merged_id_bytes": merged_file.stat().st_size if merged_file.exists() else None,
        "dedup_id_file_count": len(dedup_files) if fw_tracking.exists() else None,
        "dedup_id_bytes": _sum_existing(dedup_files) if fw_tracking.exists() else None,
        "fw_raw_id_records": raw_summary.get("raw_id_records") if raw_files else None,
        "fw_unique_pairs": merged_summary.get("unique_pairs") if merged_file.exists() else (raw_summary.get("unique_pairs") if raw_files else None),
        "deduplicated_target_pairs": target_count,
        "bw_valid_record_rows": bw_valid_record_rows,
        "bw_valid_records": bw_unique_pair_count,
        "bw_unique_pairs": bw_unique_pair_count,
        "missing_in_bw": missing_in_bw,
        "extra_in_bw": extra_in_bw,
        "target_reduction_ratio": target_reduction_ratio,
        "target_reduction_ratio_is_approximate": target_reduction_is_approx,
        "tpht_id_write_time_total_s": fw_row.get("tpht_id_write_time_total_s"),
        "tpht_id_write_count": fw_row.get("tpht_id_write_count"),
        "tpht_id_records_written_total": fw_row.get("tpht_id_records_written_total"),
        "bw_sd_selected_output_bytes": bw_row.get("sd_selected_output_bytes"),
        "bw_total_output_bytes": bw_row.get("total_output_bytes"),
        "fw_total_output_bytes": fw_row.get("total_output_bytes"),
        "bw_scientific_output_bytes": bw_row.get("scientific_output_bytes"),
        "fw_scientific_output_bytes": fw_row.get("scientific_output_bytes"),
        "bw_scientific_output_file_count": bw_row.get("scientific_output_file_count"),
        "fw_scientific_output_file_count": fw_row.get("scientific_output_file_count"),
        "bw_auxiliary_file_bytes": bw_row.get("auxiliary_file_bytes"),
        "fw_auxiliary_file_bytes": fw_row.get("auxiliary_file_bytes"),
        "coalescence_log_bytes_fw": fw_row.get("coalescence_log_bytes"),
        "coalescence_log_bytes_bw": bw_row.get("coalescence_log_bytes"),
        "tpht_handoff_bytes": handoff_bytes if fw_tracking.exists() else None,
        "tpht_total_reconstruction_bytes": total_reconstruction_bytes,
        "bytes_per_target": safe_ratio(total_reconstruction_bytes, target_count),
        "core_hours_per_target": safe_ratio(total_core, target_count),
        "estimated_full_bw_output_bytes": estimated_full_bw,
        "estimated_storage_reduction_factor": estimated_storage_factor,
        "id_epoch_sec_fw": id_epoch_fw,
        "id_epoch_sec_bw": id_epoch_bw,
        "mpi_decomposition_match": mpi_decomp_match,
        "fallback_to_sampling": fallback,
        "consistency_status": consistency_status,
        "warnings": warning_text(warnings),
    }
    return row, warnings


def build_rank_load_rows(root: Path) -> tuple[list[dict], list[str]]:
    """Build TPHT dedup-efficiency and rank-load-balance rows."""
    config = load_config()
    fw_case = {case["case_name"]: case for case in iter_cases(config, "02_tpht_3d_interest_70min")}["fw_discovery"]
    fw_tracking = case_path(root, fw_case) / "fw_tracking"
    raw_files = sorted(path for path in fw_tracking.glob("tracking_interest_ids.pe*.ids") if path.is_file())
    dedup_files = sorted(path for path in fw_tracking.glob("tracking_interest_ids_dedup.pe*.ids") if path.is_file())
    raw_summary, _raw_pairs, raw_warnings = summarize_id_files(raw_files)
    dedup_summary, _dedup_pairs, dedup_warnings = summarize_id_files(dedup_files)
    raw_per_rank = raw_summary.get("ids_per_rank") or {}
    dedup_per_rank = dedup_summary.get("ids_per_rank") or {}
    ranks = sorted(set(raw_per_rank) | set(dedup_per_rank))
    rows = []
    for rank in ranks:
        raw_count = raw_per_rank.get(rank)
        dedup_count = dedup_per_rank.get(rank)
        rows.append(
            {
                "rank": rank,
                "raw_records_per_rank": raw_count,
                "unique_pairs_per_rank": None,
                "dedup_targets_per_rank": dedup_count,
                "dedup_reduction_ratio": safe_ratio(raw_count, dedup_count),
                "rank_imbalance_raw": raw_summary.get("rank_imbalance_ratio"),
                "rank_imbalance_dedup": dedup_summary.get("rank_imbalance_ratio"),
            }
        )
    if not rows:
        rows = [
            {
                "rank": None,
                "raw_records_per_rank": None,
                "unique_pairs_per_rank": None,
                "dedup_targets_per_rank": None,
                "dedup_reduction_ratio": None,
                "rank_imbalance_raw": raw_summary.get("rank_imbalance_ratio"),
                "rank_imbalance_dedup": dedup_summary.get("rank_imbalance_ratio"),
            }
        ]
    return rows, raw_warnings + dedup_warnings


def build_feature_status_rows() -> list[dict]:
    """Report implementation status for advanced TPHT analyses."""
    return [
        {
            "feature_id": 1,
            "feature": "TPHT handoff integrity",
            "status": "implemented_if_variables_available",
            "reason": "ID-file counts and fallback/id_epoch log fields are parsed; missing/extra require BW selected-output ID variables.",
            "output_tables": "02_tpht_summary, 02_tpht_consistency",
        },
        {
            "feature_id": 2,
            "feature": "Target-reduction efficiency",
            "status": "partially_implemented",
            "reason": "bytes/core-hours per target and storage estimates are computed; full-valid-SD denominator may use FW chain count as an approximation.",
            "output_tables": "02_tpht_summary",
        },
        {
            "feature_id": 3,
            "feature": "Interest-condition decomposition",
            "status": "implemented_if_variables_available",
            "reason": "radius_only, coal_only, and both are classified once per propagated BW chain at the first record satisfying an interest condition when radius and if_coal-like variables are available.",
            "output_tables": "02_tpht_target_categories",
        },
        {
            "feature_id": 4,
            "feature": "Target discovery time",
            "status": "implemented_if_variables_available",
            "reason": "first-detected and new-target timing use TPHT_ID_DIAG records when time_s and cumulative ID counters are present; first-hour cap shading is encoded for plotting.",
            "output_tables": "02_tpht_discovery_time",
        },
        {
            "feature_id": 5,
            "feature": "Dedup efficiency and rank-load balance",
            "status": "partially_implemented",
            "reason": "raw and dedup per-rank counts and rank imbalance are computed from .ids files; unique_pairs_per_rank remains unavailable without richer per-rank merged metadata.",
            "output_tables": "02_tpht_rank_load_balance",
        },
        {
            "feature_id": 6,
            "feature": "Reconstruction validity",
            "status": "partially_implemented",
            "reason": "Adjacent-output predecessor-link fractions are computed from BW selected output; direct event-level branch validity still depends on target-linked event diagnostics.",
            "output_tables": "02_tpht_chain_validity",
        },
        {
            "feature_id": 7,
            "feature": "Target-history diagnostics",
            "status": "implemented_if_variables_available",
            "reason": "Radius/height target summaries are reconstructed by stepwise BW predecessor-link propagation when variables exist; percentile histories are intentionally not computed by loading all records into memory.",
            "output_tables": "02_tpht_target_histories",
        },
    ]


def analyze(root: Path, outdir: Path, strict: bool, options: dict | None = None) -> None:
    """Write the TPHT core summary table."""
    row, warnings = build_summary(root, options)
    rank_rows, rank_warnings = build_rank_load_rows(root)
    warnings.extend(rank_warnings)
    if strict and warnings:
        raise RuntimeError("; ".join(warnings))
    write_table_bundle([row], outdir / "tables" / "02_tpht_summary", COLUMNS)
    write_table_bundle(rank_rows, outdir / "tables" / "02_tpht_rank_load_balance", RANK_COLUMNS)
    write_table_bundle(build_feature_status_rows(), outdir / "tables" / "02_tpht_analysis_feature_status", FEATURE_STATUS_COLUMNS)


def main() -> None:
    """CLI entry point."""
    parser = build_arg_parser(__doc__ or "")
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = resolve_outdir(root, args.outdir)
    if args.dry_run:
        dry_run_message(
            Path(__file__).name,
            root,
            outdir,
            [
                "tables/02_tpht_summary.{csv,md,tex,json}",
                "tables/02_tpht_rank_load_balance.{csv,md,tex,json}",
                "tables/02_tpht_analysis_feature_status.{csv,md,tex,json}",
            ],
        )
        return
    ensure_output_dirs(outdir)
    analyze(root, outdir, args.strict, analysis_options(args))


if __name__ == "__main__":
    main()
