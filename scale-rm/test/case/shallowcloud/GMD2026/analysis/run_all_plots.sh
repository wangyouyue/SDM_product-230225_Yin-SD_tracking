#!/usr/bin/env bash
# Run all GMD2026 manuscript plotting scripts.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT}/analysis_outputs"
STRICT=()
DRY_RUN=()
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)
      ROOT="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --strict)
      STRICT=(--strict)
      shift
      ;;
    --dry-run)
      DRY_RUN=(--dry-run)
      shift
      ;;
    --quick|--metadata-only|--skip-heavy-netcdf)
      EXTRA_ARGS+=("$1")
      shift
      ;;
    --max-files|--max-records|--workers|--chunk-size)
      EXTRA_ARGS+=("$1" "$2")
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_HELPER="${ROOT}/common/load_basepy_2026_quiet.sh"
if [ -f "${ENV_HELPER}" ]; then
  # Match the GMD2026 common environment loader used by SQUID batch jobs.
  source "${ENV_HELPER}"
fi
PYTHON_BIN="${GMD2026_PYTHON:-python3}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/00_restart/plot_00_restart_timeline.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/01_benchmark/plot_01_benchmark.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/02_tpht/plot_02_tpht.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/02_tpht/plot_02_tpht_science.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/02_tpht/plot_02_tpht_predecessor_tree.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/03_sampling/plot_03_sampling.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/04_sdnc_scaling/plot_04_sdnc_scaling.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/05_outint_io/plot_05_outint_io.py" --root "${ROOT}" --outdir "${OUTDIR}" "${STRICT[@]}" "${DRY_RUN[@]}" "${EXTRA_ARGS[@]}"
