#!/usr/bin/env bash
# Run parallel heavy TPHT target-history and chain-validity analysis.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTDIR="${ROOT}/analysis_outputs"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)
      ROOT="$2"; shift 2 ;;
    --outdir)
      OUTDIR="$2"; shift 2 ;;
    *)
      EXTRA_ARGS+=("$1"); shift ;;
  esac
done

ENV_HELPER="${ROOT}/common/load_basepy_2026_quiet.sh"
if [ -f "${ENV_HELPER}" ]; then
  # Match the GMD2026 common environment loader used by SQUID batch jobs.
  source "${ENV_HELPER}"
fi
PYTHON_BIN="${GMD2026_PYTHON:-python3}"

workers="${TPHT_HEAVY_WORKERS:-8}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/analyze_02_tpht_heavy_parallel.py" \
  --root "${ROOT}" \
  --outdir "${OUTDIR}" \
  --workers "${workers}" \
  "${EXTRA_ARGS[@]}"
