#!/usr/bin/env bash
# Run lightweight TPHT summary and handoff-integrity checks.
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

"${PYTHON_BIN}" "${SCRIPT_DIR}/analyze_02_tpht.py" --root "${ROOT}" --outdir "${OUTDIR}" "${EXTRA_ARGS[@]}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/check_02_tpht_consistency.py" --root "${ROOT}" --outdir "${OUTDIR}" "${EXTRA_ARGS[@]}"
