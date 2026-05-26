#!/usr/bin/env bash
# Run full-time BW history-field sanity diagnostics for the TPHT demonstration.
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

ENV_HELPER="${GMD2026_ENV_HELPER:-${ROOT}/common/load_basepy_2026_quiet.sh}"
if [ -f "${ENV_HELPER}" ]; then
  # Match the GMD2026 common environment loader used by SQUID batch jobs.
  # shellcheck source=/dev/null
  source "${ENV_HELPER}"
fi
PYTHON_BIN="${GMD2026_PYTHON:-python3}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/analyze_02_tpht_history_sanity.py" \
  --root "${ROOT}" \
  --outdir "${OUTDIR}" \
  "${EXTRA_ARGS[@]}"
