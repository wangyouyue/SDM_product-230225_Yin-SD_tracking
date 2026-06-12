#!/usr/bin/env bash
set -euo pipefail
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/load_basepy_2026_quiet.sh"

if [ "${1:-}" = "python" ] || [ "${1:-}" = "python3" ]; then
  shift
  exec "${GMD2026_PYTHON}" "$@"
fi

exec "$@"
