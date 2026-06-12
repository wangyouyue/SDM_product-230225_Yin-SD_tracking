#!/bin/bash
# Run all static audits required for the GMD2026 experiment suite.
set -euo pipefail

cd "$(dirname "$0")"
source common/load_basepy_2026_quiet.sh

bash common/check_makefiles.sh
"${GMD2026_PYTHON}" common/audit_gmd2026_configs.py
"${GMD2026_PYTHON}" common/audit_restart_chain.py
"${GMD2026_PYTHON}" common/audit_benchmark_diagnostics.py

echo "GMD2026 static audit completed successfully."
