#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

export COLD_TPHT_VALIDATION_STAMP="${COLD_TPHT_VALIDATION_STAMP:-v13_kohler_$(date +%Y%m%dT%H%M%S)}"

# The full validation runner has been migrated to the v1.3 process-event schema:
# trigger_code is process identity, trigger_level is occurrence/significant, and
# the suite includes derived Kohler activation/deactivation occurrence smokes.
bash "$SCRIPT_DIR/run_v12_full_validation.sh"
