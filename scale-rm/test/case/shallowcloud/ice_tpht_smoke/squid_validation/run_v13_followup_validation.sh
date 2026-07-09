#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

export COLD_TPHT_VALIDATION_STAMP="${COLD_TPHT_VALIDATION_STAMP:-v13_followup_$(date +%Y%m%dT%H%M%S)}"

# Follow-up suite for v1.3+ after local Phase 0-5:
# - controlled vapor occurrence 7:1/8:1/9:1 coverage
# - significant Kohler activation/deactivation 10:2/11:2 coverage
# - spatial rank-pruning on/off comparison
# - aerosol/thermo/Kohler optional context schema checks
#
# The underlying runner name is historical; its contents are the current v1.3+
# full schema/regression suite.
bash "$SCRIPT_DIR/run_v12_full_validation.sh"
