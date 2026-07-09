#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
qsub squid_run_v12_full_validation.sh
