#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
qsub squid_run_v12_production_io_threshold.sh
