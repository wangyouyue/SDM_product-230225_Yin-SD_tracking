#!/usr/bin/env bash
# Regenerate GMD2026 manuscript-candidate figures from existing analysis outputs.
set -euo pipefail

root="$(pwd)"
input=""
outdir=""
dry_run=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --root)
      root="$2"
      shift 2
      ;;
    --input)
      input="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    --dry-run)
      dry_run="--dry-run"
      shift
      ;;
    -h|--help)
      cat <<'EOF'
Usage:
  bash analysis/regenerate_gmd_figures.sh \
    --root /path/to/GMD2026 \
    --input /path/to/GMD2026/analysis_outputs \
    --outdir /path/to/GMD2026/analysis_outputs_for_GMD

Options:
  --dry-run  Print planned outputs without reading source data.
EOF
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

case "$root" in
  */analysis)
    root="$(cd "$root/.." && pwd)"
    ;;
  */analysis/job_scripts)
    root="$(cd "$root/../.." && pwd)"
    ;;
esac

if [ -z "$input" ]; then
  input="${root}/analysis_outputs"
fi
if [ -z "$outdir" ]; then
  outdir="${root}/analysis_outputs_for_GMD"
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

env_helper="${GMD2026_ENV_HELPER:-${root}/common/load_basepy_2026_quiet.sh}"
if [ -f "${env_helper}" ]; then
  # Match the environment-loading pattern used by common SQUID analysis scripts
  # when that helper is available, while remaining portable on local machines.
  # shellcheck source=/dev/null
  source "${env_helper}"
fi

python_cmd="${GMD2026_PYTHON:-python3}"

if [ -n "$dry_run" ]; then
  "${python_cmd}" "${script_dir}/generate_gmd_candidate_outputs.py" \
    --root "$root" \
    --input "$input" \
    --outdir "$outdir" \
    --dry-run
  "${python_cmd}" "${script_dir}/check_gmd_figure_layout.py" \
    --outdir "$outdir" \
    --dry-run
  exit 0
fi

mkdir -p "${outdir}/logs"
{
  echo "GMD2026 GMD-candidate regeneration"
  echo "root=${root}"
  echo "input=${input}"
  echo "outdir=${outdir}"
  date '+start_time=%Y-%m-%dT%H:%M:%S%z'
} | tee "${outdir}/logs/regenerate_gmd_figures.log"

"${python_cmd}" "${script_dir}/generate_gmd_candidate_outputs.py" \
  --root "$root" \
  --input "$input" \
  --outdir "$outdir" 2>&1 | tee -a "${outdir}/logs/regenerate_gmd_figures.log"

"${python_cmd}" "${script_dir}/check_gmd_figure_layout.py" \
  --outdir "$outdir" 2>&1 | tee -a "${outdir}/logs/regenerate_gmd_figures.log"

date '+end_time=%Y-%m-%dT%H:%M:%S%z' | tee -a "${outdir}/logs/regenerate_gmd_figures.log"
