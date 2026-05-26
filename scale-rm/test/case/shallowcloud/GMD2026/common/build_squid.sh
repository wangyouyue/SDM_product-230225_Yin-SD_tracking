#!/usr/bin/env bash
set -euo pipefail
if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

export PAGER=cat
export LESS="${LESS:-FRX}"
export LMOD_PAGER=cat
export MODULES_PAGER=cat
export SCALE_SYS="${SCALE_SYS:-Linux64-intel-impi}"

if command -v module >/dev/null 2>&1; then
  module load BaseCPU/2024 inteloneAPI/2023.2 hdf5/1.12.3.mpi netcdf-c/4.9.2 netcdf-fortran/4.6.1
fi

make allclean
make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
make -j "${MAKE_JOBS:-8}" SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T SCALE_DYCOMS2_RF02_SDM=T
topdir=$(cd ../../../../../../.. && pwd)
for exe in scale-rm scale-rm_init; do
  if [ ! -x "${topdir}/bin/${exe}" ]; then
    echo "Missing built executable: ${topdir}/bin/${exe}" >&2
    exit 127
  fi
  cp -f "${topdir}/bin/${exe}" .
  chmod +x "./${exe}"
done
