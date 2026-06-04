# TPHT With-UV Retest

This directory is an isolated copy of the 70 min GMD2026 TPHT experiment with standard DYCOMS-II RF02 horizontal wind restored.

Purpose:

- Test whether the short-lived cloud in the current GMD2026 TPHT run is caused by the no-UV modification.
- Keep the TPHT FW/BW settings otherwise identical to `02_tpht_3d_interest_70min`.
- Avoid mixing a no-UV initialization with a with-UV large-scale forcing.

Required source state before building:

1. `scale-rm/src/preprocess/mod_mkinit.f90` must initialize DYCOMS RF02 wind as:

   ```fortran
   velx(k,i,j) =  3.0_RP + 4.3 * GRID_CZ(k)*1.E-3_RP
   vely(k,i,j) = -9.0_RP + 5.6 * GRID_CZ(k)*1.E-3_RP
   ```

2. Both TPHT case-local `code/mod_user.f90` files must use:

   ```fortran
   U_GEOS(k) =  3.0_RP + 4.3 * CZ(k)*1.E-3_RP
   V_GEOS(k) = -9.0_RP + 5.6 * CZ(k)*1.E-3_RP
   ```

On SQUID:

```bash
cd scale-rm/test/case/shallowcloud/GMD2026/02_tpht_3d_interest_70min_with_uv_retest
bash ./apply_with_uv_source_patch.sh
bash ./check_with_uv_source.sh

cd fw_discovery
../../common/build_squid.sh
cd ..

bash submit_tpht_squid.sh
```

If `check_with_uv_source.sh` fails, do not submit the run. The experiment would otherwise be physically mixed and hard to interpret.

The TPHT merge step uses the same post-processing scripts as the main `02_tpht_3d_interest_70min` case. `postprocess/serial_merge_ids.sh` requests 8 SQUID cores, defaults `TPHT_MERGE_WORKERS=8`, writes `postprocess/logs/merge_tracking_interest_ids.log`, and preserves the same `.ids` merge/dedup semantics as the no-UV TPHT workflow.
