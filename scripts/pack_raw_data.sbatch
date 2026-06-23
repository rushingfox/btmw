#!/bin/bash
# Pack btmw raw simulation data for Zenodo archiving.
#
# Usage (submit from repository root):
#   BTMW_SIM_ROOT=/path/to/sims sbatch scripts/pack_raw_data.sbatch
#
# Output: btmw_raw_data.tar.gz in the current directory.
# Users who download this archive should run:
#   tar -xzf btmw_raw_data.tar.gz
#   export BTMW_SIM_ROOT=/path/to/btmw_raw_data

#SBATCH --job-name=btmw_pack
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/pack_%j.out
#SBATCH --error=logs/pack_%j.err

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${REPO_ROOT}"

if [[ -z "${BTMW_SIM_ROOT:-}" ]]; then
    echo "ERROR: BTMW_SIM_ROOT is not set. Submit with:" >&2
    echo "  BTMW_SIM_ROOT=/path/to/sims sbatch scripts/pack_raw_data.sbatch" >&2
    exit 1
fi

# ===========================================================================
# FILES TO PACKAGE
# All paths are relative to BTMW_SIM_ROOT. Shell globs are expanded at runtime.
# Edit this list if you want to add or remove files before submitting.
# ===========================================================================
FILES=(
    # --- m12i_cdmo (fiducial, snap 0050, has VR) ---
    "m12i_cdmo_fof/soap_min_20/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_cdmo_fof/hbtplus_min_20/050/SubSnap_050.*.hdf5"
    "m12i_cdmo_fof/soap_VR_Aquarius_01_without_LRS/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_cdmo/velociraptor_0050_Aquarius_01_without_LRS/output.properties"
    "m12i_cdmo/m12i_cdmo_0050.hdf5"

    # --- m12i_btps_deep (fiducial, snap 0050, has VR) ---
    "m12i_btps_deep_fof/soap_min_20/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_btps_deep_fof/hbtplus_min_20/050/SubSnap_050.*.hdf5"
    "m12i_btps_deep_fof/soap_VR_Aquarius_01_without_LRS/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_btps_deep/velociraptor_0050_Aquarius_01_without_LRS/output.properties"
    "m12i_btps_deep/m12i_btps_deep_0050.hdf5"

    # --- m12i_btps_soft (fiducial, snap 0050, has VR) ---
    "m12i_btps_soft_fof/soap_min_20/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_btps_soft_fof/hbtplus_min_20/050/SubSnap_050.*.hdf5"
    "m12i_btps_soft_fof/soap_VR_Aquarius_01_without_LRS/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_btps_soft/velociraptor_0050_Aquarius_01_without_LRS/output.properties"
    "m12i_btps_soft/m12i_btps_soft_0050.hdf5"

    # --- m12f_cdmo (fiducial, snap 0049, no VR) ---
    "m12f_cdmo_fof/soap_min_20/SOAP_uncompressed/halo_properties_0049.hdf5"
    "m12f_cdmo_fof/hbtplus_min_20/049/SubSnap_049.*.hdf5"
    "m12f_cdmo/m12f_cdmo_0049.hdf5"

    # --- m12f_btps_deep (fiducial, snap 0049, no VR) ---
    "m12f_btps_deep_fof/soap_min_20/SOAP_uncompressed/halo_properties_0049.hdf5"
    "m12f_btps_deep_fof/hbtplus_min_20/049/SubSnap_049.*.hdf5"
    "m12f_btps_deep/m12f_btps_deep_0049.hdf5"

    # --- m12f_btps_soft (fiducial, snap 0049, no VR) ---
    "m12f_btps_soft_fof/soap_min_20/SOAP_uncompressed/halo_properties_0049.hdf5"
    "m12f_btps_soft_fof/hbtplus_min_20/049/SubSnap_049.*.hdf5"
    "m12f_btps_soft/m12f_btps_soft_0049.hdf5"

    # --- m12i_cdmo_highest (high-res, snap 0050, no VR) ---
    "m12i_cdmo_highest_fof/soap_min_20/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_cdmo_highest_fof/hbtplus_min_20/050/SubSnap_050.*.hdf5"

    # --- m12i_btps_deep_highest (high-res, snap 0050, no VR) ---
    "m12i_btps_deep_highest_fof/soap_min_20/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_btps_deep_highest_fof/hbtplus_min_20/050/SubSnap_050.*.hdf5"

    # --- m12i_btps_soft_highest (high-res, snap 0050, no VR) ---
    "m12i_btps_soft_highest_fof/soap_min_20/SOAP_uncompressed/halo_properties_0050.hdf5"
    "m12i_btps_soft_highest_fof/hbtplus_min_20/050/SubSnap_050.*.hdf5"
)
# ===========================================================================

ARCHIVE_ROOT="btmw_raw_data"
OUTFILE="${REPO_ROOT}/${ARCHIVE_ROOT}.tar.gz"
FILELIST=$(mktemp)

# Expand globs and build file list (relative paths under BTMW_SIM_ROOT)
for PATTERN in "${FILES[@]}"; do
    FULL="${BTMW_SIM_ROOT}/${PATTERN}"
    MATCHES=$(ls ${FULL} 2>/dev/null) || true
    if [[ -z "$MATCHES" ]]; then
        echo "WARNING: no files matched: ${FULL}" >&2
    else
        echo "$MATCHES" | sed "s|^${BTMW_SIM_ROOT}/||" >> "${FILELIST}"
    fi
done

echo "=== Files to archive: $(wc -l < "${FILELIST}") files ==="
echo "=== Estimated size ==="
sed "s|^|${BTMW_SIM_ROOT}/|" "${FILELIST}" | xargs du -ch 2>/dev/null | tail -1

echo ""
echo "=== Creating ${OUTFILE} ==="
# -C changes into BTMW_SIM_ROOT so paths inside the archive are relative;
# --transform then prepends btmw_raw_data/ to get btmw_raw_data/m12i_cdmo_fof/...
tar -czf "${OUTFILE}" \
    --transform "s|^|${ARCHIVE_ROOT}/|" \
    -C "${BTMW_SIM_ROOT}" \
    -T "${FILELIST}"

echo "Done: ${OUTFILE} ($(du -sh "${OUTFILE}" | cut -f1))"
echo ""
echo "Users should run:"
echo "  tar -xzf btmw_raw_data.tar.gz"
echo "  export BTMW_SIM_ROOT=/path/to/btmw_raw_data"

rm -f "${FILELIST}"
