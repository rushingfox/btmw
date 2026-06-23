#!/bin/bash
# Upload btmw_raw_data.tar.gz to the Zenodo data record.
#
# Usage (TOKEN never touches disk):
#   TOKEN=your_zenodo_token sbatch scripts/upload_zenodo.sbatch
#
# The TOKEN is passed only through environment memory and is never
# written to any log file or shell history.

#SBATCH --job-name=btmw_zenodo
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/zenodo_upload_%j.out
#SBATCH --error=logs/zenodo_upload_%j.err

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${REPO_ROOT}"

DEPOSITION_ID="20805010"
ARCHIVE="${REPO_ROOT}/btmw_raw_data.tar.gz"
ZENODO_API="https://zenodo.org/api"

# --- Checks ---
if [[ -z "${TOKEN:-}" ]]; then
    echo "ERROR: TOKEN is not set. Submit with:" >&2
    echo "  TOKEN=your_zenodo_token sbatch scripts/upload_zenodo.sbatch" >&2
    exit 1
fi

if [[ ! -f "${ARCHIVE}" ]]; then
    echo "ERROR: Archive not found: ${ARCHIVE}" >&2
    exit 1
fi

# --- Get bucket URL ---
echo "=== Fetching bucket URL for deposition ${DEPOSITION_ID} ==="
RESPONSE=$(curl -sf \
    -H "Authorization: Bearer ${TOKEN}" \
    "${ZENODO_API}/deposit/depositions/${DEPOSITION_ID}")

BUCKET_URL=$(echo "${RESPONSE}" | grep -oP '"bucket":\s*"\K[^"]+')

if [[ -z "${BUCKET_URL}" ]]; then
    echo "ERROR: Could not retrieve bucket URL. Check TOKEN and deposition ID." >&2
    exit 1
fi

echo "Bucket URL: ${BUCKET_URL}"
echo "Archive:    ${ARCHIVE} ($(du -sh "${ARCHIVE}" | cut -f1))"
echo ""

# --- Upload ---
echo "=== Uploading ... ==="
UPLOAD_RESPONSE=$(curl -# \
    -H "Authorization: Bearer ${TOKEN}" \
    --upload-file "${ARCHIVE}" \
    "${BUCKET_URL}/btmw_raw_data.tar.gz")

echo "  filename : $(echo "${UPLOAD_RESPONSE}" | grep -oP '"filename":\s*"\K[^"]+')"
echo "  size     : $(echo "${UPLOAD_RESPONSE}" | grep -oP '"filesize":\s*\K[0-9]+') bytes"
echo "  checksum : $(echo "${UPLOAD_RESPONSE}" | grep -oP '"checksum":\s*"\K[^"]+')"

echo ""
echo "=== Done. Visit https://zenodo.org/deposit/${DEPOSITION_ID} to review and publish. ==="
