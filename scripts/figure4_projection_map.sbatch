#!/bin/bash
# Render Figure 4 (matter-projection maps) for all three m12i simulations:
# m12i_cdmo, m12i_btps_deep, m12i_btps_soft.
# Uses swiftsimio to slice the snapshot, then writes a PNG per simulation.
# Figure 4 has no cache layer, so BTMW_SIM_ROOT is required.
# BTMW_REFRESH is intentionally not used here; this job always reads raw snapshots.
# Submit from the repository root.
#
# Usage:  BTMW_SIM_ROOT=/path/to/simulation_root sbatch scripts/figure4_projection_map.sbatch
# Optional mail: sbatch --mail-type=ALL --mail-user=you@example.com scripts/figure4_projection_map.sbatch

#SBATCH --job-name=btmw_proj
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/proj_%j.out
#SBATCH --error=logs/proj_%j.err

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${REPO_ROOT}"

if [[ ! -f pyproject.toml || ! -d src/btmw ]]; then
    echo "ERROR: submit this job from the BlueTilt-MWSubhalos repository root." >&2
    echo "Current directory: ${REPO_ROOT}" >&2
    exit 2
fi

module load anaconda
conda activate btmw

mkdir -p logs figures

for SIM in m12i_cdmo m12i_btps_deep m12i_btps_soft; do
    btmw projection-map --sim "${SIM}"
done
