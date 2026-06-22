#!/bin/bash
# Render Figure 5 (main-halo radial density profile) for both hosts.
# m12i uses HBT-HERONS center; m12f uses VR center (PRD version).
# Default cache-only plotting does not require BTMW_SIM_ROOT. Set it only with
# BTMW_REFRESH=1 or otherwise when reading raw simulation files.
# Submit from the repository root.
#
# Usage:  sbatch scripts/figure5_radial_density.sbatch
# Refresh usage: BTMW_REFRESH=1 BTMW_SIM_ROOT=/path/to/simulation_root sbatch scripts/figure5_radial_density.sbatch
# Optional mail: sbatch --mail-type=ALL --mail-user=you@example.com scripts/figure5_radial_density.sbatch

#SBATCH --job-name=btmw_rho
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/rho_%j.out
#SBATCH --error=logs/rho_%j.err

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

REFRESH_FLAG=()
case "${BTMW_REFRESH:-0}" in
	1|true|TRUE|True|yes|YES|Yes)
		if [[ -z "${BTMW_SIM_ROOT:-}" ]]; then
			echo "ERROR: BTMW_REFRESH=1 requires BTMW_SIM_ROOT=/path/to/simulation_root" >&2
			exit 2
		fi
		REFRESH_FLAG=(--refresh)
		;;
	0|false|FALSE|False|no|NO|No|"")
		;;
	*)
		echo "ERROR: BTMW_REFRESH must be 0/1, true/false, or yes/no; got '${BTMW_REFRESH}'" >&2
		exit 2
		;;
esac

btmw radial-density "${REFRESH_FLAG[@]}" --host m12i --center-from hbt
btmw radial-density "${REFRESH_FLAG[@]}" --host m12f --center-from vr
