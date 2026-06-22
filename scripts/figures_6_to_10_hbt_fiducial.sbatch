#!/bin/bash
# Reproduce the HBT-HERONS subhalo statistics for the fiducial-resolution
# runs (Figures 6, 7, 8, 9, 10). By default these commands render from
# checked-in caches. Set BTMW_REFRESH=1 to re-extract from raw HBT+SOAP output.
# Default cache-only plotting does not require BTMW_SIM_ROOT. Set it only with
# BTMW_REFRESH=1 or otherwise when reading raw simulation files.
# Submit from the repository root.
#
# Usage:  sbatch scripts/figures_6_to_10_hbt_fiducial.sbatch
# Refresh usage: BTMW_REFRESH=1 BTMW_SIM_ROOT=/path/to/simulation_root sbatch scripts/figures_6_to_10_hbt_fiducial.sbatch
# Optional mail: sbatch --mail-type=ALL --mail-user=you@example.com scripts/figures_6_to_10_hbt_fiducial.sbatch

#SBATCH --job-name=btmw_hbt
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/hbt_%j.out
#SBATCH --error=logs/hbt_%j.err

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

btmw hmf "${REFRESH_FLAG[@]}"
btmw hvf "${REFRESH_FLAG[@]}"
btmw hrf "${REFRESH_FLAG[@]}" --bin 6 --bin 7 --bin 8 --bin 9
btmw cmf "${REFRESH_FLAG[@]}"
btmw rvsv "${REFRESH_FLAG[@]}"
