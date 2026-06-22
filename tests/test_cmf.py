"""Regression tests: CMF extract must match the v1 archive cache files.

Reference caches live under:
  {sim}_fof/hbtplus_min_20/hbtplus_CMF_scaled{label}{snap:04d}

All tests are marked ``slow`` because they require reading the snapshot
DM-particle HDF5 (~GB) and the SubSnap HDF5. They are auto-skipped when
the snapshot file is absent (login-node safe).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures.cmf import extract, _cache_path

WORKING_ROOT = load_simulations().working_root

# All 6 fiducial sims with their snap numbers.
FIDUCIAL_SNAPS = {
    "m12i_cdmo":      50,
    "m12i_btps_deep": 50,
    "m12i_btps_soft": 50,
    "m12f_cdmo":      49,
    "m12f_btps_deep": 49,
    "m12f_btps_soft": 49,
}


def _archive_path(label: str, snap: int) -> Path:
    """Return the pre-computed CMF_scaled reference file for *label*."""
    return (
        WORKING_ROOT
        / f"{label}_fof"
        / "hbtplus_min_20"
        / f"hbtplus_CMF_scaled{label}{snap:04d}"
    )


@pytest.mark.parametrize("sim_label", list(FIDUCIAL_SNAPS.keys()))
@pytest.mark.slow
def test_cmf_extract_matches_archive(sim_label: str) -> None:
    """CMF extract output must match the v1 pre-computed archive file."""
    cfg = load_simulations()
    sim = cfg[sim_label]

    # Skip if snapshot is not available (login node).
    if not sim.snapshot_path().is_file():
        pytest.skip(f"snapshot not available: {sim.snapshot_path()}")

    snap = FIDUCIAL_SNAPS[sim_label]
    archive_file = _archive_path(sim_label, snap)
    assert archive_file.is_file(), f"archive file missing: {archive_file}"

    # Run extract with refresh to guarantee fresh computation.
    cache_file = extract(sim_label, refresh=False)

    result  = np.genfromtxt(cache_file)
    archive = np.genfromtxt(archive_file)

    assert result.shape == archive.shape, (
        f"shape mismatch: got {result.shape}, expected {archive.shape}"
    )
    np.testing.assert_allclose(result, archive, rtol=1e-12)
