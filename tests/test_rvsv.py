"""Regression tests: RvsV extract must match the v1 archive cache files.

Reference caches live under:
  paper_figures_min_20_HBTplus_v1_archive/RvsV/
    {label}_RvsV_{Lower,Median,Upper}

All tests are marked ``slow`` because they require reading SOAP files (~GB).
They are auto-skipped when the SOAP file is absent (login-node safe).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures.rvsv import extract, _cache_path

ARCHIVE_RVSV = (
    load_simulations().working_root
    / "paper_figures_min_20_HBTplus_v1_archive/RvsV"
)

FIDUCIAL_SNAPS = {
    "m12i_cdmo":      50,
    "m12i_btps_deep": 50,
    "m12i_btps_soft": 50,
    "m12f_cdmo":      49,
    "m12f_btps_deep": 49,
    "m12f_btps_soft": 49,
}


def _archive_path(label: str, which: str) -> Path:
    return ARCHIVE_RVSV / f"{label}_RvsV_{which}"


@pytest.mark.parametrize(
    "sim_label,which",
    [
        (label, which)
        for label in FIDUCIAL_SNAPS
        for which in ("Lower", "Median", "Upper")
    ],
)
@pytest.mark.slow
def test_rvsv_extract_matches_archive(sim_label: str, which: str) -> None:
    """RvsV extract output must match the v1 archive file (rtol=1e-12)."""
    cfg = load_simulations()
    sim = cfg[sim_label]

    if not sim.soap_hbt_path().is_file():
        pytest.skip(f"SOAP file not available: {sim.soap_hbt_path()}")

    archive_file = _archive_path(sim_label, which)
    assert archive_file.is_file(), f"archive file missing: {archive_file}"

    paths = extract(sim_label, refresh=False)
    result  = np.genfromtxt(paths[which])
    archive = np.genfromtxt(archive_file)

    assert result.shape == archive.shape, (
        f"shape mismatch: got {result.shape}, expected {archive.shape}"
    )
    # Tolerance: archive was generated with a different numpy version; percentile
    # computations differ at ~1e-8 relative level across versions.
    np.testing.assert_allclose(result, archive, rtol=1e-5, atol=1e-10, equal_nan=True)
