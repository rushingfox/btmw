"""Regression tests: HRF VR-SOAP extract must match the v1 archive cache files.

Reference caches live under:
  paper_figures_min_20_HBTplus_v1_archive/comparison_with_VR-SOAP/HRF_auto_bin/{bin}/
    VR-SOAP_HRF_{label}0050{bin}

All tests are marked ``slow``. Auto-skipped when VR-SOAP file is absent.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures.hrf import extract_vr, _vr_cache_path, MASS_BINS

ARCHIVE_VR_HRF_ROOT = (
    load_simulations().working_root
    / "paper_figures_min_20_HBTplus_v1_archive/comparison_with_VR-SOAP/HRF_auto_bin"
)

VR_SIMS = ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft")


@pytest.mark.parametrize("sim_label", VR_SIMS)
@pytest.mark.parametrize("mass_bin", MASS_BINS)
@pytest.mark.slow
def test_hrf_vr_extract_matches_archive(sim_label: str, mass_bin: int) -> None:
    """VR-SOAP HRF extract must match the v1 archive file (rtol=1e-5)."""
    cfg = load_simulations()
    sim = cfg[sim_label]

    if sim.soap_vr_path() is None or not sim.soap_vr_path().is_file():
        pytest.skip(f"VR-SOAP file not available: {sim.soap_vr_path()}")

    snap = sim.snapshot_num
    archive_file = ARCHIVE_VR_HRF_ROOT / str(mass_bin) / f"VR-SOAP_HRF_{sim_label}{snap:04d}{mass_bin}"
    assert archive_file.is_file(), f"archive file missing: {archive_file}"

    cache_path = extract_vr(sim_label, mass_bin, refresh=False)
    result  = np.genfromtxt(cache_path)
    archive = np.genfromtxt(archive_file)

    assert result.shape == archive.shape, (
        f"shape mismatch: got {result.shape}, expected {archive.shape}"
    )
    np.testing.assert_allclose(result, archive, rtol=1e-5, atol=1e-10, equal_nan=True)
