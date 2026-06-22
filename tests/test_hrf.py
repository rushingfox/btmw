"""Regression tests: HRF extract must match the v1 archive cache files.

Reference caches live under:
  paper_figures_min_20_HBTplus_v1_archive/HRF_auto_bin/{bin}/
    hbtplus_HRF_{label}{snap:04d}{bin}

All tests are marked ``slow`` because they require reading SOAP files (~GB).
They are auto-skipped when the SOAP file is absent (login-node safe).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures.hrf import extract, MASS_BINS, _cache_path

ARCHIVE_HRF = (
    load_simulations().working_root
    / "paper_figures_min_20_HBTplus_v1_archive/HRF_auto_bin"
)

# All 6 fiducial sims with their snap numbers.
FIDUCIAL_SNAPS = {
    "m12i_cdmo":      50,
    "m12i_btps_deep": 50,
    "m12i_btps_soft": 50,
    "m12f_cdmo":      49,
    "m12f_btps_deep": 49,
    "m12f_btps_soft": 49,
}


def _archive_path(label: str, mass_bin: int) -> Path:
    snap = FIDUCIAL_SNAPS[label]
    return ARCHIVE_HRF / str(mass_bin) / f"hbtplus_HRF_{label}{snap:04d}{mass_bin}"


def _have_soap(label: str) -> bool:
    cfg = load_simulations()
    return cfg[label].soap_hbt_path().is_file()


# ---------------------------------------------------------------------------
# Parametrize over all (label, mass_bin) combinations = 6 × 4 = 24 tests
# ---------------------------------------------------------------------------

PARAMS = [(label, mass_bin) for label in FIDUCIAL_SNAPS for mass_bin in MASS_BINS]


@pytest.mark.slow
@pytest.mark.parametrize("label,mass_bin", PARAMS)
def test_hrf_extract_matches_archive(label: str, mass_bin: int) -> None:
    """n(r)/<n> profile must match the v1 HRF archive to rtol=1e-12."""
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label} (run on compute node)")

    ref_path = _archive_path(label, mass_bin)
    if not ref_path.is_file():
        pytest.skip(f"v1 archive not available: {ref_path}")

    extract(label, mass_bin, refresh=False)

    ours = np.genfromtxt(str(_cache_path(label, mass_bin)))
    ref  = np.genfromtxt(str(ref_path))

    assert ours.shape == ref.shape, (
        f"{label} bin{mass_bin}: shape {ours.shape} != archive {ref.shape}"
    )
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label} bin{mass_bin}: values differ from v1 archive"
    )
