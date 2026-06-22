"""Regression tests: btmw HMF/HVF extracts must match the original cache files
byte-for-byte (or at minimum np.allclose).

The reference cache lives under the ``paper_figures_min_20_HBTplus_v1_archive/``
directory inside the simulation root.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures import hmf, hvf

ARCHIVE = load_simulations().working_root / "paper_figures_min_20_HBTplus_v1_archive"

FIDUCIAL = ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft",
            "m12f_cdmo", "m12f_btps_deep", "m12f_btps_soft")


def _have_soap(label: str) -> bool:
    sim = load_simulations()[label]
    return sim.soap_hbt_path().is_file()


@pytest.mark.parametrize("label", FIDUCIAL)
def test_hmf_extract_matches_archive(label):
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label}")
    sim = load_simulations()[label]
    out = hmf.extract(label, refresh=False)
    ours = np.genfromtxt(out)
    ref = np.genfromtxt(ARCHIVE / "HMF" / f"hbtplus_HMF_{label}{sim.snapshot_num:04d}")
    assert ours.shape == ref.shape
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: HMF extract differs from archive"
    )


@pytest.mark.parametrize("label", FIDUCIAL)
def test_hvf_extract_matches_archive(label):
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label}")
    sim = load_simulations()[label]
    out = hvf.extract(label, refresh=False)
    ours = np.genfromtxt(out)
    ref = np.genfromtxt(ARCHIVE / "HVF_r100" / f"hbtplus_HVF_{label}{sim.snapshot_num:04d}")
    assert ours.shape == ref.shape
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: HVF extract differs from archive"
    )
