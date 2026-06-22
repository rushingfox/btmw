"""Regression tests for figs 15-18: resolution study extracts must match the
archive cache files in resolution_study/{HMF,HVF,RvsV,HRF_auto_bin/}.

Tests are auto-skipped when the SOAP file is absent (login-node safe).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures import hmf, hvf, hrf
from btmw.figures.rvsv import extract as rvsv_extract

ARCHIVE = (
    load_simulations().working_root
    / "paper_figures_min_20_HBTplus_v1_archive/resolution_study"
)

HIGH_LABELS = (
    "m12i_cdmo_highest",
    "m12i_btps_deep_highest",
    "m12i_btps_soft_highest",
)


def _have_soap(label: str) -> bool:
    try:
        sim = load_simulations()[label]
        return sim.soap_hbt_path().is_file()
    except KeyError:
        return False


# ---------------------------------------------------------------------------
# Fig 15 — HMF high-res extract vs archive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label", HIGH_LABELS)
@pytest.mark.slow
def test_hmf_highest_extract_matches_archive(label: str) -> None:
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label}")
    sim = load_simulations()[label]
    out = hmf.extract(label, refresh=False)
    ours = np.genfromtxt(out)
    ref = np.genfromtxt(ARCHIVE / "HMF" / f"hbtplus_HMF_{label}{sim.snapshot_num:04d}")
    assert ours.shape == ref.shape, f"{label}: shape mismatch"
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), f"{label}: HMF extract differs from archive"


# ---------------------------------------------------------------------------
# Fig 16 — HVF high-res extract vs archive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label", HIGH_LABELS)
@pytest.mark.slow
def test_hvf_highest_extract_matches_archive(label: str) -> None:
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label}")
    sim = load_simulations()[label]
    out = hvf.extract(label, refresh=False)
    ours = np.genfromtxt(out)
    ref = np.genfromtxt(ARCHIVE / "HVF" / f"hbtplus_HVF_{label}{sim.snapshot_num:04d}")
    assert ours.shape == ref.shape, f"{label}: shape mismatch"
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), f"{label}: HVF extract differs from archive"


# ---------------------------------------------------------------------------
# Fig 17 — RvsV high-res extract vs archive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "label,which",
    [(lbl, w) for lbl in HIGH_LABELS for w in ("Lower", "Median", "Upper")],
)
@pytest.mark.slow
def test_rvsv_highest_extract_matches_archive(label: str, which: str) -> None:
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label}")
    archive_file = ARCHIVE / "RvsV" / f"{label}_RvsV_{which}"
    assert archive_file.is_file(), f"archive file missing: {archive_file}"
    paths = rvsv_extract(label, refresh=False)
    ours = np.genfromtxt(paths[which])
    ref = np.genfromtxt(archive_file)
    assert ours.shape == ref.shape, f"{label} {which}: shape mismatch"
    assert np.allclose(ours, ref, equal_nan=True, rtol=1e-5, atol=1e-10), (
        f"{label} {which}: RvsV extract differs from archive"
    )


# ---------------------------------------------------------------------------
# Fig 18 — HRF high-res extract vs archive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "label,mass_bin",
    [(lbl, mb) for lbl in HIGH_LABELS for mb in hrf.MASS_BINS],
)
@pytest.mark.slow
def test_hrf_highest_extract_matches_archive(label: str, mass_bin: int) -> None:
    if not _have_soap(label):
        pytest.skip(f"SOAP file not available for {label}")
    sim = load_simulations()[label]
    out = hrf.extract(label, mass_bin, refresh=False)
    ours = np.genfromtxt(str(out))
    ref_path = (
        ARCHIVE / "HRF_auto_bin" / str(mass_bin)
        / f"hbtplus_HRF_{label}{sim.snapshot_num:04d}{mass_bin}"
    )
    assert ref_path.is_file(), f"archive file missing: {ref_path}"
    ref = np.genfromtxt(str(ref_path))
    assert ours.shape == ref.shape, f"{label} bin{mass_bin}: shape mismatch"
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label} bin{mass_bin}: HRF extract differs from archive"
    )


# ---------------------------------------------------------------------------
# Smoke: plot_resolution_study functions importable and return a Path
# (no SOAP required — uses pre-existing cache)
# ---------------------------------------------------------------------------

def test_hmf_plot_resolution_study_has_cache(tmp_path):
    """plot_resolution_study must accept output kwarg and return a Path."""
    out = tmp_path / "test_hmf_res.png"
    # Skip if any required cache is missing
    cfg = load_simulations()
    for lbl in (*hmf._RES_FIDUCIAL_LABELS, *hmf._RES_HIGH_LABELS):
        sim = cfg[lbl]
        cache = hmf._cache_path(lbl, sim.snapshot_num)
        if not cache.is_file():
            pytest.skip(f"cache missing for {lbl}")
    result = hmf.plot_resolution_study(output=str(out), use_tex=False)
    assert result == out
    assert out.is_file()


def test_hvf_plot_resolution_study_has_cache(tmp_path):
    out = tmp_path / "test_hvf_res.png"
    cfg = load_simulations()
    for lbl in (*hvf._RES_FIDUCIAL_LABELS, *hvf._RES_HIGH_LABELS):
        sim = cfg[lbl]
        cache = hvf._cache_path(lbl, sim.snapshot_num)
        if not cache.is_file():
            pytest.skip(f"cache missing for {lbl}")
    result = hvf.plot_resolution_study(output=str(out), use_tex=False)
    assert result == out
    assert out.is_file()


def test_hrf_plot_resolution_study_has_cache(tmp_path):
    from btmw.figures.hrf import _cache_path as hrf_cache_path
    cfg = load_simulations()
    for lbl in (*hrf._RES_FIDUCIAL_LABELS, *hrf._RES_HIGH_LABELS):
        for mb in hrf.MASS_BINS:
            if not hrf_cache_path(lbl, mb).is_file():
                pytest.skip(f"cache missing for {lbl} bin{mb}")
    for mb in hrf.MASS_BINS:
        out = tmp_path / f"test_hrf_res_{mb}.png"
        result = hrf.plot_resolution_study(mb, output=str(out), use_tex=False)
        assert result == out
        assert out.is_file()


def test_rvsv_plot_resolution_study_has_cache(tmp_path):
    from btmw.figures import rvsv
    from btmw.figures.rvsv import _cache_path as rv_cache
    for lbl in rvsv._RES_HIGH_LABELS:
        for w in ("Lower", "Median", "Upper"):
            if not rv_cache(lbl, w).is_file():
                pytest.skip(f"cache missing for {lbl} {w}")
    out = tmp_path / "test_rvsv_res.png"
    result = rvsv.plot_resolution_study(output=str(out), use_tex=False)
    assert result == out
    assert out.is_file()
