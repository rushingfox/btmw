"""Regression tests: radial-density extract must match the archive caches.

VR tests (18 tests, login-node friendly):
  Reference: 250507 archive under
  ``working/.../paper_figures_min_20_HBTplus_version_250507/RadialDensityProfile_VR/``.

HBT tests (9 tests, SLOW — require full SWIFT snapshots, run on compute node):
  Reference: v1 cache files under
  ``working/.../m12i_{cdmo,btps_deep,btps_soft}_fof/hbtplus_min_20/``.
  Marked ``slow`` and auto-skipped when the SWIFT snapshot is absent.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from btmw.config import load_simulations
from btmw.figures.radial_density import extract_vr, extract_hbt, _vr_paths, _hbt_paths

ARCHIVE_VR = (
    load_simulations().working_root
    / "paper_figures_min_20_HBTplus_version_250507/RadialDensityProfile_VR"
)

# v1 HBT archive: cache lives inside each sim's hbtplus_min_20/ directory.
_V1_HBT_ROOT = load_simulations().working_root

def _v1_hbt_sim_path(label: str) -> Path:
    return _V1_HBT_ROOT / f"{label}_fof" / "hbtplus_min_20" / f"{label}_radial_density_simulation"

def _v1_hbt_nfw_path(label: str) -> Path:
    return _V1_HBT_ROOT / f"{label}_fof" / "hbtplus_min_20" / f"{label}_NFW_fitting"

def _v1_hbt_conv_path(label: str) -> Path:
    return _V1_HBT_ROOT / f"{label}_fof" / "hbtplus_min_20" / f"{label}_converg_radius"

# All 6 fiducial sims have VR profiles.
FIDUCIAL = (
    "m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft",
    "m12f_cdmo", "m12f_btps_deep", "m12f_btps_soft",
)

# Only m12i sims have v1 HBT cache.
M12I_LABELS = ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft")


def _have_vr(label: str) -> bool:
    sim = load_simulations()[label]
    vr_dir = sim.vr_raw_dir()
    if vr_dir is None:
        return False
    return (vr_dir / "output.profiles").is_file()


# ---------------------------------------------------------------------------
# simulation_1 cache
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label", FIDUCIAL)
def test_vr_extract_simulation_matches_archive(label):
    if not _have_vr(label):
        pytest.skip(f"VR raw files not available for {label}")

    extract_vr(label, refresh=False)

    sim_path, _, _ = _vr_paths(label)
    ours = np.genfromtxt(str(sim_path))
    ref = np.genfromtxt(str(ARCHIVE_VR / f"{label}_simulation_1"))

    assert ours.shape == ref.shape, (
        f"{label}: simulation_1 shape mismatch {ours.shape} vs {ref.shape}"
    )
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: simulation_1 values differ from 250507 archive"
    )


# ---------------------------------------------------------------------------
# NFW_1 cache
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label", FIDUCIAL)
def test_vr_extract_nfw_matches_archive(label):
    if not _have_vr(label):
        pytest.skip(f"VR raw files not available for {label}")

    extract_vr(label, refresh=False)

    _, nfw_path, _ = _vr_paths(label)
    ours = np.genfromtxt(str(nfw_path))
    ref = np.genfromtxt(str(ARCHIVE_VR / f"{label}_NFW_1"))

    assert ours.shape == ref.shape
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: NFW_1 values differ from 250507 archive"
    )


# ---------------------------------------------------------------------------
# converg_radius1 cache
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label", FIDUCIAL)
def test_vr_extract_converg_matches_archive(label):
    if not _have_vr(label):
        pytest.skip(f"VR raw files not available for {label}")

    extract_vr(label, refresh=False)

    _, _, conv_path = _vr_paths(label)
    ours = float(np.genfromtxt(str(conv_path)))
    ref = float(np.genfromtxt(str(ARCHIVE_VR / f"{label}_converg_radius1")))

    assert np.isclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: converg_radius1 {ours!r} != archive {ref!r}"
    )


# ---------------------------------------------------------------------------
# HBT regression tests (slow — require SWIFT snapshot on compute node)
# ---------------------------------------------------------------------------

def _have_hbt_snapshot(label: str) -> bool:
    """Return True only if the SWIFT snapshot for *label* exists on disk."""
    cfg = load_simulations()
    return cfg[label].snapshot_path().is_file()


@pytest.mark.slow
@pytest.mark.parametrize("label", M12I_LABELS)
def test_hbt_extract_simulation_matches_archive(label):
    """radial_density_simulation must match the v1 HBT archive file."""
    if not _have_hbt_snapshot(label):
        pytest.skip(f"SWIFT snapshot not available for {label} (run on compute node)")

    ref_path = _v1_hbt_sim_path(label)
    if not ref_path.is_file():
        pytest.skip(f"v1 HBT archive not available: {ref_path}")

    extract_hbt(label, refresh=False)

    sim_path, _, _ = _hbt_paths(label)
    ours = np.genfromtxt(str(sim_path))
    ref = np.genfromtxt(str(ref_path))

    assert ours.shape == ref.shape, (
        f"{label}: radial_density_simulation shape {ours.shape} != archive {ref.shape}"
    )
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: radial_density_simulation values differ from v1 HBT archive"
    )


@pytest.mark.slow
@pytest.mark.parametrize("label", M12I_LABELS)
def test_hbt_extract_nfw_matches_archive(label):
    """NFW_fitting must match the v1 HBT archive file."""
    if not _have_hbt_snapshot(label):
        pytest.skip(f"SWIFT snapshot not available for {label} (run on compute node)")

    ref_path = _v1_hbt_nfw_path(label)
    if not ref_path.is_file():
        pytest.skip(f"v1 HBT archive not available: {ref_path}")

    extract_hbt(label, refresh=False)

    _, nfw_path, _ = _hbt_paths(label)
    ours = np.genfromtxt(str(nfw_path))
    ref = np.genfromtxt(str(ref_path))

    assert ours.shape == ref.shape, (
        f"{label}: NFW_fitting shape {ours.shape} != archive {ref.shape}"
    )
    assert np.allclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: NFW_fitting values differ from v1 HBT archive"
    )


@pytest.mark.slow
@pytest.mark.parametrize("label", M12I_LABELS)
def test_hbt_extract_converg_matches_archive(label):
    """converg_radius must match the v1 HBT archive file."""
    if not _have_hbt_snapshot(label):
        pytest.skip(f"SWIFT snapshot not available for {label} (run on compute node)")

    ref_path = _v1_hbt_conv_path(label)
    if not ref_path.is_file():
        pytest.skip(f"v1 HBT archive not available: {ref_path}")

    extract_hbt(label, refresh=False)

    _, _, conv_path = _hbt_paths(label)
    ours = float(np.genfromtxt(str(conv_path)))
    ref = float(np.genfromtxt(str(ref_path)))

    assert np.isclose(ours, ref, rtol=1e-12, atol=0), (
        f"{label}: converg_radius {ours!r} != archive {ref!r}"
    )

