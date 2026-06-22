"""Regression tests for Fig. 14 MvsV checked-in VELOCIraptor caches."""

from __future__ import annotations

import numpy as np
import pytest

from btmw.figures.mvsv import _vr_cache_path

ALL_SIMS = [
    "m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft",
    "m12f_cdmo", "m12f_btps_deep", "m12f_btps_soft",
]
WHICH = ["Lower", "Median", "Upper"]


@pytest.mark.parametrize("sim_label", ALL_SIMS)
@pytest.mark.parametrize("which", WHICH)
def test_mvsv_vr_cache_present_and_finite(sim_label: str, which: str) -> None:
    """Fig. 14 quick path should have a checked-in cache for every host/model."""
    cache = _vr_cache_path(sim_label, which)
    assert cache.is_file(), f"missing MvsV cache: {cache}"

    data = np.genfromtxt(str(cache))
    assert data.ndim == 2
    assert data.shape[1] == 2
    assert np.all(np.isfinite(data[:, 0]))
    assert np.any(np.isfinite(data[:, 1]))
