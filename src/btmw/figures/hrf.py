"""Fig 8: Halo Radial Function (HRF) in four Mbound mass bins.
Fig 13: HBT-HERONS vs VR-SOAP HRF comparison, same four bins.

Source notebooks (Type 1 extract + Type 2 comparison):
  paper_figures_min_20_HBTplus_v1_archive/HRF_auto_bin/{6,7,8,9}/
    HRF_hbtplus.ipynb                    (Type 1 – extract HBT)
    comparison_HRF_auto_bin_{n}_HBTplus.ipynb  (Type 2 – plot fig 8)
  paper_figures_min_20_HBTplus_v1_archive/comparison_with_VR-SOAP/HRF_auto_bin/{6,7,8,9}/
    HRF_VR-SOAP.ipynb                    (Type 1 – extract VR-SOAP)
    comparison_HRF_auto_bin_{n}_HBTplus_vs_VR.ipynb  (Type 2 – plot fig 13)
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Literal

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, rc
from scipy.interpolate import interp1d

from ..config import load_simulations
from ..io_soap import read_soap_bundle, subhalo_distances_to_main
from ..paths import cache_dir, data_dir, figures_dir

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MASS_BINS: tuple[int, ...] = (6, 7, 8, 9)

# Suptitle per mass bin (TeX strings matching the reference notebooks).
_SUPTITLES: dict[int, str] = {
    6: r"$10^6<{M_{\rm sub}}/{\rm M}_{\odot}<10^{7}$",
    7: r"$10^7<{M_{\rm sub}}/{\rm M}_{\odot}<10^{8}$",
    8: r"$10^8<{M_{\rm sub}}/{\rm M}_{\odot}<10^{9}$",
    9: r"$10^9<{M_{\rm sub}}/{\rm M}_{\odot}<10^{10}$",
}

# Aquarius A1 reference data: columns are (log10(r/kpc), log10(n(r)/<n>)).
# R200c of Aquarius A1 ≈ 245 kpc → x_plot = 10^(col0) / 245.
_AQUARIUS_R200C_KPC: float = 245.0
_AQUARIUS_FILE: Path = data_dir() / "external" / "springel2008_nprofile_r200c.txt"

# All 6 fiducial sims in plot order: m12i (cdmo, deep, soft), m12f (cdmo, deep, soft).
_ALL_LABELS: tuple[str, ...] = (
    "m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft",
    "m12f_cdmo", "m12f_btps_deep", "m12f_btps_soft",
)

_LINE_STYLES: tuple[str, ...] = ("C1", "C2", "C3", "C1--", "C2--", "C3--")

_LEGEND_LABELS_TEX: tuple[str, ...] = (
    r"$\rm{m12i:PL}$",
    r"$\rm{m12i:BT\_deep}$",
    r"$\rm{m12i:BT\_soft}$",
    r"$\rm{m12f:PL}$",
    r"$\rm{m12f:BT\_deep}$",
    r"$\rm{m12f:BT\_soft}$",
    r"$\rm{Aquarius\ A1}$",
)

# ---------------------------------------------------------------------------
# Cache paths
# ---------------------------------------------------------------------------

def _hrf_cache_dir() -> Path:
    return cache_dir() / "hrf"


def _cache_path(sim_label: str, mass_bin: int) -> Path:
    return _hrf_cache_dir() / f"{sim_label}_bin{mass_bin}"


# ---------------------------------------------------------------------------
# Extract
# ---------------------------------------------------------------------------

def extract(sim_label: str, mass_bin: int, *, refresh: bool = False) -> Path:
    """Compute n(r)/<n> profile for *sim_label* in Mbound bin *mass_bin*.

    Parameters
    ----------
    sim_label : str
        One of the 6 fiducial simulation labels.
    mass_bin : int
        Log10 of the lower mass edge in solar masses (6, 7, 8, or 9).
    refresh : bool
        Force recomputation even if cache exists.

    Returns
    -------
    Path
        Path to the written cache file.
    """
    if mass_bin not in MASS_BINS:
        raise ValueError(f"mass_bin must be one of {MASS_BINS}, got {mass_bin!r}")

    out_path = _cache_path(sim_label, mass_bin)
    if out_path.is_file() and not refresh:
        return out_path

    cfg = load_simulations()
    sim = cfg[sim_label]
    soap_path = sim.soap_hbt_path()

    bundle = read_soap_bundle(soap_path, need_vmax=False)

    main_idx = bundle.main_idx
    rcut_mpc  = float(bundle.r50c[main_idx])    # cutoff  = R50c  [Mpc]
    rscale_mpc = float(bundle.r200c[main_idx])  # scale   = R200c [Mpc]

    # Mass filter: Mbound in [1e10 Msun] units.
    # bin n → 10^(n-10) < Mbound < 10^(n-9)
    m_lo = 10.0 ** (mass_bin - 10)
    m_hi = m_lo * 10.0
    mbound = bundle.mbound

    dist_mpc = subhalo_distances_to_main(bundle)  # Mpc, length N_halos

    # Mask: not main halo, within cutoff, in mass range.
    mask = (
        (np.arange(len(mbound)) != main_idx)
        & (dist_mpc < rcut_mpc)
        & (mbound > m_lo)
        & (mbound < m_hi)
    )
    dist_scaled = dist_mpc[mask] / rscale_mpc   # r / R200c
    n_within = int(mask.sum())

    if n_within == 0:
        # No halos in this bin – write empty file so cache exists.
        out_path.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(str(out_path), np.zeros((0, 2)))
        return out_path

    countdata, bins = np.histogram(dist_scaled, bins="auto")

    lhs = bins[:-1]
    rhs = bins[1:]
    center_bins = 0.5 * (lhs + rhs)

    # Shell volumes in Mpc³ (rscale converts r/R200c back to Mpc).
    bin_volumes = np.array([
        (4.0 / 3.0) * math.pi * ((rhs[i] * rscale_mpc) ** 3 - (lhs[i] * rscale_mpc) ** 3)
        for i in range(len(bins) - 1)
    ])

    n_average = n_within / ((4.0 / 3.0) * math.pi * rcut_mpc ** 3)

    density_normalized = countdata / bin_volumes / n_average

    out = np.stack((center_bins, density_normalized), axis=-1)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(str(out_path), out)
    return out_path


# ---------------------------------------------------------------------------
# Interpolation helper (same as v1 notebook's interp_data)
# ---------------------------------------------------------------------------

def _interp_data(data: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x = data[:, 0]
    y = data[:, 1]
    func = interp1d(x[y != 0], y[y != 0], kind="linear", fill_value="extrapolate")
    return np.column_stack((x_new, func(x_new)))


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(
    mass_bin: int,
    *,
    refresh: bool = False,
    output: str | Path | None = None,
    use_tex: bool = True,
) -> Path:
    """Produce one panel of Fig 8 (HRF) for the given mass bin.

    Parameters
    ----------
    mass_bin : int
        Log10 of the lower Mbound edge in solar masses: 6, 7, 8, or 9.
    refresh : bool
        Force re-extraction before plotting.
    output : path-like, optional
        Output PNG path. Defaults to ``figures/HRF_{mass_bin}_auto_bin_HBTplus.png``.
    use_tex : bool
        Enable LaTeX + Helvetica rendering (default True, matching reference).
    """
    if mass_bin not in MASS_BINS:
        raise ValueError(f"mass_bin must be one of {MASS_BINS}, got {mass_bin!r}")

    # Extract cache for all 6 sims.
    for sl in _ALL_LABELS:
        extract(sl, mass_bin, refresh=refresh)

    # TeX setup.
    if use_tex:
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
        rc("text", usetex=True)

    # Load caches.
    sim_arrays: list[np.ndarray] = []
    for sl in _ALL_LABELS:
        p = _cache_path(sl, mass_bin)
        sim_arrays.append(np.genfromtxt(str(p)))

    # Load Aquarius.
    aq_raw = np.genfromtxt(str(_AQUARIUS_FILE))
    aq_x = 10.0 ** (aq_raw[:, 0] - math.log10(_AQUARIUS_R200C_KPC))
    aq_y = 10.0 ** aq_raw[:, 1]

    # Build figure.
    fig = plt.figure(figsize=(6, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # ----- Main panel -----
    for data, label, ls in zip(sim_arrays, _LEGEND_LABELS_TEX[:-1], _LINE_STYLES):
        d = data[data[:, 1] != 0] if data.size > 0 else data
        ax0.plot(d[:, 0], d[:, 1], ls, label=label)
    ax0.plot(aq_x, aq_y, "grey", label=_LEGEND_LABELS_TEX[-1])

    ax0.set_xlim(0.1, 3)
    ax1.set_xlim(0.1, 3)
    ax0.set_ylim(0.8e-1, 60)
    ax1.set_ylim(0.35, 3.5)

    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax0.set_ylabel(r"$n\left(r\right)/\left<n\right>$", fontsize=20)
    ax0.legend(fontsize=13)

    # ----- Ratio panel -----
    # m12i: interpolate BT sims onto cdmo x-grid.
    m12i_cdmo = sim_arrays[0]
    m12i_deep = sim_arrays[1]
    m12i_soft = sim_arrays[2]
    m12f_cdmo = sim_arrays[3]
    m12f_deep = sim_arrays[4]
    m12f_soft = sim_arrays[5]

    if m12i_cdmo.size > 0:
        x_i = m12i_cdmo[:, 0]
        deep_i = _interp_data(m12i_deep, x_i)
        soft_i = _interp_data(m12i_soft, x_i)
        ratio_i_deep = np.column_stack((x_i, deep_i[:, 1] / m12i_cdmo[:, 1]))
        ratio_i_soft = np.column_stack((x_i, soft_i[:, 1] / m12i_cdmo[:, 1]))
        ax1.plot(ratio_i_deep[:, 0], ratio_i_deep[:, 1], _LINE_STYLES[1],
                 label=r"$\rm{m12i\_deep}$" if use_tex else "m12i_deep")
        ax1.plot(ratio_i_soft[:, 0], ratio_i_soft[:, 1], _LINE_STYLES[2],
                 label=r"$\rm{m12i\_soft}$" if use_tex else "m12i_soft")

    if m12f_cdmo.size > 0:
        x_f = m12f_cdmo[:, 0]
        deep_f = _interp_data(m12f_deep, x_f)
        soft_f = _interp_data(m12f_soft, x_f)
        ratio_f_deep = np.column_stack((x_f, deep_f[:, 1] / m12f_cdmo[:, 1]))
        ratio_f_soft = np.column_stack((x_f, soft_f[:, 1] / m12f_cdmo[:, 1]))
        ax1.plot(ratio_f_deep[:, 0], ratio_f_deep[:, 1], _LINE_STYLES[4],
                 label=r"$\rm{m12f\_deep}$" if use_tex else "m12f_deep")
        ax1.plot(ratio_f_soft[:, 0], ratio_f_soft[:, 1], _LINE_STYLES[5],
                 label=r"$\rm{m12f\_soft}$" if use_tex else "m12f_soft")

    ax1.set_ylabel(r"$N_{\rm BT}/N_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", label=r"$\rm{y=1}$" if use_tex else "y=1", alpha=0.3)

    # Tick styling.
    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"$r/R_{200c}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    plt.suptitle(_SUPTITLES[mass_bin], fontsize=15, y=1.02, ha="center", x=0.6)

    # Save.
    if output is None:
        out_path = figures_dir() / f"HRF_{mass_bin}_auto_bin_HBTplus.png"
    else:
        out_path = Path(output)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Fig 18 — HRF resolution study (fiducial vs high-res m12i)
# ---------------------------------------------------------------------------

_RES_FIDUCIAL_LABELS = ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft")
_RES_HIGH_LABELS     = ("m12i_cdmo_highest", "m12i_btps_deep_highest", "m12i_btps_soft_highest")
_RES_FIDUCIAL_LS     = ("C1",  "C2",  "C3")
_RES_HIGH_LS         = ("C1-.", "C2-.", "C3-.")
_RES_FIDUCIAL_TEX = (
    r"$\rm{m12i:PL-fiducial\ res}$",
    r"$\rm{m12i:BT\_deep-fiducial\ res}$",
    r"$\rm{m12i:BT\_soft-fiducial\ res}$",
)
_RES_HIGH_TEX = (
    r"$\rm{m12i:PL-high\ res}$",
    r"$\rm{m12i:BT\_deep-high\ res}$",
    r"$\rm{m12i:BT\_soft-high\ res}$",
)
_RES_RATIO_TEX = (r"$\rm{PL}$", r"$\rm{BT\_deep}$", r"$\rm{BT\_soft}$")


def plot_resolution_study(
    mass_bin: int,
    *,
    refresh: bool = False,
    output: str | Path | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 18 — HRF resolution study for one mass bin (fiducial vs high-res m12i)."""
    if mass_bin not in MASS_BINS:
        raise ValueError(f"mass_bin must be one of {MASS_BINS}, got {mass_bin!r}")

    # Ensure caches are built for all 6 sims
    for sl in (*_RES_FIDUCIAL_LABELS, *_RES_HIGH_LABELS):
        extract(sl, mass_bin, refresh=refresh)

    if use_tex:
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
        rc("text", usetex=True)

    # Load fiducial + high-res arrays
    fid_arrays  = [np.genfromtxt(str(_cache_path(sl, mass_bin))) for sl in _RES_FIDUCIAL_LABELS]
    high_arrays = [np.genfromtxt(str(_cache_path(sl, mass_bin))) for sl in _RES_HIGH_LABELS]

    # Load Aquarius
    aq_raw = np.genfromtxt(str(_AQUARIUS_FILE))
    aq_x = 10.0 ** (aq_raw[:, 0] - math.log10(_AQUARIUS_R200C_KPC))
    aq_y = 10.0 ** aq_raw[:, 1]

    fig = plt.figure(figsize=(6, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # Main panel: fiducial (solid) + high-res (dash-dot) + Aquarius (grey)
    for i, (ls, tex) in enumerate(zip(_RES_FIDUCIAL_LS, _RES_FIDUCIAL_TEX)):
        d = fid_arrays[i]
        d_nz = d[d[:, 1] != 0] if d.size > 0 else d
        ax0.plot(d_nz[:, 0], d_nz[:, 1], ls, label=tex)
    for i, (ls, tex) in enumerate(zip(_RES_HIGH_LS, _RES_HIGH_TEX)):
        d = high_arrays[i]
        d_nz = d[d[:, 1] != 0] if d.size > 0 else d
        ax0.plot(d_nz[:, 0], d_nz[:, 1], ls, label=tex)
    ax0.plot(aq_x, aq_y, "grey", label=_LEGEND_LABELS_TEX[-1])

    ax0.set_xlim(0.1, 3)
    ax1.set_xlim(0.1, 3)
    ax0.set_ylim(0.8e-1, 60)
    ax1.set_ylim(0.35, 3.5)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax0.set_ylabel(r"${n\left(r\right)}/{\left<n\right>}$", fontsize=20)
    ax0.legend(fontsize=12)

    # Ratio panel: N_high / N_fiducial (interpolate high onto fiducial grid)
    for i, (ls, tex) in enumerate(zip(_RES_FIDUCIAL_LS, _RES_RATIO_TEX)):
        fid = fid_arrays[i]
        high = high_arrays[i]
        if fid.size > 0 and high.size > 0:
            x_fid = fid[:, 0]
            high_interp = _interp_data(high, x_fid)
            ratio = np.column_stack((x_fid, high_interp[:, 1] / fid[:, 1]))
            ax1.plot(ratio[:, 0], ratio[:, 1], ls, label=tex)

    ax1.set_ylabel(r"${N_{\rm high}}/{N_{\rm fiducial}}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    ax1.legend(fontsize=9)

    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"${r}/{R_{200c}}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.suptitle(_SUPTITLES[mass_bin], fontsize=15, y=1.02, ha="center", x=0.6)

    if output is None:
        out_path = figures_dir() / f"HRF_resolution_study_{mass_bin}.png"
    else:
        out_path = Path(output)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path

# m12i sims only (archive only has m12i comparison)
_VR_LABELS: tuple[str, ...] = ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft")
_VR_LINE_STYLES: tuple[str, ...] = ("C1", "C2", "C3", "C1-.", "C2-.", "C3-.")
_VR_LEGEND_HBT: tuple[str, ...] = (
    r"$\rm{m12i:PL \sim HBT-HERONS}$",
    r"$\rm{m12i:BT\_deep \sim HBT-HERONS}$",
    r"$\rm{m12i:BT\_soft \sim HBT-HERONS}$",
)
_VR_LEGEND_VR: tuple[str, ...] = (
    r"$\rm{m12i:PL \sim VR}$",
    r"$\rm{m12i:BT\_deep \sim VR}$",
    r"$\rm{m12i:BT\_soft \sim VR}$",
)


def _vr_cache_path(sim_label: str, mass_bin: int) -> Path:
    cfg = load_simulations()
    sim = cfg[sim_label]
    snap = sim.snapshot_num
    return _hrf_cache_dir() / f"VR-SOAP_HRF_{sim_label}{snap:04d}{mass_bin}"


def extract_vr(sim_label: str, mass_bin: int, *, refresh: bool = False) -> Path:
    """Compute HRF from VR-SOAP catalog for *sim_label* in bin *mass_bin*."""
    if mass_bin not in MASS_BINS:
        raise ValueError(f"mass_bin must be one of {MASS_BINS}, got {mass_bin!r}")

    cfg = load_simulations()
    sim = cfg[sim_label]
    if not sim.has_vr:
        raise ValueError(f"{sim_label} has no VR data")

    out_path = _vr_cache_path(sim_label, mass_bin)
    if out_path.is_file() and not refresh:
        return out_path

    bundle = read_soap_bundle(sim.soap_vr_path(), need_vmax=False)

    main_idx   = bundle.main_idx
    rcut_mpc   = float(bundle.r50c[main_idx])
    rscale_mpc = float(bundle.r200c[main_idx])

    m_lo = 10.0 ** (mass_bin - 10)
    m_hi = m_lo * 10.0
    mbound = bundle.mbound
    dist_mpc = subhalo_distances_to_main(bundle)

    mask = (
        (np.arange(len(mbound)) != main_idx)
        & (dist_mpc < rcut_mpc)
        & (mbound > m_lo)
        & (mbound < m_hi)
    )
    dist_scaled = dist_mpc[mask] / rscale_mpc
    n_within = int(mask.sum())

    if n_within == 0:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(str(out_path), np.zeros((0, 2)))
        return out_path

    countdata, bins = np.histogram(dist_scaled, bins="auto")
    lhs = bins[:-1]
    rhs = bins[1:]
    center_bins = 0.5 * (lhs + rhs)

    bin_volumes = np.array([
        (4.0 / 3.0) * math.pi * ((rhs[i] * rscale_mpc) ** 3 - (lhs[i] * rscale_mpc) ** 3)
        for i in range(len(bins) - 1)
    ])
    n_average = n_within / ((4.0 / 3.0) * math.pi * rcut_mpc ** 3)
    density_normalized = countdata / bin_volumes / n_average

    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(str(out_path), np.stack((center_bins, density_normalized), axis=-1))
    return out_path


def plot_compare_vr(
    mass_bin: int,
    *,
    refresh: bool = False,
    output: str | Path | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 13: HBT-HERONS vs VR-SOAP HRF comparison for one mass bin."""
    if mass_bin not in MASS_BINS:
        raise ValueError(f"mass_bin must be one of {MASS_BINS}, got {mass_bin!r}")

    # Ensure caches exist
    for sl in _VR_LABELS:
        extract(sl, mass_bin, refresh=refresh)
        extract_vr(sl, mass_bin, refresh=refresh)

    if use_tex:
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
        rc("text", usetex=True)

    hbt_arrays = [np.genfromtxt(str(_cache_path(sl, mass_bin))) for sl in _VR_LABELS]
    vr_arrays  = [np.genfromtxt(str(_vr_cache_path(sl, mass_bin))) for sl in _VR_LABELS]

    aq_raw = np.genfromtxt(str(_AQUARIUS_FILE))
    aq_x = 10.0 ** (aq_raw[:, 0] - math.log10(_AQUARIUS_R200C_KPC))
    aq_y = 10.0 ** aq_raw[:, 1]

    fig = plt.figure(figsize=(6, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    for i, (tex, ls) in enumerate(zip(_VR_LEGEND_HBT, _VR_LINE_STYLES[:3])):
        d = hbt_arrays[i]
        d = d[d[:, 1] != 0] if d.size > 0 else d
        ax0.plot(d[:, 0], d[:, 1], ls, label=tex)

    for i, (tex, ls) in enumerate(zip(_VR_LEGEND_VR, _VR_LINE_STYLES[3:])):
        d = vr_arrays[i]
        d = d[d[:, 1] != 0] if d.size > 0 else d
        ax0.plot(d[:, 0], d[:, 1], ls, label=tex)

    ax0.plot(aq_x, aq_y, "grey", label=r"$\rm{Aquarius\ A1}$")

    ax0.set_xlim(0.1, 3)
    ax1.set_xlim(0.1, 3)
    ax0.set_ylim(0.8e-1, 60)
    ax1.set_ylim(0.35, 3.5)
    ax0.set_xscale("log"); ax0.set_yscale("log")
    ax1.set_xscale("log"); ax1.set_yscale("log")

    ax0.set_ylabel(r"${n\left(r\right)}/{\left<n\right>}$", fontsize=20)
    ax0.legend(fontsize=11)

    # Ratio panel: VR / HBT (interpolate VR onto HBT x-grid)
    for i, (ls, tex) in enumerate(zip(_VR_LINE_STYLES[:3], (
        r"$\rm{m12i:PL}$", r"$\rm{m12i:BT\_deep}$", r"$\rm{m12i:BT\_soft}$"
    ))):
        hbt = hbt_arrays[i]
        vr  = vr_arrays[i]
        if hbt.size > 0 and vr.size > 0:
            vr_interp = _interp_data(vr, hbt[:, 0])
            ratio = np.column_stack((hbt[:, 0], vr_interp[:, 1] / hbt[:, 1]))
            ax1.plot(ratio[:, 0], ratio[:, 1], ls, label=tex)

    ax1.set_ylabel(r"${N_{VR}}/{N_{HBT-HERONS}}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    ax1.legend(fontsize=10)

    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"${r}/{R_{200c}}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.suptitle(_SUPTITLES[mass_bin], fontsize=15, y=1.02, ha="center", x=0.6)

    if output is None:
        out_path = figures_dir() / f"HRF_{mass_bin}_auto_bin_HBT-HERONS_Comparison_with_VR-SOAP.png"
    else:
        out_path = Path(output)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path
