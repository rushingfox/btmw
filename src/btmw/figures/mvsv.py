"""Fig 14: M_vir vs V_max relation (MvsV).

Source notebooks:
  paper_figures_min_20_HBTplus_v1_archive/comparison_with_VR-SOAP/MvsV/
    MvsV.ipynb                         (Type 1 – extract VR from old catalog)
    MvsV_m12i.ipynb                    (Type 2 – plot VR m12i)
    MvsV_m12f.ipynb                    (Type 2 – plot VR m12f)

Algorithm (VR, Type 1):
    * Open the old VELOCIraptor ``output.properties`` catalog.
    * Use ``Vmax`` and ``Mass_BN98`` fields, excluding the main halo.
    * Fixed log-V bins: 10^arange(floor(log10(Vmin)), ceil(log10(Vmax))+0.05, 0.05).
    * For each bin: 16th / 50th / 84th percentile of M_converted (NaN if empty).
    * Output: ``{label}_VvsM_Msub_{Lower,Median,Upper}``.

Plot (Type 2):
  * Two-panel figure (height_ratios=[7, 3]), figsize=(6, 8), for m12i or m12f.
  * Main panel: 3 sim Median lines + fill_between(Lower/h, Upper/h).
    * Reference: analytic BolshoiP+MDPL/Rockstar power-law (black), hard-coded
        below rather than loaded from data/external/.
    * Unresolved axvspan from x1=0.5 to x2 (see _UNRESOLVED_X2_VR).
  * Ratio panel: BT_deep/PL and BT_soft/PL (Median only).
"""

from __future__ import annotations

import math
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, rc
from scipy.interpolate import interp1d

from ..config import load_simulations
from ..paths import cache_dir, figures_dir

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# h factor used for mass conversion (WMAP7 SWIFT sims)
_H_FACTOR = 0.702
# BolshoiP + MDPL h for the analytic reference line (Grand & White 2021 / Robert result)
# "Rockstar" here refers to the halo finder used in BolshoiP+MDPL, not the data file.
_H_BOLSHOI = 0.678

# Velocity bin width in log10 (same as HVF / RvsV)
_V_RESO = 0.1

# Percentile levels
_LOWER_PCT  = 16
_MEDIAN_PCT = 50
_UPPER_PCT  = 84

# left xlim for all MvsV plots
_X1 = 5e-1  # 0.5 km/s

# Unresolved region x2 (km/s) for fig 14 (VR-SOAP data)
# m12i: from comparison_with_VR-SOAP/MvsV/MvsV_m12i.ipynb
# m12f: from comparison_with_VR-SOAP/MvsV/MvsV_m12f.ipynb
_UNRESOLVED_X2_VR = {
    "m12i": 0.03009951542632502 * 14203994.722145732 / 1e5,   # ≈ 4.275
    "m12f": 0.024061509353586898 * 15879541.199106842 / 1e5,  # ≈ 3.820
}

# Sim groups per host
_HOST_SIMS = {
    "m12i": ["m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft"],
    "m12f": ["m12f_cdmo", "m12f_btps_deep", "m12f_btps_soft"],
}

_LABEL_MAP = {
    "m12i_cdmo":      r"$\rm{m12i:PL}$",
    "m12i_btps_deep": r"$\rm{m12i:BT\_deep}$",
    "m12i_btps_soft": r"$\rm{m12i:BT\_soft}$",
    "m12f_cdmo":      r"$\rm{m12f:PL}$",
    "m12f_btps_deep": r"$\rm{m12f:BT\_deep}$",
    "m12f_btps_soft": r"$\rm{m12f:BT\_soft}$",
}

_LINE_STYLES = ["C1", "C2", "C3"]


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _vr_cache_path(sim_label: str, which: str) -> Path:
    """Path for old-VR cache: {sim_label}_VvsM_Msub_{which} (matches archive naming)."""
    return cache_dir() / "mvsv" / f"{sim_label}_VvsM_Msub_{which}"


# ---------------------------------------------------------------------------
# Core binning logic
# ---------------------------------------------------------------------------

def _bin_mvsv(
    vmax_sel: np.ndarray,
    m_converted: np.ndarray,
    v_reso: float = _V_RESO,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Bin M-converted values by Vmax and compute 16/50/84 percentiles.

    Returns (v_center, M_lower, M_median, M_upper) — same length, NaN for
    empty bins.  Follows the notebook algorithm verbatim.
    """
    lv_min = math.floor(math.log10(vmax_sel.min()))
    lv_max = math.ceil(math.log10(vmax_sel.max()))
    bin_edges = 10.0 ** np.arange(lv_min, lv_max + v_reso, v_reso)
    v_center  = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    n_bins    = len(v_center)

    binned_m = [np.empty(0) for _ in range(n_bins)]
    for v, m in zip(vmax_sel, m_converted):
        idx = math.floor(math.log10(v) * (1.0 / v_reso)) - math.floor(lv_min * (1.0 / v_reso))
        if 0 <= idx < n_bins:
            binned_m[idx] = np.append(binned_m[idx], m)

    m_lower  = np.array([
        np.percentile(b, _LOWER_PCT)  if b.size > 0 else np.nan
        for b in binned_m
    ])
    m_median = np.array([
        np.percentile(b, _MEDIAN_PCT) if b.size > 0 else np.nan
        for b in binned_m
    ])
    m_upper  = np.array([
        np.percentile(b, _UPPER_PCT)  if b.size > 0 else np.nan
        for b in binned_m
    ])
    return v_center, m_lower, m_median, m_upper


# ---------------------------------------------------------------------------
# Extract: VR (old VelocIRaptor catalog — output.properties)
# ---------------------------------------------------------------------------

# VR bin resolution matches archive notebook (VResoInHVF = 0.05)
_VR_V_RESO = 0.05


def extract_vr(sim_label: str, *, refresh: bool = False) -> dict[str, Path]:
    """Compute MvsV percentile bands from old VelocIRaptor catalog.

    Reads ``vr_raw_dir/output.properties`` (the VELOCIraptor HDF5 output),
    using ``Vmax`` [km/s] and ``Mass_BN98`` [1e10 Msun] fields, with
    bin resolution 0.05 in log10(V).  Exactly mirrors ``MvsV.ipynb``.
    """
    paths = {w: _vr_cache_path(sim_label, w) for w in ("Lower", "Median", "Upper")}
    if all(p.is_file() for p in paths.values()) and not refresh:
        return paths

    cfg = load_simulations()
    sim = cfg[sim_label]
    if not sim.has_vr:
        raise ValueError(f"{sim_label} has no VR data")

    vr_props_path = sim.vr_raw_dir() / "output.properties"
    if not vr_props_path.is_file():
        raise FileNotFoundError(f"VR properties file not found: {vr_props_path}")

    with h5py.File(vr_props_path, "r") as f:
        halo_vmax = np.array(f["Vmax"])          # km/s
        halo_mass = np.array(f["Mass_BN98"])     # SWIFT units: 1e10 Msun
        h_raw     = f["SimulationInfo"].attrs["h_val"]
        h_factor  = float(h_raw.decode() if isinstance(h_raw, bytes) else h_raw)

    # Notebook: SWIFT_UnitMass = 1e10*Msun, Converted_UnitMass = Msun/h
    # => m_converted = m * 1e10 * h
    m_converted = halo_mass * 1e10 * h_factor  # in h^{-1} Msun

    # Exclude main halo (index 0 in VR = halo ID 1)
    # Include all halos with Vmax > 0
    mask = np.ones(len(halo_vmax), dtype=bool)
    mask[0] = False
    mask &= (halo_vmax > 0)

    vmax_sel = halo_vmax[mask]
    m_sel    = m_converted[mask]

    v_center, m_lower, m_median, m_upper = _bin_mvsv(vmax_sel, m_sel, _VR_V_RESO)

    paths["Lower"].parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(paths["Lower"],  np.stack((v_center, m_lower),  axis=-1))
    np.savetxt(paths["Median"], np.stack((v_center, m_median), axis=-1))
    np.savetxt(paths["Upper"],  np.stack((v_center, m_upper),  axis=-1))
    return paths


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _rockstar_line(v_log10: np.ndarray) -> np.ndarray:
    """Analytic M_vir-Vmax power law from BolshoiP+MDPL (Grand & White 2021, Robert's result).

    "Rockstar" refers to the halo finder used in BolshoiP+MDPL simulations.
    This is the black MvsV reference line. It is intentionally analytic and
    does not read any of the COCO M-Vmax tables under data/external/trash/.
    """
    return v_log10 * 3.2282 + 4.7355  # log10(M)


def _interp_data(data: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    mask = np.isfinite(data[:, 1]) & (data[:, 1] != 0)
    fn   = interp1d(data[mask, 0], data[mask, 1], kind="linear",
                    fill_value="extrapolate")
    return np.column_stack((x_new, fn(x_new)))


# ---------------------------------------------------------------------------
# Plot: Fig 14 — VR-SOAP data vs Rockstar reference
# ---------------------------------------------------------------------------

def plot(
    host: str,
    *,
    refresh: bool = False,
    output: str | Path | None = None,
    use_tex: bool = True,
) -> Path:
    """Render Fig 14 MvsV for *host* using SOAP-VR data ('m12i' or 'm12f').

    This reproduces the archive notebook
    ``comparison_with_VR-SOAP/MvsV/MvsV_{host}.ipynb`` which plots the
    VR-SOAP M_vir–V_max relation for the three m12i (or m12f) sims
    versus the BolshoiP+MDPL Rockstar power-law reference.
    """
    if host not in _HOST_SIMS:
        raise ValueError(f"host must be 'm12i' or 'm12f', got {host!r}")

    sim_labels = _HOST_SIMS[host]
    cfg = load_simulations()
    for sl in sim_labels:
        extract_vr(sl, refresh=refresh)

    if use_tex:
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
        rc("text", usetex=True)

    medians = [np.genfromtxt(str(_vr_cache_path(sl, "Median"))) for sl in sim_labels]
    lowers  = [np.genfromtxt(str(_vr_cache_path(sl, "Lower")))  for sl in sim_labels]
    uppers  = [np.genfromtxt(str(_vr_cache_path(sl, "Upper")))  for sl in sim_labels]

    v_rock = np.arange(0, 2 + _V_RESO, _V_RESO)
    x2 = _UNRESOLVED_X2_VR[host]

    fig = plt.figure(figsize=(6, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # Main panel: Median lines (cache stores m * h; divide by h for Msun)
    for i, sl in enumerate(sim_labels):
        d = medians[i]
        ax0.plot(d[:, 0], d[:, 1] / _H_FACTOR, _LINE_STYLES[i],
                 label=_LABEL_MAP[sl])

    # Shaded percentile bands
    for i in range(3):
        lo = lowers[i]
        up = uppers[i]
        ax0.fill_between(lo[:, 0], lo[:, 1] / _H_FACTOR, up[:, 1] / _H_FACTOR,
                         color=_LINE_STYLES[i], alpha=0.2)

    # Rockstar reference only (no coco fitting in this notebook)
    ax0.plot(10 ** v_rock, 10 ** _rockstar_line(v_rock) / _H_BOLSHOI, "black",
             label=r"$\rm{BolshoiP\&MDPL\ result}$")

    # Unresolved region
    ax0.axvspan(_X1, x2, color="grey", alpha=0.3,
                label=fr"$\rm{{unresolved\ region\ in\ {host}:PL}}$")
    ax1.axvspan(_X1, x2, color="grey", alpha=0.3)

    ax0.set_xlim(left=_X1)
    ax1.set_xlim(left=_X1)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax0.set_ylabel(r"$M_{\rm vir}\ [{\rm M}_{\odot}]$", fontsize=20)
    ax0.legend(fontsize=12)

    # Ratio panel: BT_deep/PL and BT_soft/PL (direct division, VR data)
    # Notebook divides arrays directly — interpolate onto cdmo x-grid for safety.
    base = medians[0]
    x_base = base[:, 0]  # all sims share same bin edges from same V range
    deep_interp = _interp_data(medians[1], x_base)
    soft_interp = _interp_data(medians[2], x_base)

    deep_ratio = deep_interp[:, 1] / base[:, 1]
    soft_ratio = soft_interp[:, 1] / base[:, 1]

    ax1.plot(x_base, deep_ratio, _LINE_STYLES[1])
    ax1.plot(x_base, soft_ratio, _LINE_STYLES[2])

    ax1.set_ylabel(r"$M_{\rm BT}/M_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    # ax1.legend(fontsize=8)  # commented out in archive notebook

    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"$V_{\rm max}\ \rm{[km/s]}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    if output is None:
        out_path = figures_dir() / f"MvsV_{host}.png"
    else:
        out_path = Path(output)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path

