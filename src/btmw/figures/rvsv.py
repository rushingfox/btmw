"""Fig 10: R_max vs V_max relation (RvsV).

Source notebooks:
  paper_figures_min_20_HBTplus_v1_archive/RvsV/
    RvsV_hbtplus.ipynb                 (Type 1 – extract from SOAP)
    comparison_RvsV_m12i.ipynb         (Type 2 – plot m12i)
    comparison_RvsV_m12f.ipynb         (Type 2 – plot m12f)

Algorithm (Type 1):
  * Open SOAP halo_properties (HBT-backed).
  * Main halo = argmax(BoundSubhalo/TotalMass).
  * Select subhalos within R200c of main halo (by HaloCentre distance).
  * Exclude main halo itself.
  * Rmax converted: BoundSubhalo/MaximumCircularVelocityRadiusUnsoftened [Mpc] × 1000 → [kpc].
  * Vmax: BoundSubhalo/MaximumCircularVelocityUnsoftened [km/s].
  * Fixed log-V bins: 10^arange(floor(log10(Vmin)), ceil(log10(Vmax))+0.1, 0.1).
  * For each bin: 16th, 50th, 84th percentile of Rmax.
  * Output: three files  {label}_RvsV_{Lower,Median,Upper}
    columns: (V_center [km/s], Rmax_percentile [kpc])  — NaN for empty bins.

Plot (Type 2):
  * Two panels side-by-side? No — two separate figures: m12i and m12f.
  * Each figure: two-panel (height_ratios=[7, 3]), figsize=(6, 8).
  * Main panel (ax0):
    - 3 sim Median lines + fill_between(Lower, Upper).
    - Rockstar reference (black line).
    - Grey axvspan for unresolved region.
    - Filter: data = data[(data[:,1]!=0) & ~np.isnan(data[:,1])] before plotting.
  * Ratio panel (ax1): BT_deep/PL and BT_soft/PL (Median only, interpolated).
  * xlim=(1, ∞ auto), ylim_main=(0.1, 10), ylim_ratio=(0.2, 1.5), log-log.
  * Unresolved region x2:
    m12i: 0.02995634129285346 * 14231044.216963861 / 1e5  ≈ 4.263 km/s
    m12f: 0.0247377586709585  * 15862762.136794837 / 1e5  ≈ 3.924 km/s
  * Output: figures/RvsV_m12i_HBTplus.png and figures/RvsV_m12f_HBTplus.png.
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
from ..io_soap import (
    SOAP_BOUND_TOTAL_MASS,
    SOAP_BOUND_VMAX,
    SOAP_SO200C_RADIUS,
    SOAP_HALO_CENTRE,
)
from ..paths import cache_dir, data_dir, figures_dir

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Rmax field in SOAP (not yet in io_soap constants)
_SOAP_RMAX = "/BoundSubhalo/MaximumCircularVelocityRadiusUnsoftened"

# Percentile levels (match reference notebook)
_LOWER_PCT  = 16
_MEDIAN_PCT = 50
_UPPER_PCT  = 84

# Velocity bin width in log10 (same as HVF)
_V_RESO = 0.1

# Unresolved region x2 (km/s) per host — derived from softening * particle mass
_UNRESOLVED_X2 = {
    "m12i": 0.02995634129285346 * 14231044.216963861 / 1e5,   # ≈ 4.263
    "m12f": 0.0247377586709585  * 15862762.136794837 / 1e5,   # ≈ 3.924
}

# Reference: Grand & White (2021), referred to as "Robert result" in figure legends.
# Variable named _ROCKSTAR_FILE because Rockstar is the halo finder used in that paper.
_ROCKSTAR_FILE = data_dir() / "external" / "grand2021_rmax_r200c.txt"

# Sim groups per host (plot order matches reference)
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

def _cache_dir(sim_label: str) -> Path:
    return cache_dir() / "rvsv"


def _cache_path(sim_label: str, which: str) -> Path:
    """Return path for one of Lower / Median / Upper."""
    return cache_dir() / "rvsv" / f"{sim_label}_{which}"


# ---------------------------------------------------------------------------
# Extract: Type 1  (slow — reads SOAP HDF5)
# ---------------------------------------------------------------------------

def extract(sim_label: str, *, refresh: bool = False) -> dict[str, Path]:
    """Compute RvsV percentile bands for *sim_label* and write to cache.

    Returns dict with keys 'Lower', 'Median', 'Upper' mapping to cache Paths.
    """
    paths = {w: _cache_path(sim_label, w) for w in ("Lower", "Median", "Upper")}
    if all(p.is_file() for p in paths.values()) and not refresh:
        return paths

    cfg = load_simulations()
    sim = cfg[sim_label]

    with h5py.File(sim.soap_hbt_path(), "r") as f:
        mbound      = np.array(f[SOAP_BOUND_TOTAL_MASS])
        vmax        = np.array(f[SOAP_BOUND_VMAX])
        rmax_mpc    = np.array(f[_SOAP_RMAX])
        r200c       = np.array(f[SOAP_SO200C_RADIUS])
        halo_centre = np.array(f[SOAP_HALO_CENTRE])

    main_idx   = int(np.argmax(mbound))
    r200c_main = float(r200c[main_idx])
    centre     = halo_centre[main_idx]

    # Select subhalos within R200c (exclude main halo)
    diffs = halo_centre - centre
    dists = np.linalg.norm(diffs, axis=1)
    mask  = (dists < r200c_main)
    mask[main_idx] = False

    vmax_sel = vmax[mask]
    rmax_kpc  = rmax_mpc[mask] * 1000.0   # Mpc → kpc

    # Build log-V bins
    lv_min = math.floor(math.log10(vmax_sel.min()))
    lv_max = math.ceil(math.log10(vmax_sel.max()))
    bin_edges = 10.0 ** np.arange(lv_min, lv_max + _V_RESO, _V_RESO)
    v_center  = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    n_bins    = len(v_center)

    rmax_lower  = np.full(n_bins, np.nan)
    rmax_median = np.full(n_bins, np.nan)
    rmax_upper  = np.full(n_bins, np.nan)

    for k in range(n_bins):
        in_bin = (vmax_sel >= bin_edges[k]) & (vmax_sel < bin_edges[k + 1])
        sub = rmax_kpc[in_bin]
        if sub.size > 0:
            rmax_lower[k]  = np.percentile(sub, _LOWER_PCT)
            rmax_median[k] = np.percentile(sub, _MEDIAN_PCT)
            rmax_upper[k]  = np.percentile(sub, _UPPER_PCT)

    paths["Lower"].parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(paths["Lower"],  np.stack((v_center, rmax_lower),  axis=-1))
    np.savetxt(paths["Median"], np.stack((v_center, rmax_median), axis=-1))
    np.savetxt(paths["Upper"],  np.stack((v_center, rmax_upper),  axis=-1))
    return paths


# ---------------------------------------------------------------------------
# Plot helper
# ---------------------------------------------------------------------------

def _interp_data(data: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    mask = data[:, 1] != 0
    fn   = interp1d(data[mask, 0], data[mask, 1], kind="linear",
                    fill_value="extrapolate")
    return np.column_stack((x_new, fn(x_new)))


def _plot_host(host: str, *, refresh: bool, use_tex: bool) -> Path:
    rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
    rc("text", usetex=use_tex)

    sim_labels = _HOST_SIMS[host]
    line_styles = _LINE_STYLES

    # Load caches
    medians: list[np.ndarray] = []
    lowers:  list[np.ndarray] = []
    uppers:  list[np.ndarray] = []
    for label in sim_labels:
        paths = extract(label, refresh=refresh)
        medians.append(np.genfromtxt(paths["Median"]))
        lowers.append(np.genfromtxt(paths["Lower"]))
        uppers.append(np.genfromtxt(paths["Upper"]))

    rockstar = np.genfromtxt(_ROCKSTAR_FILE)

    fig = plt.figure(figsize=(6, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    x1 = 1e0
    x2 = _UNRESOLVED_X2[host]

    # Main panel: Median lines
    for i, label in enumerate(sim_labels):
        d = medians[i]
        d_plot = d[(d[:, 1] != 0) & ~np.isnan(d[:, 1])]
        ax0.plot(d_plot[:, 0], d_plot[:, 1], line_styles[i],
                 label=_LABEL_MAP[label])

    # Rockstar reference
    ax0.plot(rockstar[:, 0], rockstar[:, 1], color="black",
             label=r"$\rm{Robert\ result\ for\ }R_{200c}$")

    # Unresolved region
    ax0.axvspan(x1, x2, color="grey", alpha=0.3,
                label=fr"$\rm{{unresolved\ region\ of\ {host}:PL}}$")
    ax1.axvspan(x1, x2, color="grey", alpha=0.3)

    # Shaded percentile bands
    for i in range(3):
        lo = lowers[i]
        up = uppers[i]
        ax0.fill_between(lo[:, 0], lo[:, 1], up[:, 1],
                         color=line_styles[i], alpha=0.2)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)
    ax0.set_ylim(0.1, 10)
    ax1.set_ylim(0.2, 1.5)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax0.set_ylabel(r"${R_{\rm max}\rm{[kpc]}}$", fontsize=20)
    ax0.legend(fontsize=13)

    # Ratio panel
    base = medians[0]
    base_filt = base[(~np.isnan(base[:, 1])) & (base[:, 1] != 0)]
    x_base = base_filt[:, 0]

    deep_interp = _interp_data(medians[1], x_base)
    soft_interp = _interp_data(medians[2], x_base)

    ratio_deep = np.column_stack((x_base, deep_interp[:, 1] / base_filt[:, 1]))
    ratio_soft = np.column_stack((x_base, soft_interp[:, 1] / base_filt[:, 1]))

    ax1.plot(ratio_deep[:, 0], ratio_deep[:, 1], line_styles[1])
    ax1.plot(ratio_soft[:, 0], ratio_soft[:, 1], line_styles[2])

    ax1.set_ylabel(r"$R_{\rm BT}/R_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="C1", linestyle="-", alpha=0.3)
    # ax1.legend(fontsize=8)  # commented out in reference

    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"${V_{\rm max}\rm{[km/s]}}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out_path = figures_dir() / f"RvsV_{host}_HBTplus.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Public plot entry point
# ---------------------------------------------------------------------------

def plot(
    host: str,
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Render RvsV figure for *host* ('m12i' or 'm12f')."""
    if host not in _HOST_SIMS:
        raise ValueError(f"host must be 'm12i' or 'm12f', got {host!r}")
    out = _plot_host(host, refresh=refresh, use_tex=use_tex)
    if output is not None:
        import shutil
        shutil.copy(out, output)
        return Path(output)
    return out


# ---------------------------------------------------------------------------
# Fig 17 — RvsV resolution study (high-res m12i only)
# ---------------------------------------------------------------------------

_RES_HIGH_LABELS = ("m12i_cdmo_highest", "m12i_btps_deep_highest", "m12i_btps_soft_highest")
_RES_HIGH_TEX = (
    r"$\rm{m12i:PL-high}$",
    r"$\rm{m12i:BT\_deep-high}$",
    r"$\rm{m12i:BT\_soft-high}$",
)
# Unresolved x2 for high-res: softening-scale Vmax * R200c_max
_RES_UNRES_HIGH_X2 = 0.01551527244439509 * 14231044.216963861 / 1e5   # ≈ 2.208 km/s


def plot_resolution_study(
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 17 — RvsV high-resolution m12i comparison."""
    rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
    rc("text", usetex=use_tex)

    medians: list[np.ndarray] = []
    lowers:  list[np.ndarray] = []
    uppers:  list[np.ndarray] = []
    for label in _RES_HIGH_LABELS:
        paths = extract(label, refresh=refresh)
        medians.append(np.genfromtxt(paths["Median"]))
        lowers.append(np.genfromtxt(paths["Lower"]))
        uppers.append(np.genfromtxt(paths["Upper"]))

    rockstar = np.genfromtxt(_ROCKSTAR_FILE)

    fig = plt.figure(figsize=(6, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    x1 = 1e0
    x2 = _RES_UNRES_HIGH_X2

    # Main panel: Median lines (high-res only)
    for i, tex in enumerate(_RES_HIGH_TEX):
        d = medians[i]
        d_plot = d[(d[:, 1] != 0) & ~np.isnan(d[:, 1])]
        ax0.plot(d_plot[:, 0], d_plot[:, 1], _LINE_STYLES[i], label=tex)

    # Rockstar reference
    ax0.plot(rockstar[:, 0], rockstar[:, 1], color="black",
             label=r"$\rm{Robert\ result\ for\ }R_{200c}$")

    # Unresolved region
    ax0.axvspan(x1, x2, color="grey", alpha=0.3,
                label=r"$\rm{unresolved\ region-high\ res}$")
    ax1.axvspan(x1, x2, color="grey", alpha=0.3)

    # Shaded percentile bands
    for i in range(3):
        lo = lowers[i]
        up = uppers[i]
        ax0.fill_between(lo[:, 0], lo[:, 1], up[:, 1],
                         color=_LINE_STYLES[i], alpha=0.2)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)
    ax0.set_ylim(0.1, 10)
    ax1.set_ylim(0.2, 1.5)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax0.set_ylabel(r"${R_{\rm max}\rm{[kpc]}}$", fontsize=20)
    ax0.legend(fontsize=13)

    # Ratio panel: BT_deep/PL and BT_soft/PL (all high-res)
    base = medians[0]
    base_filt = base[(~np.isnan(base[:, 1])) & (base[:, 1] != 0)]
    x_base = base_filt[:, 0]

    deep_interp = _interp_data(medians[1], x_base)
    soft_interp = _interp_data(medians[2], x_base)

    ratio_deep = np.column_stack((x_base, deep_interp[:, 1] / base_filt[:, 1]))
    ratio_soft = np.column_stack((x_base, soft_interp[:, 1] / base_filt[:, 1]))

    ax1.plot(ratio_deep[:, 0], ratio_deep[:, 1], _LINE_STYLES[1])
    ax1.plot(ratio_soft[:, 0], ratio_soft[:, 1], _LINE_STYLES[2])

    ax1.set_ylabel(r"$R_{\rm BT}/R_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="C1", linestyle="-", alpha=0.3)

    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"${V_{\rm max}\rm{[km/s]}}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "RvsV_resolution_study.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out
