"""Fig 9: Cumulative Mass Fraction (CMF) profile.

Source notebooks (Type 1 extract + Type 2 comparison):
  {sim}_fof/hbtplus_min_20/CMF_hbtplus_scaled.ipynb   (Type 1 – extract)
  paper_figures_min_20_HBTplus_v1_archive/CMF/
    comparison_CMF_HBTplus_scaled.ipynb                (Type 2 – plot)

Algorithm (Type 1):
  * Open SubSnap_{snap:03d}.0.hdf5 (HBTplus raw output).
  * Filter alive subhalos: SnapshotIndexOfDeath == -1.
  * Main halo = argmax(Mbound[alive]).
  * Main halo centre = ComovingMostBoundPosition[alive][main_idx]  [Mpc].
  * R200c from SOAP: Main_Halo_R200c = max(SO/200_crit/SORadius)  [Mpc].
  * Read snapshot DMParticles: ParticleIDs + Coordinates  [Mpc].
  * Distances of all DM particles from main halo centre, scaled by R200c.
  * Substructure particles = alive subhalo particles minus main halo.
  * Fixed bins: R = [0] + 10^arange(-2, 0.6, 0.1)  (26 edges).
  * Histogram both substructure and total DM particles with these bins.
  * result = cumsum(sub) / cumsum(total).
  * Output: (bins[1:], result)  [26 rows, r/R200c vs M_sub(<r)/M_tot(<r)].

Plot (Type 2):
  * Two-panel figure (height_ratios=[7, 3]), figsize=(6, 8).
  * Main panel: 6 sims + Lovell CDM reference.
    lines: ['C1','C2','C3','C1--','C2--','C3--','grey']
    Lovell CDM: x = col0 / 256.1,  y = col1.
  * Ratio panel: BT_deep / PL and BT_soft / PL, separately for m12i and m12f.
  * xlim=(1e-1, 2), ylim_main=(1e-3, 1), ylim_ratio=(2e-1, 5e1), log-log.
  * Output: figures/CMF_HBTplus_scaled.png.
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
from ..paths import cache_dir, data_dir, figures_dir

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Fixed radial bins in r/R200c units (matches reference notebook exactly).
# np.arange(-2, 0.5+0.1, 0.1) == np.arange(-2, 0.6, 0.1) → 26 values.
_BINS_LOG10 = np.arange(-2, 0.6, 0.1)
_RADIAL_BINS = np.concatenate(([0.0], 10.0 ** _BINS_LOG10))  # 27 edges

# Lovell et al. (2014) CDM reference: x = r [kpc], y = M_sub(<r)/M_tot(<r).
# Plot as x / 256.1 to put on same r/R200c axis (R200c_Lovell ≈ 256.1 kpc).
_LOVELL_R200C_KPC: float = 256.1
_LOVELL_FILE = data_dir() / "external" / "lovell2014_cmf_r200c.txt"

# 6 simulation labels in plot order (matches reference filelist).
_SIM_LABELS = [
    "m12i_cdmo",
    "m12i_btps_deep",
    "m12i_btps_soft",
    "m12f_cdmo",
    "m12f_btps_deep",
    "m12f_btps_soft",
]

_LINE_STYLES = ["C1", "C2", "C3", "C1--", "C2--", "C3--"]
_LABEL_MAP = {
    "m12i_cdmo":      r"$\rm{m12i:PL}$",
    "m12i_btps_deep": r"$\rm{m12i:BT\_deep}$",
    "m12i_btps_soft": r"$\rm{m12i:BT\_soft}$",
    "m12f_cdmo":      r"$\rm{m12f:PL}$",
    "m12f_btps_deep": r"$\rm{m12f:BT\_deep}$",
    "m12f_btps_soft": r"$\rm{m12f:BT\_soft}$",
}


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _cache_path(sim_label: str) -> Path:
    return cache_dir() / "cmf" / sim_label


# ---------------------------------------------------------------------------
# Extract: Type 1  (slow — reads snapshot + SubSnap HDF5)
# ---------------------------------------------------------------------------

def extract(sim_label: str, *, refresh: bool = False) -> Path:
    """Compute CMF for *sim_label* and write to cache.

    Reads:
      - ``{sim_fof}/hbtplus_min_20/{snap:03d}/SubSnap_{snap:03d}.0.hdf5``
      - ``{sim}/SOAP {snap}`` (R200c only)
      - ``{sim}/{label}_{snap:04d}.hdf5``  (DM particle snapshot)

    Returns the path to the cache file.
    """
    out_path = _cache_path(sim_label)
    if out_path.is_file() and not refresh:
        return out_path

    cfg = load_simulations()
    sim = cfg[sim_label]
    snap = sim.snapshot_num

    # --- SubSnap (HBTplus raw) ---
    subsnap_path = sim.hbt_subsnap_dir() / f"SubSnap_{snap:03d}.0.hdf5"
    with h5py.File(subsnap_path, "r") as f:
        subs = np.array(f["Subhalos"])               # structured array
        subhalo_particles = np.array(f["SubhaloParticles"])  # object array

    alive = subs["SnapshotIndexOfDeath"] == -1
    mbound_alive    = subs["Mbound"][alive]
    positions_alive = subs["ComovingMostBoundPosition"][alive]  # [Mpc]
    particles_alive = subhalo_particles[alive]

    main_idx  = int(np.argmax(mbound_alive))
    centre    = positions_alive[main_idx]      # (3,) [Mpc]
    main_pids = particles_alive[main_idx]

    # --- R200c from SOAP ---
    with h5py.File(sim.soap_hbt_path(), "r") as f:
        r200c_all = np.array(f["/SO/200_crit/SORadius"])
    r200c_main = float(r200c_all.max())  # Mpc (matches notebook: R200c.max())

    # --- DM snapshot ---
    with h5py.File(sim.snapshot_path(), "r") as f:
        dm = f["DMParticles"]
        hrs_pids   = np.array(dm["ParticleIDs"])
        hrs_coords = np.array(dm["Coordinates"])  # [Mpc]

    # Distances of all DM particles from main halo centre [Mpc]
    diffs = hrs_coords - centre
    halo_dist = np.linalg.norm(diffs, axis=1)

    # Substructure particles = all alive subhalo particles minus main halo
    valid_pids   = np.concatenate(particles_alive)
    mask_sub     = np.isin(hrs_pids, valid_pids)
    filt_dist    = halo_dist[mask_sub]
    filt_pids    = hrs_pids[mask_sub]
    not_main     = ~np.isin(filt_pids, main_pids)
    final_dist   = filt_dist[not_main]              # substructure particle dists

    # Scale by R200c
    final_dist_scaled = final_dist / r200c_main
    all_dist_scaled   = halo_dist  / r200c_main

    # Histogram with fixed bins
    count_sub,   _ = np.histogram(final_dist_scaled, bins=_RADIAL_BINS)
    count_total, _ = np.histogram(all_dist_scaled,   bins=_RADIAL_BINS)

    cum_sub   = np.cumsum(count_sub)
    cum_total = np.cumsum(count_total)

    result = np.divide(cum_sub, cum_total, out=np.zeros_like(cum_sub, dtype=float),
                       where=cum_total != 0)

    out = np.stack((_RADIAL_BINS[1:], result), axis=-1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_path, out)
    return out_path


# ---------------------------------------------------------------------------
# Ratio helper (matches reference notebook's `calculate_ratio`)
# ---------------------------------------------------------------------------

def _calculate_ratio(
    base_data: np.ndarray, data: np.ndarray
) -> np.ndarray:
    """Interpolate *data* onto non-zero rows of *base_data*, return ratio."""
    nz_base = np.where(base_data[:, 1] != 0)[0]
    nz_data = np.where(data[:, 1] != 0)[0]
    common  = np.intersect1d(nz_base, nz_data)

    x = base_data[common, 0]

    def _interp(d: np.ndarray) -> np.ndarray:
        mask = d[:, 1] != 0
        fn   = interp1d(d[mask, 0], d[mask, 1], kind="linear",
                        fill_value="extrapolate")
        return fn(x)

    return np.column_stack((x, _interp(data) / _interp(base_data)))


# ---------------------------------------------------------------------------
# Plot: Type 2
# ---------------------------------------------------------------------------

def plot(
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Render Fig 9 (CMF_HBTplus_scaled.png)."""

    # --- font / TeX ---
    rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
    rc("text", usetex=use_tex)

    # --- load sim caches (extract if needed) ---
    datasets: dict[str, np.ndarray] = {}
    for label in _SIM_LABELS:
        p = _cache_path(label)
        if not p.is_file() or refresh:
            extract(label, refresh=refresh)
        datasets[label] = np.genfromtxt(p)

    # --- Lovell CDM reference ---
    lovell = np.genfromtxt(_LOVELL_FILE)

    # --- figure ---
    fig = plt.figure(figsize=(6, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # Main panel
    for label, ls in zip(_SIM_LABELS, _LINE_STYLES):
        d = datasets[label]
        ax0.plot(d[:, 0], d[:, 1], ls, label=_LABEL_MAP[label])

    ax0.plot(lovell[:, 0] / _LOVELL_R200C_KPC, lovell[:, 1], "grey",
             label=r"$\rm{Lovell\ CDM}$")

    ax0.set_xlim(left=1e-1, right=2)
    ax0.set_ylim(1e-3, 1)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_ylabel(
        r"$M_{\rm substructures}\left( <r \right)/M_{\rm total}\left( <r \right)$",
        fontsize=20,
    )
    ax0.legend(fontsize=12)

    # Ratio panel
    d_m12i_cdmo      = datasets["m12i_cdmo"]
    d_m12i_btps_deep = datasets["m12i_btps_deep"]
    d_m12i_btps_soft = datasets["m12i_btps_soft"]
    d_m12f_cdmo      = datasets["m12f_cdmo"]
    d_m12f_btps_deep = datasets["m12f_btps_deep"]
    d_m12f_btps_soft = datasets["m12f_btps_soft"]

    m12i_deep_ratio = _calculate_ratio(d_m12i_cdmo, d_m12i_btps_deep)
    m12i_soft_ratio = _calculate_ratio(d_m12i_cdmo, d_m12i_btps_soft)
    m12f_deep_ratio = _calculate_ratio(d_m12f_cdmo, d_m12f_btps_deep)
    m12f_soft_ratio = _calculate_ratio(d_m12f_cdmo, d_m12f_btps_soft)

    ax1.plot(m12i_deep_ratio[:, 0], m12i_deep_ratio[:, 1], "C2",   label="m12i_deep")
    ax1.plot(m12i_soft_ratio[:, 0], m12i_soft_ratio[:, 1], "C3",   label="m12i_soft")
    ax1.plot(m12f_deep_ratio[:, 0], m12f_deep_ratio[:, 1], "C2--", label="m12f_deep")
    ax1.plot(m12f_soft_ratio[:, 0], m12f_soft_ratio[:, 1], "C3--", label="m12f_soft")

    ax1.set_xlim(left=1e-1, right=2)
    ax1.set_ylim(2e-1, 5e1)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$N_{\rm BT}/N_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="C1", linestyle="-", alpha=0.3)
    # ax1.legend(fontsize=6)   # commented out in reference notebook

    # Tick params
    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax0.xaxis.set_tick_params(labelbottom=False)

    ax1.set_xlabel(r"$r/R_{200c}$", fontsize=20)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    # --- save ---
    if output is None:
        out_path = figures_dir() / "CMF_HBTplus_scaled.png"
    else:
        out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path
