"""Fig 5 — Main-halo radial dark-matter density profile.

Two extraction backends:
  * ``'hbt'`` (paper v1 / m12i): reads HBT SubSnap + SWIFT snapshot directly.
    Cache files use the v1 naming: ``{sim}_radial_density_simulation``,
    ``{sim}_NFW_fitting``, ``{sim}_converg_radius``.
    * ``'vr'``  (PRD version / m12f): reads VR radial-profile HDF5 files
    (already binned by VelocIRaptor).
    Cache files use 250507 naming: ``{sim}_simulation_1``,
    ``{sim}_NFW_1``, ``{sim}_converg_radius1``.

Cache lives under ``data/cache/radial_density/{backend}/``.

The plot function reproduces the two-panel layout of the archive notebooks:
  * Main panel (log-log): rho vs r/R200m, with NFW fits and convergence shade.
  * Ratio panel (log-log): rho_BT / rho_PL.

Source notebooks (for numeric fidelity):
  v1 (HBT):
    working/.../paper_figures_min_20_HBTplus_v1_archive/RadialDensityProfile/
        comparison_RadialProfile_m12i.ipynb
    working/.../m12i_cdmo_fof/hbtplus_min_20/Radial_density_profile_hbtplus.ipynb
  250507 (VR):
    working/.../paper_figures_min_20_HBTplus_version_250507/RadialDensityProfile_VR/
        haloradialprofile_colossus0920.ipynb
        comparison_RadialProfile_m12i.ipynb
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Literal

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, rc
from scipy.interpolate import interp1d

from ..config import load_simulations
from ..defaults import (
    G_IN_CGS,
    KPC_IN_CGS,
    LITTLE_H,
    MPC_IN_CGS,
    SOLAR_MASS_IN_CGS,
)
from ..paths import cache_dir, figures_dir

Backend = Literal["hbt", "vr"]

# Default backend per host, matching the paper:
#   m12i: HBT particle-based (v1 archive; original centre)
#   m12f: VelocIRaptor profile-based (PRD version; better visual centre)
_DEFAULT_BACKEND: dict[str, Backend] = {"m12i": "hbt", "m12f": "vr"}

# Fiducial labels per host, in display order (PL, BT_deep, BT_soft).
_HOST_LABELS: dict[str, tuple[str, str, str]] = {
    "m12i": ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft"),
    "m12f": ("m12f_cdmo", "m12f_btps_deep", "m12f_btps_soft"),
}

_LINE_STYLES = ("C1", "C2", "C3", "C1--", "C2--", "C3--")

_LABELS_TEX = {
    "m12i": (
        r"$\rm{m12i:PL-simulation}$",
        r"$\rm{m12i:BT\_deep-simulation}$",
        r"$\rm{m12i:BT\_soft-simulation}$",
        r"$\rm{m12i:PL-NFW\ fitting}$",
        r"$\rm{m12i:BT\_deep-NFW\ fitting}$",
        r"$\rm{m12i:BT\_soft-NFW\ fitting}$",
    ),
    "m12f": (
        r"$\rm{m12f:PL-simulation}$",
        r"$\rm{m12f:BT\_deep-simulation}$",
        r"$\rm{m12f:BT\_soft-simulation}$",
        r"$\rm{m12f:PL-NFW\ fitting}$",
        r"$\rm{m12f:BT\_deep-NFW\ fitting}$",
        r"$\rm{m12f:BT\_soft-NFW\ fitting}$",
    ),
}

_LABELS_PLAIN = {
    "m12i": (
        "m12i:PL-simulation",
        "m12i:BT_deep-simulation",
        "m12i:BT_soft-simulation",
        "m12i:PL-NFW fitting",
        "m12i:BT_deep-NFW fitting",
        "m12i:BT_soft-NFW fitting",
    ),
    "m12f": (
        "m12f:PL-simulation",
        "m12f:BT_deep-simulation",
        "m12f:BT_soft-simulation",
        "m12f:PL-NFW fitting",
        "m12f:BT_deep-NFW fitting",
        "m12f:BT_soft-NFW fitting",
    ),
}


# ---------------------------------------------------------------------------
# Cache paths
# ---------------------------------------------------------------------------

def _hbt_cache_dir() -> Path:
    return cache_dir() / "radial_density" / "hbt"


def _vr_cache_dir() -> Path:
    return cache_dir() / "radial_density" / "vr"


def _hbt_paths(sim_label: str) -> tuple[Path, Path, Path]:
    """(simulation, NFW_fitting, converg_radius) for the HBT backend."""
    d = _hbt_cache_dir()
    return (
        d / f"{sim_label}_radial_density_simulation",
        d / f"{sim_label}_NFW_fitting",
        d / f"{sim_label}_converg_radius",
    )


def _vr_paths(sim_label: str) -> tuple[Path, Path, Path]:
    """(simulation_1, NFW_1, converg_radius1) for the VR backend."""
    d = _vr_cache_dir()
    return (
        d / f"{sim_label}_simulation_1",
        d / f"{sim_label}_NFW_1",
        d / f"{sim_label}_converg_radius1",
    )


# ---------------------------------------------------------------------------
# rho_crit helper (shared by both backends)
# ---------------------------------------------------------------------------

def _rho_crit_cgs(little_h: float = LITTLE_H) -> float:
    """Critical density in g cm^-3 (WMAP7 cosmology at z=0)."""
    return (
        3.0 * (little_h * 100.0) ** 2
        / (8.0 * math.pi * G_IN_CGS)
        * (1e5 / MPC_IN_CGS) ** 2
    )


# ---------------------------------------------------------------------------
# HBT-based extraction (v1 paper, m12i)
# ---------------------------------------------------------------------------

def extract_hbt(sim_label: str, *, refresh: bool = False) -> None:
    """Compute radial density from HBT particle lists + SWIFT snapshot.

    Reproduces the algorithm in
    ``m12i_cdmo_fof/hbtplus_min_20/Radial_density_profile_hbtplus.ipynb``.

    Requires the SWIFT snapshot (tens of GB) in memory; run on a compute node.
    """
    sim_path, nfw_path, conv_path = _hbt_paths(sim_label)
    if (
        not refresh
        and sim_path.exists()
        and nfw_path.exists()
        and conv_path.exists()
    ):
        return

    sim_path.parent.mkdir(parents=True, exist_ok=True)

    cfg = load_simulations()
    sim = cfg[sim_label]
    snap = sim.snapshot_num

    # ------------------------------------------------------------------
    # 1. SOAP → M200_mean, C200_mean, R200_mean of main halo
    # ------------------------------------------------------------------
    soap_path = sim.soap_hbt_path()
    with h5py.File(soap_path, "r") as f:
        m200m = np.array(f["/SO/200_mean/TotalMass"])        # 1e10 Msun
        c200m = np.array(f["/SO/200_mean/Concentration"])    # dimensionless
        r200m = np.array(f["/SO/200_mean/SORadius"])          # Mpc

    soap_main_idx = int(np.argmax(m200m))
    M200m_Msun = float(m200m[soap_main_idx]) * 1e10
    C200m = float(c200m[soap_main_idx])
    R200m_Mpc = float(r200m[soap_main_idx])
    R200m_kpc = R200m_Mpc * 1000.0

    # ------------------------------------------------------------------
    # 2. HBT SubSnap → centre + particle IDs of main halo
    # ------------------------------------------------------------------
    subsnap_dir = sim.hbt_subsnap_dir()
    subsnap_glob = sim.hbt_subsnap_glob()
    subsnap_files = sorted(subsnap_dir.glob(subsnap_glob))
    if not subsnap_files:
        raise FileNotFoundError(
            f"No SubSnap files matching {subsnap_glob!r} in {subsnap_dir}"
        )

    all_subs_list = []
    all_sp_list = []
    for sf in subsnap_files:
        with h5py.File(sf, "r") as f:
            all_subs_list.append(np.array(f["Subhalos"]))
            all_sp_list.append(np.array(f["SubhaloParticles"]))

    subs = np.concatenate(all_subs_list)
    sp = np.concatenate(all_sp_list)

    alive_mask = subs["SnapshotIndexOfDeath"] == -1
    alive_subs = subs[alive_mask]
    alive_sp = sp[alive_mask]

    main_alive_idx = int(np.argmax(alive_subs["Mbound"]))
    centre = alive_subs["ComovingMostBoundPosition"][main_alive_idx].astype(float)
    particle_ids = alive_sp[main_alive_idx]  # 1-D array of DM particle IDs

    # ------------------------------------------------------------------
    # 3. Snapshot → all DM particle IDs + coordinates
    # ------------------------------------------------------------------
    snap_path = sim.snapshot_path()
    with h5py.File(snap_path, "r") as f:
        all_ids = np.array(f["DMParticles/ParticleIDs"])
        all_coords = np.array(f["DMParticles/Coordinates"])  # comoving Mpc

    # ------------------------------------------------------------------
    # 4. Match + distances in kpc
    # ------------------------------------------------------------------
    mask = np.isin(all_ids, particle_ids)
    coords = all_coords[mask]
    distances_kpc = np.linalg.norm(coords - centre, axis=1) * 1000.0  # Mpc→kpc

    # ------------------------------------------------------------------
    # 5. Histogram (source notebook Cell 7)
    # ------------------------------------------------------------------
    bins = np.append(np.array([0.0]), 10.0 ** np.arange(-2.0, 4.1, 0.1))
    countdata, _ = np.histogram(distances_kpc, bins=bins)

    lhs_bins = bins[:-1]
    rhs_bins = bins[1:]
    centre_bins = (lhs_bins + rhs_bins) / 2.0

    volume_bins = np.array([
        4.0 / 3.0 * math.pi * (bins[i + 1] ** 3 - bins[i] ** 3)
        for i in range(len(bins) - 1)
    ])

    _pm_swift = cfg.particle_mass_in_swift_unit.get(sim.resolution)
    if _pm_swift is None:
        raise ValueError(
            f"particle_mass_in_swift_unit[{sim.resolution!r}] is not set in "
            "simulations.yaml; please fill in the value before running extract_hbt."
        )
    particle_mass_msun = _pm_swift * 1e10
    density = countdata * particle_mass_msun / volume_bins  # Msun / kpc^3

    # ------------------------------------------------------------------
    # 6. Save simulation profile (filter r < R200, source Cell 8)
    # ------------------------------------------------------------------
    out_sim = np.stack((centre_bins / R200m_kpc, density), axis=-1)
    out_sim = out_sim[out_sim[:, 0] < 1.0]
    np.savetxt(str(sim_path), out_sim)

    # ------------------------------------------------------------------
    # 7. NFW profile (source Cell 9)
    # ------------------------------------------------------------------
    rho_halo = M200m_Msun / (4.0 / 3.0 * math.pi * R200m_kpc ** 3)
    A_NFW = math.log(1.0 + C200m) - C200m / (1.0 + C200m)

    nfw_r = centre_bins[centre_bins <= R200m_kpc]
    x_nfw = nfw_r / R200m_kpc
    nfw_rho = rho_halo / (3.0 * A_NFW * x_nfw * (1.0 / C200m + x_nfw) ** 2)

    out_nfw = np.stack((nfw_r / R200m_kpc, nfw_rho), axis=-1)
    np.savetxt(str(nfw_path), out_nfw)

    # ------------------------------------------------------------------
    # 8. Convergence radius (Power+03, source Cell 10)
    # ------------------------------------------------------------------
    npart_cumsum = np.cumsum(countdata)
    mass_cumsum = npart_cumsum * particle_mass_msun  # Msun

    rho_crit = _rho_crit_cgs()
    converg_idx = 0
    for i in range(len(rhs_bins)):
        if npart_cumsum[i] <= 1:
            continue
        mean_density_cgs = (
            mass_cumsum[i] * SOLAR_MASS_IN_CGS
            / (4.0 / 3.0 * math.pi * rhs_bins[i] ** 3 * KPC_IN_CGS ** 3)
        )
        tre = (
            math.sqrt(200.0)
            / 8.0
            * npart_cumsum[i]
            / math.log(npart_cumsum[i])
            * (mean_density_cgs / rho_crit) ** (-0.5)
        )
        if tre > 0.6:
            converg_idx = i
            break
    else:
        import warnings
        warnings.warn(
            f"{sim_label}: Power+03 convergence radius not found; "
            "falling back to innermost bin (converg_idx=0).",
            RuntimeWarning,
            stacklevel=2,
        )

    np.savetxt(str(conv_path), [rhs_bins[converg_idx] / R200m_kpc])


# ---------------------------------------------------------------------------
# VR-based extraction (PRD version, m12f)
# ---------------------------------------------------------------------------

def extract_vr(sim_label: str, *, refresh: bool = False) -> None:
    """Compute radial density from VelocIRaptor profile + properties files.

    Reproduces the algorithm in
    ``paper_figures_min_20_HBTplus_version_250507/RadialDensityProfile_VR/
        haloradialprofile_colossus0920.ipynb``.

    VR profile files are small (already binned), so this can run on the
    login node.
    """
    sim_path, nfw_path, conv_path = _vr_paths(sim_label)
    if (
        not refresh
        and sim_path.exists()
        and nfw_path.exists()
        and conv_path.exists()
    ):
        return

    sim_path.parent.mkdir(parents=True, exist_ok=True)

    cfg = load_simulations()
    sim = cfg[sim_label]

    vr_dir = sim.vr_raw_dir()
    vr_base = sim.vr_raw_basename()   # "output"

    fname_profiles = str(vr_dir / vr_base) + ".profiles"
    fname_properties = str(vr_dir / vr_base) + ".properties"

    # ------------------------------------------------------------------
    # 1. VR profiles (source Cell 5-6)
    # ------------------------------------------------------------------
    with h5py.File(fname_profiles, "r") as f:
        mass_profile = np.array(f["Mass_profile"])      # shape (N_halos, N_bins), 1e10 Msun
        npart_profile = np.array(f["Npart_profile"])    # shape (N_halos, N_bins)
        bin_edges = np.array(f["Radial_bin_edges"])     # shape (N_bins+1,), Mpc

    halo_index = 0  # main halo is first in VR output

    bin_rhs = bin_edges[1:]
    bin_lhs = bin_edges[:-1]
    bin_centre = (bin_rhs + bin_lhs) / 2.0

    volume_bins = np.array([
        4.0 / 3.0 * math.pi * (bin_edges[i + 1] ** 3 - bin_edges[i] ** 3)
        for i in range(len(bin_edges) - 1)
    ])

    halo_mass_prof = mass_profile[halo_index].astype(np.float64)   # 1e10 Msun per bin
    halo_npart_prof = npart_profile[halo_index].astype(np.float64)  # N_particles per bin

    # density in 1e10 Msun / Mpc^3; multiply by 10 → Msun / kpc^3
    # (1e10 Msun/Mpc^3 = 10 Msun/kpc^3 since 1 Mpc^3 = 1e9 kpc^3)
    density = halo_mass_prof / volume_bins * 10.0   # Msun / kpc^3

    # ------------------------------------------------------------------
    # 2. VR properties (source Cell 8-9)
    # ------------------------------------------------------------------
    with h5py.File(fname_properties, "r") as f:
        R200m = float(np.array(f["R_200mean"])[halo_index])   # Mpc
        C200m = float(np.array(f["cNFW_200mean"])[halo_index])
        M200m_10 = float(np.array(f["Mass_200mean"])[halo_index])  # 1e10 Msun

    M200m_Msun = M200m_10 * 1e10

    # ------------------------------------------------------------------
    # 3. Save simulation profile (filter r < R200, source Cell 11)
    # ------------------------------------------------------------------
    out_sim = np.stack((bin_centre / R200m, density), axis=-1)
    out_sim = out_sim[out_sim[:, 0] < 1.0]
    np.savetxt(str(sim_path), out_sim)

    # ------------------------------------------------------------------
    # 4. NFW (source Cell 11, Mpc units then *1e-9 → Msun/kpc^3)
    # ------------------------------------------------------------------
    rho_halo_mpc = M200m_Msun / (4.0 / 3.0 * math.pi * R200m ** 3)  # Msun/Mpc^3
    A_NFW = math.log(1.0 + C200m) - C200m / (1.0 + C200m)

    nfw_r_mpc = bin_centre[bin_centre <= R200m]
    x_nfw = nfw_r_mpc / R200m
    nfw_rho_mpc = rho_halo_mpc / (3.0 * A_NFW * x_nfw * (1.0 / C200m + x_nfw) ** 2)
    nfw_rho = nfw_rho_mpc * 1e-9  # Msun/Mpc^3 → Msun/kpc^3

    out_nfw = np.stack((nfw_r_mpc / R200m, nfw_rho), axis=-1)
    np.savetxt(str(nfw_path), out_nfw)

    # ------------------------------------------------------------------
    # 5. Convergence radius (Power+03, source Cell 10)
    # ------------------------------------------------------------------
    halo_npart_cumsum = np.cumsum(halo_npart_prof)
    halo_mass_cumsum = np.cumsum(halo_mass_prof)   # 1e10 Msun

    rho_crit = _rho_crit_cgs()
    converg_idx = 0
    for i in range(len(bin_rhs)):
        if halo_npart_cumsum[i] <= 1:
            continue
        mean_density_cgs = (
            halo_mass_cumsum[i] * 1e10 * SOLAR_MASS_IN_CGS
            / (4.0 / 3.0 * math.pi * bin_rhs[i] ** 3 * MPC_IN_CGS ** 3)
        )
        tre = (
            math.sqrt(200.0)
            / 8.0
            * halo_npart_cumsum[i]
            / math.log(halo_npart_cumsum[i])
            * (mean_density_cgs / rho_crit) ** (-0.5)
        )
        if tre > 0.6:
            converg_idx = i
            break
    else:
        import warnings
        warnings.warn(
            f"{sim_label}: Power+03 convergence radius not found; "
            "falling back to innermost bin (converg_idx=0).",
            RuntimeWarning,
            stacklevel=2,
        )

    np.savetxt(str(conv_path), [bin_rhs[converg_idx] / R200m])


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def extract(sim_label: str, *, backend: Backend | None = None, refresh: bool = False) -> None:
    """Extract the radial density profile for one simulation.

    Parameters
    ----------
    sim_label : str
        Simulation label (e.g. ``'m12i_cdmo'``).
    backend : {'hbt', 'vr'} or None
        Which data source to use.  ``None`` (default) uses the paper default:
        ``'hbt'`` for m12i sims, ``'vr'`` for m12f sims.
        ``'hbt'`` reads HBT particle lists + snapshot (memory-heavy; run on a
        compute node).  ``'vr'`` reads VR pre-binned profile files.
    refresh : bool
        Force re-extraction even if cache files already exist.
    """
    if backend is None:
        cfg = load_simulations()
        backend = _DEFAULT_BACKEND[cfg[sim_label].host]
    if backend == "hbt":
        extract_hbt(sim_label, refresh=refresh)
    elif backend == "vr":
        extract_vr(sim_label, refresh=refresh)
    else:
        raise ValueError(f"Unknown backend {backend!r}; choose 'hbt' or 'vr'.")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _interp_data(data: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    """Interpolate (x, y) data onto x_new, skipping y=0 points."""
    x = data[:, 0]
    y = data[:, 1]
    func = interp1d(x[y != 0], y[y != 0], kind="linear", fill_value="extrapolate")
    return np.column_stack((x_new, func(x_new)))


def plot(
    host: Literal["m12i", "m12f"] = "m12i",
    *,
    backend: Backend | None = None,
    refresh: bool = False,
    output: str | Path | None = None,
    use_tex: bool | None = None,
) -> Path:
    """Produce Fig 5 (radial density profile) for the given host galaxy.

    Parameters
    ----------
    host : {'m12i', 'm12f'}
        Which host galaxy to plot.
    backend : {'hbt', 'vr'} or None
        Which cache files to read.  ``None`` (default) uses the paper default:
        ``'hbt'`` for m12i (v1 HBT particle-based centre) and ``'vr'`` for m12f
        (250507 VelocIRaptor-based centre).  Pass an explicit value to override.
    refresh : bool
        Re-extract cache before plotting.
    output : path-like, optional
        Output file path.  Defaults to
        ``figures/Radial_density_profile_{host}.png``.
    use_tex : bool or None
        Enable LaTeX rendering.  ``None`` (default) auto-selects ``True``
        for both backends, mirroring both reference notebooks which call
        ``rc('text', usetex=True)`` and ``rc('font', sans-serif=['Helvetica'])``.
        Pass ``False`` (or use ``--no-tex`` CLI flag) to disable.
    """
    if host not in _HOST_LABELS:
        raise ValueError(f"Unknown host {host!r}; choose 'm12i' or 'm12f'.")

    # Resolve backend: None → paper default for this host.
    if backend is None:
        backend = _DEFAULT_BACKEND[host]

    # Resolve use_tex: None → auto-select from backend.
    # Both reference notebooks (v1 HBT for m12i, 250507 VR for m12f) call:
    #   rc('font', family='sans-serif', sans-serif=['Helvetica'])
    #   rc('text', usetex=True)
    # So the default is True for both backends.
    if use_tex is None:
        use_tex = True

    labels = _HOST_LABELS[host]
    _paths_fn = _hbt_paths if backend == "hbt" else _vr_paths

    # Extract if needed.
    for sl in labels:
        extract(sl, backend=backend, refresh=refresh)

    # Load all 6 cache arrays.
    sim_data = []
    nfw_data = []
    converg_list = []
    for sl in labels:
        sp, np_, cp = _paths_fn(sl)
        sim_data.append(np.genfromtxt(str(sp)))
        nfw_data.append(np.genfromtxt(str(np_)))
        converg_list.append(float(np.genfromtxt(str(cp))))

    # TeX setup (mirrors v1 notebook's rc() calls).
    if use_tex:
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
        rc("text", usetex=True)

    # Choose plot style parameters based on backend / host.
    # v1 HBT (m12i): filter x>1e-4, interpolate ratio, dpi=300.
    # VR 250507 (m12f): no x>1e-4 filter in main loop, direct ratio, dpi=1200.
    use_hbt_style = (backend == "hbt")
    dpi = 300 if use_hbt_style else 1200
    legend_fontsize = 13  # both reference notebooks use ax0.legend(fontsize=13)

    all_labels = (
        _LABELS_TEX[host] if use_tex else _LABELS_PLAIN[host]
    )
    line_styles = _LINE_STYLES  # ("C1", "C2", "C3", "C1--", "C2--", "C3--")

    # Build figure (source: figsize=(6,8), height_ratios=[7,3]).
    fig = plt.figure(figsize=(6, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # ----- Main panel -----
    all_arrays = sim_data + nfw_data  # [cdmo_sim, deep_sim, soft_sim, cdmo_nfw, deep_nfw, soft_nfw]
    for i, (data, label, ls) in enumerate(zip(all_arrays, all_labels, line_styles)):
        if use_hbt_style:
            # v1: filter non-zero AND x > 1e-4
            d = data[(data[:, 1] != 0) & (data[:, 0] > 1e-4)]
        else:
            d = data[data[:, 1] != 0]
        ax0.plot(d[:, 0], d[:, 1], ls, label=label)

    # Convergence shading: v1 uses only cdmo radius; 250507 also only cdmo.
    x1 = 1e-4
    x2 = converg_list[0]  # cdmo convergence radius
    unresolved_label = (
        r"$\rm{unresolved\ region}$" if use_tex else "unresolved region"
    )
    ax0.axvspan(x1, x2, color="grey", alpha=0.3, label=unresolved_label)
    ax1.axvspan(x1, x2, color="grey", alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ylabel_rho = (
        r"$\rho\ \rm{[{\rm M}_{\odot}/(kpc)^3]}$"
        if use_tex
        else r"$\rho\ [M_{sun}/(kpc)^3]$"
    )
    ax0.set_ylabel(ylabel_rho, fontsize=20)
    ax0.legend(fontsize=legend_fontsize)

    # ----- Ratio panel -----
    cdmo_data = sim_data[0]
    deep_data = sim_data[1]
    soft_data = sim_data[2]

    if use_hbt_style:
        # v1: interpolate BT sims to cdmo x-grid (only x>1e-4, non-zero cdmo).
        cdmo_filtered = cdmo_data[(cdmo_data[:, 0] > 1e-4) & (cdmo_data[:, 1] != 0)]
        x_ref = cdmo_filtered[:, 0]
        deep_interp = _interp_data(deep_data, x_ref)
        soft_interp = _interp_data(soft_data, x_ref)
        deep_ratio = deep_interp[:, 1] / cdmo_filtered[:, 1]
        soft_ratio = soft_interp[:, 1] / cdmo_filtered[:, 1]
        ax1.plot(x_ref, deep_ratio, line_styles[1],
                 label=r"$\rm{m12i\_deep}$" if use_tex else "m12i_deep")
        ax1.plot(x_ref, soft_ratio, line_styles[2],
                 label=r"$\rm{m12i\_soft}$" if use_tex else "m12i_soft")
    else:
        # 250507: divide on the full (unfiltered) shared x-grid, then filter.
        # The simulation_1 files share the same bin edges, so arrays are same length.
        deep_ratio = np.transpose(np.array([
            cdmo_data[:, 0],
            np.divide(deep_data[:, 1], cdmo_data[:, 1]),
        ]))
        soft_ratio = np.transpose(np.array([
            cdmo_data[:, 0],
            np.divide(soft_data[:, 1], cdmo_data[:, 1]),
        ]))
        deep_ratio = deep_ratio[(deep_ratio[:, 1] != 0) & (cdmo_data[:, 1] != 0)]
        soft_ratio = soft_ratio[(soft_ratio[:, 1] != 0) & (cdmo_data[:, 1] != 0)]
        ax1.plot(deep_ratio[:, 0], deep_ratio[:, 1], line_styles[1],
                 label=f"{host}:BT_deep")
        ax1.plot(soft_ratio[:, 0], soft_ratio[:, 1], line_styles[2],
                 label=f"{host}:BT_soft")

    ratio_ylabel = r"$N_{\rm BT}/N_{\rm PL}$" if use_tex else "$N_{BT}/N_{PL}$"
    ax1.set_ylabel(ratio_ylabel, fontsize=20)
    # Both reference notebooks have ax1.legend commented out → no ax1 legend.
    ax1.axhline(y=1.0, color="black", linestyle="-", alpha=0.3)

    # Tick styling.
    for ax in (ax0, ax1):
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
        ax.yaxis.set_tick_params(labelsize=15, direction="in", which="both")

    ax0.xaxis.set_tick_params(labelbottom=False)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    xlabel = r"$r/R_{200m}$" if use_tex else r"$ r/R_{200m}$"
    ax1.set_xlabel(xlabel, fontsize=20)

    # Save.
    if output is None:
        out_path = figures_dir() / f"Radial_density_profile_{host}.png"
    else:
        out_path = Path(output)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path
