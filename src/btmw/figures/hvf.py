"""Fig 7 — Halo Velocity Function (HVF) for HBT-HERONS subhalos.
Fig 12 — HBT-HERONS vs VELOCIraptor-SOAP HVF comparison.

Source notebooks:
  HVF_r100/HVF_hbtplus.ipynb                                    (Type 1 — extract HBT)
  HVF_r100/HVF_hbtplus_m12i_m12f_comparison_fitting.ipynb       (Type 2 — plot fig 7)
  comparison_with_VR-SOAP/HVF/HVF_VR-SOAP.ipynb                 (Type 1 — extract VR-SOAP)
  comparison_with_VR-SOAP/HVF/HVF_comparison_HBTplus_vs_VR.ipynb (Type 2 — plot fig 12)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ..config import load_simulations
from ..defaults import G_IN_CGS, MPC_IN_CGS, SOLAR_MASS_IN_CGS
from ..io_soap import read_soap_bundle, subhalo_distances_to_main
from ..paths import cache_dir, external_dir, figures_dir
from ._common import (
    calculate_ratio,
    cumulative_above_log10,
    fitting_plot,
    make_hmf_hvf_figure,
    style_two_panel,
)


FIDUCIAL_LABELS = (
    "m12i_cdmo",
    "m12i_btps_deep",
    "m12i_btps_soft",
    "m12f_cdmo",
    "m12f_btps_deep",
    "m12f_btps_soft",
)
LINESTYLES = ("C1", "C2", "C3", "C1--", "C2--", "C3--")
LABELS_TEX = (
    r"$\rm{m12i:PL}$",
    r"$\rm{m12i:BT\_deep}$",
    r"$\rm{m12i:BT\_soft}$",
    r"$\rm{m12f:PL}$",
    r"$\rm{m12f:BT\_deep}$",
    r"$\rm{m12f:BT\_soft}$",
)

# Hardcoded unresolved-region bounds (Vmax/V200c) from the source plot
# notebook ``HVF_hbtplus_m12i_m12f_comparison_fitting.ipynb`` (cell 6).
# These are the 0th-percentile of the 100-DM-particle subhalo Vmax in
# the m12i and m12f PL runs, computed inside that notebook.
UNRESOLVED_BOUND = {
    "m12i": 0.02995634129285346,
    "m12f": 0.0247377586709585,
}


def _cache_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hvf" / f"hbtplus_HVF_{sim_label}{snap:04d}"


def _meta_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hvf" / f"hbtplus_HVF_{sim_label}{snap:04d}.meta.json"


def extract(sim_label: str, *, refresh: bool = False) -> Path:
    cfg = load_simulations()
    sim = cfg[sim_label]
    out = _cache_path(sim_label, sim.snapshot_num)
    meta_out = _meta_path(sim_label, sim.snapshot_num)
    if out.is_file() and meta_out.is_file() and not refresh:
        return out

    bundle = read_soap_bundle(sim.soap_hbt_path(), need_vmax=True)
    main = bundle.main_idx
    M200c_main = float(bundle.m200c[main])     # 1e10 Msun
    R200c_main = float(bundle.r200c[main])     # Mpc
    R100c_main = float(bundle.r100c[main])     # Mpc; use [main] not .max() for correctness

    V200c_in_cgs = math.sqrt(
        G_IN_CGS * 1e10 * SOLAR_MASS_IN_CGS * M200c_main
        / (R200c_main * MPC_IN_CGS)
    )

    distances = subhalo_distances_to_main(bundle)
    sel = (distances ** 2 < R100c_main ** 2)
    sel[main] = False
    # Vmax is in km/s in SOAP; *1e5 converts to cm/s
    nu = bundle.vmax[sel] * 1e5 / V200c_in_cgs

    log10_edges = np.arange(-3.0, 0.1, 0.1)    # -3 .. 0 inclusive
    cum = cumulative_above_log10(nu, log10_edges)

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, cum)
    meta_out.write_text(json.dumps({
        "sim_label": sim_label,
        "snapshot_num": sim.snapshot_num,
        "host_m200c_in_1e10Msun": M200c_main,
        "host_r200c_in_Mpc": R200c_main,
        "host_r100c_in_Mpc": R100c_main,
        "V200c_in_cgs": V200c_in_cgs,
        "n_subhalos_in_R100c": int(sel.sum()),
    }, indent=2))
    return out


def _load_or_extract(sim_label: str, *, refresh: bool) -> np.ndarray:
    cfg = load_simulations()
    sim = cfg[sim_label]
    cache = _cache_path(sim_label, sim.snapshot_num)
    if refresh or not cache.is_file():
        extract(sim_label, refresh=refresh)
    return np.genfromtxt(cache)


def plot(*, refresh: bool = False, output: str | None = None, use_tex: bool = True) -> Path:
    from .. import style as _style
    _style.apply(use_tex=use_tex)
    import matplotlib.pyplot as plt

    fig, ax0, ax1 = make_hmf_hvf_figure()

    datas = [_load_or_extract(label, refresh=refresh) for label in FIDUCIAL_LABELS]

    for i, (ls, tex) in enumerate(zip(LINESTYLES, LABELS_TEX)):
        d = datas[i]
        d_nz = d[d[:, 1] != 0]
        ax0.plot(d_nz[:, 0], d_nz[:, 1], ls, label=tex)

    coco = np.genfromtxt(external_dir() / "cautun2014_svf_r100c.txt")
    ax0.plot(10 ** coco[:, 0], 10 ** coco[:, 1], "grey",
             label=r"$\rm{Cautun\ fit\ for\ }R_{100c}$")

    x1 = 1e-2
    for host, color, tex in (
        ("m12i", "grey",
         r"$\rm{unresolved\ region\ of\ m12i:PL}$"),
        ("m12f", "brown",
         r"$\rm{unresolved\ region\ of\ m12f:PL}$"),
    ):
        x2 = UNRESOLVED_BOUND[host]
        ax0.axvspan(x1, x2, color=color, alpha=0.3, label=tex)
        ax1.axvspan(x1, x2, color=color, alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    m12i_pl, m12i_deep, m12i_soft = datas[0], datas[1], datas[2]
    m12f_pl, m12f_deep, m12f_soft = datas[3], datas[4], datas[5]

    r = calculate_ratio(m12i_pl, m12i_deep)
    ax1.plot(r[:, 0], r[:, 1], LINESTYLES[1], label=r"$\rm{m12i\_deep}$")
    r = calculate_ratio(m12i_pl, m12i_soft)
    ax1.plot(r[:, 0], r[:, 1], LINESTYLES[2], label=r"$\rm{m12i\_soft}$")
    r = calculate_ratio(m12f_pl, m12f_deep)
    ax1.plot(r[:, 0], r[:, 1], LINESTYLES[4], label=r"$\rm{m12f\_deep}$")
    r = calculate_ratio(m12f_pl, m12f_soft)
    ax1.plot(r[:, 0], r[:, 1], LINESTYLES[5], label=r"$\rm{m12f\_soft}$")

    fitting_plot(
        ax0, ax1, m12i_pl,
        center_mu=0.15, max_enhancement=3.5,
        model_index=2,
        label_name=r"$\rm{BT\_deep\ fitting:} m=3.5 ~ \nu_c = 0.15$",
        linestyle=":",
    )
    fitting_plot(
        ax0, ax1, m12i_pl,
        center_mu=0.2, max_enhancement=5,
        model_index=3,
        label_name=r"$\rm{BT\_soft\ fitting:} m=5 ~ \nu_c = 0.2$",
        linestyle=":",
    )

    style_two_panel(ax0, ax1)
    ax0.set_ylabel(r"$N(>\nu)$", fontsize=20)
    ax0.legend(fontsize=10)
    ax1.set_ylabel(r"$N_{\rm BT}/N_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="C1", linestyle="-", alpha=0.3)
    ax1.set_xlabel(r"$\nu = V_{\rm max}/V_{200c}$", fontsize=20)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "HVF_HBTplus_fitting.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig 12 — HBT-HERONS vs VR-SOAP HVF comparison (m12i only)
# ---------------------------------------------------------------------------

# m12i sims only (matches archive)
_VR_LABELS = ("m12i_cdmo", "m12i_btps_deep", "m12i_btps_soft")
_VR_LINESTYLES_HBT = ("C1",  "C2",  "C3")
_VR_LINESTYLES_VR  = ("C1:", "C2:", "C3:")
_VR_LABELS_TEX_HBT = (
    r"$\rm{m12i:PL \sim HBT-HERONS}$",
    r"$\rm{m12i:BT\_deep \sim HBT-HERONS}$",
    r"$\rm{m12i:BT\_soft \sim HBT-HERONS}$",
)
_VR_LABELS_TEX_VR = (
    r"$\rm{m12i:PL \sim VR}$",
    r"$\rm{m12i:BT\_deep \sim VR}$",
    r"$\rm{m12i:BT\_soft \sim VR}$",
)

# Hardcoded from comparison notebook cell 6 (identical to HVF UNRESOLVED_BOUND m12i)
_UNRESOLVED_HBT = 0.02995634129285346
_UNRESOLVED_VR  = 0.03009951542632502


def _vr_cache_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hvf" / f"VR-SOAP_HVF_{sim_label}{snap:04d}"


def extract_vr(sim_label: str, *, refresh: bool = False) -> Path:
    """Extract HVF from VR-SOAP catalog, cache to ``data/cache/hvf/``."""
    cfg = load_simulations()
    sim = cfg[sim_label]
    if not sim.has_vr:
        raise ValueError(f"{sim_label} has no VR data (has_vr=false)")

    out = _vr_cache_path(sim_label, sim.snapshot_num)
    if out.is_file() and not refresh:
        return out

    bundle = read_soap_bundle(sim.soap_vr_path(), need_vmax=True)
    main = bundle.main_idx
    M200c_main = float(bundle.m200c[main])
    R200c_main = float(bundle.r200c[main])
    R100c_main = float(bundle.r100c[main])     # use [main] not .max() for correctness

    V200c_in_cgs = math.sqrt(
        G_IN_CGS * 1e10 * SOLAR_MASS_IN_CGS * M200c_main
        / (R200c_main * MPC_IN_CGS)
    )

    distances = subhalo_distances_to_main(bundle)
    sel = (distances ** 2 < R100c_main ** 2)
    sel[main] = False
    nu = bundle.vmax[sel] * 1e5 / V200c_in_cgs

    log10_edges = np.arange(-3.0, 0.1, 0.1)
    cum = cumulative_above_log10(nu, log10_edges)

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, cum)
    return out


def _load_or_extract_vr(sim_label: str, *, refresh: bool) -> np.ndarray:
    cfg = load_simulations()
    sim = cfg[sim_label]
    cache = _vr_cache_path(sim_label, sim.snapshot_num)
    if refresh or not cache.is_file():
        extract_vr(sim_label, refresh=refresh)
    return np.genfromtxt(cache)


def plot_compare_vr(
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 12: HBT-HERONS vs VR-SOAP HVF comparison for m12i."""
    from .. import style as _style
    _style.apply(use_tex=use_tex)
    import matplotlib.pyplot as plt

    fig, ax0, ax1 = make_hmf_hvf_figure()

    hbt_datas = [_load_or_extract(label, refresh=refresh) for label in _VR_LABELS]
    vr_datas  = [_load_or_extract_vr(label, refresh=refresh) for label in _VR_LABELS]

    for i, tex in enumerate(_VR_LABELS_TEX_HBT):
        d = hbt_datas[i]
        ax0.plot(d[d[:, 1] != 0, 0], d[d[:, 1] != 0, 1],
                 _VR_LINESTYLES_HBT[i], label=tex)

    for i, tex in enumerate(_VR_LABELS_TEX_VR):
        d = vr_datas[i]
        ax0.plot(d[d[:, 1] != 0, 0], d[d[:, 1] != 0, 1],
                 _VR_LINESTYLES_VR[i], label=tex)

    coco = np.genfromtxt(external_dir() / "cautun2014_svf_r100c.txt")
    ax0.plot(10 ** coco[:, 0], 10 ** coco[:, 1], "grey",
             label=r"$\rm{Cautun\ fit\ for\ }R_{100c}$")

    x1 = 1e-2
    ax0.axvspan(x1, _UNRESOLVED_HBT, color="grey",  alpha=0.3,
                label=r"$\rm{unresolved\ region \sim HBT-HERONS}$")
    ax0.axvspan(x1, _UNRESOLVED_VR,  color="brown", alpha=0.3,
                label=r"$\rm{unresolved\ region \sim VR}$")
    ax1.axvspan(x1, _UNRESOLVED_HBT, color="grey",  alpha=0.3)
    ax1.axvspan(x1, _UNRESOLVED_VR,  color="brown", alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    for i, (ls, tex) in enumerate(zip(_VR_LINESTYLES_HBT, (
        r"$\rm{PL}$", r"$\rm{BT\_deep}$", r"$\rm{BT\_soft}$"
    ))):
        r = calculate_ratio(hbt_datas[i], vr_datas[i])
        ax1.plot(r[:, 0], r[:, 1], ls, label=tex)

    style_two_panel(ax0, ax1)
    ax0.set_ylabel(r"$N(>\nu)$", fontsize=20)
    ax0.legend(fontsize=10.5)
    ax1.set_ylabel(r"$N_{VR}/N_{HBT-HERONS}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_xlabel(r"$\nu = V_{\rm max}/V_{200c}$", fontsize=20)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "HVF_HBT-HERONS_comparison_with_VR-SOAP.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig 16 — HVF resolution study (fiducial vs high-res m12i)
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
# Fiducial unresolved bound: same as UNRESOLVED_BOUND["m12i"]
_RES_UNRES_FID  = 0.02995634129285346
# High-res unresolved bound (from resolution_study notebook)
_RES_UNRES_HIGH = 0.01551527244439509


def plot_resolution_study(
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 16 — HVF fiducial-resolution vs high-resolution comparison for m12i."""
    from .. import style as _style
    _style.apply(use_tex=use_tex)
    import matplotlib.pyplot as plt

    fig, ax0, ax1 = make_hmf_hvf_figure()

    fid_data  = [_load_or_extract(sl, refresh=refresh) for sl in _RES_FIDUCIAL_LABELS]
    high_data = [_load_or_extract(sl, refresh=refresh) for sl in _RES_HIGH_LABELS]

    # Main panel: fiducial (solid) + high-res (dash-dot)
    for i, (ls, tex) in enumerate(zip(_RES_FIDUCIAL_LS, _RES_FIDUCIAL_TEX)):
        d = fid_data[i]
        ax0.plot(d[d[:, 1] != 0, 0], d[d[:, 1] != 0, 1], ls, label=tex)
    for i, (ls, tex) in enumerate(zip(_RES_HIGH_LS, _RES_HIGH_TEX)):
        d = high_data[i]
        ax0.plot(d[d[:, 1] != 0, 0], d[d[:, 1] != 0, 1], ls, label=tex)

    # Cautun (2014) fit
    coco = np.genfromtxt(external_dir() / "cautun2014_svf_r100c.txt")
    ax0.plot(10 ** coco[:, 0], 10 ** coco[:, 1], "grey",
             label=r"$\rm{COCO\ fit\ for\ }R_{100c}$")

    # Unresolved regions
    x1 = 10 ** (-2.5)
    ax0.axvspan(x1, _RES_UNRES_FID,  color="grey",  alpha=0.3,
                label=r"$\rm{unresolved\ region-fiducial\ res}$")
    ax0.axvspan(x1, _RES_UNRES_HIGH, color="brown", alpha=0.3,
                label=r"$\rm{unresolved\ region-high\ res}$")
    ax1.axvspan(x1, _RES_UNRES_FID,  color="grey",  alpha=0.3)
    ax1.axvspan(x1, _RES_UNRES_HIGH, color="brown", alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    # Ratio panel: N_high / N_fiducial
    for i, (ls, tex) in enumerate(zip(_RES_FIDUCIAL_LS, _RES_RATIO_TEX)):
        r = calculate_ratio(fid_data[i], high_data[i])
        ax1.plot(r[:, 0], r[:, 1], ls, label=tex)

    style_two_panel(ax0, ax1)
    ax0.set_ylabel(r"$N(>\nu)$", fontsize=20)
    ax0.legend(fontsize=12)
    ax1.set_ylabel(r"$N_{\rm high}/N_{\rm fiducial}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    ax1.legend(fontsize=9)
    ax1.set_xlabel(r"$\nu = V_{\rm max}/V_{200c}$", fontsize=20)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "HVF_HBTplus_resolution_study.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out
