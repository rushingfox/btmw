"""Fig 6 — Halo Mass Function (HMF) for HBT-HERONS subhalos.
Fig 11 — HBT-HERONS vs VELOCIraptor-SOAP HMF comparison.

Source notebooks:
  HMF/HMF_hbtplus.ipynb                                   (Type 1 — extract HBT)
  HMF/HMF_hbtplus_comparison_with_fitting.ipynb           (Type 2 — plot fig 6)
  comparison_with_VR-SOAP/HMF/HMF_VR-SOAP.ipynb           (Type 1 — extract VR-SOAP)
  comparison_with_VR-SOAP/HMF/HMF_hbtplus_vs_VR_comparison.ipynb  (Type 2 — plot fig 11)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..config import load_simulations
from ..io_soap import read_soap_bundle, subhalo_distances_to_main
from ..paths import cache_dir, external_dir, figures_dir
from ._common import (
    calculate_ratio,
    cumulative_above_log10,
    fitting_plot,
    make_hmf_hvf_figure,
    style_two_panel,
)


# Fiducial sims (paper fig 6) in display order.
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


def _cache_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hmf" / f"hbtplus_HMF_{sim_label}{snap:04d}"


def _meta_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hmf" / f"hbtplus_HMF_{sim_label}{snap:04d}.meta.json"


def extract(sim_label: str, *, refresh: bool = False) -> Path:
    """Read SOAP, build the cumulative HMF, write to ``data/cache/hmf/``."""
    cfg = load_simulations()
    sim = cfg[sim_label]
    out = _cache_path(sim_label, sim.snapshot_num)
    meta_out = _meta_path(sim_label, sim.snapshot_num)
    if out.is_file() and meta_out.is_file() and not refresh:
        return out

    bundle = read_soap_bundle(sim.soap_hbt_path(), need_vmax=False)
    host_M200c = bundle.host_m200c           # in 1e10 Msun (SOAP native)
    host_R200c = float(bundle.r200c.max())   # source notebook uses .max()
    distances = subhalo_distances_to_main(bundle)

    # μ = M_sub / M_host, excluding the main halo, within 1.66 R200c (~ R50c)
    sel = (distances ** 2 < (1.66 * host_R200c) ** 2)
    sel[bundle.main_idx] = False
    mu = bundle.mbound[sel] / host_M200c

    log10_edges = np.arange(-8.0, -0.9, 0.1)   # -8 .. -1 inclusive, step 0.1
    cum = cumulative_above_log10(mu, log10_edges)

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, cum)
    meta_out.write_text(json.dumps({
        "sim_label": sim_label,
        "snapshot_num": sim.snapshot_num,
        "host_m200c_in_1e10Msun": host_M200c,
        "host_r200c_in_Mpc": host_R200c,
        "n_subhalos_in_R50c": int(sel.sum()),
    }, indent=2))
    return out


def _load_or_extract(sim_label: str, *, refresh: bool) -> tuple[np.ndarray, dict]:
    cfg = load_simulations()
    sim = cfg[sim_label]
    cache = _cache_path(sim_label, sim.snapshot_num)
    meta = _meta_path(sim_label, sim.snapshot_num)
    if refresh or not cache.is_file() or not meta.is_file():
        extract(sim_label, refresh=refresh)
    return np.genfromtxt(cache), json.loads(meta.read_text())


def plot(*, refresh: bool = False, output: str | None = None, use_tex: bool = True) -> Path:
    from .. import style as _style
    _style.apply(use_tex=use_tex)
    import matplotlib.pyplot as plt

    fig, ax0, ax1 = make_hmf_hvf_figure()

    datas, metas = [], []
    for label in FIDUCIAL_LABELS:
        data, meta = _load_or_extract(label, refresh=refresh)
        datas.append(data)
        metas.append(meta)

    for i, (label, ls, tex) in enumerate(zip(FIDUCIAL_LABELS, LINESTYLES, LABELS_TEX)):
        d = datas[i]
        d_nz = d[d[:, 1] != 0]
        ax0.plot(d_nz[:, 0], d_nz[:, 1], ls, label=tex)

    # COCO black dashed (fit, in log10 space)
    coco = np.genfromtxt(external_dir() / "hellwing2016_smf_r50c.txt")
    ax0.plot(10 ** coco[:, 0], 10 ** coco[:, 1], "grey",
             label=r"$\rm{COCO\ fit\ for\ }R_{50c}$")

    # Unresolved-region shading: m_DM * 20 / host_M200c (both in 1e10 Msun)
    m_dm = load_simulations().particle_mass_in_swift_unit["fiducial"]
    m12i_host = metas[0]["host_m200c_in_1e10Msun"]
    m12f_host = metas[3]["host_m200c_in_1e10Msun"]
    x1 = 1e-8
    for x2, color, tex in (
        (m_dm * 20.0 / m12i_host, "grey",
         r"$\rm{unresolved\ region\ of\ m12i:PL}$"),
        (m_dm * 20.0 / m12f_host, "brown",
         r"$\rm{unresolved\ region\ of\ m12f:PL}$"),
    ):
        ax0.axvspan(x1, x2, color=color, alpha=0.3, label=tex)
        ax1.axvspan(x1, x2, color=color, alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    # Ratios on the bottom panel
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
        center_mu=0.001, max_enhancement=2,
        model_index=4,
        label_name=r"$\rm{BT\_fitting:} m=2 ~ \mu_c = 10^{-3}$",
        linestyle="-.",
    )

    style_two_panel(ax0, ax1)
    ax0.set_ylabel(r"$N(>\mu)$", fontsize=20)
    ax0.legend(fontsize=11.5)
    ax1.set_ylabel(r"$N_{\rm BT}/N_{\rm PL}$", fontsize=20)
    ax1.axhline(y=1, color="C1", linestyle="-", alpha=0.3)
    ax1.set_xlabel(r"$\mu = M_{\rm sub}/M_{200c}$", fontsize=20)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "HMF_HBTplus_fitting.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig 11 — HBT-HERONS vs VR-SOAP HMF comparison (m12i only)
# ---------------------------------------------------------------------------

# m12i sims only (no m12f VR comparison in archive)
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


def _vr_cache_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hmf" / f"VR-SOAP_HMF_{sim_label}{snap:04d}"


def _vr_meta_path(sim_label: str, snap: int) -> Path:
    return cache_dir() / "hmf" / f"VR-SOAP_HMF_{sim_label}{snap:04d}.meta.json"


def extract_vr(sim_label: str, *, refresh: bool = False) -> Path:
    """Extract HMF from VR-SOAP catalog and cache to ``data/cache/hmf/``."""
    cfg = load_simulations()
    sim = cfg[sim_label]
    if not sim.has_vr:
        raise ValueError(f"{sim_label} has no VR data (has_vr=false in simulations.yaml)")

    out = _vr_cache_path(sim_label, sim.snapshot_num)
    meta_out = _vr_meta_path(sim_label, sim.snapshot_num)
    if out.is_file() and meta_out.is_file() and not refresh:
        return out

    # VR-SOAP uses the same SOAP HDF5 field layout as HBT-SOAP
    bundle = read_soap_bundle(sim.soap_vr_path(), need_vmax=False)
    host_M200c = bundle.host_m200c
    host_R200c = float(bundle.r200c.max())
    distances = subhalo_distances_to_main(bundle)

    sel = (distances ** 2 < (1.66 * host_R200c) ** 2)
    sel[bundle.main_idx] = False
    mu = bundle.mbound[sel] / host_M200c

    log10_edges = np.arange(-8.0, -0.9, 0.1)
    cum = cumulative_above_log10(mu, log10_edges)

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, cum)
    meta_out.write_text(json.dumps({
        "sim_label": sim_label,
        "snapshot_num": sim.snapshot_num,
        "host_m200c_in_1e10Msun": host_M200c,
        "host_r200c_in_Mpc": host_R200c,
        "n_subhalos_in_R50c": int(sel.sum()),
    }, indent=2))
    return out


def _load_or_extract_vr(sim_label: str, *, refresh: bool) -> tuple[np.ndarray, dict]:
    cfg = load_simulations()
    sim = cfg[sim_label]
    cache = _vr_cache_path(sim_label, sim.snapshot_num)
    meta = _vr_meta_path(sim_label, sim.snapshot_num)
    if refresh or not cache.is_file() or not meta.is_file():
        extract_vr(sim_label, refresh=refresh)
    return np.genfromtxt(cache), json.loads(meta.read_text())


def plot_compare_vr(
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 11: HBT-HERONS vs VR-SOAP HMF comparison for m12i."""
    from .. import style as _style
    _style.apply(use_tex=use_tex)
    import matplotlib.pyplot as plt

    fig, ax0, ax1 = make_hmf_hvf_figure()

    # Load HBT-HERONS (reuse fig-6 cache)
    hbt_datas, hbt_metas = [], []
    for label in _VR_LABELS:
        d, m = _load_or_extract(label, refresh=refresh)
        hbt_datas.append(d)
        hbt_metas.append(m)

    # Load VR-SOAP
    vr_datas, vr_metas = [], []
    for label in _VR_LABELS:
        d, m = _load_or_extract_vr(label, refresh=refresh)
        vr_datas.append(d)
        vr_metas.append(m)

    # Main panel: 3 HBT lines (solid) + 3 VR lines (dotted)
    for i, tex in enumerate(_VR_LABELS_TEX_HBT):
        d = hbt_datas[i]
        d_nz = d[d[:, 1] != 0]
        ax0.plot(d_nz[:, 0], d_nz[:, 1], _VR_LINESTYLES_HBT[i], label=tex)

    for i, tex in enumerate(_VR_LABELS_TEX_VR):
        d = vr_datas[i]
        d_nz = d[d[:, 1] != 0]
        ax0.plot(d_nz[:, 0], d_nz[:, 1], _VR_LINESTYLES_VR[i], label=tex)

    # COCO black reference
    coco = np.genfromtxt(external_dir() / "hellwing2016_smf_r50c.txt")
    ax0.plot(10 ** coco[:, 0], 10 ** coco[:, 1], "grey",
             label=r"$\rm{COCO\ fit\ for\ }R_{50c}$")

    # Unresolved-region shading
    m_dm = load_simulations().particle_mass_in_swift_unit["fiducial"]
    hbt_m200c = hbt_metas[0]["host_m200c_in_1e10Msun"]  # m12i_cdmo HBT
    vr_m200c  = vr_metas[0]["host_m200c_in_1e10Msun"]   # m12i_cdmo VR
    x1 = 1e-8
    x2_hbt = m_dm * 20.0 / hbt_m200c
    x2_vr  = m_dm * 20.0 / vr_m200c
    ax0.axvspan(x1, x2_hbt, color="grey",  alpha=0.3,
                label=r"$\rm{unresolved\ region \sim HBT-HERONS}$")
    ax0.axvspan(x1, x2_vr,  color="brown", alpha=0.3,
                label=r"$\rm{unresolved\ region \sim VR}$")
    ax1.axvspan(x1, x2_hbt, color="grey",  alpha=0.3)
    ax1.axvspan(x1, x2_vr,  color="brown", alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    # Ratio panel: VR / HBT-HERONS for each sim
    for i, (ls, tex) in enumerate(zip(_VR_LINESTYLES_HBT, (
        r"$\rm{PL}$", r"$\rm{BT\_deep}$", r"$\rm{BT\_soft}$"
    ))):
        r = calculate_ratio(hbt_datas[i], vr_datas[i])
        ax1.plot(r[:, 0], r[:, 1], ls, label=tex)

    style_two_panel(ax0, ax1)
    ax0.set_ylabel(r"$N(>\mu)$", fontsize=20)
    ax0.legend(fontsize=10.5)
    ax1.set_ylabel(r"$N_{\rm VR}/N_{\rm HBT-HERONS}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_xlabel(r"$\mu = M_{\rm sub}/M_{200c}$", fontsize=20)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "HMF_HBT-HERONS_Comparison_with_VR-SOAP.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Fig 15 — HMF resolution study (fiducial vs high-res m12i)
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
    *,
    refresh: bool = False,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Fig 15 — HMF fiducial-resolution vs high-resolution comparison for m12i."""
    from .. import style as _style
    _style.apply(use_tex=use_tex)
    import matplotlib.pyplot as plt

    fig, ax0, ax1 = make_hmf_hvf_figure()

    fid_data, fid_meta = [], []
    for sl in _RES_FIDUCIAL_LABELS:
        d, m = _load_or_extract(sl, refresh=refresh)
        fid_data.append(d)
        fid_meta.append(m)

    high_data, high_meta = [], []
    for sl in _RES_HIGH_LABELS:
        d, m = _load_or_extract(sl, refresh=refresh)
        high_data.append(d)
        high_meta.append(m)

    # Main panel: fiducial (solid) + high-res (dash-dot)
    for i, (ls, tex) in enumerate(zip(_RES_FIDUCIAL_LS, _RES_FIDUCIAL_TEX)):
        d = fid_data[i]
        ax0.plot(d[d[:, 1] != 0, 0], d[d[:, 1] != 0, 1], ls, label=tex)
    for i, (ls, tex) in enumerate(zip(_RES_HIGH_LS, _RES_HIGH_TEX)):
        d = high_data[i]
        ax0.plot(d[d[:, 1] != 0, 0], d[d[:, 1] != 0, 1], ls, label=tex)

    # COCO fit
    coco = np.genfromtxt(external_dir() / "hellwing2016_smf_r50c.txt")
    ax0.plot(10 ** coco[:, 0], 10 ** coco[:, 1], "grey",
             label=r"$\rm{COCO\ fit\ for\ }R_{50c}$")

    # Unresolved regions
    sims_cfg = load_simulations()
    m_dm_fid  = sims_cfg.particle_mass_in_swift_unit["fiducial"]   # 1e10 Msun
    m12i_fid_host  = fid_meta[0]["host_m200c_in_1e10Msun"]
    m12i_high_host = high_meta[0]["host_m200c_in_1e10Msun"]
    x1 = 1e-8
    x2_fid  = m_dm_fid        * 20.0 / m12i_fid_host
    x2_high = (m_dm_fid / 8.0) * 20.0 / m12i_high_host
    ax0.axvspan(x1, x2_fid,  color="grey",  alpha=0.3,
                label=r"$\rm{unresolved\ region-fiducial\ res}$")
    ax0.axvspan(x1, x2_high, color="brown", alpha=0.3,
                label=r"$\rm{unresolved\ region-high\ res}$")
    ax1.axvspan(x1, x2_fid,  color="grey",  alpha=0.3)
    ax1.axvspan(x1, x2_high, color="brown", alpha=0.3)

    ax0.set_xlim(left=x1)
    ax1.set_xlim(left=x1)

    # Ratio panel: N_high / N_fiducial for each model
    for i, (ls, tex) in enumerate(zip(_RES_FIDUCIAL_LS, _RES_RATIO_TEX)):
        r = calculate_ratio(fid_data[i], high_data[i])
        ax1.plot(r[:, 0], r[:, 1], ls, label=tex)

    style_two_panel(ax0, ax1)
    ax0.set_ylabel(r"$N(>\mu)$", fontsize=20)
    ax0.legend(fontsize=12)
    ax1.set_ylabel(r"$N_{\rm high}/N_{\rm fiducial}$", fontsize=20)
    ax1.axhline(y=1, color="black", linestyle="-", alpha=0.3)
    ax1.legend(fontsize=9)
    ax1.set_xlabel(r"$\mu = M_{\rm sub}/M_{200c}$", fontsize=20)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    out = Path(output) if output else figures_dir() / "HMF_HBTplus_resolution_study.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out
