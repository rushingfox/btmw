"""Fig 4 — DM projection map (200 kpc half-width) for the m12i runs.

Mirror of ``paper_figures_min_20_HBTplus_v1_archive/Projection_Map/
Projection_map0723.ipynb``.  One PNG per simulation; the published
figure stitches three panels (m12i:PL / BT_deep / BT_soft) together
externally.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..config import load_simulations
from ..io_soap import read_soap_bundle
from ..paths import figures_dir


# Mapping from sim_label -> short LaTeX title used in the published panel.
_TITLE_BY_LABEL = {
    "m12i_cdmo":      r"$\rm{m12i:PL}$",
    "m12i_btps_deep": r"$\rm{m12i:BT\_deep}$",
    "m12i_btps_soft": r"$\rm{m12i:BT\_soft}$",
    "m12f_cdmo":      r"$\rm{m12f:PL}$",
    "m12f_btps_deep": r"$\rm{m12f:BT\_deep}$",
    "m12f_btps_soft": r"$\rm{m12f:BT\_soft}$",
}


def plot(
    sim_label: str,
    *,
    width_kpc: float = 200.0,
    image_resolution: int = 600,
    output: str | None = None,
    use_tex: bool = True,
) -> Path:
    """Render the DM mass projection for one sim and write a PNG.

    Parameters
    ----------
    sim_label
        Simulation label declared in ``configs/simulations.yaml``.
    width_kpc
        Half-width of the imaged region (kpc).  The published figure uses
        200 kpc, i.e. a 400×400 kpc patch centered on the main halo.
    image_resolution
        ``project_pixel_grid`` resolution (pixels on a side).  600 in the
        source notebook.
    output
        Optional explicit output path; otherwise writes
        ``figures/<sim_label>_projection_map_<width_kpc:.0f>kpc.png``.
    """
    from .. import style as _style
    _style.apply(use_tex=use_tex)

    import matplotlib.pyplot as plt
    from swiftsimio import load
    from swiftsimio.objects import cosmo_array
    from swiftsimio.visualisation import (
        generate_smoothing_lengths,
        project_pixel_grid,
    )

    cfg = load_simulations()
    sim = cfg[sim_label]

    # --- main halo center from SOAP (matches source notebook cell 2) ---
    bundle = read_soap_bundle(sim.soap_hbt_path(), need_vmax=False)
    cx, cy, _cz = bundle.host_centre

    # --- SWIFT snapshot + DM smoothing lengths (source cell 1) ---
    data = load(str(sim.snapshot_path()))
    data.dark_matter.smoothing_length = generate_smoothing_lengths(
        data.dark_matter.coordinates,
        data.metadata.boxsize,
        kernel_gamma=1.8,
        neighbours=32,
        speedup_fac=2,
        dimension=3,
    )

    # --- region: half-width in Mpc (source uses r_around = 0.2 Mpc = 200 kpc) ---
    r_around_Mpc = width_kpc / 1000.0
    coords = data.dark_matter.coordinates
    region = cosmo_array(
        [cx - r_around_Mpc, cx + r_around_Mpc,
         cy - r_around_Mpc, cy + r_around_Mpc],
        units="Mpc",
        comoving=coords.comoving,
        cosmo_factor=coords.cosmo_factor,
    )

    dm_mass = project_pixel_grid(
        data=data.dark_matter,
        resolution=image_resolution,
        project="masses",
        parallel=True,
        region=region,
    )

    dm_mass_arr = np.asarray(getattr(dm_mass, "value", dm_mass))
    # The source notebook used the old swiftsimio API which returned surface
    # density in raw mass units (10^10 Msun / Mpc^2).  In modern swiftsimio
    # (>=12) the output is normalised to Msun / Mpc^2; rescale so the
    # vmin=1 / vmax=4 colour range from the source notebook still applies.
    dm_mass_arr = dm_mass_arr / 1e10

    fig, ax = plt.subplots()
    ax.imshow(np.log10(dm_mass_arr), cmap="inferno", vmin=1, vmax=4)
    ax.set_axis_off()
    ax.set_title(_TITLE_BY_LABEL.get(sim_label, sim_label))

    if output is None:
        out = figures_dir() / f"{sim_label}_projection_map_{int(width_kpc)}kpc.png"
    else:
        out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out, format="png", dpi=image_resolution, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return out
