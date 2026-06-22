"""h5py readers for SOAP halo catalogs (HBT-HERONS- and VR-backed).

The HMF/HVF (and most other) figures only need a small bundle of SOAP
fields plus a simple "main halo = argmax of BoundSubhalo/TotalMass"
identification. ``read_soap_bundle`` returns that bundle as a plain
dataclass so the figure code never touches h5py directly.

Conventions (verified against the original HMF/HVF source notebooks in the
paper-figure archive):

- SOAP stores TotalMass in units of 1e10 Msun.
- SOAP stores SORadius in Mpc (SWIFT internal length unit, h-free).
- SOAP stores MaximumCircularVelocityUnsoftened in km/s
  (the source notebook then multiplies by 1e5 to convert to cm/s).
- HaloCentre is in Mpc.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


# ---------------------------------------------------------------------------
# SOAP path constants (kept for figures that only need a single field).
# ---------------------------------------------------------------------------
SOAP_BOUND_TOTAL_MASS     = "/BoundSubhalo/TotalMass"
SOAP_BOUND_NPART_DM       = "/BoundSubhalo/NumberOfDarkMatterParticles"
SOAP_BOUND_VMAX           = "/BoundSubhalo/MaximumCircularVelocityUnsoftened"
SOAP_SO200C_TOTAL_MASS    = "/SO/200_crit/TotalMass"
SOAP_SO200C_RADIUS        = "/SO/200_crit/SORadius"
SOAP_SO100C_RADIUS        = "/SO/100_crit/SORadius"
SOAP_SO50C_RADIUS         = "/SO/50_crit/SORadius"
SOAP_SO200M_TOTAL_MASS    = "/SO/200_mean/TotalMass"
SOAP_SO200M_RADIUS        = "/SO/200_mean/SORadius"
SOAP_SO200M_CONCENTRATION = "/SO/200_mean/Concentration"
SOAP_HALO_CENTRE          = "/InputHalos/HaloCentre"


@dataclass(frozen=True)
class SoapBundle:
    """Bundle of SOAP arrays needed for the HMF / HVF / radial figures.

    All arrays are 1-D with length N_halos, except ``halo_centre`` which is
    (N_halos, 3). Units follow SOAP native (1e10 Msun for masses, Mpc for
    lengths, km/s for Vmax).
    """

    main_idx: int
    mbound: np.ndarray              # /BoundSubhalo/TotalMass             [1e10 Msun]
    npart_bound: np.ndarray         # /BoundSubhalo/NumberOfDarkMatterParticles
    vmax: np.ndarray                # /BoundSubhalo/MaximumCircularVelocityUnsoftened [km/s]
    m200c: np.ndarray               # /SO/200_crit/TotalMass              [1e10 Msun]
    r200c: np.ndarray               # /SO/200_crit/SORadius               [Mpc]
    r100c: np.ndarray               # /SO/100_crit/SORadius               [Mpc]
    r50c: np.ndarray                # /SO/50_crit/SORadius                [Mpc]
    halo_centre: np.ndarray         # /InputHalos/HaloCentre              [Mpc, shape (N,3)]

    @property
    def host_m200c(self) -> float:
        return float(self.m200c[self.main_idx])

    @property
    def host_r200c(self) -> float:
        return float(self.r200c[self.main_idx])

    @property
    def host_r100c(self) -> float:
        return float(self.r100c[self.main_idx])

    @property
    def host_centre(self) -> np.ndarray:
        return self.halo_centre[self.main_idx]


def read_soap_bundle(path: str | Path, *, need_vmax: bool = True) -> SoapBundle:
    """Open a SOAP halo-properties HDF5 and return the standard bundle.

    Parameters
    ----------
    path
        Absolute path to ``halo_properties_{snap:04d}.hdf5``.
    need_vmax
        If ``True`` (default), also load Vmax / Npart / R100c. Set to
        ``False`` for figures that only use mass/radius.
    """
    path = Path(path)
    with h5py.File(path, "r") as f:
        mbound      = np.array(f[SOAP_BOUND_TOTAL_MASS])
        m200c       = np.array(f[SOAP_SO200C_TOTAL_MASS])
        r200c       = np.array(f[SOAP_SO200C_RADIUS])
        r50c        = np.array(f[SOAP_SO50C_RADIUS])
        halo_centre = np.array(f[SOAP_HALO_CENTRE])
        if need_vmax:
            npart_bound = np.array(f[SOAP_BOUND_NPART_DM])
            vmax        = np.array(f[SOAP_BOUND_VMAX])
            r100c       = np.array(f[SOAP_SO100C_RADIUS])
        else:
            npart_bound = np.zeros_like(mbound, dtype=np.int64)
            vmax        = np.zeros_like(mbound)
            r100c       = np.zeros_like(mbound)

    main_idx = int(np.argmax(mbound))
    return SoapBundle(
        main_idx=main_idx,
        mbound=mbound,
        npart_bound=npart_bound,
        vmax=vmax,
        m200c=m200c,
        r200c=r200c,
        r100c=r100c,
        r50c=r50c,
        halo_centre=halo_centre,
    )


def subhalo_distances_to_main(bundle: SoapBundle) -> np.ndarray:
    """Euclidean distance of every halo to the main halo, in Mpc."""
    delta = bundle.halo_centre - bundle.host_centre
    return np.sqrt(np.einsum("ij,ij->i", delta, delta))
