"""Physical constants, unit conversions, and shared figure defaults.

Values match the originals used in the paper's Type-1 ipynbs so that
btmw reproduces the published numbers byte-for-byte. Do not "tidy" them
without first confirming the change does not affect any figure.
"""

from __future__ import annotations

# Unit conversions copied from the original notebooks.
SOLAR_MASS_IN_CGS = 1.989e33                   # g
SWIFT_UNIT_MASS_IN_CGS = 1e10 * SOLAR_MASS_IN_CGS
SWIFT_UNIT_MASS_IN_MSUN = 1e10                 # SWIFT internal mass unit, in Msun
KPC_IN_CGS = 3.08567758e21                     # cm per kpc
MPC_IN_CGS = 3.08567758e24                     # cm per Mpc
G_IN_CGS = 6.67e-8                             # cm^3 g^-1 s^-2

# Hubble constant used by the paper's WMAP7-based zoom cosmology.
LITTLE_H = 0.702

# COCO empirical relation between R50c and R200c (used by HMF / CMF / HRF).
R50C_OVER_R200C = 1.66

# Halo-finder lower limit reflected by the *_min_20 path naming.
MIN_NPART = 20

# Per-resolution DM particle mass in SWIFT internal units (1e10 Msun / h)
# Source: Type-1 ipynbs in paper_figures_min_20_HBTplus_v1_archive/HMF/.
DM_PARTICLE_MASS_FIDUCIAL_SWIFT = 4.2248334e-06
