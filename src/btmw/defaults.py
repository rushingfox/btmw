"""Physical constants and unit conversions shared by the figure modules.

These are the (rounded) values used to produce the published figures, e.g.
``G = 6.67e-8``. Changing them changes the figures.
"""

from __future__ import annotations

# Unit conversions (cgs).
SOLAR_MASS_IN_CGS = 1.989e33                   # g
KPC_IN_CGS = 3.08567758e21                     # cm per kpc
MPC_IN_CGS = 3.08567758e24                     # cm per Mpc
G_IN_CGS = 6.67e-8                             # cm^3 g^-1 s^-2

# Hubble constant used by the paper's WMAP7-based zoom cosmology.
LITTLE_H = 0.702
