"""Matplotlib style for paper figures.

Replicates the rcParams used in the original Type-2 ipynbs:

    from matplotlib import rc
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    %config InlineBackend.figure_format = 'retina'

A ``--no-tex`` switch disables LaTeX rendering for environments without a
TeX install; everything else stays the same.
"""

from __future__ import annotations

import matplotlib as mpl


def apply(*, use_tex: bool = True) -> None:
    """Apply the paper's matplotlib style. Idempotent."""
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Helvetica"]
    mpl.rcParams["text.usetex"] = bool(use_tex)
    # Retina-equivalent default DPI for figures saved with the default
    # plt.savefig settings; per-figure savefig calls also pass dpi=300.
    mpl.rcParams["figure.dpi"] = 150
    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["savefig.bbox"] = "tight"
