"""Matplotlib style for paper figures: Helvetica sans-serif with LaTeX text.

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
    # Default DPI for figures saved with the default plt.savefig settings;
    # per-figure savefig calls also pass dpi=300.
    mpl.rcParams["figure.dpi"] = 150
    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["savefig.bbox"] = "tight"
