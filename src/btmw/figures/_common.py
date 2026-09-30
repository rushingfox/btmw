"""Shared helpers for the HMF / HVF (and similar two-panel) figures.

Bin edges, the fitting form and the ratio interpolation determine the
published figures and the checked-in caches; changing them changes the
output.
"""

from __future__ import annotations

import numpy as np
from matplotlib import gridspec
from scipy.interpolate import interp1d


# ---------------------------------------------------------------------------
# Cumulative N(>x) on a fixed log10 grid
# ---------------------------------------------------------------------------

def cumulative_above_log10(values: np.ndarray, log10_edges: np.ndarray) -> np.ndarray:
    """Return ``(10**LHS_edges, N(>=LHS_edge))``.

    Histograms ``log10(values)`` on ``log10_edges``, builds
    ``N(>=RHS_of_bin) = total - cumsum`` and then adds the bin's own count
    back to convert to ``N(>=LHS_of_bin)``.
    """
    countdata, bins = np.histogram(np.log10(values), bins=log10_edges)
    LHS_bins = bins[:-1]
    data_total = np.full(len(countdata), np.sum(countdata), dtype=int)
    cumulative = np.cumsum(countdata)
    inverse_cumulative = data_total - cumulative + countdata
    return np.column_stack((10 ** LHS_bins, inverse_cumulative))


# ---------------------------------------------------------------------------
# Ratio / interpolation helpers
# ---------------------------------------------------------------------------

def interp_data(data: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x = data[:, 0]
    y = data[:, 1]
    interp_func = interp1d(x[y != 0], y[y != 0], kind="linear", fill_value="extrapolate")
    return np.column_stack((x_new[:], interp_func(x_new[:])))


def calculate_ratio(base_data: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Ratio ``data / base_data`` on the bins where both are non-zero."""
    if not np.allclose(base_data[:, 0], data[:, 0]):
        raise ValueError("calculate_ratio: the two arrays use different binning.")
    nz_base = np.where(base_data[:, 1] != 0)[0]
    nz_data = np.where(data[:, 1] != 0)[0]
    common = np.intersect1d(nz_base, nz_data)
    common_x = base_data[common, 0]
    base_interp = interp_data(base_data, common_x)
    data_interp = interp_data(data, common_x)
    return np.column_stack((common_x, np.divide(data_interp[:, 1], base_interp[:, 1])))


# ---------------------------------------------------------------------------
# Blue-tilt fitting form used by both HMF and HVF
# ---------------------------------------------------------------------------

def blue_fitting(mu: np.ndarray, center_mu: float, max_enhancement: float) -> np.ndarray:
    ratio_range = max_enhancement - 1
    scaled_x = (np.log10(mu) - np.log10(center_mu)) * 4
    return (-scaled_x / np.sqrt(1 + scaled_x ** 2) + 1) / 2 * ratio_range + 1


def fitting_plot(
    ax_main,
    ax_ratio,
    PL_data: np.ndarray,
    center_mu: float,
    max_enhancement: float,
    model_index: int,
    label_name: str,
    *,
    linestyle: str = "-.",
) -> None:
    """Overlay the analytic BT fit on both panels."""
    color_style = f"C{int(model_index)}{linestyle}"
    fitting_ratio = np.column_stack(
        (PL_data[:, 0], blue_fitting(PL_data[:, 0], center_mu, max_enhancement))
    )
    ax_ratio.plot(fitting_ratio[:, 0], fitting_ratio[:, 1], color_style)
    fitting_data = np.column_stack(
        (PL_data[:, 0], fitting_ratio[:, 1] * PL_data[:, 1])
    )
    fitting_data = fitting_data[fitting_data[:, 1] != 0]
    ax_main.plot(fitting_data[:, 0], fitting_data[:, 1], color_style, label=label_name)


# ---------------------------------------------------------------------------
# Two-panel figure scaffold (top: cumulative; bottom: ratio)
# ---------------------------------------------------------------------------

def make_hmf_hvf_figure():
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 3])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    return fig, ax0, ax1


def style_two_panel(ax0, ax1) -> None:
    for ax in (ax0, ax1):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
    ax0.xaxis.set_tick_params(labelsize=15, direction="in", which="both", labelbottom=False)
    ax0.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax1.xaxis.set_tick_params(labelsize=15, direction="in", which="both")
    ax1.yaxis.set_tick_params(labelsize=15, direction="in", which="both")
