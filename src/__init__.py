#!/usr/bin/env python
# Copyright (C) 2025 Emanuel Goncalves

import sys
import config as C
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib as mpl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Figure style
GREY, RED = "#8a949e", "#c0392b"
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    }
)
mpl.rcParams["figure.dpi"] = 300

# Figure palette
CELL_LINE_COLOUR = {"A375": "#0072B2", "HT29": "#D55E00"}
# fmt: off
ARM_COLOUR = {
    ("A375", "DO"): "#014975", ("A375", "LIN"): "#0072B2",
    ("A375", "PIC"): "#619FD1", ("A375", "SCH"): "#A7C9E7",
    ("HT29", "DO"): "#974101", ("HT29", "LIN"): "#D55E00",
    ("HT29", "PIC"): "#F79A6A", ("HT29", "SCH"): "#FFD9C6",
}
# fmt: on
PLOT_DPI = 600

# Pathway hues, for panels whose axes already carry the cell lines.
PATHWAY_COLOUR = {
    "mTOR / PI3K\npathway": "#3f8f6b",
    "MAPK/ERK\npathway": "#7b5aa6",
    "Miscellaneous\ngenes": "#9aa0a6",
}


def arm_colour(cell_line: str, arm: str) -> str:
    """The shade reserved for this cell line and arm."""
    return ARM_COLOUR[(cell_line, arm)]


# Reference lines are dashed so they never read as data.
DIAGONAL = dict(lw=0.5, color="0.55", linestyle=(0, (4, 3)), zorder=0)
ORIGIN = dict(lw=0.4, color="0.8", linestyle=(0, (3, 3)), zorder=0)


def axis_ticks(low: float, high: float, count: int = 5):
    """Round tick positions that cover the range, with a tick at each end.

    The axis limits are then set to the first and last tick, so no point is cut off
    and the axis starts and ends on a labelled tick.
    """
    step = (high - low) / max(count, 1)
    magnitude = 10 ** np.floor(np.log10(step)) if step > 0 else 1.0
    for multiple in (1, 2, 2.5, 5, 10):
        if step <= multiple * magnitude:
            step = multiple * magnitude
            break
    start = np.floor(low / step) * step
    stop = np.ceil(high / step) * step
    return np.arange(start, stop + step / 2, step)


def ticked_axis(ax, axis: str, low: float, high: float, count: int = 5):
    """Round ticks covering the range on one axis, with a tick at each end."""
    ticks = axis_ticks(low, high, count)
    if axis == "x":
        ax.set_xticks(ticks)
        ax.set_xlim(ticks[0], ticks[-1])
    else:
        ax.set_yticks(ticks)
        ax.set_ylim(ticks[0], ticks[-1])
    return ticks


def square_axes(ax, values, count: int = 5):
    """Give a scatter equal, identically ticked x and y axes covering the data."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    ticks = axis_ticks(finite.min(), finite.max(), count)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_ylim(ticks[0], ticks[-1])
    ax.set_aspect("equal")
    return ticks


def trim_spines(ax) -> None:
    """End each spine on a tick, so every axis starts and finishes at a tick mark."""
    for side, axis, limits in (
        ("bottom", ax.xaxis, ax.get_xlim()),
        ("left", ax.yaxis, ax.get_ylim()),
    ):
        if side not in ax.spines:
            continue
        low, high = sorted(limits)
        visible = [
            t for t in axis.get_majorticklocs() if low - 1e-9 <= t <= high + 1e-9
        ]
        if len(visible) >= 2:
            ax.spines[side].set_bounds(min(visible), max(visible))


def save_plot(fig, name: str, dpi: int = PLOT_DPI) -> Path:
    """Write a panel to plots/<name>.png, which is not committed."""
    for ax in fig.axes:
        trim_spines(ax)
    C.PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(C.PLOTS / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(C.PLOTS / f"{name}.svg", bbox_inches="tight", facecolor="white")
    return C.PLOTS / f"{name}.pdf"

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYVBZX*")
BLOSUM62 = pd.DataFrame([
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
    [-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
    [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1],
], index=AMINO_ACIDS, columns=AMINO_ACIDS)
assert (BLOSUM62.to_numpy() == BLOSUM62.to_numpy().T).all()


def blosum_of(reference, alternate):
    """BLOSUM62 score per substitution. dbNSFP writes a stop codon as X, not *."""
    pair = pd.MultiIndex.from_arrays(
        [reference.replace("X", "*"), alternate.replace("X", "*")]
    )
    return pd.Series(BLOSUM62.stack().reindex(pair).to_numpy(), index=reference.index)
