#!/usr/bin/env python
# Copyright (C) 2025 Emanuel Goncalves

import sys
from pathlib import Path

import matplotlib as mpl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

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
CELL_LINE_COLOUR = {"A375": "#2c6fa8", "HT29": "#c0392b"}
# fmt: off
ARM_COLOUR = {
    ("A375", "DO"): "#1f4e76", ("A375", "LIN"): "#2c6fa8",
    ("A375", "PIC"): "#5f92bd", ("A375", "SCH"): "#89aece",
    ("HT29", "DO"): "#86281e", ("HT29", "LIN"): "#c0392b",
    ("HT29", "PIC"): "#cf695e", ("HT29", "SCH"): "#dc9088",
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


def save_plot(fig, name: str, dpi: int = PLOT_DPI) -> Path:
    """Write a panel to plots/<name>.png, which is not committed."""
    C.PLOTS.mkdir(parents=True, exist_ok=True)
    path = C.PLOTS / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    return path
