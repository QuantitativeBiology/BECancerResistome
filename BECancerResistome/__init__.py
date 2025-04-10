#!/usr/bin/env python
# Copyright (C) 2025 Emanuel Goncalves

import pandas as pd
import matplotlib.pyplot as plt


# Matplotlib config
plt.rcParams["figure.figsize"] = [4, 4]
plt.rcParams["figure.dpi"] = 300

# Matplotlib set font to sans-times
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"

# Matplotlib set main axis font size
plt.rcParams["axes.titlesize"] = 7

# Matplotlib set legend font size
plt.rcParams["legend.fontsize"] = 6

# Matplotlib set legend title font size
plt.rcParams["legend.title_fontsize"] = 6

# Matplotlib set tick label font size
plt.rcParams["axes.labelsize"] = 6

# Matplotlib set tick label font size
plt.rcParams["xtick.labelsize"] = 6
plt.rcParams["ytick.labelsize"] = 6

# Matplotlib set grid line width
plt.rcParams["grid.linewidth"] = 0.5

# Matplotlib ommit top and right spines
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

# Matplotlib set grid line
plt.rcParams["axes.grid"] = True

# Matplotlib set grid line style
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams["grid.linewidth"] = 0.15

# Matplotlib set grid line color
plt.rcParams["grid.color"] = "black"

# Matplotlib set grid line alpha
plt.rcParams["grid.alpha"] = 0.5

# Matplotlib set legend frameon
plt.rcParams["legend.frameon"] = False

# Matplotlib set legend loc
plt.rcParams["legend.loc"] = "best"

# Matplotlib set axis below true
plt.rcParams["axes.axisbelow"] = True

# Matplotlib illustrator export
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def beagle2vep(r: pd.Series) -> list[str]:
    """
    Create a list of ENSEMBL VEP HGVS notations from a row of the Beagle output.

    Args:
    r: pd.Series: A row of the Beagle output. Must contain the columns "Nucleotide Edits" and "Target Transcript ID".

    Returns:
    list[str]: A list of ENSEMBL VEP HGVS notations.

    """

    r_edits = r["Nucleotide Edits"].replace(",", ";").split(";")

    return [f"{r['Target Transcript ID']}:c.{v.strip()}" for v in r_edits]
