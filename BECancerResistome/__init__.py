#!/usr/bin/env python
# Copyright (C) 2025 Emanuel Goncalves

import numpy as np
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


def convert_amino_acid_variants(variants):
    """
    Convert amino acid variants from full names to single-letter codes.
    Args:
        variants (str): A string of amino acid variants, e.g., "His307His, Ser308Ser, Pro309Ter".
    Returns:
        str: A string of converted amino acid variants, e.g., "H307H, S308S, P309*".

    # Example usage
    example_variants = "His307His, Ser308Ser, Pro309Ter, (NC)"
    converted_variants = convert_amino_acid_variants(example_variants)
    print(converted_variants)  # Output: "H307H, S308S, P309STOP"

    # Note: The function will return NaN if the input is empty or NaN.
    # The function will also filter out any variants that contain "(NC)".
    # This is useful for handling exceptions in the input data.
    """

    if not variants or pd.isna(variants):
        return np.nan

    # Dictionary to map full amino acid names to single-letter codes
    amino_acid_map = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Ter": "*",  # Added support for "Ter"
    }

    def convert_variant(variant):
        # Extract the amino acids and the position
        start_aa = amino_acid_map[variant[:3]]
        position = variant[3:-3]
        end_aa = amino_acid_map[variant[-3:]]
        return f"{start_aa}{position}{end_aa}"

    # Filter out exceptions like "(NC)"
    valid_variants = [v.strip() for v in variants.split(",") if "(NC)" not in v]
    return ", ".join(convert_variant(v) for v in valid_variants)
