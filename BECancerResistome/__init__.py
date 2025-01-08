#!/usr/bin/env python
# Copyright (C) 2025 Emanuel Goncalves

import pandas as pd


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
