"""
Project-wide settings: paths, study vocabulary and annotation conventions.

Every script imports from here, so a path or a label is changed in one place only.
"""

from pathlib import Path

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SOURCE = DATA / "source"  # original files
ANNOTATION = SOURCE / "annotation"  # BEAGLE design, Ensembl VEP output, library sheets
SCREEN = SOURCE / "screen"  # this study: A375 + HT29, CBE, LIN/PIC/SCH
COELHO = SOURCE / "coelho2024"  # independent screen: Coelho et al. 2024 Nat Genet (PMID: 39424923)


# --------------------------------------------------------------------------------------
# Studies.
# --------------------------------------------------------------------------------------
THIS_STUDY = "BECancerResistome"
COELHO_STUDY = "Coelho2024"
STUDIES = [THIS_STUDY, COELHO_STUDY]

# Editors
EDITORS = ["CBE", "ABE"]
SCREEN_EDITORS = {THIS_STUDY: ["CBE"], COELHO_STUDY: ["ABE", "CBE"]}

# Library
LIBRARY_SHEETS = {
    THIS_STUDY: "BE sgRNA library_EG",
    COELHO_STUDY: "BE sgRNA library_MC",
}

def dataset_file(study: str, kind: str, suffix: str = "csv") -> Path:
    """data/<study>_<kind>.<suffix>, so a dataset's files sort together in one block.

    Files that belong to both studies - samplesheet.xlsx, vep_input.vcf and the raw VEP
    and SpliceAI output - carry no study and are not built here.
    """
    return DATA / f"{study}_{kind}.{suffix}"


# --------------------------------------------------------------------------------------
# Nomenclature
# --------------------------------------------------------------------------------------
GUIDE_TYPES = [
    "targeting",
    "splice_essential",  # positive control: splice site of a pan-essential gene
    "splice_nonessential",  # negative control: splice site of a non-essential gene
    "intergenic",  # negative control: an intergenic locus
    "non_targeting",  # negative control: no genomic target
]
CONTROL_TYPES = ["positive", "negative"]
GENELESS_GUIDE_TYPES = ["intergenic", "non_targeting"]

OUR_CONTROL_LABELS = {
    "Intergenic control": ("non_targeting", "negative"),
    "Panlethal splice donor": ("splice_essential", "positive"),
}

COELHO_GUIDE_TYPES = {
    "exonic": ("targeting", None),
    "splice_essential": ("splice_essential", "positive"),
    "splice_nonessential": ("splice_nonessential", "negative"),
    "intergenic": ("intergenic", "negative"),
    "non_targeting": ("non_targeting", "negative"),
}


# --------------------------------------------------------------------------------------
# Annotation conventions
# --------------------------------------------------------------------------------------
# Most-severe-first ranking: Nonsense > Splice site > Missense > Intron > Silent > UTR;
# a guide with no editable base in the window gets "No edits"
SEVERITY = [
    "Nonsense",
    "Splice-donor",
    "Splice-acceptor",
    "Missense",
    "Intron",
    "Silent",
    "UTR",
    "Flank",  # outside the transcript; rare, and least severe of the lot
]
NO_EDITS = "No edits"

# BEAGLE design
SPREADSHEET_ERROR = "#"

CHROMOSOME_ACCESSIONS = {f"NC_{n:06d}": str(n) for n in range(1, 23)} | {
    "NC_000023": "X",
    "NC_000024": "Y",
    "NC_012920": "MT",
}

# Coelho 2024 report a guide's locus as the 23-nt guide+PAM window, BEAGLE as the 20-nt
GUIDE_START_SHIFT = {"left": 3, "right": 0}

# Guide genomic strand
DIRECTION_STRAND = {"left": "-", "right": "+"}

# A base editor rewrites this base anywhere in the guide's edit window.
EDITABLE_BASE = {"CBE": "C", "ABE": "A"}
EDIT_WINDOW = range(4, 9)  # 1-based positions in the 20-nt guide, BEAGLE's "4..8"

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Ensembl VEP writes "-" where a field is empty.
# Every VEP column is prefixed in the output, features can be selected with .filter(like=VEP_PREFIX).
VEP_MISSING = "-"
VEP_PREFIX = "vep_"
