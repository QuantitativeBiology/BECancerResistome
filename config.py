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
COELHO = (
    SOURCE / "coelho2024"
)  # independent screen: Coelho et al. 2024 Nat Genet (PMID: 39424923)
PLOTS = ROOT / "plots"  # figure panels


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
# This study's screen: cell lines, arms and the genes the library targets
# --------------------------------------------------------------------------------------
CELL_LINES = ["A375", "HT29"]
CELL_LINE_LABEL = {"A375": "A375", "HT29": "HT-29"}

ARMS = ["DO", "PIC", "LIN", "SCH"]
DRUG_ARMS = ["LIN", "PIC", "SCH"]  # figure order: Linsitinib, Pictilisib, SCH772984
ARM_LABEL = {
    "DO": "Untreated",
    "PIC": "Pictilisib",
    "LIN": "Linsitinib",
    "SCH": "SCH772984",
}

# Gene order and grouping follow
GENE_GROUPS = {
    "mTOR / PI3K\npathway": [
        "IGF1R",
        "EGFR",
        "ERBB2",
        "ERBB3",
        "FGFR2",
        "GAB1",
        "SHC1",
        "GRB2",
        "PIK3CA",
        "PIK3CB",
        "PIK3CD",
        "PTEN",
        "AKT1",
        "MTOR",
        "RICTOR",
    ],
    "MAPK/ERK\npathway": [
        "NF1",
        "KRAS",
        "NRAS",
        "HRAS",
        "BRAF",
        "RAF1",
        "MAP2K1",
        "MAP2K2",
        "MAPK1",
        "MAPK3",
        "MYC",
    ],
    "Miscellaneous\ngenes": ["PGD", "PKM", "PHLDA2", "ITGA3"],
}

# --------------------------------------------------------------------------------------
# External resources: cell line identifiers, drugs and their targets
# --------------------------------------------------------------------------------------
DEPMAP = DATA / "depmap26Q3"  # expression, copy number, CRISPR gene effect
STRING = DATA / "string"  # protein links v12.5, human
GDSC = DATA / "gdsc"  # fitted dose response, 24Jul22
PROTEOMICS = DATA / "proteomics"  # ProCan, indexed by Sanger model id
PROTVAR = DATA / "protvar"  # FoldX scan, ProtVar output, predicted complexes

DEPMAP_MODEL = {
    "A375": "ACH-000219",
    "HT29": "ACH-000552",
    "H23": "ACH-000900",
    "MHHES1": "ACH-000391",
    "PC9": "ACH-000779",
}
SANGER_MODEL = {
    "A375": "SIDM00795",
    "HT29": "SIDM00136",
    "H23": "SIDM00138",
    "MHHES1": "SIDM00386",
    "PC9": "SIDM00236",
}

# The genes each compound acts on, used for the network distance feature.
DRUG_TARGETS = {
    "pictilisib": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
    "linsitinib": ["IGF1R", "INSR"],
    "SCH772984": ["MAPK1", "MAPK3"],
    "sotorasib": ["KRAS"],
    "adagrasib": ["KRAS"],
    "trametinib": ["MAP2K1", "MAP2K2"],
    "dabrafenib + cetuximab": ["BRAF", "EGFR"],
    "gefitinib": ["EGFR"],
    "osimertinib": ["EGFR"],
    "olaparib": ["PARP1", "PARP2"],
    "niraparib": ["PARP1", "PARP2"],
}

# GDSC names per compound.
GDSC_DRUG_NAMES = {
    "pictilisib": ["Pictilisib"],
    "linsitinib": ["Linsitinib"],
    "SCH772984": ["SCH772984"],
    "sotorasib": [],
    "adagrasib": [],
    "trametinib": ["Trametinib"],
    "dabrafenib + cetuximab": ["Dabrafenib", "Cetuximab"],
    "gefitinib": ["Gefitinib"],
    "osimertinib": ["Osimertinib"],
    "olaparib": ["Olaparib"],
    "niraparib": ["Niraparib"],
}

STRING_CONFIDENCE = 900  # combined_score threshold, 0.9 on STRING's scale


GENE_PATHWAY = {gene: group for group, genes in GENE_GROUPS.items() for gene in genes}

GUIDE_CLASS_LABEL = {
    "non_targeting": "Non-targeting",
    "targeting": "Targeting",
    "splice_essential": "Essential\nsplice site",
}


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
VEP_MISSING = "-"
VEP_PREFIX = "vep_"
