"""
Predictive model input features beyond the screen itself.
"""

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

THREE_TO_ONE = {
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
    "Ter": "*",
}


def gene_columns(path, genes):
    """DepMap writes columns as "SYMBOL (entrez)", so pick the ones we need by symbol."""
    header = pd.read_csv(path, nrows=0).columns
    wanted = {}
    for column in header:
        symbol = column.split(" (")[0]
        if symbol in genes:
            wanted[column] = symbol
    return wanted


def depmap_matrix(path, genes, index_column=None):
    """One row per DepMap model, one column per gene symbol."""
    columns = gene_columns(path, genes)
    identifiers = [index_column] if index_column else []
    if index_column is None:
        frame = pd.read_csv(path, usecols=[0, *range(1, 1)], nrows=0)
        first = pd.read_csv(path, usecols=[0]).columns[0]
        identifiers = [first]
    extra = ["IsDefaultEntryForModel"] if index_column == "ModelID" else []
    frame = pd.read_csv(path, usecols=identifiers + extra + list(columns))
    if extra:
        frame = frame[frame.IsDefaultEntryForModel == "Yes"].drop(columns=extra)
    frame = frame.rename(columns={identifiers[0]: "ModelID", **columns})
    return frame.set_index("ModelID")


# ======================================================================================
# The gene and cell line universe
# ======================================================================================
library = pd.concat(
    [
        pd.read_csv(C.dataset_file(study, "library"), usecols=["gene", "guide_type"])
        for study in C.STUDIES
    ]
)
GENES = sorted(
    {
        gene
        for value in library.loc[library.guide_type == "targeting", "gene"].dropna()
        for gene in str(value).split(";")
    }
)
print(f"{len(GENES)} targeted genes, {len(C.DEPMAP_MODEL)} cell lines")


# ======================================================================================
# A. Cell context
# ======================================================================================
expression = depmap_matrix(
    C.DEPMAP / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", GENES, "ModelID"
)
copy_number = depmap_matrix(C.DEPMAP / "OmicsCNGeneWGS.csv", GENES, "ModelID")
crispr = depmap_matrix(C.DEPMAP / "CRISPRGeneEffect.csv", GENES)

proteomics = pd.read_csv(C.PROTEOMICS / "procan_proteomics.csv")
proteomics = proteomics[proteomics.GeneSymbol.isin(GENES)].set_index("GeneSymbol")

rows = []
for cell_line, model in C.DEPMAP_MODEL.items():
    sidm = C.SANGER_MODEL[cell_line]
    for gene in GENES:
        rows.append(
            {
                "cell_line": cell_line,
                "gene": gene,
                "expression_tpm_log2": (
                    expression.at[model, gene]
                    if model in expression.index and gene in expression.columns
                    else np.nan
                ),
                "copy_number": (
                    copy_number.at[model, gene]
                    if model in copy_number.index and gene in copy_number.columns
                    else np.nan
                ),
                "crispr_gene_effect": (
                    crispr.at[model, gene]
                    if model in crispr.index and gene in crispr.columns
                    else np.nan
                ),
                "protein_abundance": (
                    proteomics.at[gene, sidm]
                    if gene in proteomics.index and sidm in proteomics.columns
                    else np.nan
                ),
            }
        )

cell_context = pd.DataFrame(rows)
cell_context.to_csv(C.DATA / "features_cell_context.csv", index=False)

measurements = [
    "expression_tpm_log2",
    "copy_number",
    "crispr_gene_effect",
    "protein_abundance",
]
print("\nA. cell context, fraction of genes measured per cell line")
print(
    cell_context.groupby("cell_line")[measurements]
    .apply(lambda block: block.notna().mean().round(2))
    .to_string()
)


# ======================================================================================
# B. Drug context: baseline ln IC50
# ======================================================================================
gdsc = pd.read_csv(
    C.GDSC / "GDSC_fitted_dose_response_24Jul22_IC50.csv", low_memory=False
)
parts = gdsc.DRUG_INDEX.str.split(";", expand=True)
gdsc = gdsc.assign(drug_name=parts[1], release=parts[2])

samples = pd.read_excel(C.DATA / "samplesheet.xlsx", sheet_name="Samples")
arms = samples[samples.arm_kind == "drug"][
    ["study", "cell_line", "arm", "compound"]
].drop_duplicates()


def baseline_ic50(compound, cell_line):
    """GDSC2 where it exists, averaged over repeats, otherwise GDSC1."""
    column = C.SANGER_MODEL[cell_line]
    if column not in gdsc.columns:
        return np.nan, "cell line absent from GDSC"
    values, used = [], []
    for name in C.GDSC_DRUG_NAMES.get(compound, []):
        block = gdsc[gdsc.drug_name.str.casefold() == name.casefold()]
        preferred = block[block.release == "GDSC2"]
        source = "GDSC2"
        if preferred[column].notna().sum() == 0:
            preferred, source = block[block.release == "GDSC1"], "GDSC1"
        measured = preferred[column].dropna()
        if len(measured):
            values.append(measured.mean())
            used.append(f"{name}:{source}:{len(measured)}")
    if not values:
        return np.nan, "compound absent from GDSC"
    return float(np.mean(values)), "; ".join(used)


drug_context = arms.assign(
    **dict(
        zip(
            ["ln_ic50", "ic50_source"],
            zip(
                *[
                    baseline_ic50(row.compound, row.cell_line)
                    for row in arms.itertuples()
                ]
            ),
        )
    )
)
drug_context.to_csv(C.DATA / "features_drug_context.csv", index=False)
print("\nB. baseline ln IC50")
print(
    drug_context[["cell_line", "compound", "ln_ic50", "ic50_source"]].to_string(
        index=False
    )
)


# ======================================================================================
# C. Gene to drug distance in STRING
# ======================================================================================
info = pd.read_csv(C.STRING / "9606.protein.info.v12.5.txt.gz", sep="\t")
symbol = dict(zip(info["#string_protein_id"], info.preferred_name))

links = pd.read_csv(
    C.STRING / "9606.protein.links.v12.5.txt.gz",
    sep=" ",
    usecols=["protein1", "protein2", "combined_score"],
)
links = links[links.combined_score >= C.STRING_CONFIDENCE]
print(
    f"\nC. STRING v12.5: {len(links):,} edges at combined_score >= {C.STRING_CONFIDENCE}"
)

nodes = sorted(set(links.protein1) | set(links.protein2))
index_of = {node: position for position, node in enumerate(nodes)}
rows_index = links.protein1.map(index_of).to_numpy()
cols_index = links.protein2.map(index_of).to_numpy()
graph = coo_matrix(
    (np.ones(len(links)), (rows_index, cols_index)), shape=(len(nodes), len(nodes))
).tocsr()

symbol_nodes = {}
for node in nodes:
    symbol_nodes.setdefault(symbol.get(node), []).append(index_of[node])

distances = []
for compound, targets in C.DRUG_TARGETS.items():
    seeds = [i for target in targets for i in symbol_nodes.get(target, [])]
    if not seeds:
        print(f"  {compound}: no target in STRING, skipped")
        continue
    hops = dijkstra(graph, directed=False, unweighted=True, indices=seeds).min(axis=0)
    for gene in GENES:
        positions = symbol_nodes.get(gene, [])
        reach = min((hops[i] for i in positions), default=np.inf)
        distances.append(
            {
                "compound": compound,
                "gene": gene,
                "is_drug_target": gene in targets,
                "string_distance": None if np.isinf(reach) else int(reach),
            }
        )

gene_drug = pd.DataFrame(distances)
gene_drug.to_csv(C.DATA / "features_gene_drug.csv", index=False)
print(
    f"  {len(gene_drug):,} compound by gene pairs, "
    f"{gene_drug.string_distance.isna().sum()} unreachable"
)
print(gene_drug.string_distance.value_counts(dropna=False).sort_index().to_string())


# ======================================================================================
# D. Variant structure: FoldX, pockets, interfaces
# ======================================================================================
annotation = pd.concat(
    [
        pd.read_csv(
            C.dataset_file(study, "annotation"),
            low_memory=False,
            usecols=[
                "gene",
                "vep_SWISSPROT",
                "vep_Protein_position",
                "vep_Amino_acids",
                "chromosome",
                "position",
                "genomic_ref",
                "genomic_alt",
            ],
        )
        for study in C.STUDIES
    ]
)
changes = annotation[annotation.vep_Amino_acids.astype(str).str.contains("/", na=False)]
pair = changes.vep_Amino_acids.str.split("/", n=1, expand=True)
variants = pd.DataFrame(
    {
        "gene": changes.gene.to_numpy(),
        "uniprot": changes.vep_SWISSPROT.astype(str).str.split(".").str[0].to_numpy(),
        "residue": pd.to_numeric(
            changes.vep_Protein_position.astype(str).str.split("-").str[0],
            errors="coerce",
        ).astype("Int64"),
        "wild_type": pair[0].to_numpy(),
        "mutant": pair[1].to_numpy(),
        "chromosome": changes.chromosome.to_numpy(),
        "position": changes.position.to_numpy(),
        "genomic_ref": changes.genomic_ref.to_numpy(),
        "genomic_alt": changes.genomic_alt.to_numpy(),
    }
).drop_duplicates(["gene", "uniprot", "residue", "wild_type", "mutant"])
print(f"\nD. {len(variants):,} distinct protein changes over both screens")

foldx = pd.read_csv(
    C.PROTVAR / "foldx_filtered_proteins_of_interest.csv",
    usecols=[
        "uniprot_accession",
        "uniprot_position",
        "wild_type",
        "mutated_type",
        "foldx_ddg",
        "plddt",
    ],
).rename(
    columns={
        "uniprot_accession": "uniprot",
        "uniprot_position": "residue",
        "mutated_type": "mutant",
    }
)
foldx["residue"] = foldx.residue.astype("Int64")
structure = variants.merge(
    foldx, on=["uniprot", "residue", "wild_type", "mutant"], how="left"
)

# ProtVar fills what the FoldX scan misses and is the only source here for pockets.
# Every ProtVar download in data/protvar is read, so adding one needs no code change.
FOLDX_FIELD = "Foldx_prediction(foldxDdg;plddt)"
GENOMIC_BLOCK = (
    6  # the trailing fields, from Genomic_location to Position_in_structures
)


def read_protvar(path):
    """ProtVar repeats the trailing genomic block on variants that map to more than one
    position, so those rows are wider than the header. Their prefix and their last six
    fields still line up, so the repeat is dropped instead of the row."""
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        width = len(header)
        rows, repaired = [], 0
        for row in reader:
            if len(row) > width:
                row = row[: width - GENOMIC_BLOCK] + row[-GENOMIC_BLOCK:]
                repaired += 1
            if len(row) == width:
                rows.append(row)
    if repaired:
        print(
            f"    {path.name}: repaired {repaired} rows with a repeated genomic block"
        )
    return pd.DataFrame(rows, columns=header)


downloads = []
for path in sorted(C.PROTVAR.glob("*.csv")):
    columns = pd.read_csv(path, nrows=0).columns
    if FOLDX_FIELD in columns:
        downloads.append(read_protvar(path))
        print(f"  ProtVar download: {path.name} ({len(downloads[-1]):,} rows)")
assert downloads, f"no ProtVar output found in {C.PROTVAR}"
protvar = pd.concat(downloads, ignore_index=True)

# ProtVar falls back to the UniProt canonical isoform when the reference residue we ask
# about is not the one it holds at that position. Those rows describe a different
# residue, so they are discarded rather than joined.
mismatched = protvar.Notes.astype(str).str.contains(
    "does not match the UniProt sequence", na=False
)
if mismatched.any():
    affected = protvar.loc[mismatched, "Gene"].value_counts().to_dict()
    print(
        f"  discarding {int(mismatched.sum()):,} rows on a different isoform: {affected}"
    )
protvar = protvar[~mismatched]
protvar_pair = protvar.Amino_acid_change.astype(str).str.split("/", expand=True)
protvar = protvar.assign(
    gene=protvar.Gene,
    residue=pd.to_numeric(protvar.Amino_acid_position, errors="coerce").astype("Int64"),
    wild_type=protvar_pair[0].map(THREE_TO_ONE),
    mutant=protvar_pair[1].map(THREE_TO_ONE),
)
foldx_field = protvar[FOLDX_FIELD].astype(str)
protvar["protvar_ddg"] = pd.to_numeric(
    foldx_field.str.extract(r"foldxDdg:([-\d.]+)")[0], errors="coerce"
)
pocket_field = protvar["Predicted_pockets(energy;per_vol;score;resids)"].astype(str)
protvar["pocket_buriedness"] = pd.to_numeric(
    pocket_field.str.extract(r"buriedness:([\d.]+)")[0], errors="coerce"
)
protvar["pocket_score"] = pd.to_numeric(
    pocket_field.str.extract(r"score:([\d.]+)")[0], errors="coerce"
)
protvar["in_predicted_pocket"] = protvar.pocket_buriedness.notna()
protvar["protvar_annotated"] = True

supplement = protvar[
    [
        "gene",
        "residue",
        "wild_type",
        "mutant",
        "protvar_ddg",
        "pocket_buriedness",
        "pocket_score",
        "in_predicted_pocket",
        "protvar_annotated",
    ]
].drop_duplicates(["gene", "residue", "wild_type", "mutant"])
structure = structure.merge(
    supplement, on=["gene", "residue", "wild_type", "mutant"], how="left"
)
structure["foldx_ddg"] = structure.foldx_ddg.fillna(structure.protvar_ddg)
structure = structure.drop(columns="protvar_ddg")
structure["protvar_annotated"] = structure.protvar_annotated.fillna(False).astype(bool)
structure["in_predicted_pocket"] = structure.in_predicted_pocket.astype("boolean")

# Interface residues from the predicted complexes
models = pd.read_csv(
    C.PROTVAR / "suppl_ppi_models_.tsv",
    sep="\t",
    usecols=["pdockq", "uniprot_id1", "uniprot_id2", "ifresid1", "ifresid2"],
    low_memory=False,
)
accessions = set(variants.uniprot)
models = models[
    models.uniprot_id1.isin(accessions) | models.uniprot_id2.isin(accessions)
]
interface = {}
for row in models.itertuples():
    for accession, residues in (
        (row.uniprot_id1, row.ifresid1),
        (row.uniprot_id2, row.ifresid2),
    ):
        if accession not in accessions or not isinstance(residues, str):
            continue
        for token in residues.split(","):
            # residues are written as "R169", the one letter code then the position
            number = token.strip().lstrip("ABCDEFGHIKLMNPQRSTVWXYZ")
            if number.isdigit():
                key = (accession, int(number))
                counts, best = interface.get(key, (0, 0.0))
                interface[key] = (counts + 1, max(best, row.pdockq))

structure["interface_count"] = [
    interface.get((row.uniprot, row.residue), (0, 0.0))[0]
    for row in structure.itertuples()
]
structure["interface_max_pdockq"] = [
    interface.get((row.uniprot, row.residue), (0, np.nan))[1]
    for row in structure.itertuples()
]
structure["is_in_interface"] = structure.interface_count > 0

structure.to_csv(C.DATA / "features_variant_structure.csv", index=False)

# ProtVar is the only source here for pockets, and the file on disk covers the handful
# of changes an earlier query asked about. This writes the full list in ProtVar's input
# format so one submission can close the gap.
pending = structure[(structure.mutant != "*") & ~structure.protvar_annotated]
# Protein level input fails wherever UniProt numbers the canonical isoform differently
# from the MANE transcript, so the resubmission is written in genomic coordinates, which
# carry no isoform assumption. The chromosome accessions are taken from ProtVar's own
# output rather than hard coded.
located = protvar.assign(
    accession=protvar.Genomic_location.astype(str).str.split(":").str[0]
)
located = located[located.accession.str.startswith("NC_")]
accession_of = (
    located.groupby("Chromosome")
    .accession.agg(lambda values: values.mode().iat[0])
    .to_dict()
)
pending_lines = []
for row in pending.itertuples():
    accession = accession_of.get(str(row.chromosome).removesuffix(".0"))
    if accession and pd.notna(row.position):
        pending_lines.append(
            f"{accession}:g.{int(row.position)}{row.genomic_ref}>{row.genomic_alt}"
        )
(C.PROTVAR / "protvar_input_all_missense.txt").write_text(
    "\n".join(pending_lines) + "\n"
)
print(
    f"  protvar_input_all_missense.txt: {len(pending_lines):,} not returned by ProtVar"
)
missense = structure[structure.mutant != "*"]
print(
    f"  FoldX ddG: {missense.foldx_ddg.notna().mean():.1%} of {len(missense):,} missense"
    f" changes\n"
    f"  in a predicted pocket: {missense.in_predicted_pocket.fillna(False).mean():.1%}"
    f" of the {missense.protvar_annotated.mean():.1%} ProtVar annotated\n"
    f"  in a predicted interface: {missense.is_in_interface.mean():.1%}"
)

print("\nwritten:")
for name in [
    "features_cell_context",
    "features_drug_context",
    "features_gene_drug",
    "features_variant_structure",
]:
    path = C.DATA / f"{name}.csv"
    print(f"  {path.name:34s} {len(pd.read_csv(path)):>7,} rows")
