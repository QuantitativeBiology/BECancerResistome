"""
sgRNA remapping and VEP annotation. Every guide is annotated from a BEAGLE run 
and a single Ensembl VEP run, on one transcript per gene.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

DESIGN_COLUMNS = {
    "sgRNA Sequence": "guide",
    "Edit Type": "editor",
    "Edit Window": "edit_window",
    "Target Genome Sequence": "accession",
    "Target Gene Symbol": "gene",
    "Target Gene ID": "ensembl_gene",
    "Target Gene Strand": "gene_strand",
    "Target Transcript ID": "transcript_full",
    "sgRNA Context Sequence": "context_sequence",
    "sgRNA Sequence Start Pos. (global)": "target_start",
    "sgRNA Orientation": "guide_orientation",
    "Constraint Violations": "constraint_violations",
}
EDIT_COLUMNS = {
    "Nucleotide Edits (global)": "genomic_edits",
    "Nucleotide Edits": "cdna_edits",
    "Amino Acid Edits": "amino_acid_edits",
    "Mutation Category": "mutation_categories",
    "Guide Edits": "guide_edits",
}
LOCUS_COLUMNS = [
    "gene",
    "ensembl_gene",
    "gene_strand",
    "accession",
    "target_start",
    "guide_orientation",
    "context_sequence",
]
GENOMIC_EDIT = re.compile(r"^(?P<position>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
OT_COLUMNS = ["OT_0mm", "OT_1mm", "OT_2mm", "OT_3mm", "OT_4mm"]

LIBRARY_SCHEMA = [
    "guide",
    "study",
    "guide_type",
    "control_type",
    "gene",
    "n_genes",
    "ensembl_gene",
    "transcript",
    "transcript_version",
    "chromosome",
    "gene_strand",
    "genomic_strand",
    "target_start",
    "guide_orientation",
    "context_sequence",
    "edit_window",
    "constraint_violations",
    "editors_designed",
    "annotation_source",
    # Both editors appear in both files so the two share a schema; 
    # the columns of an editor a study never used stay empty.
    *[
        f"{field}_{editor.lower()}"
        for editor in C.EDITORS
        for field in ["n_variants", "mutation_category"]
    ],
    *[c.lower() for c in OT_COLUMNS],
    "ot_overflow",
]
VARIANT_SCHEMA = [
    "guide",
    "study",
    "editor",
    "variant",
    "gene",
    "transcript",
    "chromosome",
    "position",
    "genomic_ref",
    "genomic_alt",
    "cdna_edit",
    "codon_amino_acid_edit",
    "codon_mutation_category",
    "guide_edit",
    "codon_group",
]


def groups(text) -> list[str]:
    """'a, b' -> ['a', 'b']; a missing value -> []."""
    if pd.isna(text) or str(text).strip() == "":
        return []
    return [group.strip() for group in str(text).split(",")]


def bases(group: str) -> list[str]:
    """'547A>G;549A>G' -> the bases edited in one codon."""
    return [base.strip() for base in group.split(";") if base.strip()]


def chromosome_of(accession: str) -> str:
    """'NC_000012.12' -> '12'."""
    return C.CHROMOSOME_ACCESSIONS[accession.split(".")[0]]


def most_severe(categories: list[str]) -> str:
    return next((c for c in C.SEVERITY if c in set(categories)), C.NO_EDITS)


def join_unique(values: pd.Series) -> str:
    """Distinct values of a guide, for fields that do not belong to a single gene."""
    return ";".join(sorted({str(v) for v in values.dropna()}))


def join_per_gene(values: pd.Series) -> str:
    """Values in the order they are given, which the caller sorts by gene.

    Used for every field that belongs to one of a paralogue guide's genes, so that
    position i of gene, ensembl_gene, gene_strand and transcript all mean the same gene.
    """
    return ";".join(str(v) for v in values)


def split_version(transcripts: pd.Series) -> tuple[pd.Series, pd.Series]:
    """'ENST00000275493.7' -> ('ENST00000275493', '7'), semicolons kept aligned."""
    parts = transcripts.fillna("").str.split(";")
    bare = parts.map(lambda ids: ";".join(i.split(".")[0] for i in ids))
    version = parts.map(
        lambda ids: ";".join(i.split(".")[1] if "." in i else "" for i in ids)
    )
    return bare.replace("", pd.NA), version.replace("", pd.NA)


# ======================================================================================
# A. BEAGLE design -> one row per (guide, editor, edited base)
# ======================================================================================
design = pd.read_csv(C.ANNOTATION / "beagle_design.csv", low_memory=False).rename(
    columns=DESIGN_COLUMNS | EDIT_COLUMNS
)
print(
    f"BEAGLE design: {len(design):,} rows, {design.guide.nunique():,} guides, "
    f"editors {sorted(design.editor.unique())}"
)

counted = {c: design[c].map(lambda t: len(groups(t))) for c in EDIT_COLUMNS.values()}
assert pd.DataFrame(counted).nunique(axis=1).eq(1).all(), (
    "the edit columns of a design row disagree on how many codon groups it has"
)
seen_categories = set(
    design.mutation_categories.dropna().str.split(",").explode().str.strip()
)
assert seen_categories <= set(C.SEVERITY), (
    f"BEAGLE categories missing from config.SEVERITY, so they would rank as "
    f"{C.NO_EDITS!r}: {sorted(seen_categories - set(C.SEVERITY))}"
)

# One transcript per gene, applied to every guide of both libraries further down.
transcript_of = design.groupby("gene").transcript_full.agg(["nunique", "first"])
assert transcript_of["nunique"].eq(1).all(), (
    f"BEAGLE used more than one transcript for "
    f"{transcript_of.index[transcript_of['nunique'] > 1].tolist()}"
)
transcript_of = transcript_of["first"]
print(f"  one transcript per gene for {len(transcript_of)} genes")

rows, unusable = [], []
for record in design.to_dict("records"):
    genomic_groups = groups(record["genomic_edits"])
    cdna_groups = groups(record["cdna_edits"])
    aa_entries = groups(record["amino_acid_edits"])
    categories = groups(record["mutation_categories"])
    guide_edit_groups = groups(record["guide_edits"])

    for index, (genomic_group, cdna_group) in enumerate(
        zip(genomic_groups, cdna_groups)
    ):
        genomic_bases = bases(genomic_group)
        cdna_bases = bases(cdna_group)
        # "C_8;C_7" names the bases of this codon in the same order as its coordinates,
        # so each edited base gets its own notation rather than the whole group's.
        guide_bases = bases(guide_edit_groups[index])
        assert len(genomic_bases) == len(cdna_bases), (
            f"{record['guide']}: {len(genomic_bases)} genomic and {len(cdna_bases)} "
            "transcript coordinates in one codon"
        )
        for offset, (genomic, cdna) in enumerate(zip(genomic_bases, cdna_bases)):
            match = GENOMIC_EDIT.match(genomic)
            assert match is not None, f"{record['guide']}: cannot read edit {genomic!r}"
            if C.SPREADSHEET_ERROR in cdna:
                unusable.append(
                    (record["guide"], record["editor"], record["gene"], genomic, cdna)
                )
                continue
            rows.append(
                {
                    "guide": record["guide"],
                    "editor": record["editor"],
                    "variant": f"{record['transcript_full']}:c.{cdna}",
                    "gene": record["gene"],
                    "transcript": record["transcript_full"],
                    "chromosome": chromosome_of(record["accession"]),
                    "position": int(match["position"]),
                    "genomic_ref": match["ref"],
                    "genomic_alt": match["alt"],
                    "cdna_edit": cdna,
                    "codon_amino_acid_edit": aa_entries[index],
                    "codon_mutation_category": categories[index],
                    "guide_edit": guide_bases[offset]
                    if offset < len(guide_bases)
                    else guide_edit_groups[index],
                    "codon_group": index + 1,
                }
            )

variants = pd.DataFrame(rows)
print(
    f"  {len(variants):,} (guide, editor, edited base) rows, "
    f"{variants.variant.nunique():,} distinct variants; "
    f"{int(design.cdna_edits.isna().sum()):,} design rows have no editable base"
)
if unusable:
    print(f"  dropped {len(unusable)} edits whose transcript coordinate is unreadable:")
    for guide, editor, gene, genomic, cdna in unusable:
        print(f"    {guide} {editor} {gene} genomic {genomic}, transcript {cdna!r}")


# ======================================================================================
# B. Locus per guide. Gene and coordinates describe where the guide binds, so they do
# not depend on the editor. Checked per (guide, gene), which covers the paralogue guides 
# too: a guide that matches two genes has two loci.
# ======================================================================================
per_locus = design.groupby(["guide", "gene"])[
    [c for c in LOCUS_COLUMNS if c != "gene"] + ["edit_window", "constraint_violations"]
].nunique()
conflicting = per_locus.gt(1).any(axis=1)
assert not conflicting.any(), (
    "a locus is annotated differently by the two editors' design rows: "
    f"{per_locus.index[conflicting].tolist()[:3]}"
)
multi_gene = design.groupby("guide").gene.nunique() > 1
paralogue_genes = sorted(
    set(design.loc[design.guide.isin(multi_gene.index[multi_gene]), "gene"])
)
print(f"\nguides matching two paralogues: {int(multi_gene.sum())} ({paralogue_genes})")

# Sorted by gene, so join_per_gene keeps the per-gene fields aligned with each other.
# accession, target_start and context_sequence are a property of one locus, and a
# paralogue guide has two, so they describe the first gene listed.
by_gene = design.drop_duplicates(["guide", "gene"]).sort_values(["guide", "gene"])
locus = by_gene.groupby("guide").agg(
    gene=("gene", join_per_gene),
    n_genes=("gene", "nunique"),
    ensembl_gene=("ensembl_gene", join_per_gene),
    gene_strand=("gene_strand", join_per_gene),
    guide_orientation=("guide_orientation", join_per_gene),
    accession=("accession", "first"),
    target_start=("target_start", "first"),
    context_sequence=("context_sequence", "first"),
    edit_window=("edit_window", "first"),
    constraint_violations=("constraint_violations", join_unique),
)
locus["chromosome"] = locus.accession.map(chromosome_of)
# Genomic strand: a sense guide on a + strand gene reads the
# forward strand, and either one flipped reads the reverse.
locus["genomic_strand"] = np.where(
    (locus.gene_strand.str[0] == "+") == (locus.guide_orientation.str[:5] == "sense"),
    "+",
    "-",
)


# ======================================================================================
# C. Controls. Guide with no design row is a control.
# ======================================================================================
sheet = pd.read_excel(
    C.ANNOTATION / "be_library_samplesheet.xlsx",
    sheet_name=C.LIBRARY_SHEETS[C.THIS_STUDY],
)
our_labels = sheet.drop_duplicates("sgRNA sequence").set_index("sgRNA sequence")[
    "GeneID"
]
our_controls = our_labels[our_labels.isin(C.OUR_CONTROL_LABELS)]
print(f"\nthis study's control labels: {our_controls.value_counts().to_dict()}")

coelho = pd.read_csv(
    C.COELHO / "readcounts.csv",
    low_memory=False,
    usecols=[
        "guide",
        "sgRNA_type",
        "Gene",
        "Gene_ID",
        "chr_BEstimate",
        "start_BEstimate",
        "Direction",
    ],
).set_index("guide")
coelho_types = coelho["sgRNA_type"].map(C.COELHO_GUIDE_TYPES)
print(
    f"Coelho's own labels: "
    f"{coelho_types.map(lambda pair: pair[0]).value_counts().to_dict()}"
)

assert coelho["Direction"].dropna().isin(C.GUIDE_START_SHIFT).all(), (
    "Coelho report a guide Direction the coordinate shift does not cover: "
    f"{sorted(set(coelho['Direction'].dropna()) - set(C.GUIDE_START_SHIFT))}"
)
coelho_start = coelho["start_BEstimate"] + coelho["Direction"].map(C.GUIDE_START_SHIFT)
coelho_locus = pd.DataFrame(
    {
        "gene": coelho["Gene"],
        "ensembl_gene": coelho["Gene_ID"],
        "chromosome": coelho["chr_BEstimate"].astype("Int64").astype(str),
        "target_start": coelho_start.astype("Int64"),
        "genomic_strand": coelho["Direction"].map(C.DIRECTION_STRAND),
    }
).replace("<NA>", pd.NA)
COELHO_FALLBACK_COLUMNS = list(coelho_locus.columns)

both = sorted(set(design.guide) & set(coelho_locus.index[coelho_locus.gene.notna()]))
beagle_genes = design.groupby("guide").gene.agg(set).reindex(both)
agree = [
    gene in genes for gene, genes in zip(coelho_locus.gene.reindex(both), beagle_genes)
]
print(
    f"  Coelho's gene is among BEAGLE's for {sum(agree):,} of {len(both):,} guides "
    "both annotate"
)

shared_loci = design.merge(
    coelho_start.rename("coelho_start"), left_on="guide", right_index=True
).merge(coelho["Gene"], left_on="guide", right_index=True)
shared_loci = shared_loci[
    (shared_loci.gene == shared_loci.Gene) & shared_loci.coelho_start.notna()
]
assert (shared_loci.target_start == shared_loci.coelho_start).all(), (
    "the guide+PAM shift does not convert Coelho's coordinates to BEAGLE's for "
    f"{int((shared_loci.target_start != shared_loci.coelho_start).sum()):,} rows"
)
print(
    f"  their start converts to BEAGLE's for all {len(shared_loci):,} rows both annotate"
)

their_strand = coelho["Direction"].map(C.DIRECTION_STRAND).reindex(shared_loci.guide)
our_strand = np.where(
    (shared_loci.gene_strand == "+") == (shared_loci.guide_orientation == "sense"),
    "+",
    "-",
)
assert (their_strand.to_numpy() == our_strand).all(), (
    "Coelho's Direction and BEAGLE's orientation disagree on the genomic strand for "
    f"{int((their_strand.to_numpy() != our_strand).sum()):,} rows"
)
print("  their Direction agrees with BEAGLE's strand for all of them")


# ======================================================================================
# D. Off-targets (FlashFry, hg38). 'OVERFLOW' is too many sites to count, kept as a flag.
# ======================================================================================
flashfry = pd.read_excel(C.SCREEN / "offtargets_flashfry.xlsx").set_index("Guide")
off_targets = pd.DataFrame(index=flashfry.index)
off_targets["ot_overflow"] = flashfry[OT_COLUMNS].astype(str).eq("OVERFLOW").any(axis=1)
for column in OT_COLUMNS:
    off_targets[column.lower()] = pd.to_numeric(
        flashfry[column], errors="coerce"
    ).astype("Int64")
print(
    f"\nFlashFry off-targets: {len(off_targets):,} guides, "
    f"{int(off_targets.ot_overflow.sum())} with OVERFLOW"
)


# ======================================================================================
# E. Ensembl VEP. A variant is annotated on every transcript it falls in; the row kept is
# the MANE Select transcript of the guide's gene, which is also the transcript BEAGLE
# designed against. Both libraries are annotated on the same one.
# ======================================================================================
batches = sorted(C.ANNOTATION.glob("vep_results_batch_*.txt"))
assert batches, "no VEP result files found"

kept, total_rows = [], 0
for batch in batches:
    for chunk in pd.read_csv(batch, sep="\t", chunksize=200_000, low_memory=False):
        total_rows += len(chunk)
        kept.append(chunk[chunk["MANE_SELECT"].ne(C.VEP_MISSING)])
vep = pd.concat(kept, ignore_index=True).rename(
    columns={
        "#Uploaded_variation": "variant",
        "SYMBOL": "gene",
        "Feature": "transcript",
    }
)
print(
    f"\nVEP output: {total_rows:,} rows in {len(batches)} batches, "
    f"{len(vep):,} on a MANE Select transcript, {vep.variant.nunique():,} variants"
)

assert not vep.duplicated(["variant", "gene"]).any(), (
    "a variant has more than one MANE Select transcript in the same gene"
)
per_gene_mane = vep.groupby("gene").transcript.nunique()
assert per_gene_mane.eq(1).all(), (
    f"more than one MANE Select transcript for "
    f"{per_gene_mane.index[per_gene_mane > 1].tolist()}"
)
mane_of = vep.drop_duplicates("gene").set_index("gene").transcript
shared_genes = transcript_of.index.intersection(mane_of.index)
off_mane = transcript_of[shared_genes] != mane_of[shared_genes]
assert not off_mane.any(), (
    f"BEAGLE did not design against the MANE transcript for "
    f"{transcript_of[shared_genes][off_mane].to_dict()}"
)
print(f"  BEAGLE's transcript is MANE Select for all {len(shared_genes)} genes checked")

feature_columns = [c for c in vep.columns if c not in ["variant", "gene", "transcript"]]
vep[feature_columns] = vep[feature_columns].replace(C.VEP_MISSING, np.nan)
vep = vep.rename(columns={c: f"{C.VEP_PREFIX}{c}" for c in feature_columns})

annotated = variants.merge(
    vep.drop(columns="transcript"), on=["variant", "gene"], how="left", validate="m:1"
)
missing = annotated[f"{C.VEP_PREFIX}Consequence"].isna()
print(
    f"  {len(feature_columns)} VEP columns; {int(missing.sum()):,} of {len(annotated):,} "
    f"variant rows have no VEP row ({annotated.loc[missing, 'variant'].nunique():,} distinct)"
)

# Independent check of the whole parse: VEP reports the position on the genome and the
# reference base in the orientation of the HGVS it was given, so BEAGLE's genomic base
# must equal it on a + strand gene and complement it on a - strand one.
checked = annotated[annotated[f"{C.VEP_PREFIX}STRAND"].notna()]
vep_position = checked[f"{C.VEP_PREFIX}Location"].str.extract(r":(\d+)")[0].astype(int)
assert (vep_position == checked.position).all(), (
    f"{int((vep_position != checked.position).sum()):,} variants sit at a different "
    "position in the VEP output than BEAGLE's genomic coordinate"
)
# STRAND arrives as text, because VEP writes "-" for a feature that has none.
on_plus_strand = pd.to_numeric(checked[f"{C.VEP_PREFIX}STRAND"], errors="coerce").eq(1)
transcript_ref = checked.genomic_ref.where(
    on_plus_strand, checked.genomic_ref.map(C.COMPLEMENT)
)
disagree = checked[f"{C.VEP_PREFIX}REF_ALLELE"] != transcript_ref
assert not disagree.any(), (
    f"{int(disagree.sum()):,} variants disagree with VEP on the reference base"
)
print(f"  position and alleles agree with VEP for all {len(checked):,} annotated rows")


# ======================================================================================
# F. Per-study library and variant tables
# ======================================================================================
for study in C.STUDIES:
    screened = pd.read_csv(
        C.dataset_file(study, "counts"), index_col=0, usecols=[0]
    ).index
    editors = C.SCREEN_EDITORS[study]

    study_variants = annotated[
        annotated.guide.isin(screened) & annotated.editor.isin(editors)
    ].copy()
    study_variants["study"] = study
    study_variants = study_variants.reindex(
        columns=VARIANT_SCHEMA
        + [c for c in annotated.columns if c.startswith(C.VEP_PREFIX)]
    ).sort_values(["guide", "editor", "position"])

    library = pd.DataFrame(index=pd.Index(screened, name="guide"))
    library = library.join(locus.reindex(library.index)).join(
        off_targets.reindex(library.index)
    )
    library["study"] = study
    library["target_start"] = library.target_start.astype("Int64")
    library["annotation_source"] = np.where(
        library.index.isin(design.guide), "beagle", pd.NA
    )

    if study == C.COELHO_STUDY:
        fallback = coelho_locus.reindex(library.index)
        needs_fallback = library.annotation_source.isna() & fallback.gene.notna()
        for column in COELHO_FALLBACK_COLUMNS:
            library[column] = library[column].where(~needs_fallback, fallback[column])
        library["annotation_source"] = library.annotation_source.where(
            ~needs_fallback, "coelho2024"
        )

    for editor in editors:
        designed_for = pd.Series(
            library.index.isin(design.loc[design.editor == editor, "guide"]),
            index=library.index,
        )
        lost = pd.Series(
            library.index.isin([g for g, e, *_ in unusable if e == editor]),
            index=library.index,
        )
        known = designed_for & ~lost

        made = study_variants[study_variants.editor == editor]
        counted_edits = (
            made.groupby("guide")
            .size()
            .reindex(library.index)
            .fillna(0)
            .astype("Int64")
        )
        library[f"n_variants_{editor.lower()}"] = counted_edits.where(known, pd.NA)
        categories = (
            made.groupby("guide")
            .codon_mutation_category.agg(list)
            .reindex(library.index)
        )
        library[f"mutation_category_{editor.lower()}"] = [
            most_severe(cats)
            if isinstance(cats, list)
            else (C.NO_EDITS if is_known else pd.NA)
            for cats, is_known in zip(categories, known)
        ]

    study_design = design[design.editor.isin(editors)]
    library["editors_designed"] = (
        study_design.groupby("guide").editor.agg(join_unique).reindex(library.index)
    )

    types = (
        our_controls.reindex(library.index).map(C.OUR_CONTROL_LABELS)
        if study == C.THIS_STUDY
        else coelho_types.reindex(library.index)
    )
    is_control = types.notna()
    library["guide_type"] = np.where(
        library.index.isin(design.guide), "targeting", pd.NA
    )
    library["control_type"] = pd.NA
    library.loc[is_control, "guide_type"] = types[is_control].map(lambda pair: pair[0])
    library.loc[is_control, "control_type"] = types[is_control].map(
        lambda pair: pair[1]
    )

    geneless = library.guide_type.isin(C.GENELESS_GUIDE_TYPES)
    library.loc[geneless, ["gene", "ensembl_gene"]] = pd.NA
    library["n_genes"] = (
        library.gene.fillna("").str.count(";").add(1).where(library.gene.notna(), 0)
    )

    genes = library.gene.fillna("").str.split(";")

    partly_known = genes.map(
        lambda names: 0 < sum(n in transcript_of.index for n in names) < len(names)
    )
    assert not partly_known.any(), (
        f"{study}: BEAGLE has a transcript for only some genes of "
        f"{library.index[partly_known].tolist()[:3]}"
    )
    library["transcript"] = genes.map(
        lambda names: ";".join(
            transcript_of[n] for n in names if n in transcript_of.index
        )
    ).replace("", pd.NA)
    library["transcript"], library["transcript_version"] = split_version(
        library["transcript"]
    )
    library = library.reset_index().reindex(columns=LIBRARY_SCHEMA)

    assert library.guide.is_unique, f"{study}: a guide is listed twice"
    assert library.guide_type.isin(C.GUIDE_TYPES).all(), f"{study}: unlabelled guides"
    assert library.control_type.dropna().isin(C.CONTROL_TYPES).all(), (
        f"{study}: unknown control_type"
    )
    assert library.loc[library.guide_type == "targeting", "gene"].notna().all(), (
        f"{study}: a targeting guide has no gene"
    )

    annotated_guides = set(library.loc[library.annotation_source == "beagle", "guide"])
    assert set(study_variants.guide) <= annotated_guides, (
        f"{study}: a variant belongs to a guide without a BEAGLE design row"
    )
    for editor in editors:
        made = set(study_variants.loc[study_variants.editor == editor, "guide"])
        counted = library.set_index("guide")[f"n_variants_{editor.lower()}"]
        assert counted.reindex(sorted(made)).notna().all(), (
            f"{study}: a guide with {editor} variants has no {editor} variant count"
        )

    library.to_csv(C.dataset_file(study, "library"), index=False)
    study_variants.to_csv(C.dataset_file(study, "vep"), index=False)
    with_vep = study_variants[f"{C.VEP_PREFIX}Consequence"].notna()
    print(f"\n{study}_library.csv: {len(library):,} guides")
    print(f"  guide_type: {library.guide_type.value_counts().to_dict()}")
    print(
        f"  annotation_source: {library.annotation_source.value_counts(dropna=False).to_dict()}"
    )
    for editor in editors:
        counts = library[f"mutation_category_{editor.lower()}"].value_counts()
        print(f"  most severe category ({editor}): {counts.to_dict()}")
    print(
        f"{study}_vep.csv: {len(study_variants):,} (guide, editor, variant) rows, "
        f"{study_variants.guide.nunique():,} guides, {int(with_vep.sum()):,} with a VEP row"
    )
