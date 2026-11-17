"""
Annotate variants with VEP and SpliceAI scores.

VEP annotates a variant on every transcript it falls in. MANE Select transcript 
of the guide's own gene. Both studies are annotated on one transcript per gene.

Some variants fall outside their gene's MANE transcript and so have no MANE row.
Coelho's (PMID:39424923) library tiles BRAF on ENST00000644969, which reaches further than MANE's
ENST00000646891, and 3,382 BRAF variants land beyond it.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

VEP_OUTPUT = C.DATA / "vep_output.txt"
SPLICEAI_SCORES = C.DATA / "spliceai_scores.csv"
VARIANT_KEY = "variant_id"
CHUNK = 200_000


def vep_header(path: Path) -> list[str]:
    """Return the column names of the VEP output."""
    with path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            assert line.startswith("#"), f"no column header found in {path}"
            names = line.lstrip("#").rstrip("\n").split("\t")
            break
        else:
            raise AssertionError(f"{path} is empty")

    seen, unique, renamed = set(), [], []
    for name in names:
        if name in seen:
            renamed.append(name)
            name = f"{name}_dbnsfp"
        seen.add(name)
        unique.append(name)
    if renamed:
        print(f"  disambiguated duplicate VEP/dbNSFP columns: {renamed}")
    assert len(set(unique)) == len(unique), "the header is still not unique"
    return unique


def scan(path: Path, keep) -> pd.DataFrame:
    names = vep_header(path)
    kept, total = [], 0
    for chunk in pd.read_csv(
        path, sep="\t", comment="#", names=names, chunksize=CHUNK, low_memory=False
    ):
        total += len(chunk)
        kept.append(chunk[keep(chunk)])
    selected = pd.concat(kept, ignore_index=True)
    print(f"  scanned {total:,} rows, kept {len(selected):,}")
    return selected


def tidy(vep: pd.DataFrame) -> pd.DataFrame:
    basis = np.where(
        vep["MANE_SELECT"].ne(C.VEP_MISSING),
        "mane",
        np.where(vep["CANONICAL"].eq("YES"), "canonical", "other"),
    )
    vep = vep.rename(
        columns={
            "Uploaded_variation": VARIANT_KEY,
            "SYMBOL": "gene",
            "Feature": "transcript",
        }
    )
    assert not vep.duplicated([VARIANT_KEY, "gene"]).any(), (
        "a variant is kept twice for the same gene"
    )
    feature_columns = [c for c in vep.columns if c not in [VARIANT_KEY, "gene"]]
    vep[feature_columns] = vep[feature_columns].replace(C.VEP_MISSING, np.nan)
    vep = vep.rename(columns={c: f"{C.VEP_PREFIX}{c}" for c in feature_columns})
    vep[f"{C.VEP_PREFIX}transcript_basis"] = basis
    return vep


def read_spliceai() -> pd.DataFrame | None:
    if not SPLICEAI_SCORES.exists():
        print(f"{SPLICEAI_SCORES.name}: not present, skipping the splice scores")
        return None
    scores = pd.read_csv(SPLICEAI_SCORES)
    scores = scores.rename(columns={"symbol": "gene"})
    keep = [VARIANT_KEY, "gene"] + [
        c for c in scores.columns if c.startswith(("ds_", "dp_"))
    ]
    scores = scores[keep]
    # SpliceAI scores a variant once per gene it falls in, same grain as the VEP join.
    assert not scores.duplicated([VARIANT_KEY, "gene"]).any(), (
        "SpliceAI reports a variant twice for one gene"
    )
    print(
        f"{SPLICEAI_SCORES.name}: {len(scores):,} rows over "
        f"{scores[VARIANT_KEY].nunique():,} variants"
    )
    return scores.rename(
        columns={
            c: f"spliceai_{c}" for c in scores.columns if c.startswith(("ds_", "dp_"))
        }
    )


def collapse_per_transcript(annotation: pd.DataFrame) -> pd.DataFrame:
    order_column = f"{C.VEP_PREFIX}Ensembl_transcriptid"
    if order_column not in annotation.columns:
        return annotation

    listed = []
    for column in annotation.columns:
        if not column.startswith(C.VEP_PREFIX) or column == order_column:
            continue
        values = annotation[column].dropna().astype(str)
        if len(values) and values.str.contains(",").mean() > 0.5:
            listed.append(column)
    if not listed:
        return annotation

    transcripts = annotation[order_column].fillna("").str.split(";|,", regex=True)
    ours = annotation["transcript_bare"]
    index_of = [
        next((i for i, t in enumerate(row) if t == mine), None)
        for row, mine in zip(transcripts, ours)
    ]

    for column in listed:
        parts = annotation[column].fillna("").astype(str).str.split(",")
        picked = []
        for items, position in zip(parts, index_of):
            chosen = None
            if position is not None and position < len(items):
                chosen = items[position]
            if chosen in (None, "", "."):
                distinct = {i for i in items if i not in ("", ".")}
                chosen = distinct.pop() if len(distinct) == 1 else None
            picked.append(chosen)
        annotation[column] = pd.to_numeric(
            pd.Series(picked, index=annotation.index), errors="coerce"
        )
    print(
        f"  collapsed {len(listed)} per-transcript score columns to this row's transcript"
    )
    return annotation


def variant_ids(variants: pd.DataFrame) -> pd.Series:
    """The key the VCF gave each variant: chrom_pos_ref_alt."""
    return (
        variants.chromosome.astype("Int64").astype(str)
        + "_"
        + variants.position.astype("Int64").astype(str)
        + "_"
        + variants.genomic_ref
        + "_"
        + variants.genomic_alt
    )


print("MANE Select rows")
vep = tidy(scan(VEP_OUTPUT, lambda c: c["MANE_SELECT"].ne(C.VEP_MISSING)))
spliceai = read_spliceai()

wanted = set()
for study in C.STUDIES:
    variants = pd.read_csv(
        C.dataset_file(study, "variants"),
        usecols=["gene", "chromosome", "position", "genomic_ref", "genomic_alt"],
    )
    pairs = set(zip(variant_ids(variants), variants.gene))
    wanted |= pairs
unmatched = wanted - set(zip(vep[VARIANT_KEY], vep.gene))
print(
    f"{len(unmatched):,} (variant, gene) pairs have no MANE row; scanning for canonical"
)

fallback = None
if unmatched:
    rows = scan(
        VEP_OUTPUT,
        lambda c: [
            pair in unmatched for pair in zip(c["Uploaded_variation"], c["SYMBOL"])
        ],
    )
    # Deterministic pick: the gene's canonical transcript, else its first by identifier.
    rows = rows.sort_values(
        ["Uploaded_variation", "SYMBOL", "CANONICAL", "Feature"],
        ascending=[True, True, False, True],
    ).drop_duplicates(["Uploaded_variation", "SYMBOL"])
    fallback = tidy(rows)
    print(
        f"  recovered {len(fallback):,} pairs: "
        f"{fallback[f'{C.VEP_PREFIX}transcript_basis'].value_counts().to_dict()}, "
        f"{fallback[f'{C.VEP_PREFIX}Consequence'].value_counts().head(2).to_dict()}"
    )
    vep = pd.concat([vep, fallback], ignore_index=True)
    assert not vep.duplicated([VARIANT_KEY, "gene"]).any(), (
        "the fallback duplicates a pair the MANE pass already covered"
    )

for study in C.STUDIES:
    variants = pd.read_csv(C.dataset_file(study, "variants"))
    variants[VARIANT_KEY] = variant_ids(variants)

    annotation = variants.merge(
        vep, on=[VARIANT_KEY, "gene"], how="left", validate="m:1"
    )
    if spliceai is not None:
        annotation = annotation.merge(
            spliceai, on=[VARIANT_KEY, "gene"], how="left", validate="m:1"
        )

    assert len(annotation) == len(variants), "the join changed the number of edits"
    annotation["transcript_bare"] = (
        annotation[f"{C.VEP_PREFIX}transcript"].fillna("").str.split(".").str[0]
    )
    annotation = collapse_per_transcript(annotation)
    annotation = annotation.drop(columns="transcript_bare")
    annotated = annotation[f"{C.VEP_PREFIX}Consequence"].notna()
    annotation.to_csv(C.dataset_file(study, "annotation"), index=False)

    print(
        f"\n{study}_annotation.csv: {len(annotation):,} edits, {len(annotation.columns)} columns"
    )
    print(
        f"  with a VEP consequence: {int(annotated.sum()):,} "
        f"({annotated.mean():.1%}); by variant source: "
        f"{annotation.groupby('variant_source').apply(lambda g: f'{g[C.VEP_PREFIX + "Consequence"].notna().mean():.1%}', include_groups=False).to_dict()}"
    )
    for column in [
        f"{C.VEP_PREFIX}CADD_PHRED",
        f"{C.VEP_PREFIX}REVEL_score",
        "spliceai_ds_ag",
    ]:
        if column in annotation.columns:
            print(f"  {column}: {annotation[column].notna().mean():.1%} non-null")
