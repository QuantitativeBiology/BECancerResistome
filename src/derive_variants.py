"""
Derive the variants of every guide that has a locus and a strand.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

GUIDE_LENGTH = 20
VARIANT_SCHEMA = [
    "guide",
    "study",
    "editor",
    "gene",
    "chromosome",
    "position",
    "genomic_ref",
    "genomic_alt",
    "guide_position",
    "guide_edit",
    "cdna_edit",
    "genomic_strand",
    # Prediction for the codon, second source for validation.
    "codon_group",
    "codon_amino_acid_edit",
    "codon_mutation_category",
    "variant_source",
]
VCF_HEADER = [
    "##fileformat=VCFv4.2",
    '##INFO=<ID=SRC,Number=1,Type=String,Description="beagle or derived">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]


def complement(base: str) -> str:
    return C.COMPLEMENT[base]


def edits_of(guide: str, editor: str, start: int, strand: str) -> list[dict]:
    """Every edit the editor makes in this guide, in genomic coordinates."""
    editable = C.EDITABLE_BASE[editor]
    alt_on_guide = {"C": "T", "A": "G"}[editable]
    rows = []
    for position in C.EDIT_WINDOW:
        if guide[position - 1] != editable:
            continue
        if strand == "+":
            genomic_position = start + position - 1
            ref, alt = editable, alt_on_guide
        else:
            genomic_position = start + GUIDE_LENGTH - position
            ref, alt = complement(editable), complement(alt_on_guide)
        rows.append(
            {
                "editor": editor,
                "position": genomic_position,
                "genomic_ref": ref,
                "genomic_alt": alt,
                "guide_position": position,
            }
        )
    return rows


def derive(library: pd.DataFrame, editors: list[str]) -> pd.DataFrame:
    """Derive the edits of every guide that has a locus and a strand."""
    located = library[library.target_start.notna() & library.genomic_strand.notna()]
    rows = []
    for record in located.to_dict("records"):
        for editor in editors:
            for edit in edits_of(
                record["guide"],
                editor,
                int(record["target_start"]),
                record["genomic_strand"],
            ):
                rows.append(
                    {
                        "guide": record["guide"],
                        "gene": record["gene"],
                        "chromosome": record["chromosome"],
                        "genomic_strand": record["genomic_strand"],
                        **edit,
                    }
                )
    return pd.DataFrame(rows)


KEY = ["guide", "editor", "position", "genomic_ref", "genomic_alt"]
variants = {}

for study in C.STUDIES:
    library = pd.read_csv(C.dataset_file(study, "library"))
    editors = C.SCREEN_EDITORS[study]
    known = pd.read_csv(
        C.dataset_file(study, "vep"),
        low_memory=False,
        usecols=KEY
        + [
            "gene",
            "chromosome",
            "codon_group",
            "codon_amino_acid_edit",
            "codon_mutation_category",
            "guide_edit",
            "cdna_edit",
        ],
    )
    derived = derive(library, editors)
    print(
        f"{study}: BEAGLE defines {len(known):,} edits over "
        f"{known.guide.nunique():,} guides; the rule derives {len(derived):,} over "
        f"{derived.guide.nunique():,}"
    )

    single_gene = set(library.loc[library.n_genes == 1, "guide"])
    comparable = sorted(set(known.guide) & set(derived.guide) & single_gene)
    from_beagle = known[known.guide.isin(comparable)]
    from_rule = derived[derived.guide.isin(comparable)]
    merged = from_beagle.merge(from_rule, on=KEY, how="outer", indicator=True)
    unexplained = merged[merged._merge == "left_only"]
    assert unexplained.empty, (
        f"{study}: the rule does not reproduce {len(unexplained):,} BEAGLE edits, e.g.\n"
        f"{unexplained[KEY].head(3)}"
    )
    recovered = merged[merged._merge == "right_only"]
    print(
        f"  the rule reproduces all {len(from_beagle):,} BEAGLE edits of the "
        f"{len(comparable):,} single-gene guides they share"
    )
    for row in recovered.to_dict("records"):
        print(
            f"    and recovers {row['guide']} {row['editor']} "
            f"{row['position']}{row['genomic_ref']}>{row['genomic_alt']}, "
            "whose transcript coordinate was lost"
        )

    known = known.assign(variant_source="beagle")
    new = derived.merge(known[KEY], on=KEY, how="left", indicator=True)
    new = new[new._merge == "left_only"].drop(columns="_merge")
    new = new.assign(variant_source="derived")
    table = pd.concat([known, new], ignore_index=True).assign(study=study)
    # BEAGLE states the edited base as "C_7"; a derived edit knows its position but
    # not the notation. Fill each from the other so one column serves both.
    from_notation = table.guide_edit.str.extract(r"_(\d+)$")[0]
    table["guide_position"] = table.guide_position.fillna(
        pd.to_numeric(from_notation, errors="coerce")
    ).astype("Int64")
    editable = table.editor.map(C.EDITABLE_BASE)
    table["guide_edit"] = table.guide_edit.fillna(
        editable + "_" + table.guide_position.astype(str)
    )

    table = table.reindex(columns=VARIANT_SCHEMA).sort_values(
        ["chromosome", "position", "guide", "editor"]
    )
    assert not table.duplicated(KEY).any(), f"{study}: an edit is listed twice"

    table.to_csv(C.dataset_file(study, "variants"), index=False)
    variants[study] = table
    by_source = table.variant_source.value_counts().to_dict()
    print(
        f"  {study}_variants.csv: {len(table):,} edits over "
        f"{table.guide.nunique():,} guides, {by_source}"
    )

everything = pd.concat(variants.values(), ignore_index=True)
distinct = (
    everything.drop_duplicates(["chromosome", "position", "genomic_ref", "genomic_alt"])
    .sort_values(["chromosome", "position", "genomic_ref", "genomic_alt"])
    .reset_index(drop=True)
)
distinct["id"] = (
    distinct.chromosome.astype("Int64").astype(str)
    + "_"
    + distinct.position.astype("Int64").astype(str)
    + "_"
    + distinct.genomic_ref
    + "_"
    + distinct.genomic_alt
)
lines = VCF_HEADER + [
    f"{int(r['chromosome'])}\t{int(r['position'])}\t{r['id']}\t{r['genomic_ref']}\t"
    f"{r['genomic_alt']}\t.\t.\tSRC={r['variant_source']}"
    for r in distinct.to_dict("records")
]
(C.DATA / "vep_input.vcf").write_text("\n".join(lines) + "\n")
print(
    f"\nvep_input.vcf: {len(distinct):,} distinct variants "
    f"({everything.variant_source.eq('derived').sum():,} of the "
    f"{len(everything):,} guide-level edits are newly derived)"
)
