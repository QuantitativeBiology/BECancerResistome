"""
Raw counts to hit calls.

1. Log-normalise each sample on its own, log2(reads per million + 1).
2. Drop sgRNAs under-represented in the plasmid pool, z <= -3 in any pDNA sample.
3. Take the comparisons from the samplesheet.
4. Log2 fold-change per comparison, over the surviving sgRNAs.
5. Z-score every fold-change column on the non-targeting controls.
6. Collapse replicates to the smallest absolute (min abd) z, then call hits.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C

PDNA_Z_FLOOR = -3.0
FDR_THRESHOLD = 0.05
RESISTANCE, SENSITIVITY, NEUTRAL = "resistance", "sensitivity", "neutral"


def log_normalise(counts: pd.DataFrame) -> pd.DataFrame:
    per_million = counts / counts.sum() * 1e6
    return np.log2(per_million + 1)


def plasmid_passing(lognorm: pd.DataFrame, plasmids: list[str]) -> pd.Series:
    depleted = pd.DataFrame(index=lognorm.index)
    for sample in plasmids:
        column = lognorm[sample]
        z = (column - column.mean()) / column.std()
        depleted[sample] = z <= PDNA_Z_FLOOR
    return ~depleted.any(axis=1)


def null_guides(library: pd.DataFrame) -> pd.Index:
    geneless = library.guide_type.isin(C.GENELESS_GUIDE_TYPES)
    return pd.Index(library.loc[geneless, "guide"])


def condition_of(comparisons: pd.DataFrame) -> pd.Series:
    reference = comparisons.comparison_type.map(
        {"drug_vs_untreated": "untreated", "arm_vs_plasmid": "plasmid"}
    )
    return comparisons.cell_line + "_" + comparisons.arm + "_vs_" + reference


def collapse_replicates(zscores: pd.DataFrame, condition: pd.Series) -> pd.DataFrame:
    collapsed = {}
    for name, columns in condition.groupby(condition):
        block = zscores[list(columns.index)]
        least = block.abs().idxmin(axis=1)
        collapsed[name] = pd.Series(
            [
                block.at[guide, column] if isinstance(column, str) else np.nan
                for guide, column in least.items()
            ],
            index=block.index,
        )
    return pd.DataFrame(collapsed).reindex(sorted(collapsed), axis=1)


def call_hits(collapsed: pd.DataFrame, controls: pd.Index) -> pd.DataFrame:
    """Two-sided p against the condition controls, FDR correction"""
    calls = []
    for condition in collapsed.columns:
        values = collapsed[condition]
        null = values.reindex(controls).dropna()
        centre, spread = null.mean(), null.std()
        observed = values.dropna()
        p = 2 * stats.norm.sf(((observed - centre) / spread).abs())
        fdr = multipletests(p, method="fdr_bh")[1]
        significant = fdr < FDR_THRESHOLD
        calls.append(
            pd.DataFrame(
                {
                    "guide": observed.index,
                    "condition": condition,
                    "zscore": observed.to_numpy(),
                    "control_mean": centre,
                    "control_sd": spread,
                    "pvalue": p,
                    "fdr": fdr,
                    "call": np.where(
                        significant & (observed.to_numpy() > 0),
                        RESISTANCE,
                        np.where(significant, SENSITIVITY, NEUTRAL),
                    ),
                }
            )
        )
    return pd.concat(calls, ignore_index=True)


samples = pd.read_excel(C.DATA / "samplesheet.xlsx", sheet_name="Samples")
comparisons = pd.read_excel(C.DATA / "samplesheet.xlsx", sheet_name="Comparisons")

for study in C.STUDIES:
    counts = pd.read_csv(C.dataset_file(study, "counts"), index_col=0)
    library = pd.read_csv(C.dataset_file(study, "library"))
    print(
        f"\n{'=' * 70}\n{study}: {counts.shape[0]:,} sgRNAs x {counts.shape[1]} samples"
    )

    # ---- 1. log-normalise ------------------------------------------------------------
    lognorm = log_normalise(counts)

    # ---- 2. plasmid filter -----------------------------------------------------------
    plasmids = samples.loc[
        (samples.study == study) & (samples.arm_kind == "plasmid"), "sample_id"
    ].tolist()
    passing = plasmid_passing(lognorm, plasmids)
    lognorm = lognorm[passing]
    print(
        f"  pDNA filter over {len(plasmids)} plasmid sample(s): kept "
        f"{int(passing.sum()):,}, dropped {int((~passing).sum()):,}"
    )
    lognorm.to_csv(C.dataset_file(study, "lognorm"))

    # ---- 3. comparisons --------------------------------------------------------------
    study_comparisons = comparisons[comparisons.study == study].copy()
    study_comparisons["condition"] = condition_of(study_comparisons)
    usable = study_comparisons[
        study_comparisons.numerator_sample.isin(lognorm.columns)
        & study_comparisons.denominator_sample.isin(lognorm.columns)
    ]
    assert len(usable) == len(study_comparisons), (
        f"{study}: a comparison names a sample the counts do not have"
    )

    # ---- 4. log2 fold-changes --------------------------------------------------------
    lfc = pd.DataFrame(
        {
            row.comparison_id: lognorm[row.numerator_sample]
            - lognorm[row.denominator_sample]
            for row in usable.itertuples()
        }
    )
    lfc.index.name = counts.index.name
    lfc.to_csv(C.dataset_file(study, "lfc"))
    print(f"  {len(lfc.columns)} replicate-level comparisons")

    # ---- 5. z-score on the controls --------------------------------------------------
    controls = null_guides(library).intersection(lfc.index)
    assert len(controls) > 10, f"{study}: only {len(controls)} control guides survive"
    zscores = (lfc - lfc.loc[controls].mean()) / lfc.loc[controls].std()
    zscores.to_csv(C.dataset_file(study, "zscores"))
    print(
        f"  z-scored on {len(controls):,} non-targeting controls "
        f"({', '.join(sorted(set(library.set_index('guide').loc[controls, 'guide_type'])))})"
    )

    # ---- 6. collapse replicates and call hits ----------------------------------------
    condition = usable.set_index("comparison_id").condition
    collapsed = collapse_replicates(zscores, condition)
    hits = call_hits(collapsed, controls)
    hits.to_csv(C.dataset_file(study, "hits"), index=False)

    per_condition = hits.groupby("condition").call.value_counts().unstack(fill_value=0)
    print(f"  {len(collapsed.columns)} conditions after collapsing replicates")
    print(f"  calls at FDR < {FDR_THRESHOLD}: {hits.call.value_counts().to_dict()}")
    drug = per_condition[per_condition.index.str.endswith("_vs_untreated")]
    if len(drug):
        print("  drug versus untreated, by condition:")
        print(drug.to_string().replace("\n", "\n    "))
