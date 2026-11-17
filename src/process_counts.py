"""
Raw counts to hit calls.

1. Log-normalise each sample on its own, log2(reads per million + 1).
2. Drop sgRNAs under-represented in the plasmid pool, z <= -3 in any pDNA sample.
3. Take the comparisons from the samplesheet.
4. Log2 fold-change per comparison, over the surviving sgRNAs.
5. Z-score every fold-change column on the non-targeting controls.
6. Collapse replicates to the smallest absolute (min abd) z, then call hits.

ABE and CBE edit the same guide differently, so they are separate perturbations: a
condition keeps them apart wherever a study screened both, and omits the editor from the
name where there is only one. Collapsing by min |z| bounds the magnitude but does not
require the replicates to agree in direction, so `replicates_agree` and `n_replicates`
record that instead of changing the published call.
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
COMPARISON_REFERENCE = {"drug_vs_untreated": "untreated", "arm_vs_plasmid": "plasmid"}


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


def condition_of(comparisons: pd.DataFrame, editors: list[str]) -> pd.Series:
    reference = comparisons.comparison_type.map(COMPARISON_REFERENCE)
    unmapped = sorted(set(comparisons.comparison_type[reference.isna()]))
    assert not unmapped, f"comparison_type with no reference arm: {unmapped}"

    pieces = [comparisons.cell_line]
    if len(editors) > 1:
        pieces.append(comparisons.editor)
    pieces += [comparisons.arm, reference.radd("vs_")]
    return pieces[0].str.cat(pieces[1:], sep="_")


def collapse_replicates(zscores: pd.DataFrame, condition: pd.Series) -> pd.DataFrame:
    frames = []
    for name, columns in condition.groupby(condition, sort=True):
        block = zscores[list(columns.index)].to_numpy(dtype=float)
        present = ~np.isnan(block)
        # inf parks the missing entries where argmin will never choose them
        magnitude = np.where(present, np.abs(block), np.inf)
        zscore = block[np.arange(len(block)), magnitude.argmin(axis=1)]
        zscore[~present.any(axis=1)] = np.nan
        frames.append(
            pd.DataFrame(
                {
                    "guide": zscores.index,
                    "condition": name,
                    "zscore": zscore,
                    "n_replicates": present.sum(axis=1),
                    "replicates_agree": ~(
                        (block > 0).any(axis=1) & (block < 0).any(axis=1)
                    ),
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def call_hits(collapsed: pd.DataFrame, controls: pd.Index) -> pd.DataFrame:
    """Two-sided p against the condition controls, FDR correction"""
    calls = []
    for _, block in collapsed.groupby("condition", sort=True):
        block = block.dropna(subset="zscore")
        null = block.loc[block.guide.isin(controls), "zscore"]
        centre, spread = null.mean(), null.std()
        observed = block.zscore.to_numpy()
        p = 2 * stats.norm.sf(np.abs(observed - centre) / spread)
        fdr = multipletests(p, method="fdr_bh")[1]
        significant = fdr < FDR_THRESHOLD
        calls.append(
            block.assign(
                control_mean=centre,
                control_sd=spread,
                pvalue=p,
                fdr=fdr,
                call=np.where(
                    significant & (observed > centre),
                    RESISTANCE,
                    np.where(significant, SENSITIVITY, NEUTRAL),
                ),
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
    study_comparisons["condition"] = condition_of(
        study_comparisons, C.SCREEN_EDITORS[study]
    )
    assert study_comparisons.comparison_id.is_unique, (
        f"{study}: a comparison_id is listed twice, so one would overwrite the other"
    )
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
    assert len(lfc.columns) == len(usable), (
        f"{study}: {len(usable)} comparisons gave {len(lfc.columns)} fold-change columns"
    )
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
    kept_apart = ", kept apart by editor" if len(C.SCREEN_EDITORS[study]) > 1 else ""
    print(
        f"  {collapsed.condition.nunique()} conditions after collapsing "
        f"replicates{kept_apart}"
    )
    print(f"  calls at FDR < {FDR_THRESHOLD}: {hits.call.value_counts().to_dict()}")
    significant = hits[hits.call != NEUTRAL]
    discordant = int((~significant.replicates_agree).sum())
    print(
        f"  of those, {discordant:,} rest on replicates that disagree in direction "
        f"({discordant / max(len(significant), 1):.1%})"
    )
    drug = per_condition[per_condition.index.str.endswith("_vs_untreated")]
    if len(drug):
        print("  drug versus untreated, by condition:")
        print(drug.to_string().replace("\n", "\n    "))
