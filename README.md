# BECancerResistome

BECancerResistome is a project leveraging CRISPR base editing to systematically identify mutations associated with drug resistance in cancer.
This project aims to unravel the genetic basis of resistance mechanisms by screening cancer cell lines under various drug treatments.

## Contributions

Instituto Superior Técnico (U Lisboa) and The Broad Insitute of MIT and Harvard.

## License

This project is licensed under the BSD 3-Clause "New" or "Revised" License. See the LICENSE file for details.


# Data Processing 

## Base-Editing Screen Analysis Pipeline (CoelhoNatGen & dt studies)

Notebook Processing_final_20260927.ipynb processes raw per-guide read counts from base-editing CRISPR
screens (the `CoelhoNatGen` and `dt` studies), computes log2 fold-changes and
Z-scores against Intergenic-control guides, classifies significant hits
(Plasmid- and Control-relative), runs replicate-correlation QC, and filters
guides by off-target predictions.

## Requirements

- Python 3.10+
- `poola` (installed via `!pip install poola` in the first cell)
- `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`

## Required input files

Place these **6 files in the same directory as the notebook** before running:

1. `counts_with_mismatches_CBE_activity_A375_DO_PIC.csv`
2. `counts_with_mismatches_CBE_activity_A375_LIN_SCH.csv`
3. `counts_with_mismatches_CBE_activity_HT29_DO_PIC.csv`
4. `counts_with_mismatches_CBE_activity_HT29_LIN_SCH.csv`
5. `MC_readcouns.csv`
6. `Base Editing Screens Samplesheet.xlsx`

Two additional files are required partway through the notebook:

- `EG_off-targets.xlsx` — FlashFry off-target predictions (used in the
  **Off-target filtering** section)
- Intermediate CSVs the notebook writes and then re-reads from disk (see
  **Important: run cells in order** below)

> **Naming note:** internally the pipeline no longer uses "MC"/"EG" — those
> studies are now called `CoelhoNatGen` and `dt` respectively (see the table
> in **Outputs**). The three items above still carry their original `MC`/`EG`.

## Pipeline overview

### 1. Setup
Installs `poola`.

### 2. Core processing module (single large code cell)
Defines the reusable pipeline. `process_df` is the only function most callers
need — everything else is a building block it calls internally:

- `poola_normalize_and_filter` — log2-normalizes read counts via `poola` and
  drops guides whose plasmid (pDNA) abundance is >3 SD below the mean. **From
  this point on every count column holds LOG2 values, not raw counts.**
- `process_CoelhoNatGen` / `process_dt` — reshape each study's raw counts
  into one wide table (Guide/Gene/Editor + one column per condition/replicate).
- `compute_log2fc_and_zscore` — for each comparison in the samplesheet's
  `Comparisons` sheet, computes `log2(numerator) - log2(denominator)` (a
  subtraction, since inputs are already log2), then Z-scores each fold-change
  column using only `Intergenic control` guides as the null (no-effect)
  distribution.
- `collapse_replicates_min` — collapses `_RepA`, `_RepB`, ... columns per
  condition down to one value per guide, keeping whichever replicate's value
  is closest to zero (the more conservative call).
- `compute_lfc_zscore_and_save` — runs the above and writes three CSVs per
  study/comparison (see Outputs below).
- `process_df(file1, study, file2=None, comparison=None)` — dispatches to the
  `CoelhoNatGen` or `dt` pipeline based on `study` (`'CoelhoNatGen'` or
  `'dt'`).

### 3. Run the pipeline for each study
```python
df_CoelhoNatGen, lfc_data_CoelhoNatGen, zscores_CoelhoNatGen, zscores_min_CoelhoNatGen = process_df(
    'MC_readcouns.csv', 'CoelhoNatGen', comparison='CoelhoNatGen_processed'
)

df, lfc_data, zscores, zscores_min = process_df(
    'counts_with_mismatches_CBE_activity_A375_DO_PIC.csv', 'dt',
    'counts_with_mismatches_CBE_activity_A375_LIN_SCH.csv',
    comparison='dt'
)
```
(HT29's input files are hardcoded inside the module, not passed as arguments.)

### 4. Merge Plasmid- and Control-relative datasets
Re-reads `dt_zscores_min.csv` and `CoelhoNatGen_zscores_min.csv` from disk,
splits each into "Plasmid" vs "Control" condition columns, renames columns
for consistency between the two sources, and outer-merges them on
Guide/Gene/Editor into `plasmid_data` and `control_data`.

### 5. Hit assignment
Run **twice** — once on `plasmid_data`, once on `control_data`:
- Z-scores are converted to two-sided p-values against the Intergenic-control
  null distribution, then FDR-corrected (Benjamini–Hochberg, α = 0.05).
- `classify_hit` labels each guide `"positive"`, `"negative"`, or `"non-hit"`
  based on FDR < 0.05 and the sign of the Z-score.
- Results: `df_plasmid` and `df_control` (displayed inline, not written to
  disk by the notebook as-is).

### 6. Screen QC — replicate correlation plots
Reads `dt_zscores_rep.csv`, plots per-drug/cell-line scatterplots of
Replicate A vs. Replicate B Z-scores (density-colored, with Pearson r),
arranged in a 2×3 grid (cell lines × drugs).

### 7. Off-target filtering
Loads `EG_off-targets.xlsx` (FlashFry predictions against hg38, all possible
SpG-Cas9 NGN PAMs) and `dt_zscores_min.csv`. Keeps only guides with exactly
one perfect genomic match (`OT_0mm` in `{0, 1}`), then inner-merges that
filter onto the Z-score table.

## Outputs

| File | Produced by |
|---|---|
| `CoelhoNatGen_processed.csv` | `CoelhoNatGen` pipeline run |
| `CoelhoNatGen_LFC_rep.csv`, `CoelhoNatGen_zscores_rep.csv`, `CoelhoNatGen_zscores_min.csv` | `CoelhoNatGen` pipeline run |
| `dt_processed.csv` | `dt` pipeline run |
| `dt_LFC_rep.csv`, `dt_zscores_rep.csv`, `dt_zscores_min.csv` | `dt` pipeline run |
| `replicate_correlations.png` / `.pdf` | Screen QC section |
| `dt_zscores_min_off-targets-filtered.csv` | Off-target filtering section |

## Important: run cells in order

Several sections **re-read CSVs from disk** rather than reusing the in-memory
DataFrames from earlier cells (e.g. the merge step re-reads
`dt_zscores_min.csv` / `CoelhoNatGen_zscores_min.csv`, and the QC plot reads
`dt_zscores_rep.csv`). This means the notebook must be run top-to-bottom on a
first pass — skipping or reordering cells will raise `FileNotFoundError` or
silently operate on stale data.

## Single-cell Perturb-seq analysis (Supplementary Figure 4)

Notebooks in `Perturb-seq_SuppFig4/` use the ABE Perturb-seq data from Coelho et al. (2024) to test whether
single-cell gene expression can predict which base-editing guides cause drug resistance.

### Requirements
Python 3.10+
numpy, pandas, scipy, scikit-learn, matplotlib (processing)
anndata, mofapy2, xgboost, optuna, openpyxl (machine learning)

### Required input files
Place these in a `data/` folder next to the notebooks:

logcounts_ABE.csv and meta_ABE.csv — download from Zenodo (doi:10.5281/zenodo.16951444)
protein-coding_genes_HGNC.txt and hgnc-symbol-check.csv — included in this repository
zscores-unambiguous-VEPannotated-processed-control.csv.gz — VEP annotations per guide, included in this repository
41588_2024_1948_MOESM4_ESM.xlsx — Supplementary Table 2 of Coelho et al. (2024), from the paper's Supplementary Information

### 1. Perturb-seq_processing_final.ipynb
Prepares the expression matrix for modelling:
- Keeps single-guide cells, and double-guide cells whose two guides have different variant classes. For these,
  the guide matching the cell's highest-impact class is recorded as responsible for the phenotype.
- Merges canonical drug resistance and drug addiction into a binary resistance vs control label.
- Selects highly variable genes with a two-component Gaussian mixture on per-gene standard deviation.
- Keeps protein-coding genes (HGNC), renaming genes with outdated symbols using the HGNC Multi-Symbol Checker.

### 2. ML_single_cell.ipynb
Classifies cells as resistance or control:
- Reduces gene expression to MOFA+ latent factors (5–25 factors).
- Optionally adds VEP annotations for each cell's responsible guide.
- Trains XGBoost over 10 repeated train/test splits. Controls are downsampled to balance each training set,
  hyperparameters are tuned with Optuna on the training sets only, and models are evaluated on held-out test cells.

Run the processing notebook first: the ML notebook reads its outputs from `output/`.

### Outputs
File	Produced by
output/data_checkpoint.csv	Processing notebook
output/double_guide_cells.csv	Processing notebook
output/symbols_for_hgnc_check.txt	Processing notebook (input for the HGNC Multi-Symbol Checker)
output/mofa/	ML notebook (MOFA+ models, factor matrices, variance explained)
output/xgb_summary.csv	ML notebook
output/figures/f1_and_gap_comparison.png	ML notebook
