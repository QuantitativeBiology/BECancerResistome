# BECancerResistome

BECancerResistome is a project leveraging CRISPR base editing to systematically identify
mutations associated with drug resistance in cancer. This project aims to unravel the
genetic basis of resistance mechanisms by screening cancer cell lines under various drug
treatments.

## Analysis

Preprocessing runs from `src/`: sgRNA mapping and predicted edits (`map_guides.py`,
`derive_variants.py`), variant annotation with Ensembl VEP (`run_vep.sh`,
`annotate_variants.py`), read count normalisation and hit calling (`process_counts.py`),
and the feature tables (`build_features.py`).

The notebooks in `notebooks/` then run in order, each reading what the previous ones
write:

**01_screen_qc** — assesses replicate agreement, separation of essential and non-targeting guides, fold changes by guide class, and hit counts per gene and condition.

**02_variant_analysis** — collapses guides into amino acid variants, compares the two cell lines across drugs, maps variant effects along proteins, and assesses agreement with the independent screen by Coelho et al.

**03_feature_exploration** — assembles and selects annotation features based on coverage, variability and variant-level specificity. Redundant features (Spearman correlation ≥0.9 with a previously selected feature) are removed. The selected features are saved to `data/features_selected.csv` and reused in notebooks 04 and 05 to ensure a consistent feature set.

**04_xgboost_model** — trains a three-class classifier (sensitivity, neutral and resistance). Neutral rows are downsampled to match the size of the second-largest class. Performance is assessed using average precision over 30 repeats, with folds grouped by genomic edit. The notebook compares feature subsets, validates predictions against the external screen, reports SHAP feature attributions, and exports design matrices for reuse in other notebooks.

**05_perturbseq_prediction** — tests whether the same variant annotations improve the prediction of drug resistance from single-cell expression data using the Perturb-seq screen by Coelho et al. It also distinguishes drug addiction from canonical resistance.

**06_held_out_validation** — evaluates model generalisation by withholding one gene, cell line or compound at a time, retraining the model and predicting the withheld group.

Figures are written to a gitignored `plots/` folder.

## Environment

```
conda env create -f environment.yml
conda activate becr
```

## Contributions

INESC-ID & Instituto Superior Técnico - Universidade de Lisboa and The Broad Insitute of MIT and Harvard.

## License

This project is licensed under the BSD 3-Clause "New" or "Revised" License. See the
LICENSE file for details.
