# Machine-learning result tables

This folder contains the final FoldX-output machine-learning result tables.

## ML tasks

| Task | Target |
|---|---|
| WT regression | FoldX wild-type energy per residue |
| Mutation regression | FoldX mutation ΔΔG |
| Mutation classification | Strict stabilising/destabilising FoldX mutation class |

For mutation classification, neutral mutations were excluded:

- stabilising: ΔΔG ≤ -0.5 kcal/mol
- destabilising: ΔΔG ≥ +0.5 kcal/mol
- neutral: -0.5 < ΔΔG < +0.5 kcal/mol, excluded

## Included tables

| File | Description |
|---|---|
| `ml_datasets_summary.csv` | Dataset sizes, unique accessions, targets, and feature lists |
| `wt_regression_performance.csv` | WT FoldX energy-per-residue regression performance |
| `mutation_ddg_regression_performance.csv` | Mutation ΔΔG regression performance |
| `mutation_classification_performance.csv` | Strict mutation-classification performance |
| `wt_regression_cv_predictions.csv` | WT regression cross-validated predictions |
| `mutation_ddg_regression_cv_predictions.csv` | Mutation ΔΔG cross-validated predictions |
| `mutation_classification_cv_predictions.csv` | Mutation-classification cross-validated predictions |
| `wt_regression_permutation_importance.csv` | WT regression permutation importance |
| `mutation_ddg_regression_permutation_importance.csv` | Mutation ΔΔG permutation importance |
| `mutation_classification_permutation_importance.csv` | Mutation-classification permutation importance |
| `MAIN_FOLDX_ML_REPORT.md` | Full generated ML report |
