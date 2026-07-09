# Main FoldX machine-learning report

## Scope

The main machine-learning analyses were restricted to FoldX-derived stability targets. Docking-derived substrate-binding retention was not used as a major supervised ML target and was retained only as a secondary candidate-filtering layer.

## Datasets

- WT stability regression rows: 858
- WT unique accessions: 679
- Mutation ΔΔG regression rows: 2830
- Mutation unique accessions: 198
- Strict binary mutation-classification rows: 1682
- Strict binary mutation-classification accessions: 197

Strict mutation classification used:

- stabilising: ΔΔG <= -0.5 kcal/mol
- destabilising: ΔΔG >= +0.5 kcal/mol
- neutral: -0.5 < ΔΔG < +0.5 kcal/mol, excluded from binary classification

Class counts used in the classifier:

| class         |   rows |
|:--------------|-------:|
| destabilising |   1232 |
| stabilising   |    450 |

## Best WT energy-per-residue regression model

- Model: random_forest
- MAE: 0.3038 ± 0.0095
- RMSE: 0.4720 ± 0.0181
- R²: 0.8341 ± 0.0148

## Best mutation ΔΔG regression model

- Model: random_forest
- MAE: 1.1325 ± 0.1305
- RMSE: 2.1447 ± 0.4987
- R²: 0.4536 ± 0.0967

## Best strict mutation-stability classifier

- Model: random_forest
- Balanced accuracy: 0.6738 ± 0.0198
- Macro F1: 0.6927 ± 0.0184
- PR-AUC for stabilising class: 0.6624 ± 0.0715

## Interpretation note

These models predict computational FoldX-derived stability proxies, not experimentally measured thermostability or catalytic activity. Model performance was evaluated using accession-grouped cross-validation to reduce leakage between related mutation rows from the same protein.
