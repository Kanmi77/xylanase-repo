# Full Snakemake test-workflow status

The full Snakemake workflow has been executed successfully in test mode.

## Workflow coverage

The workflow runs the following stages:

1. UniProt enzyme record fetching
2. GH10/GH11 bacterial/fungal curation
3. Sequence-feature calculation
4. MAFFT alignment
5. FastTree phylogeny
6. Conservation summary
7. PDB structure inventory and structural-feature extraction
8. FoldX WT stability calculation
9. FoldX mutation screening
10. Receptor cleaning and PDBQT preparation
11. AutoDock Vina docking with xylobiose and xylotriose
12. Machine-learning workflow execution
13. Integrated candidate ranking
14. Final workflow report generation

## Final test-mode output counts

- Curated enzyme records: 577
- Sequence-feature rows: 577
- Structure-feature rows: 452
- FoldX WT stability rows: 4
- FoldX mutation rows: 12
- Docking rows: 24
- Successful docking rows: 24
- Failed docking rows: 0
- ML summary rows: 6
- Integrated candidate rows: 12

## Test-mode group coverage

- bacterial GH10: 3 mutation rows
- bacterial GH11: 3 mutation rows
- fungal GH10: 3 mutation rows
- fungal GH11: 3 mutation rows

## Docking summary

AutoDock Vina docking completed successfully for both ligands.

- Ligands: xylobiose and xylotriose
- Total docking jobs: 24
- Successful jobs: 24
- Failed jobs: 0
- Mean Vina affinity: -7.441917 kcal/mol
- Minimum Vina affinity: -9.165 kcal/mol
- Maximum Vina affinity: -6.004 kcal/mol

## Important note

This is a test-mode reproducibility run. FoldX and docking were intentionally limited to a small balanced enzyme subset to verify that the workflow executes correctly. Final thesis-scale interpretation requires switching test mode off in the configuration and running the larger dataset.
