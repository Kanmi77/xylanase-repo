# Final workflow report

This report was generated automatically from the Snakemake workflow.

## Workflow mode

The current run is a test-mode reproducibility run. FoldX and docking were executed on a small balanced subset of bacterial GH10, bacterial GH11, fungal GH10, and fungal GH11 enzymes.

## Output summary

- Curated enzyme records: 577
- Sequence-feature rows: 577
- Structure-feature rows: 452
- FoldX WT stability rows: 4
- FoldX mutation rows: 12
- Docking rows: 24
- Successful docking rows: 24
- ML summary rows: 6
- Integrated candidate rows: 12

## Group coverage in integration table

| organism_type   | gh_family   |   rows |
|:----------------|:------------|-------:|
| bacterial       | GH10        |      3 |
| bacterial       | GH11        |      3 |
| fungal          | GH10        |      3 |
| fungal          | GH11        |      3 |

## Docking score summary

|       |   vina_affinity |
|:------|----------------:|
| count |       24        |
| mean  |       -7.44192  |
| std   |        0.950457 |
| min   |       -9.165    |
| 25%   |       -7.82675  |
| 50%   |       -7.4905   |
| 75%   |       -6.51875  |
| max   |       -6.004    |

## Top integrated candidates

|   rank | uniprot_accession   | organism_type   | gh_family   | pdb_id   | foldx_mutation_code   |   foldx_ddg |   mean_vina_affinity |   best_vina_affinity | foldx_class   |   integrated_score | candidate_category              |
|-------:|:--------------------|:----------------|:------------|:---------|:----------------------|------------:|---------------------:|---------------------:|:--------------|-------------------:|:--------------------------------|
|      1 | A0A0M9BNX9          | bacterial       | GH11        | 7KV0     | WB211V                |  -0.0172052 |              -8.3045 |               -9.165 | neutral       |           0.988399 | neutral_with_docking_support    |
|      2 | A0A0S2I9J6          | bacterial       | GH10        | 9UPA     | AA353G                |  -0.0871928 |              -7.7715 |               -8.537 | neutral       |           0.860756 | neutral_with_docking_support    |
|      3 | A0A0S2I9J6          | bacterial       | GH10        | 9UPA     | AB353V                |  -0.111008  |              -7.5065 |               -8.528 | neutral       |           0.812113 | neutral_with_docking_support    |
|      4 | A0A0M9BNX9          | bacterial       | GH11        | 7KV0     | HB26G                 |   2.07149   |              -8.29   |               -9.031 | destabilising |           0.724429 | lower_priority_or_destabilising |
|      5 | A0A2T4BZZ5          | fungal          | GH11        | 8YJI     | FA93G                 |   0.515756  |              -6.8345 |               -7.497 | destabilising |           0.551475 | lower_priority_or_destabilising |
|      6 | A0A1L9WG58          | fungal          | GH10        | 6Q8M     | VB343S                |   1.11226   |              -6.9975 |               -7.514 | destabilising |           0.511238 | lower_priority_or_destabilising |
|      7 | A0A0M9BNX9          | bacterial       | GH11        | 7KV0     | HA26A                 |   4         |              -8.3055 |               -9.018 | destabilising |           0.492168 | lower_priority_or_destabilising |
|      8 | A0A1L9WG58          | fungal          | GH10        | 6Q8M     | VA343G                |   1.04355   |              -6.904  |               -7.288 | destabilising |           0.489539 | lower_priority_or_destabilising |
|      9 | A0A0S2I9J6          | bacterial       | GH10        | 9UPA     | PA1A                  |   4         |              -7.7955 |               -8.53  | destabilising |           0.367904 | lower_priority_or_destabilising |
|     10 | A0A2T4BZZ5          | fungal          | GH11        | 8YJI     | SA190V                |   1.91406   |              -6.7485 |               -7.493 | destabilising |           0.364623 | lower_priority_or_destabilising |

## ML branch note

The ML branch executed successfully, but in test mode the dataset is too small for thesis-scale interpretation. Final ML interpretation should use the full FoldX and docking tables.
