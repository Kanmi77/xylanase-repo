# Previous PDB-reference MD MDP provenance

The generic MDP files in `scripts/08_md/mdp/` are tracked in git, but they should not be treated as the exact parameter files used for all previous PDB-reference MD runs.

The completed PDB-reference MD folders were inspected directly to recover the actual run-level MDP files. These were found under:

- `md/_archive/20260223_123717/systems/`
- `md_10ns/systems/`

A comparison between the tracked repo-level MDP files and the completed-run MDP files showed that many actual run files differed from `scripts/08_md/mdp/*.mdp`, especially for production `md.mdp`.

Therefore, for methodology and reproducibility of the previous PDB-reference MD workflow, the authoritative MDP files are stored under:

- `scripts/08_md/legacy_pdb_reference_md/mdp_actual_previous/`

The comparison reports are stored as:

- `results/reports/previous_pdb_md_actual_mdp_comparison.csv`
- `results/reports/previous_pdb_md_actual_vs_repo_mdp_comparison.csv`

The workflow scripts remain under:

- `scripts/08_md/run_stage1_top15_3temps.sh`
- `scripts/08_md/run_md_one.sh`
- `scripts/08_md/02_run_stage1_md_top15_3temps.sh`
- `scripts/08_md/09_analyze_stage1_top15_3temps.sh`
- `scripts/08_md/10_summarize_stage1_md_metrics.py`

For thesis writing, cite the recovered run-level MDPs, not the generic template MDP folder, as the source of the previous PDB-reference MD parameters.
