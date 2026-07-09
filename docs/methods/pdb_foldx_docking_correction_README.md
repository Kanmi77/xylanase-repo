# PDB-aware FoldX and docking correction layer

This workflow was added to correct the earlier mutation/docking layer by including PDB-derived xylanase structures in addition to MODELLER-derived structures.

## Purpose

The original mutation-screening layer was dominated by MODELLER-derived structures and did not adequately include experimentally resolved PDB-derived representatives. This correction layer introduced a PDB-aware subset, performed FoldX BuildModel mutation screening, prepared mutant and wild-type receptors for AutoDock Vina, calculated paired mutant-vs-WT docking differences, and generated a final non-tiered ranked candidate table.

## Main correction steps

1. Reconciled PDB-derived xylanase structures with the curated master dataset.
2. Selected representative PDB structures across GH10/GH11 and bacterial/fungal groups.
3. Extracted chain sequences and residue mappings.
4. Generated mutation candidates.
5. Ran FoldX BuildModel for PDB-derived structures.
6. Parsed FoldX ΔΔG values and selected stabilizing mutants.
7. Prepared mutant receptors for docking.
8. Prepared matching wild-type receptors using the same receptor preparation logic.
9. Ran paired AutoDock Vina docking with xylobiose and xylotriose.
10. Calculated mutant-minus-WT docking deltas.
11. Generated a final non-tiered PDB-aware ranked candidate table.

## Interpretation

The PDB-aware correction provides an experimentally anchored structural correction layer. FoldX ΔΔG values were used to identify stabilizing mutations, while paired docking deltas were used only as comparative docking support. Because fallback receptor PDBQT preparation used simplified zero-charge receptor formatting, docking scores were interpreted comparatively rather than as absolute binding free energies.

## Files committed

The Git repository stores:
- correction scripts
- parsed FoldX summaries
- docking score summaries
- mutant-vs-WT paired delta tables
- final ranked candidate tables
- reports and README files

The repository does not store:
- raw FoldX work directories
- generated mutant PDB structures
- receptor PDBQT files
- Vina pose output files
- large intermediate docking folders
