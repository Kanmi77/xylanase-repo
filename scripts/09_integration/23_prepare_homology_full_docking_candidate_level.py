#!/usr/bin/env python3

from pathlib import Path
import argparse
import re
import pandas as pd
import numpy as np

ROOT = Path(".").resolve()

DEFAULT_INPUT = "results/correlation/08_mutation_ddg_vs_docking_change_input_used.csv"
DEFAULT_OUT = "results/integration/homology_full_docking_candidate_level.csv"
DEFAULT_REPORT = "results/integration/HOMOLOGY_FULL_DOCKING_CANDIDATE_LEVEL_REPORT.md"


def find_col(df, regexes):
    for rgx in regexes:
        for c in df.columns:
            if re.search(rgx, c, flags=re.I):
                return c
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--report", default=DEFAULT_REPORT)
    args = parser.parse_args()

    inp = ROOT / args.input
    out = ROOT / args.out
    report = ROOT / args.report

    out.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)

    if not inp.exists():
        raise SystemExit(f"Input file not found: {inp}")

python scripts/10_integration/23_prepare_homology_full_docking_candidate_level.pypyd worst mutant-minus-WT binding differences.\n\n")
Saved candidate-level homology docking table: /home/ubuntu/xylanase-thesis/results/integration/homology_full_docking_candidate_level.csv
Saved report: /home/ubuntu/xylanase-thesis/results/integration/HOMOLOGY_FULL_DOCKING_CANDIDATE_LEVEL_REPORT.md
Input rows: 300
Candidate-level rows: 150

protein uniprot_accession mutation candidate_id  n_ligand_rows              ligands  foldx_ddg  foldx_ddg_first  delta_binding_mut_minus_wt_mean  delta_binding_mut_minus_wt_best  delta_binding_mut_minus_wt_worst  docking_improved_or_retained_all_ligands  docking_improved_or_retained_any_ligand  delta_top3_mut_minus_wt_mean  delta_top3_mut_minus_wt_best  delta_top3_mut_minus_wt_worst  wt_binding_energy_mean  mut_binding_energy_mean  wt_mean_top3_binding_mean  mut_mean_top3_binding_mean  foldx_predicted_stabilising  mean_docking_improved_vs_wt
 D5UGW9            D5UGW9    LA23A D5UGW9_LA23A              2 xylobiose;xylotriose     -5.482           -5.482                          -0.0065                           -0.042                             0.029                                     False                                     True                     -0.049833                     -0.083667                      -0.016000                 -6.3215                  -6.3280                  -6.215833                   -6.265667                         True                         True
 C7NP66            C7NP66    NA12F C7NP66_NA12F              2 xylobiose;xylotriose     -4.269           -4.269                           0.0910                            0.014                             0.168                                     False                                    False                     -0.003667                     -0.022333                       0.015000                 -6.0500                  -5.9590                  -5.844833                   -5.848500                         True                        False
 Q9HEN6            Q9HEN6    YA14P Q9HEN6_YA14P              2 xylobiose;xylotriose     -4.185           -4.185                           0.1780                           -0.018                             0.374                                     False                                     True                      0.049667                     -0.143000                       0.242333                 -6.2620                  -6.0840                  -6.004333                   -5.954667                         True                        False
 Q9HGX1            Q9HGX1    DA11K Q9HGX1_DA11K              2 xylobiose;xylotriose     -3.863           -3.863                           0.1915                            0.118                             0.265                                     False                                    False                      0.056000                     -0.041667                       0.153667                 -6.3160                  -6.1245                  -6.154333                   -6.098333                         True                        False
 B0ZSE5            B0ZSE5    NA65G B0ZSE5_NA65G              2 xylobiose;xylotriose     -3.423           -3.423                           0.5310                            0.282                             0.780                                     False                                    False                      0.563667                      0.344333                       0.783000                 -7.0295                  -6.4985                  -6.907667                   -6.344000                         True                        False
 E3WF08            E3WF08    YA25F E3WF08_YA25F              2 xylobiose;xylotriose     -3.339           -3.339                           0.8515                            0.640                             1.063                                     False                                    False                      0.727333                      0.611000                       0.843667                 -6.3525                  -5.5010                  -6.175667                   -5.448333                         True                        False
 Q9X584            Q9X584    VA60T Q9X584_VA60T              2 xylobiose;xylotriose     -3.285           -3.285                          -0.1335                           -0.242                            -0.025                                      True                                     True                     -0.274500                     -0.355000                      -0.194000                 -6.9230                  -7.0565                  -6.428833                   -6.703333                         True                         True
 Q5XQ46            Q5XQ46    DA35G Q5XQ46_DA35G              2 xylobiose;xylotriose     -3.245           -3.245                           0.0070                            0.001                             0.013                                     False                                    False                     -0.013000                     -0.020000                      -0.006000                 -6.6290                  -6.6220                  -6.511333                   -6.524333                         True                        False
 Q8J0T4            Q8J0T4    AA53G Q8J0T4_AA53G              2 xylobiose;xylotriose     -3.231           -3.231                           0.0415                            0.036                             0.047                                     False                                    False                      0.008667                     -0.028667                       0.046000                 -6.9895                  -6.9480                  -6.795000                   -6.786333                         True                        False
 Q6U892            Q6U892    EA63G Q6U892_EA63G              2 xylobiose;xylotriose     -3.150           -3.150                          -0.0565                           -0.060                            -0.053                                      True                                     True                      0.013167                     -0.020000                       0.046333                 -6.3460                  -6.4025                  -6.297667                   -6.284500                         True                         True
 D2XV89            D2XV89    AA55G D2XV89_AA55G              2 xylobiose;xylotriose     -3.118           -3.118                           0.5325                           -0.174                             1.239                                     False                                     True                      0.574000                     -0.013000                       1.161000                 -6.7120                  -6.1795                  -6.645000                   -6.071000                         True                        False
 A5FIE5            A5FIE5    TA27S A5FIE5_TA27S              2 xylobiose;xylotriose     -3.074           -3.074                           0.3165                           -0.001                             0.634                                     False                                     True                      0.317333                      0.015333                       0.619333                 -7.2125                  -6.8960                  -7.142667                   -6.825333                         True                        False
 Q2U7D0            Q2U7D0    EA47D Q2U7D0_EA47D              2 xylobiose;xylotriose     -3.059           -3.059                           0.2010                            0.152                             0.250                                     False                                    False                      0.222000                      0.200000                       0.244000                 -6.4735                  -6.2725                  -6.417167                   -6.195167                         True                        False
 B6VF01            B6VF01    FA11K B6VF01_FA11K              2 xylobiose;xylotriose     -3.046           -3.046                           0.4810                            0.106                             0.856                                     False                                    False                      0.364333                      0.097667                       0.631000                 -5.3680                  -4.8870                  -5.159833                   -4.795500                         True                        False
 Q9HEN4            Q9HEN4    YA14P Q9HEN4_YA14P              2 xylobiose;xylotriose     -2.903           -2.903                          -0.1095                           -0.288                             0.069                                     False                                     True                     -0.030667                     -0.154667                       0.093333                 -5.8005                  -5.9100                  -5.724833                   -5.755500                         True                         True
 C9Z2V1            C9Z2V1    AA64G C9Z2V1_AA64G              2 xylobiose;xylotriose     -2.686           -2.686                          -0.2385                           -0.570                             0.093                                     False                                     True                     -0.237833                     -0.616333                       0.140667                 -6.8895                  -7.1280                  -6.820667                   -7.058500                         True                         True
 Q3YAW6            Q3YAW6    EA50A Q3YAW6_EA50A              2 xylobiose;xylotriose     -2.679           -2.679                          -2.2590                           -4.048                            -0.470                                      True                                     True                     -2.251667                     -4.265000                      -0.238333                 -2.5440                  -4.8030                  -2.396500                   -4.648167                         True                         True
 Q14RS0            Q14RS0    HA45P Q14RS0_HA45P              2 xylobiose;xylotriose     -2.623           -2.623                          -0.1485                           -0.993                             0.696                                     False                                     True                      0.105333                     -0.416667                       0.627333                 -7.0010                  -7.1495                  -6.652333                   -6.547000                         True                         True
 Q60041            Q60041    NA33G Q60041_NA33G              2 xylobiose;xylotriose     -2.612           -2.612                          -0.2425                           -0.507                             0.022                                     False                                     True                     -0.249167                     -0.327333                      -0.171000                 -7.3820                  -7.6245                  -7.085000                   -7.334167                         True                         True
 Q6TDT4            Q6TDT4    YA15F Q6TDT4_YA15F              2 xylobiose;xylotriose     -2.584           -2.584                           0.0760                            0.006                             0.146                                     False                                    False                      0.135667                     -0.040667                       0.312000                 -5.4000                  -5.3240                  -5.310500                   -5.174833                         True                        False
(base) ubuntu@lsi:~/xylanase-thesis$ cd ~/xylanase-thesis

python - <<'PY'
import pandas as pd

p = "results/integration/homology_full_docking_candidate_level.csv"
df = pd.read_csv(p)

print("Rows:", len(df))
print("Columns:", len(df.columns))
print(df.columns.tolist())
print(df.head(20).to_string(index=False))
