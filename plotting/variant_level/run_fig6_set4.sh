#!/usr/bin/env bash
# Figure 6 / Figure 7 composition + phase-separation panels, design: set4
# Variant-level: ClinVar P/LP vs gnomAD benign NMD-escape PTCs WITHIN the same
# disease genes, so no gene-level matching ratio applies.
# MIN_TAIL = 20 aa (the published setting; the script default of 2 is NOT it).
set -euo pipefail
python composition_lmm_run.py master_table.csv 20
python plot_composition_run.py .
