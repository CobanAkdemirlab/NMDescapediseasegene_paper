#!/usr/bin/env python3
"""
aa20_features.py  --  20-amino-acid-fraction feature panel ONLY.

Drop-in replacement for the Mensah localCIDER/BioPython panel: computes just the
fraction of each of the 20 standard amino acids in a protein/region sequence.
No localCIDER, no BioPython, no metapredict, no network -- pure Python.

Use it two ways:
  1. As a library:   from aa20_features import compute_aa20, AA20_NAMES
  2. As a driver on the wt_vs_mutant pipeline output (see analyze_master_table.py)
"""
import numpy as np

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA20_SET = set(AA20)
AA20_NAMES = [f"frac_{a}" for a in AA20]   # frac_A ... frac_Y

def _clean(seq):
    return "".join(c for c in str(seq).upper() if c in AA20_SET)

def compute_aa20(seq):
    """Return {frac_A:.., frac_C:.., ...} for the 20 standard amino acids.
    Non-standard characters (X, *, gaps, lowercase) are ignored; the denominator
    is the number of standard residues. Empty/too-short sequences -> all NaN."""
    s = _clean(seq)
    n = len(s)
    if n < 1:
        return {f"frac_{a}": np.nan for a in AA20}
    return {f"frac_{a}": s.count(a) / n for a in AA20}
