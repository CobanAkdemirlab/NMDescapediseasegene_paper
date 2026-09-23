#!/usr/bin/env python3
"""Remove the rows listed in rows_to_rebuild.csv from the master table so that a
--resume run of wt_vs_mutant_genomic_fixed.py rebuilds only those variants (plus the
10 variants that are still missing).
    python drop_rows_for_rebuild.py master_table_all.csv rows_to_rebuild.csv
The original file is kept as <name>.before_rebuild.csv"""
import sys, shutil
import pandas as pd
master, drop = sys.argv[1], sys.argv[2]
shutil.copy(master, master.replace('.csv', '.before_rebuild.csv'))
m = pd.read_csv(master, low_memory=False)
keys = set(pd.read_csv(drop).key)
keep = m[~m.key.isin(keys)]
keep = keep[[c for c in keep.columns if not c.startswith('DELTA_')]]   # recomputed at the end
keep.to_csv(master, index=False)
print('removed %d rows; %d remain' % (len(m) - len(keep), len(keep)))
