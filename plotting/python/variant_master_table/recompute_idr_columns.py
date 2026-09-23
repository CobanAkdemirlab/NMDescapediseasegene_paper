#!/usr/bin/env python3
"""
Recompute ONLY the IDR columns of an existing master table from a real IDR
coordinate file, without rebuilding any protein.

The master table was built without --idr-coords, so every IDR column (IDRseq,
IDR<feature>, WTIDR<feature>, IDRLength, TrueIDRstart/end, IDR_intervals, ...) was
taken from a "C-terminal third" placeholder (IDR_source = fallback_C_terminal_third).
This script replaces them using real IDR intervals and recomputes the matching
DELTA_IDR* columns. Sequences and all Fullpep/PFS columns are left untouched.

Usage (same folder as wt_vs_mutant_genomic_fixed.py, same Python environment):
    python recompute_idr_columns.py master_table_all.csv IDR_FILE.csv master_table_all_idr.csv

IDR_FILE: columns ID, Start, End (1-based, inclusive). ID may be
  * an ENSP (HG38_CanonicalPEP_IDRs_coord.csv), matched via the table's protein_id, or
  * a UniProt accession (e.g. from fetch_uniprot_idrs.py), matched via uniprotswissprot.
UniProt is tried first, then ENSP (same priority as the main script). A protein that
is not in the file gets no IDR (IDR columns NaN, IDRLength 0), not a placeholder.
Optional 4th argument: minimum IDR length in aa (default 0 = keep all intervals;
use 20 to match the Figure 4 definition).
Optional 5th argument: gene=ACC pairs for rows with no uniprotswissprot, e.g.
  GLI2=P10070,CELF2=O95319

Isoform remapping: if IDR_FILE's companion FASTA (same name, .fasta; written by
fetch_uniprot_idrs_nodeps.py) is present and the UniProt canonical sequence differs
from the table's WTSequence, UniProt coordinates are mapped onto the Ensembl protein
by global alignment (IDR_source = catalog_uniprot_remapped); residues with no
aligned counterpart are dropped and segments shorter than the minimum are removed.
"""
import sys
import numpy as np
import pandas as pd
import wt_vs_mutant_genomic_fixed as w

SRC, IDR_FILE, OUT = sys.argv[1:4]
MIN_LEN = int(sys.argv[4]) if len(sys.argv) > 4 else 0

import os
from Bio import Align
FILL = dict(p.split('=') for p in sys.argv[5].split(',')) if len(sys.argv) > 5 else {}
coords = w.load_idr_coords(IDR_FILE)
fasta = IDR_FILE.rsplit('.', 1)[0] + '.fasta'
useqs = {}
if os.path.exists(fasta):
    name = None
    for line in open(fasta):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]; useqs[name] = ''
        elif name:
            useqs[name] += line
    print('UniProt sequences for remapping: %d' % len(useqs))
else:
    print('WARNING: %s not found - UniProt coordinates are used as-is (no isoform check)' % fasta)
_aligner = Align.PairwiseAligner()
_aligner.mode = 'global'
_aligner.match_score, _aligner.mismatch_score = 2, -1
_aligner.open_gap_score, _aligner.extend_gap_score = -5, -0.5
_map_cache = {}

def remap(acc, wt, iv):
    """Map 0-based [s,e) UniProt intervals onto the Ensembl protein wt."""
    u = useqs.get(acc)
    if not u or u == wt:
        return iv, False
    key = (acc, wt)
    if key not in _map_cache:
        aln = _aligner.align(u, wt)[0]
        pos = {}
        for (us, ue), (ws, we) in zip(*aln.aligned):
            for k in range(ue - us):
                pos[us + k] = ws + k
        _map_cache[key] = pos
    pos = _map_cache[key]
    out = []
    for s0, e0 in iv:
        hit = [pos[p] for p in range(s0, e0) if p in pos]
        if hit and (max(hit) - min(hit) + 1) >= max(MIN_LEN, 1):
            out.append((min(hit), max(hit) + 1))
    return out, True

if MIN_LEN:
    coords = {k: [(s, e) for s, e in v if e - s >= MIN_LEN] for k, v in coords.items()}
    coords = {k: v for k, v in coords.items() if v}
print('IDR file: %d proteins with >=1 interval' % len(coords))

m = pd.read_csv(SRC, low_memory=False).reset_index(drop=True)
m = m[[c for c in m.columns if not c.startswith('DELTA_')]]      # recomputed below
idr_cols = [c for c in m.columns if (c.startswith('IDR') or c.startswith('WTIDR'))
            and c not in ('IDR_intervals', 'IDR_n_segments', 'IDR_source')]

def clean(x):
    return '' if pd.isna(x) else str(x).split('.')[0].strip()

cache, new_rows, src_count = {}, [], {}
for i, r in enumerate(m.itertuples(index=False)):
    acc, ensp = clean(getattr(r, 'uniprotswissprot', '')), clean(getattr(r, 'protein_id', ''))
    if not acc and r.gene_symbol in FILL:
        acc = FILL[r.gene_symbol]
        m.at[m.index[i], 'uniprotswissprot'] = acc
    wt, mut = str(r.WTSequence), str(r.Sequence)
    if acc and acc in coords:
        iv, moved = remap(acc, wt, coords[acc])
        src = 'catalog_uniprot_remapped' if moved else 'catalog_uniprot'
    elif ensp and ensp in coords:
        iv, src = coords[ensp], 'catalog_ensp'
    else:
        iv, src = [], 'catalog_no_idr'
    src_count[src] = src_count.get(src, 0) + 1
    wt_idr, mut_idr = w.slice_idrs(wt, iv), w.slice_idrs(mut, iv)
    rec = {}
    for seq, prefix in ((mut_idr, 'IDR'), (wt_idr, 'WTIDR')):
        k = (prefix, seq)
        if seq not in cache:
            cache[seq] = w.compute_features(seq) if seq else {n: np.nan for n in w.ALL_FEATURE_NAMES}
        rec.update({prefix + n: v for n, v in cache[seq].items()})
    rec.update({'IDRseq': mut_idr, 'WTIDRseq': wt_idr, 'IDRLength': len(mut_idr),
                'WTIDRLength': len(wt_idr),
                'TrueIDRstart': iv[0][0] if iv else np.nan,
                'TrueIDRend': iv[-1][1] if iv else np.nan,
                'IDR_intervals': ';'.join('%d..%d' % (s + 1, e) for s, e in iv),
                'IDR_n_segments': len(iv), 'IDR_source': src})
    new_rows.append(rec)
    if (i + 1) % 500 == 0:
        print('  %d/%d' % (i + 1, len(m)))

new = pd.DataFrame(new_rows, index=m.index)
remapped = sorted({g for g, rr in zip(m.gene_symbol, new_rows) if rr['IDR_source'] == 'catalog_uniprot_remapped'})
print('isoform-remapped genes (%d): %s' % (len(remapped), ', '.join(remapped)))
for c in new.columns:
    m[c] = new[c]
m = w.add_delta_columns(m)
m.to_csv(OUT, index=False)
print('IDR source:', src_count)
print('rows with an IDR: %d of %d; mutant IDR shorter than WT IDR (truncated/altered): %d' % (
    (m.WTIDRLength > 0).sum(), len(m), (m.IDRLength < m.WTIDRLength).sum()))
print('wrote', OUT)
