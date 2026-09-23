#!/usr/bin/env python3
"""
Figure 4 (v9): gene-level flag effects, NMD-escape disease genes vs 1:1 length-matched
controls. Python port of gene_flag_effectsize.R, with the Protein and Domain flags removed.

    flag ~ is_case + (1 | pair_id)      linear-probability mixed model, matched design
    BH correction across the flags shown, within each class (FS, SNV)

Seven flags: PPI, Pfam, SLiM, PTM, NLS, LCS and IDR. The first six come from
gene_flags_per_gene.csv (design == "matched"); IDR is built exactly as in the R script:
the NMD-escape region (CDS nt -> aa, floor((nt-1)/3)+1) overlaps a UniProt IDR of >= 20 aa,
0 if the accession has no such IDR. MoRF is excluded (n < 5), as before.

Usage:
  python fig4_gene_flags_v9.py gene_flags_per_gene.csv OUT_PREFIX \\
         [uniprot_human_idrs.csv] [gene_all_withflags_0826.csv]

Without the IDR file the IDR panel is left out and the other six flags are fitted and
BH-corrected among themselves.
Outputs: OUT.png, OUT.pdf, OUT.csv (per-flag statistics), OUT_gene_flags.csv (flags used).
"""
import sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SRC = sys.argv[1] if len(sys.argv) > 1 else 'gene_flags_per_gene.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'Figure4'
IDR_FILE = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ('', '-') else None
COORD_FILE = sys.argv[4] if len(sys.argv) > 4 else 'gene_all_withflags_0826.csv'
MIN_IDR_AA = 20

FLAGS = [('PPI', 'PPI_flag'), ('Pfam', 'Pfam_flag'), ('SLiM', 'SLiM_flag'),
         ('PTM', 'PTM_flag'), ('NLS', 'NLS_flag'), ('LCS', 'LCS_flag'), ('IDR', 'IDR_flag')]
STRATA = ('FS', 'SNV')

d = pd.read_csv(SRC, low_memory=False)
d = d[d.design == 'matched'].copy()
d['is_case'] = d.is_disease.astype(int)

# ------------------------------------------------------------------ IDR flag
if IDR_FILE:
    idr = pd.read_csv(IDR_FILE)
    idr = idr[(idr.End - idr.Start + 1) >= MIN_IDR_AA].copy()
    idr['acc'] = idr.ID.astype(str).str.split('-').str[0].str.split('.').str[0]
    by = {a: x[['Start', 'End']].values for a, x in idr.groupby('acc')}

    w = pd.read_csv(COORD_FILE, low_memory=False).dropna(subset=['NMD_region_start', 'NMD_region_end'])
    tx = w.drop_duplicates('ensembl_transcript_id').set_index('ensembl_transcript_id')
    up = w.drop_duplicates('uniprot').set_index('uniprot')
    d['ns'] = d.ensembl_transcript_id.map(tx.NMD_region_start)
    d['ne'] = d.ensembl_transcript_id.map(tx.NMD_region_end)
    m = d.ns.isna()
    d.loc[m, 'ns'] = d.loc[m, 'uniprot'].map(up.NMD_region_start)
    d.loc[m, 'ne'] = d.loc[m, 'uniprot'].map(up.NMD_region_end)

    def hit(acc, ns, ne):
        if pd.isna(ns) or pd.isna(ne):
            return 0
        a = max(1, int((ns - 1) // 3) + 1)
        b = int((ne - 1) // 3) + 1
        a, b = min(a, b), max(a, b)
        iv = by.get(str(acc).split('-')[0].split('.')[0])
        if iv is None:
            return 0
        return int(np.any(np.minimum(b, iv[:, 1]) - np.maximum(a, iv[:, 0]) + 1 > 0))

    d['IDR_flag'] = [hit(a, s, e) for a, s, e in zip(d['uniprot'], d['ns'], d['ne'])]
    print('IDR flag (escape region overlaps an IDR >= %d aa): %d/%d genes positive'
          % (MIN_IDR_AA, int(d['IDR_flag'].sum()), len(d)))
    print('  positives by stratum/cohort:\n%s'
          % d.groupby(['stratum', 'is_case']).IDR_flag.mean().round(3).to_string())
else:
    print('no IDR file given - fitting the six non-IDR flags only')

flags = [(l, c) for l, c in FLAGS if c in d.columns]
d.to_csv(OUT + '_gene_flags.csv', index=False,
         columns=['pair_id', 'stratum', 'is_case', 'hgnc_symbol', 'uniprot',
                  'ensembl_transcript_id'] + [c for _, c in flags])

# ------------------------------------------------------------------ models
rows = []
for st in STRATA:
    for lab, col in flags:
        x = d[d.stratum == st][['pair_id', 'is_case', col]].rename(columns={col: 'y'}).dropna(subset=['y'])
        rec = dict(stratum=st, flag=lab, n_pairs=x.pair_id.nunique(),
                   n_case=int(x.is_case.sum()), n_ctrl=int((1 - x.is_case).sum()),
                   rate_case=x[x.is_case == 1].y.mean(), rate_ctrl=x[x.is_case == 0].y.mean())
        if x.y.nunique() < 2 or x.pair_id.nunique() < 3:
            rec.update(beta=np.nan, se=np.nan, p=np.nan)
        else:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    f = smf.mixedlm('y ~ is_case', x, groups=x['pair_id']).fit(reml=True)
                rec.update(beta=f.params['is_case'], se=f.bse['is_case'], p=f.pvalues['is_case'])
            except Exception as e:
                print('fit failed', st, lab, e)
                rec.update(beta=np.nan, se=np.nan, p=np.nan)
        rows.append(rec)
res = pd.DataFrame(rows)
res['q'] = np.nan
for st in STRATA:
    m = (res.stratum == st) & res.p.notna()
    res.loc[m, 'q'] = multipletests(res.loc[m, 'p'], method='fdr_bh')[1]
res['lo'], res['hi'] = res.beta - 1.96 * res.se, res.beta + 1.96 * res.se
res.to_csv(OUT + '.csv', index=False)
print(res.round(4).to_string(index=False))

# ------------------------------------------------------------------ plot (no title)
FS_C, SG_C, FS_L, SG_L = '#1F5FA6', '#B2182B', '#BBD3E8', '#E6B7B2'
labels = [l for l, _ in flags]
fig, ax = plt.subplots(figsize=(10.5, 5.2))
ymin = min(-0.05, res.lo.min() - 0.07)
ystar = ymin + 0.02
for i, lab in enumerate(labels):
    for st, dx, mk, c, cl in (('FS', -0.16, 'o', FS_C, FS_L), ('SNV', 0.16, 'D', SG_C, SG_L)):
        r = res[(res.stratum == st) & (res.flag == lab)]
        if r.empty or pd.isna(r.beta.iloc[0]):
            continue
        r = r.iloc[0]
        sig = r.q < 0.05
        col = c if sig else cl
        ax.plot([i + dx] * 2, [r.lo, r.hi], color=col, lw=1.8, zorder=2)
        ax.scatter(i + dx, r.beta, marker=mk, s=55, color=col, zorder=3)
        if sig:
            ax.text(i + dx, r.hi + 0.012, f'{r.beta:+.3f}', ha='center', va='bottom',
                    fontsize=8.5, fontweight='bold', color=c)
            star = '***' if r.q < 0.001 else '**' if r.q < 0.01 else '*'
            ax.text(i + dx, ystar, star, ha='center', va='center', fontsize=11,
                    fontweight='bold', color=c)
ax.axhline(0, color='grey', ls='--', lw=0.9)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=10.5, fontweight='bold')
ax.set_ylabel('Effect size (β)', fontsize=11)
ax.set_xlim(-0.6, len(labels) - 0.4)
ax.set_ylim(ymin, max(0.15, res.hi.max() + 0.20))
ax.grid(axis='y', color='#e5e5e5', lw=0.7)
ax.set_axisbelow(True)
h = [Line2D([], [], marker='o', ls='', color=FS_C, ms=7, label='Frameshift - significant'),
     Line2D([], [], marker='D', ls='', color=SG_C, ms=6.5, label='Stopgain - significant'),
     Line2D([], [], marker='o', ls='', color=FS_L, ms=7, label='Not significant after BH')]
ax.legend(handles=h, loc='upper right', frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(OUT + '.png', dpi=200, facecolor='white')
plt.savefig(OUT + '.pdf', facecolor='white')
print('saved', OUT + '.png')
