#!/usr/bin/env python3
"""
Gene-level flag effects, NMD-escape disease variants vs gene-matched controls.

For each flag and each variant class (frameshift, stopgain):
    flag ~ is_case + (1 + is_case | gene)   linear mixed model (statsmodels mixedlm)
restricted to genes carrying >=1 disease AND >=1 control variant of that class.
The random slope lets the disease-control difference vary between genes, so
genes with hundreds of variants cannot drive the p-value on their own. The
random-intercept-only model (1 | gene) is reported alongside as beta_ri / p_ri;
it gives similar betas but far too narrow CIs here (anti-conservative).
BH adjustment across the 9 flags within each class.

Flags (1/TRUE -> 1, 0/FALSE -> 0, missing -> excluded):
  PPI     variant_ppi_overlap
  Protein variant_protein_flag
  Pfam    ptc_before_max_pfam_end   (PTC upstream of the last Pfam domain end)
  Domain  variant_domain_flag
  SLiM    variant_slim_flag
  PTM     variant_ptm_flag
  NLS     variant_nls_flag
  LCS     variant_LCS_flag
  IDR     region from the variant to the C-terminus overlaps a predicted IDR
          by >= 20 aa (metapredict v3 disorder domains on the WT protein
          translated from the `coding` column)

Usage: python flag_lme.py variant_all5_0918_1.csv [out_prefix] [idr_coords.csv]
  With idr_coords.csv (ID,Start,End; e.g. from fetch_uniprot_idrs.py) the IDR flag
  uses those curated intervals instead of metapredict predictions.
"""
import sys, warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else 'variant_all5_0918_1.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'flag_effects_lme'
MIN_IDR_OVERLAP = 20

d = pd.read_csv(SRC, low_memory=False)
d['gene'] = d['hgnc_symbol']
d['cls'] = np.where(d.group.str.startswith('fs'), 'fs', 'sg')
d['is_case'] = d.group.str.endswith('disease').astype(int)


def to01(s):
    m = {'TRUE': 1, 'FALSE': 0, '1': 1, '0': 0, '1.0': 1, '0.0': 0,
         True: 1, False: 0, 1: 1, 0: 0, 1.0: 1, 0.0: 0}
    return s.map(lambda x: m.get(x if not isinstance(x, str) else x.strip().upper(), np.nan))


FLAGS = [('PPI', 'variant_ppi_overlap'), ('Protein', 'variant_protein_flag'),
         ('Pfam', 'ptc_before_max_pfam_end'), ('Domain', 'variant_domain_flag'),
         ('SLiM', 'variant_slim_flag'), ('PTM', 'variant_ptm_flag'),
         ('NLS', 'variant_nls_flag'), ('LCS', 'variant_LCS_flag'), ('IDR', 'idr_flag')]
for lab, col in FLAGS[:-1]:
    d['f_' + lab] = to01(d[col])

# ------------------------------------------------------------------ IDR flag
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}


def translate(cds):
    cds = str(cds).upper()
    p = []
    for i in range(0, len(cds) - 2, 3):
        a = CODON.get(cds[i:i + 3], 'X')
        if a == '*':
            break
        p.append(a)
    return ''.join(p)


IDR_FILE = sys.argv[3] if len(sys.argv) > 3 else None
idr = {}
if IDR_FILE and 'IDR_intervals' in pd.read_csv(IDR_FILE, nrows=0).columns:
    # master table: UniProt IDRs already remapped onto the Ensembl protein (IDR_intervals,
    # 1-based inclusive "s..e;s..e"); one transcript per gene, WT protein identical
    mt = pd.read_csv(IDR_FILE, usecols=['gene_symbol', 'IDR_intervals'], low_memory=False)
    g2iv = {}
    for g, s in mt.drop_duplicates('gene_symbol').values:
        iv = [] if pd.isna(s) or not str(s).strip() else [
            (int(p.split('..')[0]) - 1, int(p.split('..')[1])) for p in str(s).split(';')]
        g2iv[g] = [(a, b) for a, b in iv if b - a >= MIN_IDR_OVERLAP]
    for tid, g in d.drop_duplicates('ensembl_transcript_id')[['ensembl_transcript_id', 'gene']].values:
        idr[tid] = g2iv.get(g)
    IDR_SOURCE = 'master table ' + IDR_FILE
    print('IDR intervals from %s: %d of %d proteins have >=1 IDR' % (
        IDR_FILE, sum(bool(v) for v in idr.values()), len(idr)))
elif IDR_FILE:
    # curated IDR coordinates (ID = UniProt accession or ENSP; Start/End 1-based inclusive)
    cat = pd.read_csv(IDR_FILE)
    cat['ID'] = cat['ID'].astype(str).str.split('.').str[0]
    byid = {k: [(int(s) - 1, int(e)) for s, e in zip(g.Start, g.End) if e - s + 1 >= MIN_IDR_OVERLAP]
            for k, g in cat.groupby('ID')}
    ucol = 'uniprotswissprot' if 'uniprotswissprot' in d else 'uniprot'
    for tid, acc in d.drop_duplicates('ensembl_transcript_id')[['ensembl_transcript_id', ucol]].values:
        acc = str(acc).split('.')[0] if pd.notna(acc) else ''
        idr[tid] = byid.get(acc, [])        # protein in catalogue with no IDR -> no overlap
    IDR_SOURCE = 'catalogue ' + IDR_FILE
    print('IDR intervals from %s: %d of %d proteins have >=1 IDR' % (
        IDR_FILE, sum(bool(v) for v in idr.values()), len(idr)))
else:
    import metapredict as meta
    for tid, cds in d.drop_duplicates('ensembl_transcript_id')[['ensembl_transcript_id', 'coding']].values:
        prot = translate(cds)
        try:
            dom = meta.predict_disorder_domains(prot)
            idr[tid] = [(int(s), int(e)) for s, e in dom.disordered_domain_boundaries]  # 0-based, end-excl
        except Exception as e:
            warnings.warn(f'metapredict failed for {tid}: {e}')
            idr[tid] = None
    IDR_SOURCE = 'metapredict'
    print('IDR predictions:', sum(v is not None for v in idr.values()), 'proteins')


def idr_flag(r):
    ints = idr.get(r.ensembl_transcript_id)
    if ints is None or pd.isna(r.cds_mutation_loc):
        return np.nan
    start = int((r.cds_mutation_loc - 1) // 3)          # 0-based aa of the variant
    ov = sum(max(0, e - max(s, start)) for s, e in ints)
    return int(ov >= MIN_IDR_OVERLAP)


d['f_IDR'] = d.apply(idr_flag, axis=1)
pd.DataFrame([(t, ';'.join(f'{s + 1}..{e}' for s, e in v) if v else '') for t, v in idr.items()],
             columns=['transcript', 'IDR_intervals_metapredict']).to_csv(OUT + '_idr_intervals.csv', index=False)

# ------------------------------------------------------------------ models
rows = []
for cls in ('fs', 'sg'):
    for lab, _ in FLAGS:
        x = d[(d.cls == cls)].dropna(subset=['f_' + lab])[['gene', 'is_case', 'f_' + lab]]
        x = x.rename(columns={'f_' + lab: 'y'})
        both = x.groupby('gene').is_case.agg(['min', 'max'])
        x = x[x.gene.isin(both[(both['min'] == 0) & (both['max'] == 1)].index)]
        rec = dict(cls=cls, flag=lab, n_case=int(x.is_case.sum()),
                   n_ctrl=int((1 - x.is_case).sum()), n_gene=x.gene.nunique(),
                   rate_case=x[x.is_case == 1].y.mean(), rate_ctrl=x[x.is_case == 0].y.mean())
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                f = smf.mixedlm('y ~ is_case', x, groups=x['gene'],
                                re_formula='~is_case').fit(reml=True, method='lbfgs')
                f0 = smf.mixedlm('y ~ is_case', x, groups=x['gene']).fit(reml=True, method='lbfgs')
            ci = f.conf_int().loc['is_case']
            rec.update(beta=f.params['is_case'], lo=ci[0], hi=ci[1], p=f.pvalues['is_case'],
                       converged=f.converged,
                       beta_ri=f0.params['is_case'], p_ri=f0.pvalues['is_case'])
        except Exception as e:
            rec.update(beta=np.nan, lo=np.nan, hi=np.nan, p=np.nan)
            print('fit failed', cls, lab, e)
        rows.append(rec)
res = pd.DataFrame(rows)
res['q'] = np.nan
for cls in ('fs', 'sg'):
    m = (res.cls == cls) & res.p.notna()
    res.loc[m, 'q'] = multipletests(res.loc[m, 'p'], method='fdr_bh')[1]
res.to_csv(OUT + '.csv', index=False)
print(res.round(4).to_string())

# ------------------------------------------------------------------ plot
FS, SG = '#1f5aa6', '#b2182b'
FS_L, SG_L = '#b9cfe8', '#e8b4b4'
fig, ax = plt.subplots(figsize=(10.5, 5.2))
ymin = min(res.lo.min(), 0) - 0.07
ystar = ymin + 0.025
labels = [l for l, _ in FLAGS]
for i, lab in enumerate(labels):
    for cls, dx, mk, c, cl in (('fs', -0.16, 'o', FS, FS_L), ('sg', 0.16, 'D', SG, SG_L)):
        r = res[(res.cls == cls) & (res.flag == lab)].iloc[0]
        if pd.isna(r.beta):
            continue
        sig = r.q < 0.05
        col = c if sig else cl
        ax.plot([i + dx] * 2, [r.lo, r.hi], color=col, lw=1.8, zorder=2)
        ax.scatter(i + dx, r.beta, marker=mk, s=55, color=col, zorder=3)
        if sig:
            ax.text(i + dx, r.hi + 0.012, f'{r.beta:+.3f}', ha='center', va='bottom',
                    fontsize=8.5, fontweight='bold', color=c)
            st = '***' if r.q < 0.001 else '**' if r.q < 0.01 else '*'
            ax.text(i + dx, ystar, st, ha='center', va='center', fontsize=11,
                    fontweight='bold', color=c)
ax.axhline(0, color='grey', ls='--', lw=0.9)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=10.5, fontweight='bold')
ax.set_ylabel('Effect size (β)', fontsize=11)
ax.set_xlim(-0.6, len(labels) - 0.4)
top = res.hi.max() + 0.08
ax.set_ylim(ymin, top)
ax.grid(axis='y', color='#e5e5e5', lw=0.7)
ax.set_axisbelow(True)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker='o', ls='', color=FS, ms=7, label='Frameshift - significant'),
     Line2D([], [], marker='D', ls='', color=SG, ms=6.5, label='Stopgain - significant'),
     Line2D([], [], marker='o', ls='', color=FS_L, ms=7, label='Not significant after BH')]
ax.legend(handles=h, loc='upper right', frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(OUT + '.png', dpi=200, facecolor='white')
plt.savefig(OUT + '.pdf', facecolor='white')
print('saved', OUT + '.png')
