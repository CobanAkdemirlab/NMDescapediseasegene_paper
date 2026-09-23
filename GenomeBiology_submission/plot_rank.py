#!/usr/bin/env python3
"""Per-gene ranking of the within-gene pathogenic−benign difference in Δ(mutant − WT),
from comp_variants.csv (composition_lmm.py output).  Genes with >= MIN variants per cohort.
Usage: python plot_rank.py comp_variants.csv fs|sg PS_proteins.txt out_prefix"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSV, CLS, PSF, OUT = sys.argv[1:5]
MIN = 2
CLASSES = {'Aromatic (FWY)': 'FWY', 'Negatively charged (DE)': 'DE',
           'Nonpolar aliphatic (AGILMV)': 'AGILMV', 'Polar uncharged (STCNQ)': 'STCNQ',
           'Positively charged (RHK)': 'RHK', 'Special cases (P)': 'P'}
TOP = {'fs': [('Kappa (charge patterning)', 'f_kappa'), ('Frac. charged (FCR)', 'f_FracCharged'),
              ('Frac. positive', 'f_FracPos'), ('Mean net charge', 'f_MeanNetCharge')],
       'sg': [('Frac. negative', 'f_FracNeg'), ('Arginine (R)', 'aa_R'),
              ('Disorder (metapredict)', 'f_disorder_score'), ('pLDDT (metapredict)', 'f_pLDDT')]}[CLS]

d = pd.read_csv(CSV)
d = d[d.cls == CLS].copy()
for lab, aas in CLASSES.items():
    d['cls_' + lab] = d[['aa_' + a for a in aas]].sum(axis=1)
FEATS = TOP + [(lab, 'cls_' + lab) for lab in CLASSES]
d['cohort'] = np.where(d.is_case == 1, 'disease', 'control')
tb = d.groupby(['gene', 'cohort']).size().unstack(fill_value=0)
keep = tb[(tb.disease >= MIN) & (tb.control >= MIN)].index
d = d[d.gene.isin(keep)]
ps = {l.strip().upper() for l in open(PSF) if l.strip() and l.strip().upper() != 'PROTEIN'}
print('%s: %d genes with >=%d variants per cohort; %d P/LP, %d gnomAD' % (
    CLS, d.gene.nunique(), MIN, (d.cohort == 'disease').sum(), (d.cohort == 'control').sum()))

fig, axes = plt.subplots(4, 3, figsize=(17, 20))
axes = axes.flatten()
summary = []
for ax, (title, col) in zip(axes, FEATS):
    g = d.groupby(['gene', 'cohort'])[col].mean().unstack()
    n = d.groupby(['gene', 'cohort']).size().unstack(fill_value=0)
    eff = (g['disease'] - g['control']).dropna()
    order = eff.sort_values()
    sel = pd.concat([order.head(10), order.tail(10)])
    sel = sel[~sel.index.duplicated()].sort_values()
    summary.append((title, ', '.join(order.tail(5).index[::-1]), ', '.join(order.head(5).index)))
    span = max(abs(sel).max(), 1e-9)
    ax.axvline(0, color='grey', lw=0.6)
    for i, (gn, v) in enumerate(sel.items()):
        ax.barh(i, v, color='#c0392b' if v < 0 else '#2166ac', height=0.72)
        ax.text(v + (span * 0.02 if v >= 0 else -span * 0.02), i,
                'n=%d/%d' % (n.loc[gn, 'disease'], n.loc[gn, 'control']),
                va='center', ha='left' if v >= 0 else 'right', fontsize=6.3, color='#333')
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels([('★ ' + g_ if g_.upper() in ps else g_) for g_ in sel.index], fontsize=7.5)
    for tk, g_ in zip(ax.get_yticklabels(), sel.index):
        if g_.upper() in ps:
            tk.set_color('#b8860b'); tk.set_fontweight('bold')
    ax.set_xlabel('%s\nwithin-gene Δ(mutant−WT), P/LP − gnomAD' % title, fontsize=8)
    ax.margins(x=0.24)
for ax in axes[len(FEATS):]:
    ax.axis('off')
plt.tight_layout()
plt.savefig(OUT + '.png', dpi=165, facecolor='white')
plt.savefig(OUT + '.pdf', facecolor='white')
pd.DataFrame(summary, columns=['feature', 'top5_higher_in_PLP', 'top5_lower_in_PLP']).to_csv(OUT + '_top.csv', index=False)
print(pd.DataFrame(summary, columns=['feature', 'higher', 'lower']).to_string(index=False))
