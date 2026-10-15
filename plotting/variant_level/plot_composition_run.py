#!/usr/bin/env python3
"""Plot Figure 6A-D (frameshift, tails >= 20 aa) and Figure 7A-C (stopgain) from the
CSV output of composition_lmm.py. Panel letters only, no panel titles."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

D = sys.argv[1] if len(sys.argv) > 1 else '.'
CLS_COL = {'Aromatic': '#7b4fa0', 'Negatively charged': '#c0392b', 'Nonpolar aliphatic': '#43a047',
           'Polar uncharged': '#f0896a', 'Positively charged': '#2166ac', 'Special cases': '#8c6d31'}
FAM_COL = {'Hydrophobic clusters': '#c0392b', 'Aromatic / cation-π': '#1a9c8c',
           'Charge patterning': '#3b6ea8', 'Sequence complexity': '#e8a33d',
           'Sticker / prion-like': '#7b4fa0'}
PS_FAM = {'Prion-like (global)': 'Sticker / prion-like', 'Prion-like (max window)': 'Sticker / prion-like',
          'Sticker fraction': 'Sticker / prion-like', 'Pi-pi score': 'Aromatic / cation-π',
          'SCD (charge decoration)': 'Charge patterning', 'Low-complexity fraction': 'Sequence complexity',
          'Hydrophobic cluster frac': 'Hydrophobic clusters'}
UP, DOWN, TREND_UP, TREND_DOWN, NS = '#b2182b', '#2166ac', '#f4b6a6', '#a8c8e4', '#b5b5b5'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9})


def star(q):
    return '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'


def ptxt(p, q):
    return 'p=%s %s (BH=%s)' % (('%.2g' % p), star(p) if p < 0.05 else 'ns', ('%.2g' % q))


def letter(ax, s, dx=-0.02, dy=1.01):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=15, fontweight='bold', va='bottom', ha='right')


def panel_aa(ax, A, xlabel):
    A = A.sort_values('beta')
    y = np.arange(len(A))
    ax.barh(y, A.beta, color=[CLS_COL[c] for c in A.cls_name], edgecolor='#333', lw=0.5, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(['%s  %s' % (f, n) for f, n in zip(A.feature, A.name)], fontsize=9)
    ax.axvline(0, color='#555', lw=0.7)
    lim = np.abs(A.beta).max() * 1.15
    ax.set_xlim(-lim, lim)
    for yi, (p, q) in zip(y, zip(A.p, A.q)):
        ax.text(1.02, yi, ptxt(p, q), transform=ax.get_yaxis_transform(), fontsize=7,
                va='center', color='#222' if q < 0.05 else '#666',
                fontweight='bold' if q < 0.05 else 'normal')
    ax.set_xlabel(xlabel)
    ax.legend(handles=[Patch(color=c, label=k) for k, c in CLS_COL.items()], title='Property',
              loc='lower right', fontsize=7.5, title_fontsize=8, frameon=True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def panel_class(ax, B):
    order = ['Aromatic', 'Negatively charged', 'Nonpolar aliphatic', 'Polar uncharged',
             'Positively charged', 'Special cases']
    B = B.set_index('cls_name').loc[order]
    y = np.arange(len(B))
    ax.barh(y, B.sum_beta, color=[CLS_COL[c] for c in order], edgecolor='#333', lw=0.5, height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.axvline(0, color='#555', lw=0.7)
    lim = np.abs(B.sum_beta).max() * 2.3
    ax.set_xlim(-lim, lim)
    for yi, (v, n) in zip(y, zip(B.sum_beta, B.n_aa)):
        ax.text(v + (lim * 0.03 if v >= 0 else -lim * 0.03), yi, 'Σβ=%+.4f (%d aa)' % (v, n),
                va='center', ha='left' if v >= 0 else 'right', fontsize=7.5)
    ax.set_xlabel('Σ β over amino acids in class')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def panel_forest(ax, C, xlabel):
    C = C.sort_values('beta')
    y = np.arange(len(C))
    for yi, r in zip(y, C.itertuples()):
        if r.q < 0.05:
            col = UP if r.beta > 0 else DOWN
        elif r.q < 0.1:
            col = TREND_UP if r.beta > 0 else TREND_DOWN
        else:
            col = NS
        ax.plot([r.lo, r.hi], [yi, yi], color=col, lw=2, solid_capstyle='round')
        ax.plot(r.beta, yi, 'o', ms=8, color=col, mec='#333', mew=0.6, zorder=3)
        lab = ('▲' if r.beta > 0 else '▼') + ' q=%.3f %s' % (r.q, star(r.q) if r.q < 0.05 else
                                                        ('trend' if r.q < 0.1 else 'n.s.'))
        ax.text(r.hi + 0.02, yi, lab, va='center', fontsize=7.5)
    ax.set_yticks(y); ax.set_yticklabels(C.feature)
    ax.axvline(0, color='#999', lw=0.8)
    m = max(np.nanmax(np.abs(C[['lo', 'hi']].values)), 0.2)
    ax.set_xlim(-m * 1.15, m * 1.75)
    ax.set_xlabel(xlabel)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def panel_ps(ax, Dp):
    Dp = Dp.sort_values('beta')
    y = np.arange(len(Dp))
    ax.barh(y, Dp.beta, color=[FAM_COL[PS_FAM[f]] for f in Dp.feature], edgecolor='#333',
            lw=0.5, height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(Dp.feature)
    ax.axvline(0, color='#555', lw=0.7)
    lim = np.abs(Dp.beta).max() * 1.2
    ax.set_xlim(-lim, lim)
    for yi, (p, q) in zip(y, zip(Dp.p, Dp.q)):
        ax.text(1.02, yi, ptxt(p, q), transform=ax.get_yaxis_transform(), fontsize=7.5,
                va='center', color='#222' if q < 0.05 else '#666',
                fontweight='bold' if q < 0.05 else 'normal')
    ax.set_xlabel('Mixed-model β (P/LP − gnomAD) for Δ propensity metric')
    ax.legend(handles=[Patch(color=c, label=k) for k, c in FAM_COL.items()], title='Mechanism family',
              loc='lower right', fontsize=7.5, title_fontsize=8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def load(tag):
    A = pd.read_csv('%s/comp_%s_A_aa.csv' % (D, tag))
    B = pd.read_csv('%s/comp_%s_B_class.csv' % (D, tag))
    C = pd.read_csv('%s/comp_%s_C_feat.csv' % (D, tag))
    return A, B, C


# ---------------------------------------------------------------- Figure 6 A-D
A, B, C = load('fs')
Dp = pd.read_csv('%s/comp_fs_D_ps.csv' % D)
fig = plt.figure(figsize=(14, 13.5))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.15, 1.05], width_ratios=[1.05, 1],
                      hspace=0.42, wspace=0.85, left=0.13, right=0.97, top=0.97, bottom=0.05)
axA = fig.add_subplot(gs[0:2, 0]); panel_aa(axA, A, 'Mixed-model β (P/LP − gnomAD) for ΔPFS amino-acid fraction'); letter(axA, 'A', -0.14)
axB = fig.add_subplot(gs[0, 1]); panel_class(axB, B); letter(axB, 'B', -0.42, 1.02)
axC = fig.add_subplot(gs[1, 1]); panel_forest(axC, C, 'Mixed-model β (P/LP − gnomAD), SD units'); letter(axC, 'C', -0.52, 1.02)
sub = gs[2, :].subgridspec(1, 3, width_ratios=[0.22, 1, 0.45])
axD = fig.add_subplot(sub[0, 1]); panel_ps(axD, Dp); letter(axD, 'D', -0.33)
fig.savefig('fig6_ABCD.png', dpi=200, facecolor='white')
plt.close(fig)

# ---------------------------------------------------------------- Figure 7 A-C
A, B, C = load('sg')
fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.2], width_ratios=[1.05, 1],
                      hspace=0.35, wspace=0.85, left=0.13, right=0.97, top=0.96, bottom=0.07)
axA = fig.add_subplot(gs[:, 0]); panel_aa(axA, A, 'Mixed-model β (P/LP − gnomAD) for ΔFull amino-acid fraction'); letter(axA, 'A', -0.14)
axB = fig.add_subplot(gs[0, 1]); panel_class(axB, B); letter(axB, 'B', -0.42, 1.02)
axC = fig.add_subplot(gs[1, 1]); panel_forest(axC, C, 'Mixed-model β (P/LP − gnomAD), SD units'); letter(axC, 'C', -0.52, 1.02)
fig.savefig('fig7_new.png', dpi=200, facecolor='white')
fig.savefig('fig7_new.pdf', facecolor='white')
print('saved')
