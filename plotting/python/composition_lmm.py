#!/usr/bin/env python3
"""
Regenerate the escape-region composition analyses (Figures 6A-D and 7) from the
new master table built from the corrected variant list.

Frameshift (Figure 6A-D): change relative to wild-type in the post-frameshift tail
    Δ = feature(PFSseq) − feature(WTPFSseq)
Stopgain (Figure 7): change between truncated mutant and full-length wild type
    Δ = feature(Sequence) − feature(WTSequence)

Model per feature:  Δ ~ is_case + (1 | transcript)   (statsmodels MixedLM, REML)
restricted to transcripts carrying both P/LP (disease) and gnomAD (control) variants.
Panel C / D features are standardised (SD units) before fitting; amino-acid
fractions are left on their natural scale. BH correction within each panel.
A random-slope sensitivity model  Δ ~ is_case + (1 + is_case | transcript)  is
reported alongside (columns *_rs).

Usage: python composition_lmm.py master_table_all_1.csv [min_tail_aa]
Outputs: comp_fs_*.csv, comp_sg_*.csv, per-variant feature cache comp_variants.csv
"""
import sys, warnings, math
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from collections import Counter

SRC = sys.argv[1] if len(sys.argv) > 1 else 'master_table_all_1.csv'
MIN_TAIL = int(sys.argv[2]) if len(sys.argv) > 2 else 2
AAS = 'ACDEFGHIKLMNPQRSTVWY'
AA_NAME = dict(A='Alanine', C='Cysteine', D='Aspartate', E='Glutamate', F='Phenylalanine',
               G='Glycine', H='Histidine', I='Isoleucine', K='Lysine', L='Leucine',
               M='Methionine', N='Asparagine', P='Proline', Q='Glutamine', R='Arginine',
               S='Serine', T='Threonine', V='Valine', W='Tryptophan', Y='Tyrosine')
CLASSES = {'Aromatic': 'FWY', 'Negatively charged': 'DE', 'Nonpolar aliphatic': 'AGILMV',
           'Polar uncharged': 'STCNQ', 'Positively charged': 'RHK', 'Special cases': 'P'}
AA_CLASS = {a: c for c, s in CLASSES.items() for a in s}
FEATS = [('Kappa (charge patterning)', 'kappa'), ('Disorder (metapredict)', 'disorder_score'),
         ('Frac. disorder-promoting', 'FracDisoPromoting'), ('pLDDT (metapredict)', 'pLDDT'),
         ('Frac. aromatic', 'frac_aromatic'), ('Mean net charge', 'MeanNetCharge'),
         ('Frac. acidic', 'frac_acidic'), ('Frac. negative', 'FracNeg'),
         ('Frac. basic', 'frac_basic'), ('Helix propensity', 'sec_struc_helix'),
         ('Frac. R+K', 'frac_RK'), ('Frac. positive', 'FracPos'),
         ('Frac. charged (FCR)', 'FracCharged')]

# ---- phase-separation metrics (same definitions as phase_sep_propensity.py)
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}


def _shannon(seq):
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in Counter(seq).values()) if n else 0.0


def _cover(flags, win, n):
    out = np.zeros(n, bool)
    for s in np.flatnonzero(flags):
        out[s:s + win] = True
    return out.mean()


def ps_metrics(seq):
    n = len(seq)
    if n == 0:
        return None
    c = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    inset = lambda s: np.isin(c, np.frombuffer(s.encode('ascii'), dtype=np.uint8))
    pm = inset('QNGSY')
    prion = pm.mean()
    if n < 100:
        prion_w = prion
    else:
        cs = np.concatenate([[0], np.cumsum(pm)])
        prion_w = (cs[100:] - cs[:-100]).max() / 100.0
    aro, arg = inset('YFW').mean(), (c == ord('R')).mean()
    stick = inset('YFWRQN').mean()
    pos = np.flatnonzero(inset('KRDE'))
    if n < 2 or len(pos) < 2:
        scd = 0.0
    else:
        q = np.where(np.isin(c[pos], np.frombuffer(b'KR', dtype=np.uint8)), 1.0, -1.0)
        d = np.sqrt(np.abs(pos[None, :] - pos[:, None]).astype(float))
        scd = np.triu((q[:, None] * q[None, :]) * d, 1).sum() / (n * n)
    if n < 40:
        lc = 0.0 if _shannon(seq) >= 2.5 else 1.0
    else:
        ind = np.array([(c == ord(a)) for a in AAS], dtype=np.int32)
        cs = np.concatenate([np.zeros((20, 1), np.int32), np.cumsum(ind, axis=1)], axis=1)
        p = (cs[:, 40:] - cs[:, :-40]) / 40.0
        with np.errstate(divide='ignore', invalid='ignore'):
            ent = -np.nansum(np.where(p > 0, p * np.log2(p), 0.0), axis=0)
        lc = _cover(ent < 2.5, 40, n)
    if n < 10:
        hyd = 0.0
    else:
        kd = np.array([KD.get(a, 0.0) for a in seq])
        avg = np.lib.stride_tricks.sliding_window_view(kd, 10).mean(axis=1)
        hyd = _cover(avg > 0.45, 10, n)
    return {'Prion-like (global)': prion, 'Prion-like (max window)': prion_w,
            'Pi-pi score': aro * arg * 100, 'Sticker fraction': stick,
            'SCD (charge decoration)': scd, 'Low-complexity fraction': lc,
            'Hydrophobic cluster frac': hyd}


def aa_frac(seq):
    n = len(seq)
    cnt = Counter(seq)
    return {a: cnt.get(a, 0) / n for a in AAS}


# ---------------------------------------------------------------- data
m = pd.read_csv(SRC, low_memory=False)
m['cls'] = np.where(m.group.str.startswith('fs'), 'fs', 'sg')
m['is_case'] = m.group.str.endswith('disease').astype(int)

rows = []
for r in m.itertuples(index=False):
    if r.cls == 'fs':
        mut, wt, mp, wp = r.PFSseq, r.WTPFSseq, 'PFS', 'WTPFS'
    else:
        mut, wt, mp, wp = r.Sequence, r.WTSequence, 'Fullpep', 'WTFullpep'
    mut = mut if isinstance(mut, str) else ''
    wt = wt if isinstance(wt, str) else ''
    if len(mut) < MIN_TAIL or len(wt) < MIN_TAIL:
        continue
    rec = {'cls': r.cls, 'is_case': r.is_case, 'transcript': r.transcript,
           'gene': r.gene_symbol, 'key': r.key, 'mut_len': len(mut), 'wt_len': len(wt)}
    fm, fw = aa_frac(mut), aa_frac(wt)
    for a in AAS:
        rec['aa_' + a] = fm[a] - fw[a]
    for lab, col in FEATS:
        rec['f_' + col] = getattr(r, mp + col) - getattr(r, wp + col)
    if r.cls == 'fs':
        pm, pw = ps_metrics(mut), ps_metrics(wt)
        for k in pm:
            rec['ps_' + k] = pm[k] - pw[k]
    rows.append(rec)
V = pd.DataFrame(rows)
V.to_csv('comp_variants.csv', index=False)


def shared(df):
    tb = df.groupby('transcript').is_case.agg(['min', 'max'])
    return df[df.transcript.isin(tb[(tb['min'] == 0) & (tb['max'] == 1)].index)]


def fit(df, col, standardise):
    x = df[['transcript', 'is_case', col]].dropna().rename(columns={col: 'y'})
    if standardise:
        sd = x.y.std()
        if not sd or not np.isfinite(sd):
            return dict(beta=np.nan, lo=np.nan, hi=np.nan, p=np.nan, n=len(x))
        x['y'] = (x.y - x.y.mean()) / sd
    out = dict(n=len(x), n_case=int(x.is_case.sum()), n_ctrl=int((1 - x.is_case).sum()))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            f = smf.mixedlm('y ~ is_case', x, groups=x['transcript']).fit(reml=True, method='lbfgs')
            ci = f.conf_int().loc['is_case']
            out.update(beta=f.params['is_case'], lo=ci[0], hi=ci[1], p=f.pvalues['is_case'])
        except Exception as e:
            out.update(beta=np.nan, lo=np.nan, hi=np.nan, p=np.nan)
        try:
            g = smf.mixedlm('y ~ is_case', x, groups=x['transcript'],
                            re_formula='~is_case').fit(reml=True, method='lbfgs')
            out.update(beta_rs=g.params['is_case'], p_rs=g.pvalues['is_case'])
        except Exception:
            out.update(beta_rs=np.nan, p_rs=np.nan)
    return out


def panel(df, cols, labels, standardise):
    res = []
    for c, lab in zip(cols, labels):
        r = fit(df, c, standardise); r['feature'] = lab; r['col'] = c
        res.append(r)
    res = pd.DataFrame(res)
    ok = res.p.notna()
    res['q'] = np.nan
    res.loc[ok, 'q'] = multipletests(res.loc[ok, 'p'], method='fdr_bh')[1]
    ok = res.p_rs.notna()
    res['q_rs'] = np.nan
    res.loc[ok, 'q_rs'] = multipletests(res.loc[ok, 'p_rs'], method='fdr_bh')[1]
    return res


for cls in ('fs', 'sg'):
    d = shared(V[V.cls == cls])
    print('\n==== %s: %d transcripts, %d P/LP vs %d gnomAD (tail >= %d aa)' % (
        cls, d.transcript.nunique(), d.is_case.sum(), (1 - d.is_case).sum(), MIN_TAIL))
    A = panel(d, ['aa_' + a for a in AAS], [a for a in AAS], False)
    A['name'] = A.feature.map(AA_NAME); A['cls_name'] = A.feature.map(AA_CLASS)
    A.to_csv('comp_%s_A_aa.csv' % cls, index=False)
    B = A.groupby('cls_name').agg(sum_beta=('beta', 'sum'), n_aa=('beta', 'size')).reset_index()
    B.to_csv('comp_%s_B_class.csv' % cls, index=False)
    C = panel(d, ['f_' + c for _, c in FEATS], [l for l, _ in FEATS], True)
    C.to_csv('comp_%s_C_feat.csv' % cls, index=False)
    print(A[['feature', 'beta', 'p', 'q', 'q_rs']].round(4).to_string(index=False))
    print(B.round(4).to_string(index=False))
    print(C[['feature', 'n', 'beta', 'p', 'q', 'beta_rs', 'q_rs']].round(4).to_string(index=False))
    if cls == 'fs':
        pcols = [c for c in V.columns if c.startswith('ps_')]
        D = panel(d, pcols, [c[3:] for c in pcols], False)
        D.to_csv('comp_fs_D_ps.csv', index=False)
        print(D[['feature', 'beta', 'p', 'q', 'beta_rs', 'q_rs']].round(4).to_string(index=False))
