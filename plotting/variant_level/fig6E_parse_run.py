#!/usr/bin/env python3
"""
Figure 6E -- escape-tail phase-separation (ParSe v2), recomputed from the
master-table peptides.

Variant level: ClinVar P/LP vs gnomAD benign NMD-escape frameshift variants
WITHIN the same disease genes, so no gene-level matching ratio applies.

Pipeline (the convention of the 2026-10-05 set 4 recompute, reproduced here so
both designs sit on the same footing):
  1. frameshift rows only; mutant tail PFSseq vs its wild-type counterpart
     WTPFSseq, both required to be >= MIN_TAIL (20 aa).
  2. ParSe v2 per peptide (parse_v2_port, a line-for-line port of
     stevewhitten/Parse_v2_FASTA at WINDOW_SIZE=25, REGION_LENGTH=20,
     CUTOFF=0.9) -> p_sum ("sum of classifier distance over P residues") and
     longest_ps_idr (longest PS IDR, reported only when >= 20 aa).
  3. Length-normalised deltas, mutant minus wild type:
       cd_delta = var_psum / var_len - wt_psum / wt_len
       ps_delta = var_idr  / var_len - wt_idr  / wt_len
  4. Fit on genes that carry BOTH arms (a gene with only disease or only gnomAD
     variants cannot inform a within-gene contrast):
       y ~ cohort_n + (1 | gene),  y = delta / sd(delta over the FULL fs frame)
     ML (reml=False), lbfgs. BH over the two metrics.

Usage: python fig6E_parse_run.py <master_table.csv> <design_label> <out_png>
"""
import sys, os, warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_v2_port import parse_v2

SRC = sys.argv[1]
DESIGN = sys.argv[2]
OUT_PNG = sys.argv[3]
OUT_PDF = OUT_PNG.replace(".png", ".pdf")
OUT_STATS = OUT_PNG.replace(".png", "_stats.csv")
OUT_VARIANTS = OUT_PNG.replace(".png", "_per_variant.csv")
MIN_TAIL = 20

m = pd.read_csv(SRC, low_memory=False)
m = m[m.group.astype(str).str.startswith("fs")].copy()
m["cohort_n"] = m.group.astype(str).str.endswith("disease").astype(int)
m["gene"] = m.gene_symbol.astype(str)

rows = []
for r in m.itertuples(index=False):
    mut = r.PFSseq if isinstance(r.PFSseq, str) else ""
    wt = r.WTPFSseq if isinstance(r.WTPFSseq, str) else ""
    mut = "".join(ch for ch in mut.upper() if ch.isalpha())
    wt = "".join(ch for ch in wt.upper() if ch.isalpha())
    if len(mut) < MIN_TAIL or len(wt) < MIN_TAIL:
        continue
    vp, vi = parse_v2(mut)
    wp, wi = parse_v2(wt)
    if vp is None or wp is None:
        continue
    rows.append({"key": r.key, "gene": r.gene, "transcript": r.transcript,
                 "cohort_n": r.cohort_n, "var_len": len(mut), "wt_len": len(wt),
                 "var_psum": vp, "var_idr": vi, "wt_psum": wp, "wt_idr": wi,
                 "cd_delta": vp / len(mut) - wp / len(wt),
                 "ps_delta": vi / len(mut) - wi / len(wt)})
V = pd.DataFrame(rows)
V.to_csv(OUT_VARIANTS, index=False)

arms = V.groupby("gene").cohort_n.agg(["min", "max"])
both = set(arms[(arms["min"] == 0) & (arms["max"] == 1)].index)
F = V[V.gene.isin(both)].copy()

print("%s: %d fs variants with both tails >= %d aa, %d genes; "
      "fitted on %d variants in %d genes carrying both arms (%d P/LP, %d gnomAD)"
      % (DESIGN, len(V), MIN_TAIL, V.gene.nunique(), len(F), F.gene.nunique(),
         int(F.cohort_n.sum()), int((1 - F.cohort_n).sum())))


def fit(col):
    sd = V[col].std() or 1.0          # sd from the FULL fs frame, not the fitted subset
    x = F.dropna(subset=[col]).copy()
    x["y"] = x[col] / sd
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = smf.mixedlm("y ~ cohort_n", x, groups=x["gene"]).fit(reml=False, method="lbfgs")
    b, se, p = f.params["cohort_n"], f.bse["cohort_n"], f.pvalues["cohort_n"]
    return dict(beta=b, se=se, lo=b - 1.96 * se, hi=b + 1.96 * se, p=p,
                n=len(x), genes=x.gene.nunique())


LABELS = [("Classifier / residue", "cd_delta"), ("PS-IDR fraction", "ps_delta")]
S = pd.DataFrame([dict(design=DESIGN, feature=lab, col=col, **fit(col)) for lab, col in LABELS])
S["q"] = multipletests(S.p, method="fdr_bh")[1]
S.to_csv(OUT_STATS, index=False)
print(S[["feature", "beta", "lo", "hi", "p", "q", "n", "genes"]].to_string(index=False))

# ---- forest -----------------------------------------------------------------
INK = "#1a1a1a"
fig, ax = plt.subplots(figsize=(9, 3.4))
y = np.arange(len(S))[::-1]
sig = S.q < 0.05
ax.axvline(0, color="grey", lw=0.8, ls="--", zorder=1)
ax.hlines(y, S.lo, S.hi, color=INK, lw=1.6, zorder=2)
ax.scatter(S.beta, y, s=110, zorder=3, edgecolor=INK, linewidth=0.8,
           c=np.where(sig, "#2166AC", "white"))
for yy, (_, r) in zip(y, S.iterrows()):
    qtxt = "  q = %.3f" % r.q if r.q >= 1e-3 else "  q = %.1e" % r.q
    ax.text(ax.get_xlim()[1], yy, qtxt, va="center", fontsize=11,
            color=INK, fontweight="bold" if r.q < 0.05 else "normal")
ax.set_yticks(y)
ax.set_yticklabels(S.feature, fontsize=12, color=INK)
ax.set_xlabel("Standardised effect of disease vs gnomAD (SD units)", fontsize=12,
              fontweight="bold", color=INK)
ax.set_title("Escape-tail phase separation (ParSe v2) \u2014 %s" % DESIGN,
             fontsize=13, fontweight="bold", color=INK)
ax.text(0.5, -0.38, "Frameshift tails, mutant minus wild type \u00b7 "
        "value ~ cohort + (1 | gene) \u00b7 n = %d variants in %d genes (%d P/LP vs %d gnomAD)"
        % (S.n.iloc[0], S.genes.iloc[0], int(F.cohort_n.sum()), int((1 - F.cohort_n).sum())),
        transform=ax.transAxes, ha="center", va="top", fontsize=10, color="#444444")
ax.tick_params(colors=INK, labelsize=11)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
print("Figure -> %s / %s\nStats  -> %s" % (OUT_PNG, OUT_PDF, OUT_STATS))
