"""Figure S5 - overlap of NMD-escape regions with UniProt IDRs (>=20 aa).
usage: python figS5_idr_overlap.py ROWS_CSV STATS_CSV OUT_STEM "SUBTITLE"
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

ROWS, STATS, STEM, SUB = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
rows = pd.read_csv(ROWS); st = pd.read_csv(STATS)

METRICS = [("n_idr", "N IDRs\noverlapping", "n_idr"),
           ("tot_aa", "Total IDR\noverlap (aa)", "tot_aa"),
           ("frac", "IDR overlap\nfraction", "frac"),
           ("any_idr", "Any IDR\noverlap", "any_idr"),
           ("longest", "Longest IDR\noverlap (aa)", "longest")]
STRATA = [("FS", "Frameshift", "FS Control", "#4472A8", "#C6DBEF"),
          ("SNV", "Stopgain", "Stopgain Control", "#C2405A", "#F7C6AD")]

def stars(q):
    if not np.isfinite(q): return "n/a"
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < .05 else "ns"

def violin(ax, x, v, face, edge="black"):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) == 0: return
    lo, hi = v.min(), v.max()
    nu = len(np.unique(v))
    if hi > lo and nu >= 2:
        g = np.linspace(lo, hi, 256)
        bw = 0.30 if nu <= 3 else None          # narrow bw keeps binary lobes separated
        try:
            w = gaussian_kde(v, bw_method=bw)(g)
        except Exception:
            w = np.ones_like(g)
    else:
        g = np.array([lo - 1e-9, lo + 1e-9]); w = np.ones_like(g)
    w = w / w.max() * 0.34
    ax.fill_betweenx(g, x - w, x + w, facecolor=face, edgecolor=edge, linewidth=.7, zorder=2)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    if q3 > q1:
        ax.add_patch(plt.Rectangle((x - .105, q1), .21, q3 - q1, facecolor=face,
                                   edgecolor="black", linewidth=.9, zorder=3))
    ax.plot([x - .155, x + .155], [med, med], color="black", lw=2.0, zorder=4,
            solid_capstyle="butt")
    ax.scatter(np.full(len(v), x) + np.random.default_rng(0).normal(0, .045, len(v)), v,
               s=2.2, color="#3a3a3a", alpha=.35, linewidths=0, zorder=5)

fig, axes = plt.subplots(2, 5, figsize=(16.4, 7.3))
ylims = {}
for m, _, _ in METRICS:
    v = rows[m].dropna().values
    hi = 1.0 if m in ("frac", "any_idr") else np.percentile(v, 99.9) * 1.12 + 1
    ylims[m] = (-0.04 * hi, hi)

for ri, (stt, lab_c, lab_k, col_c, col_k) in enumerate(STRATA):
    d = rows[rows.stratum == stt]
    for ci, (m, title, _) in enumerate(METRICS):
        ax = axes[ri, ci]
        violin(ax, 0, d.loc[d.role == "case", m], col_c)
        violin(ax, 1, d.loc[d.role == "control", m], col_k)
        s = st[(st.stratum == stt) & (st.metric == m)].iloc[0]
        lo, hi = ylims[m]; ax.set_ylim(lo, hi)
        yb = hi * 0.845
        ax.plot([0, 1], [yb, yb], color="#333333", lw=.9, zorder=6)
        sg = stars(s["q"])
        sig = sg != "ns"
        ax.text(.5, yb + hi * .012, sg, ha="center", va="bottom", zorder=6,
                fontsize=13 if sig else 9.5, fontweight="bold" if sig else "normal",
                color="#B3202C" if sig else "#555555")
        txt = (f"n = {int(s['n_pairs'])} pairs\n"
               f"$\\beta$ = {s['beta']:+.3f}\n"
               f"p = {'<0.001' if s['p'] < 1e-3 else f'{s[chr(112)]:.3f}'}\n"
               f"q = {'<0.001' if s['q'] < 1e-3 else f'{s[chr(113)]:.3f}'}")
        ax.text(.035, .975, txt, transform=ax.transAxes, ha="left", va="top", fontsize=6.6,
                color="#B3202C" if sig else "#606060", linespacing=1.28,
                bbox=dict(boxstyle="round,pad=.28", fc="white",
                          ec="#B3202C" if sig else "#B8B8B8", lw=.8, alpha=.92), zorder=7)
        ax.set_xticks([0, 1])
        n_c = int((d.role == "case")[d[m].notna()].sum())
        n_k = int((d.role == "control")[d[m].notna()].sum())
        ax.set_xticklabels([f"{lab_c}\n(n={n_c})", f"{lab_k}\n(n={n_k})"], fontsize=8.2)
        ax.set_xlim(-.62, 1.62)
        ax.set_ylabel(m if ci else m, fontsize=8.5)
        if ri == 0: ax.set_title(title, fontsize=10, fontweight="bold", pad=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="y", color="#E8E8E8", lw=.6); ax.set_axisbelow(True)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)

fig.suptitle("IDR Overlap Features \u2014 NMD-Escape Region (UniProt IDRs \u2265 20 aa)",
             fontsize=14, fontweight="bold", y=.982)
fig.text(.5, .944, "LME: outcome ~ is_case + (1 | pair_id) \u00b7 BH-adjusted per mutation group "
         "(5 metrics each) \u00b7 all 416 accessions retrieved from UniProt; the 113 with no Region >= 20 aa contribute a true 0 \u00b7 "
         "top row: frameshift pairs \u00b7 bottom row: stopgain pairs",
         ha="center", fontsize=7.6, color="#444444")
fig.text(.5, .924, SUB, ha="center", fontsize=8.2, color="#222222")
fig.legend(handles=[Patch(fc=c, ec="black", lw=.6, label=l) for l, c in
                    [("Frameshift (case)", "#4472A8"), ("FS Control", "#C6DBEF"),
                     ("Stopgain (case)", "#C2405A"), ("Stopgain Control", "#F7C6AD")]],
           loc="lower center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(.5, .003))
fig.text(.5, -.004, "Stars from BH-adjusted q within mutation group; "
         "*** q<0.001, ** q<0.01, * q<0.05. Unadjusted p also shown.",
         ha="center", fontsize=7, color="#666666")
fig.tight_layout(rect=[0, .045, 1, .915])
fig.savefig(f"{STEM}.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{STEM}.pdf", bbox_inches="tight")
print("wrote", STEM)
