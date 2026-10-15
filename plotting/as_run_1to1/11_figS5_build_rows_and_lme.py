# Figure S5: per-gene overlap between the NMD-escape region and UniProt disordered Regions
# (>= 20 aa), then outcome ~ is_case + (1 | pair_id) per stratum, BH across the 5 metrics.
# The catalogue is the SYMMETRIC rebuild: every accession in all four arms was queried and
# returned, so a protein with no qualifying Region contributes a true 0 - it is not imputed.
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

idr = pd.read_csv("uniprot_idrs.csv")            # columns: ID, Start, End (protein coords)
reg = {}
for acc, s, e in idr.itertuples(index=False):
    reg.setdefault(acc, []).append((int(s), int(e)))

g = pd.read_csv("gene_flags_per_gene.csv")
g = g[g.design == "matched"].copy()

def metrics(row):
    if pd.isna(row.uniprot) or pd.isna(row.NMD_region_start) or pd.isna(row.NMD_region_end):
        return pd.Series({k: np.nan for k in ("n_idr", "tot_aa", "frac", "any_idr", "longest")})
    a0 = int(np.ceil(row.NMD_region_start / 3))  # CDS nt -> protein aa
    a1 = int(np.floor(row.NMD_region_end / 3))
    L = max(a1 - a0 + 1, 0)
    ov = [o for s, e in reg.get(row.uniprot, []) if (o := min(e, a1) - max(s, a0) + 1) > 0]
    return pd.Series({"n_idr": len(ov), "tot_aa": sum(ov), "frac": (sum(ov) / L if L else np.nan),
                      "any_idr": int(len(ov) > 0), "longest": max(ov) if ov else 0})

rows = pd.concat([g, g.apply(metrics, axis=1)], axis=1)
rows["is_case"] = (rows.role == "case").astype(int)
rows.to_csv("FigureS5_rows_1to1.csv", index=False)

METRICS = [("n_idr", "N IDRs\\noverlapping"), ("tot_aa", "Total IDR\\noverlap (aa)"),
           ("frac", "IDR overlap\\nfraction"), ("any_idr", "Any IDR\\noverlap"),
           ("longest", "Longest IDR\\noverlap (aa)")]
out = []
for strat in ("FS", "SNV"):
    sub = rows[rows.stratum == strat].dropna(subset=["n_idr"])
    recs = []
    for m, lab in METRICS:
        fit = smf.mixedlm(f"{m} ~ is_case", sub, groups=sub.pair_id).fit(reml=True)
        recs.append(dict(stratum=strat, metric=m, label=lab, n_pairs=sub.pair_id.nunique(),
                         n_rows=len(sub), beta=fit.params["is_case"], p=fit.pvalues["is_case"],
                         mean_case=sub.loc[sub.is_case == 1, m].mean(),
                         mean_ctrl=sub.loc[sub.is_case == 0, m].mean()))
    r = pd.DataFrame(recs); r["q"] = multipletests(r.p, method="fdr_bh")[1]; out.append(r)
stats = pd.concat(out, ignore_index=True)
stats.insert(0, "design", "set4")
stats.to_csv("FigureS5_stats_1to1.csv", index=False)
print(stats.round(4).to_string(index=False))
