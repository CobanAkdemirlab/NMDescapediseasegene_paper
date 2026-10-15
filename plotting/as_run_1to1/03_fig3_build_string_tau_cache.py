# Figure 3 panels F-G: build ppi_tau_4groups.csv for updated_builder_ppitau.R.
# degree_centrality in the gene table is all-NA (the repo's STRING download fallback),
# so degree is counted from the local STRING v11.5 edge list at combined_score >= 400
# and tau is read from the symbol-keyed GTEx V8 derivative.
import pandas as pd
from collections import Counter

F3 = "/Users/qkelly/Desktop/NMDesc/plotting_run_20260923/reproduced/Figure3"
pairs = pd.read_csv("matched_pairs_full.csv")
genes = sorted(set(pairs.case_hgnc.dropna()) | set(pairs.control_hgnc.dropna()))

info = pd.read_csv(f"{F3}/9606.protein.info.v11.5.txt.gz", sep="\t", usecols=[0, 1],
                   names=["string_id", "preferred_name"], header=0)
need = {g: dict(zip(info.preferred_name, info.string_id)).get(g) for g in genes}
unmapped = [g for g, v in need.items() if v is None]
if unmapped:                      # alias fallback; empty on the set1+set4 union (416/416)
    al = pd.read_csv(f"{F3}/9606.protein.aliases.v11.5.txt.gz", sep="\t", usecols=[0, 1],
                     names=["string_id", "alias"], header=0)
    a2s = dict(zip(al.alias, al.string_id))
    for g in unmapped:
        if g in a2s: need[g] = a2s[g]

sids = {v for v in need.values() if v}
deg = Counter()
for ch in pd.read_csv(f"{F3}/9606.protein.links.v11.5.txt.gz", sep=" ", chunksize=2_000_000):
    ch = ch[ch.combined_score >= 400]
    deg.update(ch.loc[ch.protein1.isin(sids), "protein1"])

tau = pd.read_csv(f"{F3}/Tau_gene_V8_derived.csv")
tau_lookup = dict(zip(tau.gene, tau.tau))

rows = []
for grp, col, sub in [("Frameshift", "case_hgnc", pairs[pairs.case_group == "fs"]),
                      ("FS Control", "control_hgnc", pairs[pairs.case_group == "fs"]),
                      ("Stopgain", "case_hgnc", pairs[pairs.case_group == "snv"]),
                      ("Stopgain Control", "control_hgnc", pairs[pairs.case_group == "snv"])]:
    for g in sorted(sub[col].dropna().unique()):
        rows.append(dict(gene=g, Group=grp, degree=deg.get(need.get(g)), tau=tau_lookup.get(g)))
cache = pd.DataFrame(rows)
cache.to_csv("ppi_tau_4groups.csv", index=False)
# report per-arm missingness next to the panels; it was 0 in all four arms here
print(cache.groupby("Group").agg(n=("gene", "size"), degree_na=("degree", lambda s: s.isna().sum()),
                                 tau_na=("tau", lambda s: s.isna().sum()),
                                 degree_median=("degree", "median")).to_string())
