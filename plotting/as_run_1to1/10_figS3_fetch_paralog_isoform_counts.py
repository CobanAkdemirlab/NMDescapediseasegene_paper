# Figure S3 panels A-B: paralogue and protein-coding isoform counts for EVERY gene in the
# 1:1 design, from one source (Ensembl REST), so the control arm is not topped up from a
# second source. Writes paralog_counts.csv / isoform_counts.csv for the two comparison
# scripts, which drop a pair when either member has no count (they do not score it 0).
import json, pandas as pd, numpy as np, requests
from concurrent.futures import ThreadPoolExecutor

H = {"Content-Type": "application/json", "Accept": "application/json"}
pairs = pd.read_csv("matched_pairs_full.csv")
genes = sorted(set(pairs.case_hgnc.dropna()) | set(pairs.control_hgnc.dropna()))

def lookup_batch(b):                     # batches of 20; 100 with expand=1 times out
    for _ in range(4):
        try:
            r = requests.post("https://rest.ensembl.org/lookup/symbol/homo_sapiens",
                              headers=H, data=json.dumps({"symbols": b, "expand": 1}), timeout=120)
            if r.status_code == 200: return r.json()
        except Exception: pass
    return {}

batches = [genes[i:i + 20] for i in range(0, len(genes), 20)]
iso = {}
with ThreadPoolExecutor(4) as ex:
    for b, r in zip(batches, ex.map(lookup_batch, batches)):
        for g in b:
            e = r.get(g)
            iso[g] = sum(t["biotype"] == "protein_coding" for t in e["Transcript"]) if e and e.get("Transcript") else np.nan

def para(g):
    for _ in range(4):
        try:
            rr = requests.get(f"https://rest.ensembl.org/homology/symbol/human/{g}",
                              headers={"Accept": "application/json"},
                              params={"type": "paralogues", "format": "condensed", "target_taxon": 9606},
                              timeout=120)
            if rr.status_code == 200:
                d = rr.json()["data"]
                return g, (len(d[0]["homologies"]) if d else 0)
            if rr.status_code == 404: return g, 0
        except Exception: pass
    return g, np.nan

with ThreadPoolExecutor(5) as ex:
    pl = dict(ex.map(para, genes))

ens = pd.DataFrame({"hgnc_symbol": genes, "n_isoform": [iso[g] for g in genes],
                    "n_paralog": [pl[g] for g in genes]})
ens.to_csv("ensembl_paralog_isoform_counts.csv", index=False)
# the comparison scripts read one column each; a gene with no count is simply absent
ens.dropna(subset=["n_paralog"])[["hgnc_symbol", "n_paralog"]].to_csv("paralog_counts.csv", index=False)
ens.dropna(subset=["n_isoform"])[["hgnc_symbol", "n_isoform"]].to_csv("isoform_counts.csv", index=False)
print("paralog missing:", int(ens.n_paralog.isna().sum()), "| isoform missing:", int(ens.n_isoform.isna().sum()))
