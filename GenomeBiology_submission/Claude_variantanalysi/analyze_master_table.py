#!/usr/bin/env python3
"""
analyze_master_table.py  --  20-amino-acid composition analysis on the output of
wt_vs_mutant_genomic.py, gene-matched LMM.  ONLY the 20 amino-acid fractions.

FRAMESHIFT comparison (region=PFS, the correct/default one):
    mutant PFS  = PFSseq   = mut_prot[fs_start:]  (novel out-of-frame tail -> new PTC)
    WT at PTC   = WTPFSseq  = wt_prot[fs_start:]   (WT protein from the same position)
    Both are anchored at the frameshift position; we model  delta = mutant - WT.

STEP 1 (needs Ensembl, run once on ALL variants):
    python wt_vs_mutant_genomic.py \
        --variants variants_all0901.csv --output master_table_all.csv --cache-dir ./vep_cache

STEP 2 (this script, no network):
    python analyze_master_table.py \
        --master master_table_all.csv --variants variants_all0901.csv \
        --region PFS --out-prefix fs_aa20_ALL
"""
import argparse, sys, numpy as np, pandas as pd, warnings
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from aa20_features import compute_aa20, AA20
warnings.simplefilter("ignore")

# region -> (mutant seq col, WT seq col or None, mutant len col, WT len col or None)
REGION_COLS = {
    "PFS":  ("PFSseq",   "WTPFSseq",   "PFSseqLength", "WTPFSseqLength"),  # frameshift novel tail vs WT-at-PTC
    "IDR":  ("IDRseq",   "WTIDRseq",   "IDRLength",    "WTIDRLength"),
    "Full": ("Sequence", "WTSequence", "FullLength",   "WTFullLength"),
    "MUT":  ("PFSseq",   None,          "PFSseqLength", None),             # absolute composition (e.g. stopgain lost tail)
}

def _len(seq): return len("".join(c for c in str(seq).upper() if c.isalpha()))

def load_join(master, variants):
    m = pd.read_csv(master, low_memory=False)
    v = pd.read_csv(variants)
    if "key" not in m.columns or "key" not in v.columns:
        sys.exit("FATAL: both --master and --variants need a 'key' column (chrom:pos|REF|ALT)")
    keep = ["key"] + [c for c in ("source","group","uniprot","uniprotswissprot") if c in v.columns]
    d = m.merge(v[keep], on="key", how="left", suffixes=("","_var"))
    src = (d["source"] if "source" in d.columns else d.get("group", pd.Series([""]*len(d)))).astype(str)
    d["cohort"] = np.where(src.isin(["fs","snv","fs_disease","snv_disease"]), "disease",
                    np.where(src.str.contains("control"), "control", ""))
    d.loc[d["cohort"]=="", "cohort"] = np.nan
    gene = None
    for cand in ("uniprot","uniprotswissprot","uniprotswissprot_var","gene_symbol"):
        if cand in d.columns: gene = d[cand]; break
    d["gene"] = gene.astype(str).str.split("-").str[0] if gene is not None else d.get("transcript","").astype(str)
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--region", default="PFS", choices=list(REGION_COLS))
    ap.add_argument("--min-len", type=int, default=20,
                    help="min length (aa) of the MUTANT region (paper filter: PFSseqLength>=20)")
    ap.add_argument("--both-min-len", action="store_true",
                    help="also require the WT region >= --min-len (default: only the mutant region)")
    ap.add_argument("--length-covariate", action="store_true",
                    help="add mutant & WT region length as fixed-effect covariates in the LMM")
    ap.add_argument("--out-prefix", default="aa20")
    args = ap.parse_args()

    d = load_join(args.master, args.variants)
    mut_col, wt_col, mlen_col, wlen_col = REGION_COLS[args.region]
    if mut_col not in d.columns:
        sys.exit(f"FATAL: '{mut_col}' not in master table. Seq cols present: "
                 f"{[c for c in d.columns if 'seq' in c.lower() or c in ('Sequence','WTSequence')]}")

    # explicit, self-documenting statement of the comparison
    print("="*70)
    print(f"REGION = {args.region}")
    if wt_col:
        print(f"  comparison : delta = mutant('{mut_col}') - WT('{wt_col}')  [both anchored at the same position]")
    else:
        print(f"  comparison : absolute composition of '{mut_col}' (no WT subtraction)")
    print("="*70)

    d = d.dropna(subset=["cohort","gene"]).copy()
    d["_mlen"] = d[mut_col].map(_len)
    d = d[d["_mlen"] >= args.min_len]
    if wt_col:
        d["_wlen"] = d[wt_col].map(_len)
        n_before = len(d)
        if args.both_min_len:
            d = d[d["_wlen"] >= args.min_len]
        print(f"  variants with mutant {mut_col} >= {args.min_len} aa: {n_before}")
        if args.both_min_len:
            print(f"  ... and WT {wt_col} >= {args.min_len} aa (both): {len(d)}")
        print(f"  (WT region runs to the natural stop, so it is usually longer than the mutant)")
    # --- cohort sanity check: fail loudly instead of drawing an empty plot ---
    nd0=int((d["cohort"]=="disease").sum()); nc0=int((d["cohort"]=="control").sum())
    print(f"  cohort after join+length filter: disease={nd0}  control={nc0}")
    src_counts = (d.get("source", pd.Series(dtype=str))).value_counts(dropna=False).to_dict()
    print(f"  source breakdown of used rows: {src_counts}")
    if nd0==0 or nc0==0:
        sys.exit("FATAL: one cohort is empty (disease=%d, control=%d).\n"
                 "  The master table is missing all variants of one class.\n"
                 "  For a frameshift PFS run you need source=fs (disease) AND source=fs_control (benign);\n"
                 "  check that both were in the pipeline input and that they produced a non-empty PFSseq."
                 % (nd0,nc0))
    d = d.copy(); d["cohort_n"] = (d["cohort"]=="disease").astype(int)

    MF = pd.DataFrame([compute_aa20(s) for s in d[mut_col]]); MF.index = d.index
    if wt_col:
        WF = pd.DataFrame([compute_aa20(s) for s in d[wt_col]]); WF.index = d.index
        for a in AA20: d["y_"+a] = MF["frac_"+a] - WF["frac_"+a]
        xl = "within-gene pathogenic\u2212benign \u0394(mutant\u2212WT) frac \u00b7 SD units"
    else:
        for a in AA20: d["y_"+a] = MF["frac_"+a]
        xl = "within-gene pathogenic\u2212benign composition \u00b7 SD units"

    import statsmodels.formula.api as smf, statsmodels.stats.multitest as mt
    cov = ""
    if args.length_covariate:
        d["_mlz"] = (d["_mlen"]-d["_mlen"].mean())/(d["_mlen"].std() or 1)
        cov = " + _mlz"
        if wt_col:
            d["_wlz"] = (d["_wlen"]-d["_wlen"].mean())/(d["_wlen"].std() or 1); cov += " + _wlz"
    def lmm(col):
        x = d.dropna(subset=[col]).copy(); sd = x[col].std() or 1; x["yy"]=x[col]/sd
        tb = x.groupby("gene")["cohort_n"].nunique(); x = x[x.gene.isin(tb[tb==2].index)]
        if x.gene.nunique() < 3: return np.nan,np.nan,np.nan
        try:
            mm = smf.mixedlm(f"yy ~ cohort_n{cov}", x, groups=x["gene"]).fit(reml=False, method="lbfgs")
            return mm.params["cohort_n"], mm.bse["cohort_n"], mm.pvalues["cohort_n"]
        except Exception: return np.nan,np.nan,np.nan
    R = pd.DataFrame([dict(aa=a, **dict(zip(["b","se","p"], lmm("y_"+a)))) for a in AA20])
    R["q"] = mt.multipletests(R.p.fillna(1), method="fdr_bh")[1]
    R = R.sort_values("b").reset_index(drop=True)
    R.to_csv(f"{args.out_prefix}_stats.csv", index=False)

    def sig(q): return "***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else "trend" if q<0.10 else "n.s."
    def colr(b,q):
        if pd.isna(q) or q>=0.10: return "#c0c0c0"
        hi=b>0; base=("#b2182b","#ef8a62","#fddbc7") if hi else ("#2166ac","#67a9cf","#d1e5f0")
        return base[0] if q<0.01 else base[1] if q<0.05 else base[2]
    fig, ax = plt.subplots(figsize=(8,9))
    for i,r in R.iterrows():
        se = r.se if pd.notna(r.se) else 0
        ax.errorbar(r.b, i, xerr=1.96*se, fmt="o", color=colr(r.b,r.q), ms=9, lw=2.3, mec="black", mew=0.5)
        lbl = ("\u25b2 " if r.b>0 else "\u25bc ")+(f"q={r.q:.3f} {sig(r.q)}" if pd.notna(r.q) else "n/a")
        ax.text(r.b+1.96*se+0.03, i, lbl, va="center", fontsize=9)
    ax.axvline(0, color="grey", lw=0.8); ax.set_yticks(range(len(R))); ax.set_yticklabels([f"frac_{a}" for a in R.aa], fontsize=10)
    ax.set_ylim(-0.6, len(R)-0.4); ax.margins(x=0.34); ax.set_xlabel(xl, fontsize=10)
    ng=d.gene.nunique(); nd=int((d.cohort=='disease').sum()); nc=int((d.cohort=='control').sum())
    ttl=f"20-AA composition, gene-matched LMM \u2014 region={args.region}"+(" (mutant\u2212WT)" if wt_col else " (absolute)")
    ax.set_title(f"{ttl}\n{len(d)} variants, {ng} genes, {nd} disease / {nc} control"+(", length-adjusted" if args.length_covariate else ""), fontsize=11, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.png", dpi=180, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{args.out_prefix}.pdf", bbox_inches="tight", facecolor="white")
    print(f"\nvariants used: {len(d)}  genes {ng}  disease {nd} / control {nc}")
    print("significant AAs (q<0.05):", list(R[R.q<0.05].aa))
    print(f"wrote {args.out_prefix}.png/.pdf and {args.out_prefix}_stats.csv")

if __name__ == "__main__":
    main()
