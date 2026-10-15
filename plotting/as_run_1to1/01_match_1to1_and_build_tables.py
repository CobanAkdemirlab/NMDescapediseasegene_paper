# Recovered verbatim from the artifact lineage of matched_pairs_full_set4_1to1.csv.
# 1:1 length matching (CDS length + NMDesc/CDS ratio, rectangular assignment) and the
# three canonical tables: gene_flags_per_gene.csv, matched_pairs_full.csv, gene_all_AD.csv.
# Gene-level Pfam uses the house rule max_pfam_end * 3 >= NMD_region_start, NA preserved.

import os
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import ranksums

# ── skill: draw-4-plots helpers ──────────────────────────────────────────────
REPO_TO_PLOT_FLAGS = {
    "pfam_overlap_flag": "Pfam_flag",
    "ppi_overlap": "PPI_flag",
    "gene_protein_flag": "Protein_feature_flag",
    "gene_domains_flag": "Structured_domain_flag",
    "gene_slim_flag": "SLiM_flag",
    "gene_morf_flag": "MoRF_flag",
    "gene_ptm_flag": "PTM_flag",
    "gene_nls_flag": "NLS_flag",
    "gene_LCS_flag": "LCS_flag",
}

PLOT_FLAGCOLS = ["Pfam_flag", "Protein_feature_flag", "Structured_domain_flag",
                 "SLiM_flag", "MoRF_flag", "PTM_flag", "NLS_flag", "LCS_flag", "PPI_flag"]

EXTRA_GENE_COLS = ["cds_length", "NMDesc_region_length", "nmdesc_gc_content",
                   "pfam_overlap_length", "pfam_overlap_fraction"]

GROUP_TO_STRATUM = {"snv": "SNV", "snv_control": "SNV", "fs": "FS", "fs_control": "FS"}


def add_plot_flag_aliases(gene_all):
    out = gene_all.copy()
    for src, dst in REPO_TO_PLOT_FLAGS.items():
        if src in out.columns:
            out[dst] = out[src]
    present = [c for c in PLOT_FLAGCOLS if c in out.columns]
    out["n_flags_set"] = out[present].astype("float").ge(1).sum(axis=1)
    out["n_flags_available"] = out[present].notna().sum(axis=1)
    return out


def gene_flags_sets(matched_pairs, gene_all, flag_cols=None, extra_cols=None):
    if flag_cols is None:
        flag_cols = PLOT_FLAGCOLS
    if extra_cols is None:
        extra_cols = EXTRA_GENE_COLS
    g = gene_all.copy()
    if "stratum" not in g.columns:
        g["stratum"] = g["group"].map(GROUP_TO_STRATUM)
    gu = g.drop_duplicates(subset=["ensembl_transcript_id", "stratum"]) \
          .set_index(["ensembl_transcript_id", "stratum"])
    take = [c for c in list(flag_cols) + list(extra_cols) if c in g.columns]
    rows = []
    for stratum, m in matched_pairs.groupby("stratum"):
        for i, (ctx, sub) in enumerate(m.groupby("case_transcript", sort=True), start=1):
            pid = "%s_%03d" % (sub["case_group"].iloc[0], i)
            members = [("case", sub["case_hgnc"].iloc[0], ctx, 1)]
            members += [("control", r.control_hgnc, r.control_transcript, 0)
                        for r in sub.itertuples()]
            for role, sym, tx, isd in members:
                s = gu.loc[(tx, stratum)] if (tx, stratum) in gu.index else None
                rec = dict(design="matched", pair_id=pid, stratum=stratum, role=role,
                           is_disease=isd, hgnc_symbol=sym, ensembl_transcript_id=tx,
                           uniprot=(s.get("uniprot") if s is not None else None),
                           group=(s.get("group") if s is not None else None))
                rec.update({c: (s.get(c) if s is not None else None) for c in take})
                rows.append(rec)
    return pd.DataFrame(rows)


# ── skill: nmdesc-all-figures helpers ────────────────────────────────────────
MATCHED_PAIRS_COLUMNS = (
    "case_group", "case_hgnc", "case_transcript", "case_cds_length",
    "case_NMD_region_start", "case_NMD_region_end", "case_NMDesc_region_length",
    "case_nmd_cds_ratio",
    "control_group", "control_hgnc", "control_transcript", "control_cds_length",
    "control_NMD_region_start", "control_NMD_region_end",
    "control_NMDesc_region_length", "control_nmd_cds_ratio",
    "cds_log10_diff", "ratio_abs_diff",
)


def nmdesc_pairs_to_wide(pairs, gene_all, case_col="case_transcript",
                         control_col="control_transcript", group_col="case_group"):
    need = ("ensembl_transcript_id", "group", "hgnc_symbol", "cds_length",
            "NMD_region_start", "NMD_region_end", "NMDesc_region_length")
    absent = [c for c in need if c not in gene_all.columns]
    if absent:
        raise ValueError("gene_all is missing %s" % absent)
    ga = gene_all.drop_duplicates(["ensembl_transcript_id", "group"])

    def side(tx_col, grp_col, prefix):
        m = pairs[[tx_col, grp_col]].rename(
            columns={tx_col: "ensembl_transcript_id", grp_col: "group"})
        j = m.merge(ga[list(need)], on=["ensembl_transcript_id", "group"], how="left")
        j = j.rename(columns={"group": prefix + "_group", "hgnc_symbol": prefix + "_hgnc",
                               "ensembl_transcript_id": prefix + "_transcript",
                               "cds_length": prefix + "_cds_length",
                               "NMD_region_start": prefix + "_NMD_region_start",
                               "NMD_region_end": prefix + "_NMD_region_end",
                               "NMDesc_region_length": prefix + "_NMDesc_region_length"})
        j[prefix + "_nmd_cds_ratio"] = (j[prefix + "_NMDesc_region_length"]
                                        / j[prefix + "_cds_length"])
        return j.reset_index(drop=True)

    ctrl_group = "control_group" if "control_group" in pairs.columns else group_col
    wide = pd.concat([side(case_col, group_col, "case"),
                      side(control_col, ctrl_group, "control")], axis=1)
    wide["cds_log10_diff"] = np.log10(wide["control_cds_length"] / wide["case_cds_length"])
    wide["ratio_abs_diff"] = (wide["control_nmd_cds_ratio"]
                              - wide["case_nmd_cds_ratio"]).abs()
    return wide[list(MATCHED_PAIRS_COLUMNS)]


# ── house Pfam rule ──────────────────────────────────────────────────────────
def house_pfam_flag(gene_all):
    out = gene_all.copy()
    ok = out.max_pfam_end.notna() & out.NMD_region_start.notna()
    out["Pfam_flag"] = np.where(ok, (out.max_pfam_end * 3 >= out.NMD_region_start).astype(float),
                                np.nan)
    return out


# ── matching helper ──────────────────────────────────────────────────────────
def match_k_relaxed(cases, pool, k, cds_col="cds_length",
                    region_col="NMDesc_region_length", label="case"):
    def feats(df):
        cds = np.log10(df[cds_col].astype(float).values)
        ratio = df[region_col].astype(float).values / df[cds_col].astype(float).values
        return np.column_stack([cds, ratio])

    fc, fp = feats(cases), feats(pool)
    both = np.vstack([fc, fp])
    mu, sd = np.nanmean(both, 0), np.nanstd(both, 0)
    sd[sd == 0] = 1.0
    zc, zp = (fc - mu) / sd, (fp - mu) / sd
    cost = ((zc[:, None, :] - zp[None, :, :]) ** 2).sum(axis=2)
    ri, ci = linear_sum_assignment(np.repeat(cost, k, axis=0))
    out = []
    for r, c in zip(ri, ci):
        cs, ct = cases.iloc[r // k], pool.iloc[c]
        out.append(dict(
            case_group=label,
            case_hgnc=cs["hgnc_symbol"], case_transcript=cs["ensembl_transcript_id"],
            case_cds_length=cs[cds_col], case_NMDesc_region_length=cs[region_col],
            case_nmd_cds_ratio=float(cs[region_col]) / float(cs[cds_col]),
            control_hgnc=ct["hgnc_symbol"], control_transcript=ct["ensembl_transcript_id"],
            control_cds_length=ct[cds_col], control_NMDesc_region_length=ct[region_col],
            control_nmd_cds_ratio=float(ct[region_col]) / float(ct[cds_col]),
            cds_log10_diff=abs(np.log10(float(cs[cds_col])) - np.log10(float(ct[cds_col]))),
            ratio_abs_diff=abs(float(cs[region_col]) / float(cs[cds_col])
                               - float(ct[region_col]) / float(ct[cds_col]))))
    df = pd.DataFrame(out)
    df["requested_k"] = k
    df["feasible"] = len(pool) >= k * len(cases)
    return df


def rekey_on_group(g, gene_all, flag_cols, extra_cols):
    ga = gene_all.drop_duplicates(["ensembl_transcript_id", "group"]).set_index(
        ["ensembl_transcript_id", "group"])
    take = [c for c in list(flag_cols) + list(extra_cols) if c in gene_all.columns]
    out = g.copy()
    exp = out.stratum.str.lower().map({"fs": "fs", "snv": "snv"}) + np.where(out.role.eq("case"), "", "_control")
    out["group"] = exp
    idx = pd.MultiIndex.from_arrays([out.ensembl_transcript_id, out.group])
    miss = ~idx.isin(ga.index)
    if miss.any():
        raise ValueError("no (transcript, group) feature row for %d design rows" % int(miss.sum()))
    pulled = ga.loc[idx, take + ["uniprot"]].reset_index(drop=True)
    for c in take + ["uniprot"]:
        out[c] = pulled[c].values
    return out


# ── constants ────────────────────────────────────────────────────────────────
STRAT_LAB = {"SNV": "snv", "FS": "fs"}
FLAGCOLS = list(PLOT_FLAGCOLS) + ["IDR_flag"]
EXTRA = [c for c in EXTRA_GENE_COLS] + ["NMD_region_start", "NMD_region_end"]
PLOTF = PLOT_FLAGCOLS

# ── load data ────────────────────────────────────────────────────────────────
u = pd.read_csv('/Users/qkelly/.claude-science/orgs/1fc7f1be-05bf-47d6-99ab-3e536d468e52/artifacts/proj_7b35821b8355/11ff0bb3-5040-45a6-8bc7-3eedb4bbce4d/vec38bf7b_union_all_gene_withfeatures_fixed.csv', low_memory=False)
u_house = house_pfam_flag(add_plot_flag_aliases(u))

gf = pd.read_csv('/Users/qkelly/.claude-science/orgs/1fc7f1be-05bf-47d6-99ab-3e536d468e52/artifacts/proj_7b35821b8355/d60c36c4-55b6-4567-bdc5-4d9ff57edb2b/v8267f3f6_gene_flags_per_gene_set4_1to2.csv')

# ── gene set 4 case symbols ──────────────────────────────────────────────────
case_sym = {st: sorted(gf[(gf.stratum == st) & (gf.is_disease == 1)].hgnc_symbol.unique())
            for st in ("SNV", "FS")}

allcase = set(case_sym["SNV"]) | set(case_sym["FS"])

cases = {st: u_house[(u_house.group == g) & (u_house.hgnc_symbol.isin(case_sym[st]))].copy()
         for st, g in STRAT_LAB.items()}
pools = {st: u_house[(u_house.group == g + "_control") & (~u_house.hgnc_symbol.isin(allcase))].copy()
         for st, g in STRAT_LAB.items()}

# ── 1:1 matching for set 4 ──────────────────────────────────────────────────
ps = []
for st in ("FS", "SNV"):
    c, p = cases[st], pools[st]
    pr = match_k_relaxed(c, p, k=1, label=STRAT_LAB[st])
    pr["stratum"] = st
    ps.append(pr)
pairs_set4 = pd.concat(ps, ignore_index=True)
pairs_set4["control_group"] = pairs_set4["case_group"] + "_control"

# ── build gene_flags and rekey on group ─────────────────────────────────────
g = gene_flags_sets(pairs_set4, u_house, flag_cols=FLAGCOLS, extra_cols=EXTRA)
g["n_flags_set"] = g[PLOTF].astype("float").ge(1).sum(axis=1)
g["n_flags_available"] = g[PLOTF].notna().sum(axis=1)
g = rekey_on_group(g, u_house, FLAGCOLS, EXTRA)
g["n_flags_set"] = g[PLOTF].astype("float").ge(1).sum(axis=1)
g["n_flags_available"] = g[PLOTF].notna().sum(axis=1)

# ── write wide matched pairs ─────────────────────────────────────────────────
os.makedirs("one2one", exist_ok=True)
w = nmdesc_pairs_to_wide(pairs_set4, u_house)
w.to_csv("matched_pairs_full_set4_1to1.csv", index=False)