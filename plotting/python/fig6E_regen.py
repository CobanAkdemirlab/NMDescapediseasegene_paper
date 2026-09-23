#!/usr/bin/env python3
"""
fig11_regen.py -- regenerate the escape-tail two-panel figure:
  * Panel A: ParSe PS, ONLY the Delta-vs-WT rows (direct rows removed)
  * Panel B: composition/charge/disorder; disorder = REAL metapredict (from per_variant_disorder.csv)
  * No text label on n.s. (grey) rows.  Gene-matched LMM value ~ cohort + (1|gene).
Inputs: variants_all0901.csv, NMD_region_DivergentPos_CIDER.csv, nmd_length_DivergentPos_PARSE_v2.csv,
        gene_all_withflags_0826.csv, per_variant_disorder.csv (real metapredict)
"""
import os, math, numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf, statsmodels.stats.multitest as mt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
warnings.simplefilter("ignore")
def pick(f, bases=(".","/mnt/user-data/uploads")):
    stem,ext=os.path.splitext(f); cands=[f,stem+"__1_"+ext,stem.replace("__1_","")+ext]
    for b in bases:
        for c in cands:
            if os.path.exists(os.path.join(b,c)): return os.path.join(b,c)
    raise FileNotFoundError(f)
AA=set("ACDEFGHIKLMNPQRSTVWY")
PKA={'Nterm':9.69,'Cterm':2.34,'D':3.65,'E':4.25,'C':8.18,'Y':10.07,'H':6.00,'K':10.54,'R':12.48}; DPROM=set("SPEKQGAD")
def cpH(s,pH):
    c=1.0/(1.0+10**(pH-PKA['Nterm']))-1.0/(1.0+10**(PKA['Cterm']-pH))
    for a in s:
        if a in('D','E','C','Y'): c+=-1.0/(1.0+10**(PKA[a]-pH))
        elif a in('H','K','R'):   c+= 1.0/(1.0+10**(pH-PKA[a]))
    return c
def pI(s):
    lo,hi=0.0,14.0
    while hi-lo>0.01:
        mid=(lo+hi)/2; lo,hi=(mid,hi) if cpH(s,mid)>0 else (lo,mid)
    return (lo+hi)/2
def scd(s):
    if len(s)<2: return 0.0
    q=[1.0 if a in('K','R') else -1.0 if a in('D','E') else 0.0 for a in s]; N=len(s)
    ci=[(i,q[i]) for i in range(N) if q[i]]; x=0.0
    for a in range(len(ci)):
        for b in range(a+1,len(ci)): x+=ci[a][1]*ci[b][1]*math.sqrt(abs(ci[b][0]-ci[a][0]))
    return x/(N*N)
def comp(seq):
    seq="".join(x for x in str(seq).upper() if x.isalpha()); n=len(seq)
    if n==0: return {}
    fp=sum(a in('R','K') for a in seq)/n; fn=sum(a in('D','E') for a in seq)/n
    return dict(frac_disorder_promoting=sum(a in DPROM for a in seq)/n,prion_like=sum(a in set('QNGSY') for a in seq)/n,
                SCD=scd(seq),FCR=fp+fn,NCPR=fp-fn,pI=pI(seq))
# common variant set
m=pd.read_csv(pick("variants_all0901.csv")); m["key_us"]=m["key"].str.replace(":","_").str.replace("|","_"); MASTER=set(m.key_us)
w=pd.read_csv(pick("gene_all_withflags_0826.csv")); w=w[w.group=="fs"]; w["up"]=w.uniprot.map(lambda x:str(x).split("-")[0]); fsup=set(w.up)
cid=pd.read_csv(pick("NMD_region_DivergentPos_CIDER.csv")); cid=cid[cid.source_folder.str.startswith("fs_")].copy(); cid["key"]=cid.file_name
par=pd.read_csv(pick("nmd_length_DivergentPos_PARSE_v2.csv")); par=par[par['id'].str.split("_").str[0]=="fs"].copy(); par["key"]=par['id'].str.split("-",n=1).str[1]
cdc=[c for c in par.columns if "classifier" in c]; psc=[c for c in par.columns if "PS IDR" in c or "longest" in c]
for c in cdc+psc+["var_NMD_length","wt_NMD_length"]: par[c]=pd.to_numeric(par[c],errors="coerce")
d=cid.merge(par[["key"]+cdc+psc+["var_NMD_length","wt_NMD_length"]],on="key",how="inner")
d=d[d.key.isin(MASTER)].copy(); d["cohort"]=d.source_folder.str.split("_").str[1]; d["up"]=d.uniprot_id.map(lambda x:str(x).split("-")[0])
d=d[d.up.isin(fsup)].copy(); d["gene"]=d.up
d["vlen"]=d.var_sequence_after_NMD.map(lambda s:sum(x in AA for x in str(s).upper())); d=d[(d.vlen>=20)&(d.var_NMD_length>0)&(d.wt_NMD_length>0)]
d["cohort_n"]=(d.cohort=="disease").astype(int); d=d.reset_index(drop=True)
d["cd_delta"]=d[cdc[0]]/d.var_NMD_length-d[cdc[1]]/d.wt_NMD_length
d["ps_delta"]=d[psc[0]]/d.var_NMD_length-d[psc[1]]/d.wt_NMD_length
vm=pd.DataFrame([comp(s) for s in d.var_sequence_after_NMD]).add_prefix("var_")
wm=pd.DataFrame([comp(s) for s in d.wt_sequence_after_NMD ]).add_prefix("wt_")
d=pd.concat([d,vm,wm],axis=1)
for k in ["frac_disorder_promoting","prion_like","SCD","FCR","NCPR","pI"]: d["d_"+k]=d["var_"+k]-d["wt_"+k]
# real metapredict disorder (same row order & n)
pv=pd.read_csv(pick("per_variant_disorder.csv"))
if len(pv)==len(d):
    d["d_mean_disorder"]=(pv["var_mean_disorder"]-pv["wt_mean_disorder"]).values
    d["d_frac_disordered"]=(pv["var_frac_disordered"]-pv["wt_frac_disordered"]).values
    DIS_SRC="metapredict"
else:
    raise SystemExit(f"per_variant_disorder.csv has {len(pv)} rows but analysis has {len(d)} -- re-run fig9 to regenerate it")
def lmm(col):
    x=d.dropna(subset=[col]).copy(); sd=x[col].std() or 1; x["y"]=x[col]/sd
    tb=x.groupby("gene")["cohort_n"].nunique(); x=x[x.gene.isin(tb[tb==2].index)]
    if x.gene.nunique()<3: return np.nan,np.nan,np.nan
    try:
        mm=smf.mixedlm("y ~ cohort_n",x,groups=x["gene"]).fit(reml=False,method="lbfgs")
        return mm.params["cohort_n"],mm.bse["cohort_n"],mm.pvalues["cohort_n"]
    except Exception: return np.nan,np.nan,np.nan
def run(items):
    R=pd.DataFrame([dict(lab=l,**dict(zip(["b","se","p"],lmm(c)))) for l,c in items]); R["q"]=mt.multipletests(R.p.fillna(1),method="fdr_bh")[1]; return R
A=run([("Classifier / residue","cd_delta"),("PS-IDR fraction","ps_delta")])   # Delta vs WT only, direct removed
Bitems=[("disorder (mean)","d_mean_disorder"),("disorder (frac)","d_frac_disordered"),
        ("disorder-promoting frac","d_frac_disorder_promoting"),("prion-like","d_prion_like"),
        ("SCD (charge patterning)","d_SCD"),("FCR","d_FCR"),("NCPR (net charge)","d_NCPR"),("pI","d_pI")]
B=run(Bitems)
def sig(q): return "***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else "trend" if q<0.10 else "n.s."
def colr(b,q):
    if pd.isna(q) or q>=0.10: return "#c0c0c0"
    hi=b>0; base=("#b2182b","#ef8a62","#fddbc7") if hi else ("#2166ac","#67a9cf","#d1e5f0")
    return base[0] if q<0.01 else base[1] if q<0.05 else base[2]
def panel(ax,R,title,green=None):
    labs=list(R.lab)[::-1]
    for i,lab in enumerate(labs):
        r=R[R.lab==lab].iloc[0]; b=r.b; se=r.se if pd.notna(r.se) else 0; q=r.q
        if green and lab in green: ax.axhspan(i-0.45,i+0.45,color="#eef6ee",zorder=0)
        ax.errorbar(b,i,xerr=1.96*se,fmt="o",color=colr(b,q),ms=9,lw=2.3,mec="black",mew=0.5)
        if pd.notna(q) and q<0.10:                       # NO label on grey (n.s.) rows
            ax.text(b+1.96*se+0.03,i,("\u25b2 " if b>0 else "\u25bc ")+f"q={q:.3f} {sig(q)}",va="center",fontsize=9,color="#222")
    ax.axvline(0,color="grey",lw=0.8); ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs,fontsize=10)
    ax.set_ylim(-0.6,len(labs)-0.4); ax.margins(x=0.36); ax.set_xlabel("pathogenic\u2212benign \u00b7 SD",fontsize=10)
    ax.set_title(title,fontsize=12,fontweight="bold",loc="left")
ng=d.gene.nunique(); dn=int((d.cohort=="disease").sum()); cn=int((d.cohort=="control").sum())
fig=plt.figure(figsize=(15,6)); gs=fig.add_gridspec(1,2,width_ratios=[1,1.05],wspace=0.5)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1])
panel(axA,A,"A. Escape-tail ParSe phase separation")
panel(axB,B,"B. Escape-tail composition / charge / disorder",green={"disorder (mean)","disorder (frac)"})
fig.legend(handles=[Line2D([0],[0],marker='o',color='#b2182b',ls='',label='higher in pathogenic (q<0.05)'),
                    Line2D([0],[0],marker='o',color='#2166ac',ls='',label='lower in pathogenic (q<0.05)'),
                    Line2D([0],[0],marker='o',color='#c0c0c0',ls='',label='n.s.')],
           loc="lower center",ncol=3,fontsize=9,frameon=False,bbox_to_anchor=(0.5,-0.02))
fig.suptitle(f"Escape-tail, frameshift, gene-matched LMM (value ~ cohort + (1|gene)) \u2014 {ng} genes, {dn} P/LP / {cn} benign",fontsize=13,fontweight="bold",y=1.0)
plt.figtext(0.5,-0.06,f"Divergent-position tail only.  Variants in variants_all0901 \u2229 CIDER-div \u2229 ParSe-div, fs, \u226520aa (n={len(d)}).  ParSe: \u0394-vs-WT only.  Disorder = {DIS_SRC}.",ha="center",fontsize=8,color="#555")
plt.tight_layout(rect=[0,0.02,1,0.99]); plt.savefig("escape_tail_LMM_metapredict.png",dpi=180,bbox_inches="tight",facecolor="white")
plt.savefig("escape_tail_LMM_metapredict.pdf",bbox_inches="tight",facecolor="white")
A.to_csv("escape_tail_ParSe_deltaWT_LMM.csv",index=False); B.to_csv("escape_tail_composition_metapredict_LMM.csv",index=False)
print("A:\n",A.round(3).to_string(index=False)); print("B:\n",B.round(3).to_string(index=False)); print("saved fig11")
