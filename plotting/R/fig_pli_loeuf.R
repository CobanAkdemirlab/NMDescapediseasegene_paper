#!/usr/bin/env Rscript
# =============================================================================
# FIGURE — pLI / LOEUF category distributions, paired (case vs matched control).
#   A pLI  Frameshift   B pLI  Stopgain
#   C LOEUF Frameshift   D LOEUF Stopgain
# Per-category paired McNemar's exact test. gnomAD v4.1 (MANE/canonical pref).
#
# Two-stage: if CACHE (pli_loeuf_genes.csv, written by build stage) exists, plot
# straight from it (no gnomAD download needed); else build from gnomAD + pairs.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(forcats); library(patchwork); library(scales)
})

CACHE       <- "pli_loeuf_genes.csv"
PAIRS_FILE  <- "matching_log_v2.csv"
GNOMAD_FILE <- Sys.getenv("GNOMAD_FILE", "gnomad.v4.1.constraint_metrics.tsv")

grp_cols <- c("Frameshift"="#2166AC","FS Control"="#92C5DE",
              "Stopgain"="#B2182B","Stopgain Control"="#F4A582")
edge_cols<- c("Frameshift"="#15406F","FS Control"="#5E97BE",
              "Stopgain"="#7E1019","Stopgain Control"="#D17A55")

# ── category functions (published thresholds) ────────────────────────────────
cat_pli <- function(x) factor(dplyr::case_when(
  is.na(x) ~ NA_character_,
  x <= 0.35            ~ "Low\n(pLI \u2264 0.35)",
  x > 0.35 & x < 0.66  ~ "Medium\n(0.35 < pLI < 0.66)",
  x >= 0.66            ~ "High\n(pLI \u2265 0.66)"),
  levels=c("Low\n(pLI \u2264 0.35)","Medium\n(0.35 < pLI < 0.66)","High\n(pLI \u2265 0.66)"))
cat_loeuf <- function(x) factor(dplyr::case_when(
  is.na(x) ~ NA_character_,
  x <= 0.2            ~ "Low\n(LOEUF \u2264 0.2)",
  x > 0.2 & x <= 0.6  ~ "Medium\n(0.2 < LOEUF \u2264 0.6)",
  x > 0.6             ~ "High\n(LOEUF > 0.6)"),
  levels=c("Low\n(LOEUF \u2264 0.2)","Medium\n(0.2 < LOEUF \u2264 0.6)","High\n(LOEUF > 0.6)"))

annotate_num <- function(df) df %>% mutate(
  study_pliCat=cat_pli(study_pLI),   ctrl_pliCat=cat_pli(ctrl_pLI),
  study_loeufCat=cat_loeuf(study_LOEUF), ctrl_loeufCat=cat_loeuf(ctrl_LOEUF))

# ── data: cache or build ─────────────────────────────────────────────────────
if (file.exists(CACHE)) {
  message("Cache found -> plotting from ", CACHE)
  allann <- read_csv(CACHE, show_col_types=FALSE)
} else {
  message("No cache -> building from gnomAD + pairs ...")
  if (!file.exists(GNOMAD_FILE))
    stop("gnomAD constraint file not found. Download v4.1 constraint metrics and set GNOMAD_FILE.")
  pairs <- read_csv(PAIRS_FILE, show_col_types=FALSE)
  gr <- read_tsv(GNOMAD_FILE, show_col_types=FALSE) %>%
    dplyr::select(gene, transcript, canonical, mane_select,
                  pLI=`lof.pLI`, LOEUF=`lof.oe_ci.upper`) %>% filter(!is.na(gene)) %>%
    mutate(prio=as.integer(mane_select)*2L+as.integer(canonical)) %>%
    arrange(gene, desc(prio), desc(pLI)) %>% distinct(gene, .keep_all=TRUE)
  pli<-setNames(gr$pLI,gr$gene); loe<-setNames(gr$LOEUF,gr$gene)
  allann <- pairs %>% transmute(study_gene, ctrl_gene, group,
    study_pLI=pli[study_gene], ctrl_pLI=pli[ctrl_gene],
    study_LOEUF=loe[study_gene], ctrl_LOEUF=loe[ctrl_gene])
  write_csv(allann, CACHE)
}
allann <- annotate_num(allann)
ann_fs  <- allann %>% filter(group=="FS")
ann_snv <- allann %>% filter(group=="Stopgain")

# ── McNemar per category ─────────────────────────────────────────────────────
mcnemar_cat <- function(sc, cc, target){
  ok<-!is.na(sc)&!is.na(cc); s<-sc[ok]==target; c<-cc[ok]==target
  b<-sum(s&!c); d<-sum(!s&c)
  if((b+d)==0) return(NA_real_)
  min(2*pbinom(min(b,d), b+d, 0.5), 1)        # exact (mid-binomial) McNemar
}
sig <- function(p) if(is.na(p)) "" else
  if(p>=.05)"ns" else if(p>=.01)"*" else if(p>=.001)"**" else "***"

# ── per-panel data ───────────────────────────────────────────────────────────
panel_df <- function(ann, study_lbl, ctrl_lbl, scol, ccol, levs){
  ok<-!is.na(ann[[scol]])&!is.na(ann[[ccol]]); a<-ann[ok,]; n<-nrow(a)
  ts<-prop.table(table(factor(a[[scol]],levels=levs)))*100
  tc<-prop.table(table(factor(a[[ccol]],levels=levs)))*100
  d<-tibble(Group=factor(rep(c(study_lbl,ctrl_lbl),each=length(levs)),
                         levels=names(grp_cols)),
            Category=factor(rep(levs,2),levels=levs),
            pct=c(as.numeric(ts),as.numeric(tc)))
  ps<-sapply(levs,function(L) mcnemar_cat(a[[scol]],a[[ccol]],L))
  list(d=d, p=ps, n=n, study=study_lbl, ctrl=ctrl_lbl)
}

theme_big <- theme_classic(base_size=17) + theme(
  plot.title  = element_text(face="bold", size=19, hjust=0.5, margin=margin(b=10)),
  axis.title  = element_text(size=16, face="bold"),
  axis.text.x = element_text(size=13, colour="#1a1a1a", lineheight=0.9),
  axis.text.y = element_text(size=14, colour="#1a1a1a"),
  axis.line   = element_line(linewidth=0.9, colour="#2a2a2a"),
  axis.ticks  = element_line(linewidth=0.9, colour="#2a2a2a"),
  panel.grid.major.y = element_line(linewidth=0.45, colour="#ececec"),
  plot.tag    = element_text(face="bold", size=26),
  legend.position="bottom", legend.title=element_blank(),
  legend.text=element_text(size=13))

make_panel <- function(P, xlab, title, tag, nlab){
  d<-P$d; dodge<-0.66
  d$xnum <- as.integer(d$Category) + ifelse(d$Group==P$study, -0.18, 0.18)
  cats<-levels(d$Category)
  top<-sapply(seq_along(cats),function(i) max(d$pct[as.integer(d$Category)==i]))
  br<-tibble(i=seq_along(cats),
             xs=seq_along(cats)-0.18, xc=seq_along(cats)+0.18,
             yb=top+9, p=P$p, lab=paste0("p = ",ifelse(is.na(P$p),"NA",sprintf("%.3f",P$p)),
                                          "  ", sapply(P$p,sig)))
  ggplot(d, aes(Category, pct, fill=Group)) +
    geom_col(aes(colour=Group), position=position_dodge(width=0.72),
             width=dodge, linewidth=0.8, alpha=0.92) +
    geom_text(aes(label=sprintf("%d%%",round(pct)), group=Group),
              position=position_dodge(width=0.72), vjust=-0.4,
              size=5.0, fontface="bold", colour="#1a1a1a") +
    # significance brackets
    geom_segment(data=br, inherit.aes=FALSE, aes(x=xs,xend=xc,y=yb,yend=yb),
                 colour="#444444", linewidth=0.7) +
    geom_segment(data=br, inherit.aes=FALSE, aes(x=xs,xend=xs,y=yb,yend=yb-2.5),
                 colour="#444444", linewidth=0.7) +
    geom_segment(data=br, inherit.aes=FALSE, aes(x=xc,xend=xc,y=yb,yend=yb-2.5),
                 colour="#444444", linewidth=0.7) +
    geom_text(data=br, inherit.aes=FALSE, aes(x=i, y=yb+3.5, label=lab,
              fontface=ifelse(!is.na(p)&p<0.05,"bold","plain")),
              size=4.5, colour="#1a1a1a") +
    scale_fill_manual(values=grp_cols, labels=nlab, drop=FALSE, limits=names(grp_cols)) +
    scale_colour_manual(values=edge_cols, drop=FALSE, limits=names(grp_cols), guide="none") +
    scale_y_continuous(limits=c(0,120), breaks=seq(0,100,20),
                       labels=function(x) paste0(x,"%"), expand=expansion(mult=c(0,0))) +
    labs(title=title, x=xlab, y="Percentage (%)", tag=tag) +
    guides(fill=guide_legend(override.aes=list(colour=NA))) +
    theme_big
}

pli_levs   <- levels(cat_pli(0));  loeuf_levs <- levels(cat_loeuf(0))
Pa<-panel_df(ann_fs ,"Frameshift","FS Control","study_pliCat","ctrl_pliCat",pli_levs)
Pb<-panel_df(ann_snv,"Stopgain","Stopgain Control","study_pliCat","ctrl_pliCat",pli_levs)
Pc<-panel_df(ann_fs ,"Frameshift","FS Control","study_loeufCat","ctrl_loeufCat",loeuf_levs)
Pd<-panel_df(ann_snv,"Stopgain","Stopgain Control","study_loeufCat","ctrl_loeufCat",loeuf_levs)

nlab <- setNames(c(sprintf("Frameshift (n=%d)",Pa$n), sprintf("FS Control (n=%d)",Pa$n),
                   sprintf("Stopgain (n=%d)",Pb$n),   sprintf("Stopgain Control (n=%d)",Pb$n)),
                 names(grp_cols))

pA<-make_panel(Pa,"pLI category","pLI Category Distribution \u2014 Frameshift","A",nlab)
pB<-make_panel(Pb,"pLI category","pLI Category Distribution \u2014 Stopgain","B",nlab)
pC<-make_panel(Pc,"LOEUF category","LOEUF Category Distribution \u2014 Frameshift","C",nlab)
pD<-make_panel(Pd,"LOEUF category","LOEUF Category Distribution \u2014 Stopgain","D",nlab)

fig <- (pA | pB) / (pC | pD) +
  plot_layout(guides="collect") +
  plot_annotation(title="pLI and LOEUF category distributions",
    theme=theme(plot.title=element_text(face="bold", size=23, hjust=0.5, margin=margin(b=6)))) &
  theme(legend.position="bottom", plot.tag.position=c(0,1))

suppressWarnings({
  ggsave("Figure_pli_loeuf.png", fig, width=17, height=12.5, dpi=300, bg="white")
  ggsave("Figure_pli_loeuf.pdf", fig, width=17, height=12.5, bg="white")
})
cat("saved Figure_pli_loeuf\n")
