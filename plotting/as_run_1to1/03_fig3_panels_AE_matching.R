#!/usr/bin/env Rscript
# Combined matching figure: A-C length distributions (violins) + D-E balance Love plot.
# Inputs : matched_pairs_full.csv, gene_all_AD.csv
# Outputs: fig_matching_combined.png / .pdf
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(scales)
})

## ---- palette ----------------------------------------------------------
COL  <- c("Frameshift"="#2E5E8C","FS Control"="#A9C6DE",
          "Stopgain"="#A8423F","Stopgain Control"="#E6B6AC")
EDGE <- c("Frameshift"="#1C3F61","FS Control"="#6E97B8",
          "Stopgain"="#7A2B29","Stopgain Control"="#C58A7E")
GORDER <- c("Frameshift","FS Control","Stopgain","Stopgain Control")
XPOS   <- c("Frameshift"=1,"FS Control"=2.2,"Stopgain"=3.8,"Stopgain Control"=5.0)

base_theme <- theme_classic(base_size = 15) +
  theme(plot.title   = element_text(face="bold", size=17, hjust=0.5,
                                     margin=margin(b=10)),
        axis.title   = element_text(size=15),
        axis.text    = element_text(size=12.5, colour="#3a3a3a"),
        axis.line    = element_line(linewidth=0.8, colour="#3a3a3a"),
        axis.ticks   = element_line(linewidth=0.8, colour="#3a3a3a"),
        panel.grid.major.y = element_line(linewidth=0.5, colour="#e8e8e8"),
        plot.tag     = element_text(face="bold", size=22))

## ---- data -------------------------------------------------------------
mp <- read.csv("matched_pairs_full.csv")
g  <- read.csv("gene_all_AD.csv")
fs <- mp %>% filter(case_group=="fs")
sg <- mp %>% filter(case_group=="snv")

make_long <- function(cc, kc){
  rbind(
    data.frame(group="Frameshift",       value=fs[[cc]]),
    data.frame(group="FS Control",       value=fs[[kc]]),
    data.frame(group="Stopgain",         value=sg[[cc]]),
    data.frame(group="Stopgain Control", value=sg[[kc]])
  ) |> mutate(group = factor(group, levels=GORDER),
              xpos  = XPOS[as.character(group)])
}
pv <- function(a,b) wilcox.test(a,b, exact=FALSE)$p.value
pfmt <- function(p){
  txt <- if(p>=.10) sprintf("P = %.2f",p) else if(p>=.01) sprintf("P = %.3f",p) else "P < 0.001"
  sig <- if(p>=.05) "ns" else if(p>=.01) "*" else if(p>=.001) "**" else "***"
  paste0(txt,"  (",sig,")")
}

## ---- violin panel builder --------------------------------------------
violin_panel <- function(cc, kc, title, ylab, logy){
  d <- make_long(cc, kc)
  if(logy) d$y <- log10(d$value) else d$y <- d$value
  p_fs <- pv(d$y[d$group=="Frameshift"], d$y[d$group=="FS Control"])
  p_sg <- pv(d$y[d$group=="Stopgain"],   d$y[d$group=="Stopgain Control"])

  set.seed(7)
  pl <- ggplot(d, aes(xpos, y, group=group)) +
    geom_violin(aes(fill=group, colour=group), width=0.78, alpha=0.60,
                linewidth=1.6, position=position_identity()) +
    geom_jitter(aes(colour=group), width=0.075, height=0, size=1.3,
                alpha=0.40, show.legend=FALSE) +
    geom_boxplot(aes(colour=group), fill="white", width=0.18, alpha=0.95,
                 linewidth=1.4, outlier.shape=NA, fatten=2.4,
                 position=position_identity(), show.legend=FALSE) +
    scale_fill_manual(values=COL) + scale_colour_manual(values=EDGE) +
    scale_x_continuous(breaks=XPOS,
                       labels=c("Frameshift","FS\ncontrol","Stopgain","Stopgain\ncontrol"),
                       limits=c(0.3,5.7)) +
    labs(title=title, x=NULL, y=ylab) +
    guides(fill="none", colour="none") + base_theme

  rng  <- range(d$y); span <- diff(rng)
  br_y <- function(x1,x2,top,p){
    y <- top + 0.10*span; h <- 0.022*span
    list(annotate("segment", x=c(x1,x1,x2), xend=c(x1,x2,x2),
                  y=c(y,y+h,y+h), yend=c(y+h,y+h,y), colour="#4a4a4a", linewidth=0.8),
         annotate("text", x=(x1+x2)/2, y=y+h+0.075*span, label=pfmt(p),
                  size=4.6, colour="#2a2a2a"))
  }
  top_fs <- max(d$y[d$group %in% c("Frameshift","FS Control")])
  top_sg <- max(d$y[d$group %in% c("Stopgain","Stopgain Control")])
  pl <- pl + br_y(XPOS[["Frameshift"]],XPOS[["FS Control"]],top_fs,p_fs) +
             br_y(XPOS[["Stopgain"]],XPOS[["Stopgain Control"]],top_sg,p_sg) +
    coord_cartesian(ylim=c(rng[1]-0.05*span, rng[2]+0.40*span), clip="off")

  if(logy){
    ints <- seq(floor(rng[1]), ceiling(rng[2]))
    pl <- pl + scale_y_continuous(breaks=ints,
                 labels=function(b) formatC(10^b, format="d", big.mark=",")) +
      annotate("text", x=0.35, y=rng[2]+0.38*span, label="log scale",
               hjust=0, fontface="italic", size=3.6, colour="#8a8a8a")
  }
  pl
}

pA <- violin_panel("case_cds_length","control_cds_length",
                   "CDS length","CDS length (bp)", TRUE)
pB <- violin_panel("case_NMDesc_region_length","control_NMDesc_region_length",
                   "NMD escape region length","NMD escape region length (bp)", TRUE)
pC <- violin_panel("case_nmd_cds_ratio","control_nmd_cds_ratio",
                   "NMD / CDS length ratio","NMD escape / CDS length", FALSE)

## ---- Love plot panels -------------------------------------------------
g$ratio   <- g$NMDesc_region_length / g$cds_length
g$log_cds <- log10(g$cds_length)
g$log_nmd <- log10(g$NMDesc_region_length)
smd <- function(t,c,s) (mean(t)-mean(c))/s
VLEV <- c("NMD region length","NMD / CDS ratio","CDS length")  # bottom->top

love_df <- function(case_grp, ctrl_grp, label){
  cs <- g[g$group==case_grp,]; pool <- g[g$group==ctrl_grp,]
  used <- mp$control_transcript[mp$case_group==case_grp]
  mt <- g[g$group==ctrl_grp & g$ensembl_transcript_id %in% used,]
  specs <- list(c("log_cds","CDS length","*"),
                c("ratio","NMD / CDS ratio","*"),
                c("log_nmd","NMD region length",""))
  do.call(rbind, lapply(specs, function(v){
    s <- sd(cs[[v[1]]])
    data.frame(class=label, variable=v[2], star=v[3],
               before=smd(cs[[v[1]]], pool[[v[1]]], s),
               after =smd(cs[[v[1]]], mt[[v[1]]], s))
  }))
}

love_panel <- function(case_grp, ctrl_grp, label, ncase, npool){
  df <- love_df(case_grp, ctrl_grp, label)
  df$variable <- factor(df$variable, levels=VLEV)
  df$ylab <- paste0(as.character(df$variable), ifelse(df$star=="*"," *",""))
  lab_map <- setNames(df$ylab, as.character(df$variable))
  ggplot(df) +
    annotate("rect", xmin=-0.1, xmax=0.1, ymin=0.4, ymax=3.6,
             fill="#e7f0e9") +
    geom_vline(xintercept=0, colour="#9a9a9a", linewidth=0.7) +
    geom_vline(xintercept=c(-0.1,0.1), colour="#7bb38a", linetype="22", linewidth=0.7) +
    geom_segment(aes(x=before, xend=after, y=variable, yend=variable),
                 colour="#bdbdbd", linewidth=1.4) +
    geom_point(aes(x=before, y=variable), shape=21, fill="white",
               colour="#9a9a9a", stroke=1.5, size=4.6) +
    geom_point(aes(x=after, y=variable), shape=21, fill=COL[[label]],
               colour="white", stroke=1.0, size=5.0) +
    scale_y_discrete(labels=lab_map, expand=expansion(add=0.6)) +
    scale_x_continuous(limits=c(-0.15,0.62)) +
    labs(title=sprintf("%s  (n = %d vs %d)", label, ncase, npool),
         x="Standardized mean difference (case - control)", y=NULL) +
    base_theme +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(size=16))
}

n_fs_case <- sum(g$group=="fs");  n_fs_pool <- sum(g$group=="fs_control")
n_sg_case <- sum(g$group=="snv"); n_sg_pool <- sum(g$group=="snv_control")
pD <- love_panel("fs","fs_control","Frameshift", n_fs_case, n_fs_pool)
pE <- love_panel("snv","snv_control","Stopgain",  n_sg_case, n_sg_pool)

## ---- assemble ---------------------------------------------------------
fig <- (pA | pB | pC) / (pD | pE) +
  plot_layout(heights=c(1, 0.78)) +
  plot_annotation(
    title="Length-matched controls reproduce the size distribution of PTC variants",
    tag_levels="A",
    theme=theme(plot.title=element_text(face="bold", size=20, hjust=0.5,
                                         margin=margin(b=6)))) &
  theme(plot.tag.position=c(0,1))

ggsave("fig_matching_combined.png", fig, width=16, height=12, dpi=300, bg="white")
ggsave("fig_matching_combined.pdf", fig, width=16, height=12, bg="white")
cat("rank-sum P  CDS:", pv(log10(fs$case_cds_length),log10(fs$control_cds_length)),
    pv(log10(sg$case_cds_length),log10(sg$control_cds_length)), "\nsaved fig_matching_combined\n")
