#!/usr/bin/env Rscript
# =============================================================================
# FIGURE 2 — Paired analysis: PPI degree centrality (STRING) & tissue
# specificity (Tau), length-matched pairs.
#   A. STRING v11.5 degree centrality (score >= 400), log10
#   B. Tissue specificity Tau (gene-level, GTEx-derived Tau_gene_V8.csv)
# Test: paired Wilcoxon signed-rank, BH-adjusted within each metric.
#
# Two-stage design:
#   * If CACHE_FILE (ppi_tau_4groups.csv) exists -> skip the heavy STRING/Bioc
#     build and plot directly. Otherwise build it (needs internet + Bioconductor),
#     write the cache, then plot.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(forcats); library(patchwork); library(scales)
})

CACHE_FILE <- "ppi_tau_4groups.csv"        # per-gene degree + tau (build-stage output)
# Pairs come from the matcher. Prefer the unified output of match_and_plot.R
# (matched_pairs_full.csv: case_hgnc / control_hgnc / case_group in {fs,snv});
# fall back to the legacy matching_log_v2.csv (study_gene/ctrl_gene/group {FS,Stopgain}).
PAIRS_FILE_NEW <- "matched_pairs_full.csv"
PAIRS_FILE_OLD <- "matching_log_v2.csv"

grp_levels <- c("Frameshift","FS Control","Stopgain","Stopgain Control")
grp_cols   <- c("Frameshift"="#2166AC","FS Control"="#92C5DE",
                "Stopgain"="#B2182B","Stopgain Control"="#F4A582")
edge_cols  <- c("Frameshift"="#15406F","FS Control"="#5E97BE",
                "Stopgain"="#7E1019","Stopgain Control"="#D17A55")

# ── pairs (used for paired tests + connecting lines) ─────────────────────────
load_pairs <- function() {
  if (file.exists(PAIRS_FILE_NEW)) {
    r <- read_csv(PAIRS_FILE_NEW, show_col_types = FALSE)
    stopifnot(all(c("case_hgnc","control_hgnc","case_group") %in% names(r)))
    message("Pairs <- ", PAIRS_FILE_NEW, " (new matched list, n=", nrow(r), ")")
    data.frame(study_gene = r$case_hgnc, ctrl_gene = r$control_hgnc,
               group = dplyr::recode(as.character(r$case_group),
                        fs="FS", snv="Stopgain", Frameshift="FS", Stopgain="Stopgain"),
               stringsAsFactors = FALSE)
  } else if (file.exists(PAIRS_FILE_OLD)) {
    message("Pairs <- ", PAIRS_FILE_OLD, " (legacy)")
    as.data.frame(read_csv(PAIRS_FILE_OLD, show_col_types = FALSE))
  } else stop("No pairs file: need matched_pairs_full.csv or matching_log_v2.csv")
}
pairs_log <- load_pairs()
pairs_fs  <- pairs_log |> filter(group == "FS")
pairs_snv <- pairs_log |> filter(group == "Stopgain")

# =============================================================================
# STAGE 1 — build per-gene degree + Tau  (skipped if cache present)
# =============================================================================
if (!file.exists(CACHE_FILE)) {
  message("No cache found -> building degree + Tau (needs internet + Bioconductor)...")
  ## --- packages for the build stage only ---
  for (p in c("org.Hs.eg.db","AnnotationDbi","STRINGdb","igraph")) {
    if (!requireNamespace(p, quietly=TRUE)) {
      if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
      BiocManager::install(p, ask=FALSE)
    }
  }
  suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi)
                                   library(STRINGdb); library(igraph) })
  all_genes <- unique(c(pairs_fs$study_gene, pairs_snv$study_gene,
                        pairs_fs$ctrl_gene,  pairs_snv$ctrl_gene))

  ## --- STRING degree (v11.5, score>=400) ---
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400,
                            network_type="full", input_directory=".")
  gmap <- string_db$map(data.frame(gene=all_genes), "gene", removeUnmappedRows=FALSE)
  deg_all <- igraph::degree(string_db$get_graph())
  ppi_df <- gmap |> filter(!is.na(STRING_id)) |>
    mutate(degree = deg_all[STRING_id]) |> filter(!is.na(degree)) |>
    group_by(gene) |> slice_max(degree, n=1, with_ties=FALSE) |> ungroup() |>
    dplyr::select(gene, degree)
  ppi_lookup <- setNames(ppi_df$degree, ppi_df$gene)

  ## --- Tau (gene-level, Ensembl -> HGNC) — robust mapping ---
  tau_raw <- read_csv("Tau_gene_V8.csv", show_col_types=FALSE)
  gid_col <- intersect(c("gene_id","gene","ensembl_gene_id","ensembl","Gene","ENSEMBL"),
                       names(tau_raw))[1]
  tau_col <- intersect(c("tau","Tau","TAU"), names(tau_raw))[1]
  if (is.na(gid_col) || is.na(tau_col))
    stop("Tau_gene_V8.csv: need an Ensembl-ID column and a 'tau' column; found: ",
         paste(names(tau_raw), collapse=", "))
  ens <- sub("\\..*$", "", as.character(tau_raw[[gid_col]]))   # strip version suffix ENSG...N
  # mapIds returns a NAMED VECTOR with the SYMBOL explicitly -> no missing-column pitfalls
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys=ens,
                          column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  tau_df <- data.frame(gene = unname(sym), tau = tau_raw[[tau_col]],
                       stringsAsFactors = FALSE) |>
    filter(!is.na(gene), gene != "", !is.na(tau)) |>
    group_by(gene) |> summarise(tau = dplyr::first(tau), .groups = "drop")
  tau_lookup <- setNames(tau_df$tau, tau_df$gene)

  df <- bind_rows(
    data.frame(gene=unique(pairs_fs$study_gene),  Group="Frameshift"),
    data.frame(gene=unique(pairs_fs$ctrl_gene),   Group="FS Control"),
    data.frame(gene=unique(pairs_snv$study_gene), Group="Stopgain"),
    data.frame(gene=unique(pairs_snv$ctrl_gene),  Group="Stopgain Control")) |>
    mutate(Group=factor(Group, levels=grp_levels),
           degree=ppi_lookup[gene], tau=tau_lookup[gene])
  write_csv(df, CACHE_FILE)
} else {
  message("Cache found -> plotting from ", CACHE_FILE)
  df <- read_csv(CACHE_FILE, show_col_types=FALSE) |>
    mutate(Group = factor(Group, levels=grp_levels))
}
ppi_lookup <- setNames(df$degree, df$gene)
tau_lookup <- setNames(df$tau,    df$gene)

# =============================================================================
# STAGE 2 — paired stats + the FIGURE  (A & B only, large fonts, no grey labels)
# =============================================================================
XPOS <- c("Frameshift"=1,"FS Control"=2.2,"Stopgain"=3.8,"Stopgain Control"=5.0)

paired_p <- function(pairs_df, lkp) {
  sv <- unname(lkp[pairs_df$study_gene]); cv <- unname(lkp[pairs_df$ctrl_gene])
  ok <- is.finite(sv) & is.finite(cv)
  if (sum(ok) < 3) return(c(p=NA, n=sum(ok)))
  p <- suppressWarnings(wilcox.test(sv[ok], cv[ok], paired=TRUE, exact=FALSE)$p.value)
  c(p=p, n=sum(ok))
}
sig <- function(p) if (is.na(p)) "" else
  if (p>=.05) "ns" else if (p>=.01) "*" else if (p>=.001) "**" else "***"

# BH within each metric (2 comparisons), as in the original
ppi_raw <- rbind(paired_p(pairs_fs, ppi_lookup), paired_p(pairs_snv, ppi_lookup))
tau_raw <- rbind(paired_p(pairs_fs, tau_lookup), paired_p(pairs_snv, tau_lookup))
ppi_adj <- p.adjust(ppi_raw[,"p"], "BH"); tau_adj <- p.adjust(tau_raw[,"p"], "BH")

big <- theme_classic(base_size=18) + theme(
  plot.title = element_text(face="bold", size=22, hjust=0.5, margin=margin(b=14)),
  axis.title.y = element_text(size=18, face="bold"),
  axis.text  = element_text(size=16, colour="#1a1a1a"),
  axis.line  = element_line(linewidth=0.9, colour="#2a2a2a"),
  axis.ticks = element_line(linewidth=0.9, colour="#2a2a2a"),
  plot.tag   = element_text(face="bold", size=28),
  plot.margin= margin(8, 10, 8, 40),
  panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(linewidth=0.4, colour="#ececec"))

panel <- function(metric, title, ylab, logy, padj, nvec, tag) {
  d <- df |> filter(is.finite(.data[[metric]]))
  if (logy) d <- d |> filter(.data[[metric]] > 0)
  d$x <- XPOS[as.character(d$Group)]
  d$yv <- if (logy) log10(d[[metric]]) else d[[metric]]

  # connecting lines (study<->ctrl per pair)
  mkln <- function(pdf, sg, cg) {
    L <- data.frame(pid=seq_len(nrow(pdf)),
                    sv=unname((if(metric=="degree") ppi_lookup else tau_lookup)[pdf$study_gene]),
                    cv=unname((if(metric=="degree") ppi_lookup else tau_lookup)[pdf$ctrl_gene]),
                    sg=sg, cg=cg)
    L <- L[is.finite(L$sv) & is.finite(L$cv), ]
    if (logy) { L <- L[L$sv>0 & L$cv>0, ]; L$sv<-log10(L$sv); L$cv<-log10(L$cv) }
    rbind(data.frame(pid=L$pid, grp=sg, x=unname(XPOS[sg]), y=L$sv, key=paste(L$pid,sg)),
          data.frame(pid=L$pid, grp=cg, x=unname(XPOS[cg]), y=L$cv, key=paste(L$pid,sg)))
  }
  ln <- rbind(mkln(pairs_fs,"Frameshift","FS Control"),
              mkln(pairs_snv,"Stopgain","Stopgain Control"))

  rng <- range(d$yv); span <- diff(rng)
  brk <- function(x1,x2,top,p,n){y<-top+0.10*span;h<-0.03*span
    list(annotate("segment",x=c(x1,x1,x2),xend=c(x1,x2,x2),y=c(y,y+h,y+h),
                  yend=c(y+h,y+h,y),colour="#333333",linewidth=0.9),
         annotate("text",x=(x1+x2)/2,y=y+h+0.085*span,
                  label=paste0(sig(p),"\n(n=",n," pairs)"),
                  size=5.4,lineheight=0.9,fontface=if(sig(p)=="ns")"plain" else "bold",
                  colour="#1a1a1a"))}

  p <- ggplot(d, aes(x, yv)) +
    geom_line(data=ln, aes(x=x, y=y, group=key),
              colour="grey70", linewidth=0.35, alpha=0.45) +
    geom_boxplot(aes(group=Group, fill=Group, colour=Group), width=0.5,
                 alpha=0.80, linewidth=1.1, fatten=2.6, outlier.shape=21,
                 outlier.size=2, outlier.fill=NA, outlier.stroke=0.7,
                 position=position_identity()) +
    scale_fill_manual(values=grp_cols) + scale_colour_manual(values=edge_cols) +
    scale_x_continuous(breaks=XPOS,
                       labels=c("Frameshift","FS\nControl","Stopgain","Stopgain\nControl"),
                       limits=c(0.4, 5.6)) +
    labs(title=title, x=NULL, y=ylab, tag=tag) +
    guides(fill="none", colour="none") + big +
    brk(XPOS[["Frameshift"]],XPOS[["FS Control"]],
        max(d$yv[d$Group %in% c("Frameshift","FS Control")]), padj[1], nvec[1]) +
    brk(XPOS[["Stopgain"]],XPOS[["Stopgain Control"]],
        max(d$yv[d$Group %in% c("Stopgain","Stopgain Control")]), padj[2], nvec[2]) +
    coord_cartesian(ylim=c(rng[1]-0.05*span, rng[2]+0.42*span), clip="off") +
    theme(plot.tag.position=c(0,1))
  if (logy) { ints <- seq(floor(rng[1]), ceiling(rng[2]))
    p <- p + scale_y_continuous(breaks=ints,
                labels=function(b) formatC(10^b, format="d", big.mark=",")) }
  p
}

pA <- panel("degree","PPI Degree Centrality (STRING)",
            "PPI degree (partners)", TRUE,
            ppi_adj, ppi_raw[,"n"], "F")
pB <- panel("tau","Tissue Specificity (Tau)",
            "Tau (tissue specificity)", FALSE,
            tau_adj, tau_raw[,"n"], "G")

RESULT <- (pA | pB) +
  plot_annotation(title="Protein-network and tissue-expression properties",
    theme=theme(plot.title=element_text(face="bold", size=22, hjust=0.5, margin=margin(t=2,b=16)))) &
  theme(plot.tag.position=c(0,1))

# ---- save (was missing in the original) + report ---------------------------
ggsave("fig_ppi_tau.png", RESULT, width=14, height=7, dpi=300, bg="white")
ggsave("fig_ppi_tau.pdf", RESULT, width=14, height=7, bg="white")
cat(sprintf("PPI  paired BH-P: FS=%.3g (n=%d)   Stopgain=%.3g (n=%d)\n",
            ppi_adj[1], ppi_raw[1,"n"], ppi_adj[2], ppi_raw[2,"n"]))
cat(sprintf("Tau  paired BH-P: FS=%.3g (n=%d)   Stopgain=%.3g (n=%d)\n",
            tau_adj[1], tau_raw[1,"n"], tau_adj[2], tau_raw[2,"n"]))
cat("saved fig_ppi_tau.png / .pdf\n")
