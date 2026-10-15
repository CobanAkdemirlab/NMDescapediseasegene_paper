#!/usr/bin/env Rscript
# =============================================================================
# go_robustness_permutation.R
# Reviewer rebuttal for "GO analysis robustness": AD disease genes as a class
# are enriched for developmental terms, so ANY AD subset looks enriched against
# a whole-genome background. This script tests whether the frameshift / stopgain
# NMD-escape sets are enriched BEYOND random AD gene sets, two ways:
#
#   (A) AD-BACKGROUND enrichment: hypergeometric test with universe = the AD
#       gene list (not the genome). Terms significant here are enriched relative
#       to AD genes generally, not just relative to the genome.
#
#   (B) PERMUTATION NULL (the experiment the reviewer asked for): draw N_PERM
#       random subsets of the same size from the AD list, run the SAME
#       whole-genome enrichment used originally, and build the null distribution
#       of (i) number of FDR<5% terms and (ii) strength of the top term.
#       The observed set's percentile / empirical p tells you whether it is more
#       enriched than random AD subsets.
#
# Inputs (same dir or /mnt/user-data/uploads):
#   GMT_FILE  GO:BP gene sets, gene SYMBOLS  (e.g. MSigDB c5.go.bp.*.symbols.gmt)
#   AD_FILE   AD disease-gene universe, one HGNC symbol per line
#             (e.g. derived from OMIM genemap2.txt, AD phenotypes)
#   gene_flags_per_gene.csv   -> the frameshift / stopgain case gene sets
# Outputs: go_robustness_permutation.png/.pdf, go_robustness_stats.csv
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
set.seed(1)

N_PERM   <- 2000
FDR      <- 0.05
MIN_SET  <- 10      # ignore tiny GO terms
MAX_SET  <- 500     # ignore huge unspecific GO terms

pick <- function(f, required = TRUE) {
  cc <- c(f, file.path("/mnt/user-data/uploads", f)); hit <- cc[file.exists(cc)]
  if (!length(hit)) { if (required) stop("missing input: ", f, call. = FALSE) else return(NA) }
  hit[1] }

read_gmt <- function(path) {
  ln <- readLines(path); sets <- list()
  for (l in ln) { p <- strsplit(l, "\t")[[1]]; if (length(p) < 3) next
    sets[[p[1]]] <- unique(toupper(p[-c(1,2)])) }
  sets }

# ---- build gene x term incidence once (fast permutations) ------------------
build_incidence <- function(gmt) {
  genes_all <- sort(unique(unlist(gmt)))
  gi <- setNames(seq_along(genes_all), genes_all)
  # sparse-ish: store, per term, the integer indices of its genes
  term_idx <- lapply(gmt, function(gs) gi[intersect(gs, genes_all)])
  list(genes_all = genes_all, gi = gi, term_idx = term_idx)
}
# membership counts k for a set of gene-indices, per term (fast)
count_k <- function(set_idx, term_idx) {
  hit <- logical(length(INC$genes_all)); hit[set_idx] <- TRUE
  vapply(term_idx, function(ti) sum(hit[ti]), integer(1))
}
enrich_stats <- function(set_genes, K, Nuniv, keep) {
  si <- INC$gi[intersect(set_genes, names(INC$gi))]
  n  <- length(si)
  k  <- count_k(si, INC$term_idx)
  p  <- phyper(k - 1, K, Nuniv - K, n, lower.tail = FALSE)
  padj <- rep(NA_real_, length(p)); padj[keep] <- p.adjust(p[keep], "BH")
  c(n_sig = sum(padj < FDR, na.rm = TRUE),
    top   = suppressWarnings(max(-log10(padj), na.rm = TRUE)))
}

# ---- load inputs -----------------------------------------------------------
gmt_path <- NA
for (cand in c("GO_BP.gmt", "c5.all.v2026.1.Hs.symbols.gmt", "c5.go.bp.symbols.gmt",
               "/mnt/user-data/uploads/GO_BP.gmt",
               "/mnt/user-data/uploads/c5.go.bp.v2023.2.Hs.symbols.gmt"))
  if (file.exists(cand)) { gmt_path <- cand; break }
if (is.na(gmt_path))
  stop("Need a GO:BP gene-set file (.gmt, gene symbols), e.g. MSigDB c5.go.bp.*.symbols.gmt", call. = FALSE)
gmt <- read_gmt(gmt_path)
# ---- restrict to GO Biological Process terms ONLY --------------------------
# (c5.all mixes GOBP_/GOCC_/GOMF_/HP_ ; we keep only GOBP_)
is_bp <- grepl("^GOBP_", names(gmt))
if (any(is_bp)) {
  message(sprintf("Restricted to GO:BP: kept %d GOBP_ terms, dropped %d non-BP (GOCC_/GOMF_/HP_/...)",
                  sum(is_bp), sum(!is_bp)))
  gmt <- gmt[is_bp]
} else {
  message("No 'GOBP_' terms found -- assuming the supplied .gmt is already GO:BP-only (",
          length(gmt), " terms)")
}
ad  <- unique(toupper(trimws(readLines(pick("AD_gene_universe.txt")))))
ad  <- ad[nzchar(ad) & ad != "HGNC_SYMBOL"]
g   <- read.csv(pick("gene_flags_per_gene.csv"), stringsAsFactors = FALSE)
g   <- g[g$design == "matched", ]
CASE <- list(
  Frameshift = unique(toupper(g$hgnc_symbol[g$role == "case" & g$stratum == "FS"])),
  Stopgain   = unique(toupper(g$hgnc_symbol[g$role == "case" & g$stratum == "SNV"])))
ad <- union(ad, unlist(CASE))

INC <- build_incidence(gmt)
Kall <- vapply(INC$term_idx, length, integer(1))          # term sizes in genome
keep <- Kall >= MIN_SET & Kall <= MAX_SET
N_genome <- length(INC$genes_all)
ad_in <- intersect(ad, INC$genes_all)                     # AD genes that are annotated
K_ad  <- count_k(INC$gi[ad_in], INC$term_idx); N_ad <- length(ad_in)
cat(sprintf("AD universe: %d (annotated %d) | GO:BP sets used: %d | genome universe: %d\n",
            length(ad), N_ad, sum(keep), N_genome))

# ---- analysis per group ----------------------------------------------------
plots <- list(); rows <- list()
for (grp in names(CASE)) {
  genes <- intersect(CASE[[grp]], ad_in); n <- length(genes)

  # (A) AD-background enrichment: universe = AD genes
  si <- INC$gi[genes]; k <- count_k(si, INC$term_idx)
  p_ad <- phyper(k - 1, K_ad, N_ad - K_ad, n, lower.tail = FALSE)
  padj_ad <- rep(NA_real_, length(p_ad)); padj_ad[keep] <- p.adjust(p_ad[keep], "BH")
  n_sig_adbg <- sum(padj_ad < FDR, na.rm = TRUE)

  # observed whole-genome enrichment (matches the original analysis)
  obs <- enrich_stats(genes, Kall, N_genome, keep)

  # (B) permutation: random AD subsets, whole-genome enrichment
  null <- t(replicate(N_PERM, {
    samp <- INC$genes_all[sample(INC$gi[ad_in], n)]
    enrich_stats(samp, Kall, N_genome, keep) }))
  p_nsig <- (1 + sum(null[, "n_sig"] >= obs["n_sig"])) / (N_PERM + 1)
  p_top  <- (1 + sum(null[, "top"]   >= obs["top"],   na.rm = TRUE)) / (N_PERM + 1)

  rows[[grp]] <- data.frame(group = grp, n_genes = n,
      obs_n_sig_genome = obs["n_sig"], perm_p_n_sig = p_nsig,
      obs_top_neglog10 = round(obs["top"], 2), perm_p_top = p_top,
      n_sig_AD_background = n_sig_adbg)

  nd <- data.frame(n_sig = null[, "n_sig"])
  plots[[grp]] <- ggplot(nd, aes(n_sig)) +
    geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    geom_vline(xintercept = obs["n_sig"], colour = "#B2182B", linewidth = 1.2) +
    annotate("text", x = Inf, y = Inf, vjust = 1.4, hjust = 1.03,
             label = sprintf("observed = %d\nperm p = %.3f\nAD-background sig = %d",
                             as.integer(obs["n_sig"]), p_nsig, n_sig_adbg),
             colour = "#B2182B", fontface = "bold", size = 4) +
    labs(title = sprintf("%s (n=%d)", grp, n),
         subtitle = "Null = random AD subsets, whole-genome enrichment",
         x = "# FDR<5% GO:BP terms", y = "random AD subsets") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 10))
}
stats <- bind_rows(rows); write.csv(stats, "go_robustness_stats.csv", row.names = FALSE)
print(stats, row.names = FALSE)

fig <- wrap_plots(plots, ncol = length(plots)) +
  plot_annotation(
    title = "GO enrichment robustness: observed vs random AD-gene subsets",
    subtitle = sprintf("%d permutations \u00b7 red line = observed # significant GO:BP terms \u00b7 empirical p and AD-background counts annotated", N_PERM),
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = .5),
                  plot.subtitle = element_text(size = 11, hjust = .5, colour = "grey30")))
ggsave("go_robustness_permutation.png", fig, width = 7.5 * length(plots), height = 6, dpi = 200, bg = "white")
ggsave("go_robustness_permutation.pdf", fig, width = 7.5 * length(plots), height = 6, bg = "white")
cat("\nsaved go_robustness_permutation.png / .pdf and go_robustness_stats.csv\n")
