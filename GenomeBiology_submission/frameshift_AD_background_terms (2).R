#!/usr/bin/env Rscript
# =============================================================================
# frameshift_AD_background_terms.R
# Lists the GO:BP terms enriched in the FRAMESHIFT NMD-escape set when the
# background is the AD gene list (not the genome) -- i.e. the terms that survive
# the reviewer-proof test. Writes a clean table for the rebuttal.
#
# Inputs (same dir or /mnt/user-data/uploads):
#   GO:BP gene-set .gmt (symbols)  -> GO_BP.gmt or c5.go.bp.*.symbols.gmt
#   AD_gene_universe.txt           -> AD disease-gene background, one symbol/line
#   gene_flags_per_gene.csv        -> frameshift case set
# Output: frameshift_AD_background_GO_terms.csv (+ printed table)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(stringr) })
FDR_CUT <- 0.05; MIN_SET <- 10; MAX_SET <- 500

pick <- function(f) { cc <- c(f, file.path("/mnt/user-data/uploads", f))
  hit <- cc[file.exists(cc)]; if (!length(hit)) stop("missing input: ", f, call. = FALSE); hit[1] }
read_gmt <- function(path) { ln <- readLines(path); s <- list()
  for (l in ln) { p <- strsplit(l, "\t")[[1]]; if (length(p) < 3) next
    s[[p[1]]] <- unique(toupper(p[-c(1,2)])) }; s }
clean_term <- function(x) tolower(gsub("_", " ", sub("^GO[A-Z]*_", "", x)))

gmt_path <- NA
for (cand in c("GO_BP.gmt", "c5.go.bp.v2023.2.Hs.symbols.gmt", "c5.go.bp.symbols.gmt", "c5.all.v2026.1.Hs.symbols.gmt", "c5.all.v2026.1.Hs.symbols.gmt",
               "/mnt/user-data/uploads/GO_BP.gmt",
               "/mnt/user-data/uploads/c5.go.bp.v2023.2.Hs.symbols.gmt",
               "/mnt/user-data/uploads/c5.all.v2026.1.Hs.symbols.gmt"))
  if (file.exists(cand)) { gmt_path <- cand; break }
if (is.na(gmt_path)) stop("Need a GO:BP .gmt (gene symbols), e.g. MSigDB c5.go.bp.*.symbols.gmt", call. = FALSE)
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
fs  <- unique(toupper(g$hgnc_symbol[g$role == "case" & g$stratum == "FS"]))
ad  <- union(ad, fs)                                   # case genes must be in the universe

# ---- AD-background hypergeometric enrichment -------------------------------
universe <- intersect(ad, unique(unlist(gmt)))         # AD genes that carry GO annotation
genes    <- intersect(fs, universe)
N <- length(universe); n <- length(genes)

rows <- lapply(names(gmt), function(t) {
  tg <- intersect(gmt[[t]], universe); K <- length(tg)
  if (K < MIN_SET || K > MAX_SET) return(NULL)
  hit <- intersect(genes, tg); k <- length(hit)
  data.frame(ID = t, Description = clean_term(t),
             case_hits = k, term_size = K,
             fold_enrichment = (k / n) / (K / N),
             p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
             genes = paste(sort(hit), collapse = ", "), stringsAsFactors = FALSE)
})
res <- bind_rows(rows)
res$FDR <- p.adjust(res$p, "BH")

sig <- res %>% filter(FDR < FDR_CUT) %>% arrange(FDR) %>%
  transmute(ID, Description,
            overlap = sprintf("%d/%d", case_hits, term_size),
            fold = round(fold_enrichment, 2),
            p = signif(p, 3), FDR = signif(FDR, 3), genes)
write.csv(sig, "frameshift_AD_background_GO_terms.csv", row.names = FALSE)

cat(sprintf("\nFrameshift set: %d annotated case genes | AD-annotated background: %d genes | GO:BP terms tested: %d\n",
            n, N, sum(!is.na(res$FDR))))
cat(sprintf("GO:BP terms significant at FDR < %.0f%% against the AD background: %d\n\n",
            100 * FDR_CUT, nrow(sig)))
print(as.data.frame(sig[, c("Description", "overlap", "fold", "FDR")]), row.names = FALSE, right = FALSE)
cat("\nFull table (with gene lists) written to frameshift_AD_background_GO_terms.csv\n")

# ---- figure: lollipop of AD-background-significant terms --------------------
TOP_SHOW <- 30   # cap rows for legibility; set to nrow(sig) to show all
if (nrow(sig) >= 1) {
  d <- sig %>% arrange(FDR) %>% head(TOP_SHOW) %>%
    mutate(k = as.integer(sub("/.*", "", overlap)),
           neglogFDR = -log10(as.numeric(FDR)),
           lab = str_to_sentence(str_wrap(Description, 42)))
  d$lab <- factor(d$lab, levels = rev(d$lab))   # most significant at top
  ttl_n <- nrow(sig)
  p <- ggplot(d, aes(neglogFDR, lab)) +
    geom_segment(aes(x = 0, xend = neglogFDR, yend = lab),
                 colour = "grey85", linewidth = 0.7) +
    geom_point(aes(size = k, colour = fold)) +
    scale_colour_gradient(low = "#F6B8B2", high = "#8B0A1A",
                          name = "Fold\nenrichment\n(vs AD)") +
    scale_size_continuous(range = c(3, 9), name = "Frameshift\ngenes") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Frameshift NMD-escape: GO:BP terms enriched above AD background",
         subtitle = sprintf("%d Biological-Process terms at FDR < 5%% with the AD gene list as background  ·  %d frameshift genes",
                            ttl_n, n),
         x = expression(-log[10](FDR)), y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13.5),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey35", size = 10),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10, colour = "black"))
  h <- max(4, 0.42 * nrow(d) + 1.6)
  ggsave("frameshift_AD_background_GO.png", p, width = 11.5, height = h, dpi = 200, bg = "white")
  ggsave("frameshift_AD_background_GO.pdf", p, width = 11.5, height = h, bg = "white")
  cat(sprintf("Figure written: frameshift_AD_background_GO.png / .pdf (%d terms shown)\n",
              nrow(d)))
}
