#!/usr/bin/env Rscript
# =============================================================================
# stopgain_AD_background_terms.R
# GO:BP enrichment of the STOP-GAIN (SNV) NMD-escape set against the AD gene
# background. Even when NO term passes FDR<5%, this ALWAYS draws a plot of the
# top terms by raw p, with a dashed line marking the FDR<5% threshold, so the
# stop-gain landscape can be shown alongside the frameshift figure (panel B).
#
# Inputs (same dir or /mnt/user-data/uploads):
#   GO:BP gene-set .gmt (symbols)  -> GO_BP.gmt or c5.go.bp.*.symbols.gmt
#   AD_gene_universe.txt           -> AD disease-gene background, one symbol/line
#   gene_flags_per_gene.csv        -> stop-gain case set (stratum == "SNV")
# Output: stopgain_AD_background_GO_terms.csv, stopgain_AD_background_GO.png/.pdf
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(stringr) })
FDR_CUT <- 0.05; MIN_SET <- 10; MAX_SET <- 500
TOP_SHOW <- 20            # how many top terms to display when none are significant

pick <- function(f) { cc <- c(f, file.path("/mnt/user-data/uploads", f))
  hit <- cc[file.exists(cc)]; if (!length(hit)) stop("missing input: ", f, call. = FALSE); hit[1] }
read_gmt <- function(path) { ln <- readLines(path); s <- list()
  for (l in ln) { p <- strsplit(l, "\t")[[1]]; if (length(p) < 3) next
    s[[p[1]]] <- unique(toupper(p[-c(1,2)])) }; s }
clean_term <- function(x) tolower(gsub("_", " ", sub("^GO[A-Z]*_", "", x)))

gmt_path <- NA
for (cand in c("GO_BP.gmt", "c5.all.v2026.1.Hs.symbols.gmt", "c5.go.bp.symbols.gmt", "c5.all.v2026.1.Hs.symbols.gmt",
               "/mnt/user-data/uploads/GO_BP.gmt",
               "/mnt/user-data/uploads/c5.go.bp.v2023.2.Hs.symbols.gmt",
               "/mnt/user-data/uploads/c5.all.v2026.1.Hs.symbols.gmt"))
  if (file.exists(cand)) { gmt_path <- cand; break }
if (is.na(gmt_path)) stop("Need a GO:BP .gmt (gene symbols), e.g. MSigDB c5.go.bp.*.symbols.gmt", call. = FALSE)
gmt <- read_gmt(gmt_path)
is_bp <- grepl("^GOBP_", names(gmt))
if (any(is_bp)) { message(sprintf("Restricted to GO:BP: kept %d, dropped %d non-BP", sum(is_bp), sum(!is_bp))); gmt <- gmt[is_bp] } else
  message("No 'GOBP_' terms found -- assuming .gmt is already GO:BP-only (", length(gmt), " terms)")

ad  <- unique(toupper(trimws(readLines(pick("AD_gene_universe.txt")))))
ad  <- ad[nzchar(ad) & ad != "HGNC_SYMBOL"]
g   <- read.csv(pick("gene_flags_per_gene.csv"), stringsAsFactors = FALSE)
g   <- g[g$design == "matched", ]
sg  <- unique(toupper(g$hgnc_symbol[g$role == "case" & g$stratum == "SNV"]))
ad  <- union(ad, sg)

universe <- intersect(ad, unique(unlist(gmt)))
genes    <- intersect(sg, universe)
N <- length(universe); n <- length(genes)

rows <- lapply(names(gmt), function(t) {
  tg <- intersect(gmt[[t]], universe); K <- length(tg)
  if (K < MIN_SET || K > MAX_SET) return(NULL)
  hit <- intersect(genes, tg); k <- length(hit)
  data.frame(ID = t, Description = clean_term(t), case_hits = k, term_size = K,
             fold_enrichment = (k / n) / (K / N),
             p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
             genes = paste(sort(hit), collapse = ", "), stringsAsFactors = FALSE)
})
res <- bind_rows(rows); res$FDR <- p.adjust(res$p, "BH")
n_sig <- sum(res$FDR < FDR_CUT, na.rm = TRUE)

# full table (all tested terms, sorted) + a significant-only convenience column
out <- res %>% arrange(p) %>%
  transmute(ID, Description, overlap = sprintf("%d/%d", case_hits, term_size),
            fold = round(fold_enrichment, 2), p = signif(p, 3), FDR = signif(FDR, 3),
            passes_FDR05 = FDR < FDR_CUT, genes)
write.csv(out, "stopgain_AD_background_GO_terms.csv", row.names = FALSE)
cat(sprintf("\nStop-gain set: %d annotated case genes | AD-annotated background: %d | GO:BP tested: %d\n", n, N, sum(!is.na(res$FDR))))
cat(sprintf("Terms at FDR < %.0f%% vs AD background: %d\n\n", 100*FDR_CUT, n_sig))

# ---- figure: ALWAYS drawn. Show significant terms if any; else top TOP_SHOW by raw p ----
if (n_sig >= 1) { d0 <- res %>% filter(FDR < FDR_CUT) %>% arrange(FDR); shown <- nrow(d0); mode_lab <- sprintf("%d term(s) at FDR < 5%%", n_sig) } else {
  d0 <- res %>% arrange(p) %>% head(TOP_SHOW); shown <- nrow(d0); mode_lab <- sprintf("no term reaches FDR < 5%%; top %d terms by raw p shown", shown) }

# dashed line = the -log10(FDR) position of the FDR<5% cutoff, mapped back to a raw-p scale
# (the smallest raw p whose BH-FDR would be 0.05 given the number of tests)
m_tests <- sum(!is.na(res$FDR))
p_at_fdr05 <- FDR_CUT * rank(res$p)[which.min(abs(res$FDR - FDR_CUT))] / m_tests   # approx BH boundary
thr_x <- -log10(min(p_at_fdr05, 0.05))

d <- d0 %>% mutate(k = case_hits, neglogP = -log10(p),
                   lab = str_to_sentence(str_wrap(Description, 42)))
d$lab <- factor(d$lab, levels = rev(d$lab))
p <- ggplot(d, aes(neglogP, lab)) +
  geom_segment(aes(x = 0, xend = neglogP, yend = lab), colour = "grey85", linewidth = 0.7) +
  geom_point(aes(size = k, colour = fold_enrichment)) +
  geom_vline(xintercept = thr_x, linetype = "dashed", colour = "grey40") +
  annotate("text", x = thr_x, y = 0.6, label = "FDR = 0.05", hjust = -0.05, vjust = 0,
           size = 3.2, colour = "grey40") +
  scale_colour_gradient(low = "#F6B8B2", high = "#8B0A1A", name = "Fold\nenrichment\n(vs AD)") +
  scale_size_continuous(range = c(3, 9), name = "Stop-gain\ngenes") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Stop-gain NMD-escape: GO:BP terms vs AD background",
       subtitle = sprintf("%s  \u00b7  %d stop-gain genes  \u00b7  AD gene list as background", mode_lab, n),
       x = expression(-log[10](italic(p))), y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13.5), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 10),
        panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10, colour = "black"))
h <- max(4, 0.42 * nrow(d) + 1.7)
ggsave("stopgain_AD_background_GO.png", p, width = 11.5, height = h, dpi = 200, bg = "white")
ggsave("stopgain_AD_background_GO.pdf", p, width = 11.5, height = h, bg = "white")
cat(sprintf("Figure written: stopgain_AD_background_GO.png / .pdf (%d terms shown; %s)\n", shown, mode_lab))
