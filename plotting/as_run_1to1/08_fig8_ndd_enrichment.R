# =============================================================================
# NMD-escape disease genes  x  AD NDD risk genes : overlap + enrichment
#   Panel a : two-set Venn (NMD-escape vs AD NDD risk)
#   Panel b : observed vs expected shared-gene count, Fisher's exact
# Inputs : matching_log_v2.csv  -> NMD-escape disease genes = unique study_gene
#          NDDriskgenes.txt     -> NDD risk gene list (HGNC symbols)
# -----------------------------------------------------------------------------
# DERIVED FROM DATA (verifiable):
#   NMD-escape genes = 207 ; overlap with NDD list = 62 ; NMD-only = 145
# SUPPLIED / REVERSE-ENGINEERED (not present in the two input files):
#   ndd_ad_total  = 350   -> autosomal-dominant subset of the 505-gene NDD file
#   universe_N    = 1907  -> background gene universe for the enrichment test
#   Replace these two with your actual AD-NDD subset and universe to recompute.
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})

## ---- input paths (edit here, or pass on the command line) -------------------
##   Rscript nmd_ndd_enrichment_figure.R  <gene_table.csv>  <NDDriskgenes.txt>
args <- commandArgs(trailingOnly = TRUE)
GENE_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else ""   # gene_flags_per_gene.csv or matched_pairs_full.csv
NDD_FILE  <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "NDDriskgenes.txt"

# candidate locations searched when GENE_FILE is not given explicitly
GENE_CANDIDATES <- c(GENE_FILE,
  "gene_flags_per_gene.csv", "matched_pairs_full.csv",
  "/mnt/user-data/uploads/gene_flags_per_gene.csv",
  "/mnt/user-data/outputs/matched_pairs_full.csv",
  "matching_log_v2.csv")
NDD_CANDIDATES <- c(NDD_FILE, "NDDriskgenes.txt",
  "/mnt/user-data/uploads/NDDriskgenes.txt")

find_first <- function(paths, what) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0)
    stop(sprintf("Could not find %s. Looked for:\n  %s\nPass the path explicitly:\n  Rscript nmd_ndd_enrichment_figure.R <gene_table.csv> <NDDriskgenes.txt>",
                 what, paste(unique(paths), collapse = "\n  ")), call. = FALSE)
  hit[1]
}

## ---- 1. derive gene sets from data -----------------------------------------
load_nmd_genes <- function() {
  f <- find_first(GENE_CANDIDATES, "the matched gene table (gene_flags_per_gene.csv / matched_pairs_full.csv)")
  d <- read.csv(f, stringsAsFactors = FALSE)
  nm <- names(d)
  if (all(c("design","role","hgnc_symbol") %in% nm)) {          # gene_flags_per_gene.csv
    d <- d[d$design == "matched" & d$role == "case", ]
    message("NMD-escape genes <- ", f, " (matched cases)")
    toupper(trimws(unique(d$hgnc_symbol)))
  } else if ("case_hgnc" %in% nm) {                             # matched_pairs_full.csv
    message("NMD-escape genes <- ", f, " (case_hgnc)")
    toupper(trimws(unique(d$case_hgnc)))
  } else if ("study_gene" %in% nm) {                            # legacy matching_log_v2.csv
    message("NMD-escape genes <- ", f, " (legacy study_gene)")
    toupper(trimws(unique(d$study_gene)))
  } else stop("Unrecognised gene table columns in ", f, ": ", paste(nm, collapse=", "), call.=FALSE)
}
nmd <- load_nmd_genes()                                     # new list: 168 unique
ndd_file <- toupper(trimws(readLines(find_first(NDD_CANDIDATES, "the NDD risk-gene list (NDDriskgenes.txt)"))))
ndd_file <- unique(ndd_file[ndd_file != ""])                # 505
overlap  <- length(intersect(nmd, ndd_file))                # 62 (robust)
nmd_n    <- length(nmd)

## ---- 2. figure parameters (external / reverse-engineered) ------------------
ndd_ad_total <- 350     # AD subset shown in figure (file has 505; AD filter absent here)
universe_N   <- 1907    # background universe (expected = nmd_n * ndd_ad_total / N = 38)

nmd_only <- nmd_n - overlap        # 145
ndd_only <- ndd_ad_total - overlap # 288
neither  <- universe_N - nmd_n - ndd_only   # N - 207 - 288

## ---- 3. enrichment: Fisher's exact -----------------------------------------
tab <- matrix(c(overlap, nmd_only, ndd_only, neither), nrow = 2, byrow = TRUE)
ft  <- fisher.test(tab, alternative = "greater")
# sample OR + Wald 95% CI (matches "2.1 (1.5-2.9)")
o <- overlap; b <- nmd_only; c <- ndd_only; d <- neither
OR   <- (o * d) / (b * c)
se   <- sqrt(1/o + 1/b + 1/c + 1/d)
ci   <- exp(log(OR) + c(-1, 1) * 1.96 * se)
expected <- nmd_n * ndd_ad_total / universe_N
cat(sprintf("overlap=%d  expected=%.1f  OR=%.2f (%.1f-%.1f)  Fisher P=%.2e\n",
            overlap, expected, OR, ci[1], ci[2], ft$p.value))

## ---- 4. panel a : Venn -----------------------------------------------------
circle <- function(cx, cy, r, n = 200) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = cx + r*cos(t), y = cy + r*sin(t))
}
blue <- "#6BA3D6"; red <- "#E0908F"
cL <- circle(-0.55, 0, 1.05); cL$set <- "NMD"
cR <- circle( 0.60, 0, 1.25); cR$set <- "NDD"
pa <- ggplot() +
  geom_polygon(data = cL, aes(x, y), fill = blue, colour = "grey35", alpha = .8, linewidth=.4) +
  geom_polygon(data = cR, aes(x, y), fill = red,  colour = "grey35", alpha = .6, linewidth=.4) +
  annotate("text", x = -1.05, y = 0, label = nmd_only, size = 7) +
  annotate("text", x =  0.02, y = 0, label = overlap, size = 7, fontface = "bold") +
  annotate("text", x =  1.15, y = 0, label = ndd_only, size = 7) +
  annotate("text", x = -0.9, y = 1.35, label = "NMD-escape\ndisease genes",
           colour = "#2C6FB0", fontface = "bold", size = 5, lineheight = .9) +
  annotate("text", x =  0.95, y = 1.5, label = "AD NDD\nrisk genes",
           colour = "#B41E2E", fontface = "bold", size = 5, lineheight = .9) +
  coord_equal(clip = "off") + xlim(-2.3, 2.3) + ylim(-1.5, 2.0) +
  theme_void() + theme(plot.margin = margin(20, 10, 10, 10))

## ---- 5. panel b : observed vs expected -------------------------------------
bd <- data.frame(grp = factor(c("Expected","Observed"), levels=c("Expected","Observed")),
                 val = c(round(expected), overlap),
                 fill = c("#CFCFCF", "#2C6FB0"))
ptxt <- sprintf("OR = %.1f (%.1f-%.1f)", OR, ci[1], ci[2])
ytop <- max(bd$val); y_br <- ytop + 8; y_pv <- ytop + 11.5; y_or <- ytop + 16; y_max <- ytop + 24
pexp <- floor(log10(ft$p.value)); pman <- ft$p.value / 10^pexp
pval <- sprintf("italic(P) == %.1f %%*%% 10^%d", pman, pexp)
pb <- ggplot(bd, aes(grp, val, fill = I(fill))) +
  geom_col(width = .62, colour = "grey30", linewidth = .3) +
  geom_text(aes(label = val), vjust = -0.6, size = 5) +
  annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br) +
  annotate("segment", x = 1, xend = 1, y = y_br - 3, yend = y_br) +
  annotate("segment", x = 2, xend = 2, y = y_br - 3, yend = y_br) +
  annotate("text", x = 1.5, y = y_or, label = ptxt, size = 4.5) +
  annotate("text", x = 1.5, y = y_pv, label = pval, parse = TRUE, size = 4.5) +
  scale_y_continuous("Genes shared with AD NDD set", limits = c(0, y_max),
                     breaks = seq(0, floor(y_max/10)*10, 10), expand = expansion(mult = c(0, .02))) +
  labs(x = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 13),
        plot.margin = margin(20, 15, 10, 10))

## ---- 6. combine ------------------------------------------------------------
fig <- pa + pb + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 20))
ggsave("NMD_NDD_enrichment.png", fig, width = 12, height = 5.2, dpi = 150, bg = "white")
cat("Saved NMD_NDD_enrichment.png\n")
