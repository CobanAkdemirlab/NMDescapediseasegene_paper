# =============================================================================
# complex_comparison.R
# Protein-complex membership across the four categories
#   Frameshift / FS Control / Stopgain / Stopgain Control
# Mirrors paralog_comparison.R / isoform_comparison.R: paired within
# length-matched pairs, same styling (black labels, larger fonts, RdBu palette,
# shaded backgrounds, paired lines, paired Wilcoxon brackets, sqrt y-axis).
#
# Metric: number of distinct CORUM complexes each protein is a subunit of.
#   0 = not a subunit of any known multiprotein complex ("or not")
#   higher = participates in more complexes
#
# CORUM is curated and not exhaustive, so absence = "no KNOWN complex" (a real 0,
# not missing data). Matching uses gene symbol AND UniProt to maximise recall.
#
# >>> ONE-TIME SETUP: download the CORUM core-complexes table <<<
#   https://mips.helmholtz-muenchen.de/corum/  ->  Download  ->  coreComplexes.txt
#   (tab-separated; place it next to this script, or set CORUM_PATH below)
#
# Outputs: complex_comparison.pdf / .png, complex_comparison_stats.csv,
#          complex_counts.csv (per-gene table, for transparency)
# =============================================================================

CORUM_PATH        <- "corum_humanComplexes.txt"  # path to the CORUM download
MULTIPROTEIN_ONLY <- TRUE                 # keep complexes with >= 2 distinct subunits
HUMAN_ONLY        <- TRUE                 # restrict to human complexes

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- macOS-safe figure saver (native quartz; no XQuartz/cairo needed) --------
save_fig <- function(path, plot, width, height, dpi = 200) {
  is_pdf <- grepl("\\.pdf$", path, ignore.case = TRUE)
  if (isTRUE(capabilities("aqua"))) {
    if (is_pdf) grDevices::quartz(file = path, type = "pdf", width = width, height = height)
    else        grDevices::quartz(file = path, type = "png", width = width, height = height, dpi = dpi)
    print(plot); grDevices::dev.off()
  } else if (is_pdf) {
    ggsave(path, plot, width = width, height = height, device = cairo_pdf)
  } else {
    ggsave(path, plot, width = width, height = height, dpi = dpi, type = "cairo")
  }
  invisible(path)
}

# ---- shared styling ----------------------------------------------------------
grp_cols <- c("Frameshift"       = "#2166AC",
              "FS Control"       = "#92C5DE",
              "Stopgain"         = "#B2182B",
              "Stopgain Control" = "#F4A582")
bg_fs <- "#2166AC"; bg_sv <- "#B2182B"; bg_alpha <- 0.06

theme_big <- theme_bw(base_size = 18) +
  theme(plot.title    = element_text(face = "bold", size = 20, hjust = 0.5, color = "black"),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
        plot.tag      = element_text(face = "bold", size = 23),
        plot.tag.position = c(0.04, 0.98),
        axis.title = element_text(size = 17, face = "bold", color = "black"),
        axis.text  = element_text(size = 15, color = "black"),
        axis.line  = element_line(linewidth = 1.0),
        axis.ticks = element_line(linewidth = 1.0),
        panel.border = element_rect(linewidth = 1.0, color = "black"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

# ---- build CORUM membership lookups (by gene symbol and by UniProt) ----------
build_corum <- function(path, human_only = TRUE, multiprotein_only = TRUE) {
  if (!file.exists(path))
    stop("CORUM file not found at '", path,
         "'. Download coreComplexes.txt from https://mips.helmholtz-muenchen.de/corum/")
  cx <- read.delim(path, sep = "\t", quote = "", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
  nm <- names(cx)
  pick <- function(pat) { i <- grep(pat, nm, ignore.case = TRUE); if (length(i)) i[1] else NA }
  i_id  <- pick("^ComplexID$"); if (is.na(i_id)) i_id <- pick("complex.?id")
  i_org <- pick("^Organism$")
  i_sym <- pick("gene.?name")
  i_uni <- pick("uniprot")
  if (is.na(i_id) || is.na(i_sym))
    stop("Could not locate ComplexID / gene-name columns in CORUM file; columns are: ",
         paste(nm, collapse = ", "))
  if (human_only && !is.na(i_org))
    cx <- cx[grepl("human", cx[[i_org]], ignore.case = TRUE), , drop = FALSE]

  splitcell <- function(s) trimws(unlist(strsplit(ifelse(is.na(s), "", s), "[;,]")))
  sym_rows <- list(); uni_rows <- list()
  for (r in seq_len(nrow(cx))) {
    id   <- cx[[i_id]][r]
    syms <- unique(splitcell(cx[[i_sym]][r])); syms <- syms[nzchar(syms)]
    if (multiprotein_only && length(syms) < 2) next
    if (length(syms)) sym_rows[[length(sym_rows) + 1]] <-
        data.frame(key = syms, cx = id, stringsAsFactors = FALSE)
    if (!is.na(i_uni)) {
      unis <- unique(sub("-.*", "", splitcell(cx[[i_uni]][r]))); unis <- unis[nzchar(unis)]
      if (length(unis)) uni_rows[[length(uni_rows) + 1]] <-
          data.frame(key = unis, cx = id, stringsAsFactors = FALSE)
    }
  }
  sym_long <- if (length(sym_rows)) unique(do.call(rbind, sym_rows)) else data.frame(key=character(), cx=character())
  uni_long <- if (length(uni_rows)) unique(do.call(rbind, uni_rows)) else data.frame(key=character(), cx=character())
  list(sym = split(sym_long$cx, sym_long$key),
       uni = split(uni_long$cx, uni_long$key))
}

complex_count <- function(symv, univ, maps) {
  vapply(seq_along(symv), function(i) {
    s <- symv[i]; u <- sub("-.*", "", univ[i])
    ids <- c(if (!is.na(s) && nzchar(s)) maps$sym[[s]],
             if (!is.na(u) && nzchar(u)) maps$uni[[u]])
    length(unique(ids))
  }, integer(1))
}

# ---- load matched pairs + count complexes ------------------------------------
pairs <- read.csv("matched_pairs_full.csv", stringsAsFactors = FALSE)
stopifnot(all(c("case_group", "case_hgnc", "control_hgnc") %in% names(pairs)))
if (!"case_uniprot"    %in% names(pairs)) pairs$case_uniprot    <- NA_character_
if (!"control_uniprot" %in% names(pairs)) pairs$control_uniprot <- NA_character_

maps <- build_corum(CORUM_PATH, HUMAN_ONLY, MULTIPROTEIN_ONLY)

pairs$case_cx <- complex_count(pairs$case_hgnc,    pairs$case_uniprot,    maps)
pairs$ctrl_cx <- complex_count(pairs$control_hgnc, pairs$control_uniprot, maps)
pairs$pair_id <- seq_len(nrow(pairs))

# per-gene transparency table
write.csv(unique(rbind(
  data.frame(hgnc = pairs$case_hgnc,    n_complex = pairs$case_cx),
  data.frame(hgnc = pairs$control_hgnc, n_complex = pairs$ctrl_cx))),
  "complex_counts.csv", row.names = FALSE)

P  <- pairs
fs <- P[P$case_group == "fs",  , drop = FALSE]
sv <- P[P$case_group == "snv", , drop = FALSE]

# ---- binary: member of any human complex? + Wilcoxon rank-sum (unpaired) -----
P$case_in <- as.integer(P$case_cx > 0)
P$ctrl_in <- as.integer(P$ctrl_cx > 0)
fs <- P[P$case_group == "fs",  , drop = FALSE]
sv <- P[P$case_group == "snv", , drop = FALSE]

# Wilcoxon rank-sum (Mann-Whitney) on binary 0/1 membership; ties -> normal approx
wilcox_grp <- function(case_in, ctrl_in)
  suppressWarnings(wilcox.test(case_in, ctrl_in, paired = FALSE, exact = FALSE)$p.value)
p_fs <- wilcox_grp(fs$case_in, fs$ctrl_in)
p_sv <- wilcox_grp(sv$case_in, sv$ctrl_in)
padj <- p.adjust(c(p_fs, p_sv), "BH")
praw <- c(p_fs, p_sv)
plab <- paste0("P = ", formatC(praw, format = "g", digits = 2))  # raw p shown on figure
psig <- praw < 0.05

stats <- data.frame(
  comparison        = c("Frameshift vs FS Control", "Stopgain vs Stopgain Control"),
  n_case            = c(nrow(fs), nrow(sv)),
  n_ctrl            = c(nrow(fs), nrow(sv)),
  frac_case_in_cplx = c(mean(fs$case_in), mean(sv$case_in)),
  frac_ctrl_in_cplx = c(mean(fs$ctrl_in), mean(sv$ctrl_in)),
  p_value_wilcoxon  = c(p_fs, p_sv),
  p_BH              = padj)
write.csv(stats, "complex_comparison_stats.csv", row.names = FALSE)
print(stats)

# ---- Panel: fraction in a complex, per category (Wilson 95% CI) --------------
wilson <- function(k, n) {
  if (n == 0) return(c(NA, NA, NA))
  p <- k / n; z <- 1.96; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(p, max(0, ctr - hw), min(1, ctr + hw))
}
summ <- data.frame(
  grp = factor(c("Frameshift", "FS Control", "Stopgain", "Stopgain Control"),
               levels = c("Frameshift", "FS Control", "Stopgain", "Stopgain Control")),
  k = c(sum(fs$case_in), sum(fs$ctrl_in), sum(sv$case_in), sum(sv$ctrl_in)),
  n = c(nrow(fs), nrow(fs), nrow(sv), nrow(sv)))
ci <- t(mapply(wilson, summ$k, summ$n))
summ$prop <- ci[, 1]; summ$lo <- ci[, 2]; summ$hi <- ci[, 3]
summ$pct  <- sprintf("%.0f%%", 100 * summ$prop)

brkA <- data.frame(x1 = c(1, 3), x2 = c(2, 4),
                   y = c(max(summ$hi[1:2]), max(summ$hi[3:4])) + 0.07,
                   lab = plab, sig = psig)

pA <- ggplot(summ, aes(grp, prop, fill = grp)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = bg_fs, alpha = bg_alpha) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = bg_sv, alpha = bg_alpha) +
  geom_col(width = 0.68, alpha = 0.9, colour = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.7, colour = "black") +
  geom_text(aes(y = prop, label = pct), vjust = 1.6, size = 5.2,
            colour = "white", fontface = "bold") +
  geom_segment(data = brkA, inherit.aes = FALSE,
               aes(x = x1, xend = x2, y = y, yend = y), colour = "grey20", linewidth = 0.6) +
  geom_text(data = brkA, inherit.aes = FALSE,
            aes(x = (x1 + x2) / 2, y = y, label = lab,
                fontface = ifelse(sig, "bold", "plain")),
            vjust = -0.4, size = 5.3, colour = "black") +
  scale_fill_manual(values = grp_cols) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                     limits = c(0, 1.16), expand = expansion(mult = c(0, 0))) +
  labs(title = "Protein-Complex Membership by Category",
       subtitle = "Fraction of proteins that are a subunit of \u22651 human CORUM complex (Wilson 95% CI)",
       x = NULL, y = "In a protein complex (%)") +
  theme_big +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, color = "black"))

# ---- assemble + save ---------------------------------------------------------
fig <- pA +
  plot_annotation(
    title = "Protein-Complex Membership \u2014 NMD-Escape Disease Genes vs Length-Matched Controls",
    subtitle = paste0("Frameshift vs FS Control (n=", nrow(fs),
                      " each) \u00b7 Stopgain vs Stopgain Control (n=", nrow(sv),
                      " each) \u00b7 Wilcoxon rank-sum (unpaired)"),
    caption = "P-values shown are unadjusted; BH-adjusted values in complex_comparison_stats.csv",
    theme = theme(plot.title    = element_text(face = "bold", size = 26, hjust = 0.5, color = "black"),
                  plot.subtitle = element_text(size = 16, hjust = 0.5, color = "black"),
                  plot.caption  = element_text(size = 13, hjust = 0.5, color = "black")))

save_fig("complex_comparison.pdf", fig, 10, 8)
save_fig("complex_comparison.png", fig, 10, 8)
message("Done: complex_comparison.pdf / .png and complex_comparison_stats.csv")
