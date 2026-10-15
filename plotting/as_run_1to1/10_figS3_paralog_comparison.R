# =============================================================================
# paralog_comparison.R
# Paralog-count comparison across the four categories
#   Frameshift / FS Control / Stopgain / Stopgain Control
# Styled to match repeat_homopolymer.R (black labels, larger fonts, RdBu palette,
# shaded per-class backgrounds, paired connecting lines, paired Wilcoxon brackets).
#
# Paired design (consistent with the other figures): uses matched_pairs_full.csv,
# i.e. each NMD-escape disease gene (case) vs its length-matched control gene.
#
# Paralog counts come from Ensembl (biomaRt) and are CACHED to paralog_counts.csv,
# so the slow query only runs once. Run this locally (needs internet + biomaRt).
#
# Outputs: paralog_comparison.pdf / .png, paralog_comparison_stats.csv,
#          paralog_counts.csv (cache)
# =============================================================================

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

# ---- shared styling (matches repeat_homopolymer.R) ---------------------------
grp_cols <- c("Frameshift"       = "#2166AC",
              "FS Control"       = "#92C5DE",
              "Stopgain"         = "#B2182B",
              "Stopgain Control" = "#F4A582")
bg_fs <- "#2166AC"; bg_sv <- "#B2182B"; bg_alpha <- 0.06
MINUS <- "\u2212"; DELTA <- "\u0394"

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

# ---- paralog counts via biomaRt (cached) -------------------------------------
# hgnc_symbol lives on the "feature" attribute page and paralog attributes on the
# "homolog" page; biomaRt forbids mixing pages in one query, so we go in 2 steps.
get_paralog_counts <- function(genes, cache = "paralog_counts.csv",
                               paralog_table = "paralog_pairs.tsv") {
  genes <- unique(genes[!is.na(genes) & genes != ""])

  ## 1) cached counts ----------------------------------------------------------
  if (file.exists(cache)) {
    message("Loading cached paralog counts from ", cache)
    return(read.csv(cache, stringsAsFactors = FALSE))
  }

  ## 2) OFFLINE: count from a local paralog-pairs table (no internet) ----------
  ##    Any TSV/CSV with a gene-symbol column and a paralog column
  ##    (one row per gene-paralog relationship).
  if (file.exists(paralog_table)) {
    message("Computing paralog counts from local table ", paralog_table)
    tt <- if (grepl("\\.tsv$|\\.txt$", paralog_table))
            read.delim(paralog_table, stringsAsFactors = FALSE)
          else read.csv(paralog_table, stringsAsFactors = FALSE)
    nm <- names(tt)
    gcol <- nm[grepl("hgnc|gene.?name|symbol|^gene$", nm, ignore.case = TRUE)][1]
    pcol <- nm[grepl("paralog", nm, ignore.case = TRUE)][1]
    if (is.na(gcol) || is.na(pcol))
      stop("Could not find gene-symbol and paralog columns in ", paralog_table,
           "; columns: ", paste(nm, collapse = ", "), call. = FALSE)
    tt  <- tt[!is.na(tt[[pcol]]) & tt[[pcol]] != "" & tt[[gcol]] %in% genes, , drop = FALSE]
    agg <- aggregate(tt[[pcol]], by = list(hgnc_symbol = tt[[gcol]]),
                     FUN = function(x) length(unique(x)))
    names(agg)[2] <- "n_paralog"
    # a paralog-relationship table lists pairs, so genes absent from it have 0 paralogs
    unmapped <- setdiff(genes, agg$hgnc_symbol)
    if (length(unmapped)) agg <- rbind(agg, data.frame(hgnc_symbol = unmapped, n_paralog = 0))
    write.csv(agg, cache, row.names = FALSE)
    message("Wrote ", cache, " (", nrow(agg), " genes; absent-from-table treated as 0 paralogs)")
    return(agg)
  }

  ## 3) biomaRt (needs internet + Bioconductor) -- wrapped for a clear error ---
  out <- tryCatch({
    if (!requireNamespace("biomaRt", quietly = TRUE))
      stop("package 'biomaRt' is not installed", call. = FALSE)
    suppressPackageStartupMessages(library(biomaRt))
    mart <- tryCatch(
      useEnsembl("genes", dataset = "hsapiens_gene_ensembl"),
      error = function(e)
        useEnsembl("genes", dataset = "hsapiens_gene_ensembl", mirror = "useast"))
    map <- getBM(c("hgnc_symbol", "ensembl_gene_id"),
                 filters = "hgnc_symbol", values = genes, mart = mart)
    map <- map[map$ensembl_gene_id != "" & map$hgnc_symbol != "", , drop = FALSE]
    egids <- unique(map$ensembl_gene_id)
    para <- tryCatch(
      getBM(c("ensembl_gene_id", "hsapiens_paralog_associated_gene_name"),
            filters = "ensembl_gene_id", values = egids, mart = mart),
      error = function(e)
        getBM(c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene"),
              filters = "ensembl_gene_id", values = egids, mart = mart))
    pcol <- setdiff(names(para), "ensembl_gene_id")[1]
    para$has <- para[[pcol]] != "" & !is.na(para[[pcol]])
    cnt <- aggregate(has ~ ensembl_gene_id, data = para, FUN = sum)
    names(cnt)[2] <- "n_paralog"
    miss <- setdiff(egids, cnt$ensembl_gene_id)
    if (length(miss)) cnt <- rbind(cnt, data.frame(ensembl_gene_id = miss, n_paralog = 0))
    m <- merge(map, cnt, by = "ensembl_gene_id", all.x = TRUE)
    m$n_paralog[is.na(m$n_paralog)] <- 0
    aggregate(n_paralog ~ hgnc_symbol, data = m, FUN = max)
  }, error = function(e) {
    stop("Could not obtain paralog counts.\n",
         "  biomaRt/Ensembl query failed: ", conditionMessage(e), "\n",
         "  Provide ONE of these next to the script and re-run:\n",
         "    (a) paralog_counts.csv  columns: hgnc_symbol, n_paralog\n",
         "    (b) paralog_pairs.tsv   gene-symbol + paralog columns (one row per relationship)\n",
         call. = FALSE)
  })
  unmapped <- setdiff(genes, out$hgnc_symbol)
  if (length(unmapped)) out <- rbind(out, data.frame(hgnc_symbol = unmapped, n_paralog = NA))
  write.csv(out, cache, row.names = FALSE)
  message("Wrote ", cache, " (", sum(is.na(out$n_paralog)), " unmapped of ", nrow(out), " genes)")
  out
}

# ---- load matched pairs ------------------------------------------------------
# ---- load matched pairs: new list via gene_flags_per_gene.csv, else matched_pairs_full.csv
load_pairs <- function() {
  if (file.exists("matched_pairs_full.csv")) {
    message("Pairs <- matched_pairs_full.csv")
    return(read.csv("matched_pairs_full.csv", stringsAsFactors = FALSE))
  }
  if (file.exists("gene_flags_per_gene.csv")) {
    message("Pairs <- gene_flags_per_gene.csv (matched)")
    g  <- read.csv("gene_flags_per_gene.csv", stringsAsFactors = FALSE)
    g  <- g[g$design == "matched", ]
    ca <- g[g$role == "case", ]; co <- g[g$role == "control", ]
    m  <- merge(
      data.frame(pair_id = ca$pair_id, stratum = ca$stratum,
                 case_hgnc = ca$hgnc_symbol, case_uniprot = ca$uniprot,
                 stringsAsFactors = FALSE),
      data.frame(pair_id = co$pair_id,
                 control_hgnc = co$hgnc_symbol, control_uniprot = co$uniprot,
                 stringsAsFactors = FALSE),
      by = "pair_id")
    m$case_group <- ifelse(m$stratum == "FS", "fs", "snv")
    return(m[, c("case_group","case_hgnc","control_hgnc","case_uniprot","control_uniprot")])
  }
  stop("No pairs found: need matched_pairs_full.csv or gene_flags_per_gene.csv", call. = FALSE)
}
pairs <- load_pairs()
stopifnot(all(c("case_group", "case_hgnc", "control_hgnc") %in% names(pairs)))

genes <- unique(c(pairs$case_hgnc, pairs$control_hgnc))
pc <- get_paralog_counts(genes)
lk <- setNames(pc$n_paralog, pc$hgnc_symbol)

pairs$case_par <- lk[pairs$case_hgnc]
pairs$ctrl_par <- lk[pairs$control_hgnc]
pairs$pair_id  <- seq_len(nrow(pairs))

ok <- !is.na(pairs$case_par) & !is.na(pairs$ctrl_par)
if (sum(!ok)) message(sum(!ok), " pair(s) dropped (paralog count unmapped for case or control)")
P <- pairs[ok, , drop = FALSE]
P$delta <- P$case_par - P$ctrl_par

fs <- P[P$case_group == "fs",  , drop = FALSE]
sv <- P[P$case_group == "snv", , drop = FALSE]

# ---- paired Wilcoxon (within length-matched pairs), BH-adjusted --------------
p_fs <- suppressWarnings(wilcox.test(fs$case_par, fs$ctrl_par, paired = TRUE)$p.value)
p_sv <- suppressWarnings(wilcox.test(sv$case_par, sv$ctrl_par, paired = TRUE)$p.value)
padj <- p.adjust(c(p_fs, p_sv), "BH")
praw <- c(p_fs, p_sv)
plab <- paste0("P = ", formatC(praw, format = "g", digits = 2))  # raw p shown on figure
psig <- praw < 0.05

stats <- data.frame(
  comparison   = c("Frameshift vs FS Control", "Stopgain vs Stopgain Control"),
  n_pairs      = c(nrow(fs), nrow(sv)),
  median_case  = c(median(fs$case_par), median(sv$case_par)),
  median_ctrl  = c(median(fs$ctrl_par), median(sv$ctrl_par)),
  median_delta = c(median(fs$delta),    median(sv$delta)),
  p_value      = c(p_fs, p_sv),
  p_BH         = padj)
write.csv(stats, "paralog_comparison_stats.csv", row.names = FALSE)
print(stats)

# ---- Panel A: paralog count by category (paired) -----------------------------
case_lab <- ifelse(P$case_group == "fs", "Frameshift", "Stopgain")
ctrl_lab <- ifelse(P$case_group == "fs", "FS Control", "Stopgain Control")
longA <- rbind(
  data.frame(pair_id = P$pair_id, grp = case_lab, value = P$case_par),
  data.frame(pair_id = P$pair_id, grp = ctrl_lab, value = P$ctrl_par))
longA$grp <- factor(longA$grp,
                    levels = c("Frameshift", "FS Control", "Stopgain", "Stopgain Control"))

yA <- max(longA$value, na.rm = TRUE)
brkA <- data.frame(x1 = c(1, 3), x2 = c(2, 4), y = yA * 1.02, lab = plab, sig = psig)

longA$xn <- as.numeric(longA$grp)
pA <- ggplot(longA, aes(xn, value, fill = grp, colour = grp)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = bg_fs, alpha = bg_alpha) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = bg_sv, alpha = bg_alpha) +
  geom_line(aes(group = pair_id), colour = "grey60", linewidth = 0.3, alpha = 0.45) +
  geom_boxplot(aes(group = grp), width = 0.55, outlier.shape = NA, alpha = 0.35, colour = "black") +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.6) +
  geom_segment(data = brkA, inherit.aes = FALSE,
               aes(x = x1, xend = x2, y = y, yend = y), colour = "grey20", linewidth = 0.6) +
  geom_text(data = brkA, inherit.aes = FALSE,
            aes(x = (x1 + x2) / 2, y = y, label = lab,
                fontface = ifelse(sig, "bold", "plain")),
            vjust = -0.4, size = 5.3, colour = "black") +
  scale_fill_manual(values = grp_cols) +
  scale_colour_manual(values = grp_cols) +
  scale_x_continuous(breaks = 1:4, labels = levels(longA$grp), limits = c(0.4, 4.6)) +
  scale_y_sqrt(expand = expansion(mult = c(0.02, 0.14))) +
  labs(title = "Paralog Count by Category",
       subtitle = "Paralogs per gene (Ensembl), paired within length-matched pairs",
       x = NULL, y = "Number of paralogs (sqrt scale)") +
  theme_big +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, color = "black"))

# ---- assemble + save ---------------------------------------------------------
fig <- pA +
  plot_annotation(
    title = "Paralog Count \u2014 NMD-Escape Disease Genes vs Length-Matched Controls",
    subtitle = paste0("Frameshift vs FS Control (n=", nrow(fs),
                      " pairs) \u00b7 Stopgain vs Stopgain Control (n=", nrow(sv),
                      " pairs) \u00b7 Paired Wilcoxon signed-rank"),
    caption = "P-values shown are unadjusted; BH-adjusted values in paralog_comparison_stats.csv",
    theme = theme(plot.title    = element_text(face = "bold", size = 26, hjust = 0.5, color = "black"),
                  plot.subtitle = element_text(size = 16, hjust = 0.5, color = "black"),
                  plot.caption  = element_text(size = 13, hjust = 0.5, color = "black")))

save_fig("paralog_comparison.pdf", fig, 10, 8)
save_fig("paralog_comparison.png", fig, 10, 8)
message("Done: paralog_comparison.pdf / .png and paralog_comparison_stats.csv")
