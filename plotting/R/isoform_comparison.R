# =============================================================================
# isoform_comparison.R
# Isoform (transcript) count comparison across the four categories
#   Frameshift / FS Control / Stopgain / Stopgain Control
# Mirrors paralog_comparison.R: paired within length-matched pairs, same styling
# (black labels, larger fonts, RdBu palette, shaded backgrounds, paired lines,
# paired Wilcoxon brackets, sqrt y-axis).
#
# Isoform counts come from Ensembl (biomaRt) and are CACHED to isoform_counts.csv,
# so the slow query only runs once. Run locally (needs internet + biomaRt).
#
# CODING_ONLY = TRUE counts protein-coding transcripts only (recommended for
# protein-isoform biology); set FALSE to count all annotated transcripts.
#
# Outputs: isoform_comparison.pdf / .png, isoform_comparison_stats.csv,
#          isoform_counts.csv (cache)
# =============================================================================

CODING_ONLY <- TRUE   # protein-coding transcripts only; FALSE = all transcripts

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

# ---- shared styling (matches paralog_comparison.R) ---------------------------
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

# ---- isoform (transcript) counts via biomaRt (cached) ------------------------
# hgnc_symbol, ensembl_transcript_id and transcript_biotype are all on the same
# attribute page, so a single query suffices; count distinct transcripts per gene.
get_isoform_counts <- function(genes, cache = "isoform_counts.csv",
                               coding_only = TRUE,
                               transcript_table = "ensembl_transcripts.tsv") {
  genes <- unique(genes[!is.na(genes) & genes != ""])

  ## 1) cached counts (fastest; produced by any earlier run) -------------------
  if (file.exists(cache)) {
    message("Loading cached isoform counts from ", cache)
    return(read.csv(cache, stringsAsFactors = FALSE))
  }

  ## 2) OFFLINE: count transcripts from a local Ensembl table (no internet) -----
  ##    Accepts any TSV/CSV with a gene-symbol column, a transcript-id column,
  ##    and (optionally) a transcript-biotype column. e.g. a biomaRt/GTF export.
  if (file.exists(transcript_table)) {
    message("Computing isoform counts from local table ", transcript_table)
    tt <- if (grepl("\\.tsv$|\\.txt$", transcript_table))
            read.delim(transcript_table, stringsAsFactors = FALSE)
          else read.csv(transcript_table, stringsAsFactors = FALSE)
    nm <- names(tt)
    gcol <- nm[grepl("hgnc|gene.?name|symbol", nm, ignore.case = TRUE)][1]
    tcol <- nm[grepl("transcript.?id|ensembl_transcript", nm, ignore.case = TRUE)][1]
    bcol <- nm[grepl("biotype", nm, ignore.case = TRUE)][1]
    if (is.na(gcol) || is.na(tcol))
      stop("Could not find gene-symbol and transcript-id columns in ",
           transcript_table, "; columns: ", paste(nm, collapse = ", "), call. = FALSE)
    if (coding_only && !is.na(bcol))
      tt <- tt[tt[[bcol]] == "protein_coding", , drop = FALSE]
    tt  <- tt[tt[[gcol]] %in% genes, , drop = FALSE]
    agg <- aggregate(tt[[tcol]], by = list(hgnc_symbol = tt[[gcol]]),
                     FUN = function(x) length(unique(x)))
    names(agg)[2] <- "n_isoform"
    unmapped <- setdiff(genes, agg$hgnc_symbol)
    if (length(unmapped)) agg <- rbind(agg, data.frame(hgnc_symbol = unmapped, n_isoform = NA))
    write.csv(agg, cache, row.names = FALSE)
    message("Wrote ", cache, " (", sum(is.na(agg$n_isoform)), " unmapped of ", nrow(agg), ")")
    return(agg)
  }

  ## 3) biomaRt (needs internet + Bioconductor) -- wrapped so failures are clear
  agg <- tryCatch({
    if (!requireNamespace("biomaRt", quietly = TRUE))
      stop("package 'biomaRt' is not installed", call. = FALSE)
    suppressPackageStartupMessages(library(biomaRt))
    mart <- tryCatch(
      useEnsembl("genes", dataset = "hsapiens_gene_ensembl"),
      error = function(e)
        useEnsembl("genes", dataset = "hsapiens_gene_ensembl", mirror = "useast"))
    bm <- getBM(c("hgnc_symbol", "ensembl_transcript_id", "transcript_biotype"),
                filters = "hgnc_symbol", values = genes, mart = mart)
    bm <- bm[bm$hgnc_symbol != "" & bm$ensembl_transcript_id != "", , drop = FALSE]
    if (coding_only) bm <- bm[bm$transcript_biotype == "protein_coding", , drop = FALSE]
    a <- aggregate(ensembl_transcript_id ~ hgnc_symbol, data = bm,
                   FUN = function(x) length(unique(x)))
    names(a)[2] <- "n_isoform"; a
  }, error = function(e) {
    stop("Could not obtain isoform counts.\n",
         "  biomaRt/Ensembl query failed: ", conditionMessage(e), "\n",
         "  Provide ONE of these next to the script and re-run:\n",
         "    (a) isoform_counts.csv       columns: hgnc_symbol, n_isoform\n",
         "    (b) ensembl_transcripts.tsv  gene-symbol, transcript-id (+ biotype) columns\n",
         call. = FALSE)
  })
  unmapped <- setdiff(genes, agg$hgnc_symbol)
  if (length(unmapped)) agg <- rbind(agg, data.frame(hgnc_symbol = unmapped, n_isoform = NA))
  write.csv(agg, cache, row.names = FALSE)
  message("Wrote ", cache, " (", sum(is.na(agg$n_isoform)), " unmapped of ", nrow(agg),
          " genes; coding_only = ", coding_only, ")")
  agg
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
ic <- get_isoform_counts(genes, coding_only = CODING_ONLY)
lk <- setNames(ic$n_isoform, ic$hgnc_symbol)

pairs$case_iso <- lk[pairs$case_hgnc]
pairs$ctrl_iso <- lk[pairs$control_hgnc]
pairs$pair_id  <- seq_len(nrow(pairs))

ok <- !is.na(pairs$case_iso) & !is.na(pairs$ctrl_iso)
if (sum(!ok)) message(sum(!ok), " pair(s) dropped (isoform count unmapped for case or control)")
P <- pairs[ok, , drop = FALSE]

fs <- P[P$case_group == "fs",  , drop = FALSE]
sv <- P[P$case_group == "snv", , drop = FALSE]

# ---- paired Wilcoxon (within length-matched pairs), BH-adjusted --------------
p_fs <- suppressWarnings(wilcox.test(fs$case_iso, fs$ctrl_iso, paired = TRUE)$p.value)
p_sv <- suppressWarnings(wilcox.test(sv$case_iso, sv$ctrl_iso, paired = TRUE)$p.value)
padj <- p.adjust(c(p_fs, p_sv), "BH")
praw <- c(p_fs, p_sv)
plab <- paste0("P = ", formatC(praw, format = "g", digits = 2))  # raw p shown on figure
psig <- praw < 0.05

stats <- data.frame(
  comparison   = c("Frameshift vs FS Control", "Stopgain vs Stopgain Control"),
  n_pairs      = c(nrow(fs), nrow(sv)),
  median_case  = c(median(fs$case_iso), median(sv$case_iso)),
  median_ctrl  = c(median(fs$ctrl_iso), median(sv$ctrl_iso)),
  median_delta = c(median(fs$case_iso - fs$ctrl_iso), median(sv$case_iso - sv$ctrl_iso)),
  p_value      = c(p_fs, p_sv),
  p_BH         = padj)
write.csv(stats, "isoform_comparison_stats.csv", row.names = FALSE)
print(stats)

# ---- Panel: isoform count by category (paired) -------------------------------
case_lab <- ifelse(P$case_group == "fs", "Frameshift", "Stopgain")
ctrl_lab <- ifelse(P$case_group == "fs", "FS Control", "Stopgain Control")
longA <- rbind(
  data.frame(pair_id = P$pair_id, grp = case_lab, value = P$case_iso),
  data.frame(pair_id = P$pair_id, grp = ctrl_lab, value = P$ctrl_iso))
longA$grp <- factor(longA$grp,
                    levels = c("Frameshift", "FS Control", "Stopgain", "Stopgain Control"))

yA <- max(longA$value, na.rm = TRUE)
brkA <- data.frame(x1 = c(1, 3), x2 = c(2, 4), y = yA * 1.02, lab = plab, sig = psig)

iso_kind <- if (CODING_ONLY) "protein-coding" else "all"
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
  labs(title = "Isoform Count by Category",
       subtitle = paste0("Transcript isoforms per gene (Ensembl, ", iso_kind,
                          "), paired within length-matched pairs"),
       x = NULL, y = "Number of isoforms (sqrt scale)") +
  theme_big +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, color = "black"))

# ---- assemble + save ---------------------------------------------------------
fig <- pA +
  plot_annotation(
    title = "Isoform Count \u2014 NMD-Escape Disease Genes vs Length-Matched Controls",
    subtitle = paste0("Frameshift vs FS Control (n=", nrow(fs),
                      " pairs) \u00b7 Stopgain vs Stopgain Control (n=", nrow(sv),
                      " pairs) \u00b7 Paired Wilcoxon signed-rank"),
    caption = "P-values shown are unadjusted; BH-adjusted values in isoform_comparison_stats.csv",
    theme = theme(plot.title    = element_text(face = "bold", size = 26, hjust = 0.5, color = "black"),
                  plot.subtitle = element_text(size = 16, hjust = 0.5, color = "black"),
                  plot.caption  = element_text(size = 13, hjust = 0.5, color = "black")))

save_fig("isoform_comparison.pdf", fig, 10, 8)
save_fig("isoform_comparison.png", fig, 10, 8)
message("Done: isoform_comparison.pdf / .png and isoform_comparison_stats.csv")
