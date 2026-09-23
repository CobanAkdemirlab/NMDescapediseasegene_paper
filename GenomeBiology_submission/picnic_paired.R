#!/usr/bin/env Rscript
# =============================================================================
# PICNIC Score >= 0.5  --  Paired Analysis Across 4 Gene Groups
# R port of preview_picnic.py (faithful reproduction of pipeline + 4-panel figure)
#
# Inputs:
#   PICNIC-9606-data.csv   columns: ID, "PICNIC score", Genes, "Organism ID"
#   matching_log_v2.csv    columns: study_gene, ctrl_gene, group ('FS' / 'Stopgain')
#
# Pipeline:
#   1. Keep human PICNIC entries (Organism ID == 9606); map score -> HGNC symbol
#      using the first space-separated gene token; per gene take the highest score.
#   2. Binary outcome PICNIC >= 0.5.
#   3. Paired McNemar's exact (binomial on discordant pairs) for the binary outcome;
#      paired Wilcoxon signed-rank for the continuous score. BH-adjust within test.
#   4. Plot: (A) bar %, (B) violin+box, (C) ECDF, (D) paired scatter w/ quadrants.
#
# Usage:
#   Rscript picnic_paired.R [picnic_csv] [pairs_csv] [out_png]
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork)
})

# ---- Config & paths ---------------------------------------------------------
THRESH <- 0.5
args   <- commandArgs(trailingOnly = TRUE)
PICNIC_CSV <- if (length(args) >= 1) args[1] else "PICNIC-9606-data.csv"
PAIRS_CSV  <- if (length(args) >= 2) args[2] else "matched_pairs_full.csv"
OUT_PNG    <- if (length(args) >= 3) args[3] else "picnic_paired.png"
OUT_STATS  <- file.path(dirname(OUT_PNG), "picnic_stats.csv")

group_order <- c("Frameshift", "FS Control", "Stopgain", "Stopgain Control")
palette <- c("Frameshift" = "#2166AC", "FS Control" = "#92C5DE",
             "Stopgain"   = "#B2182B", "Stopgain Control" = "#F4A582")
# scatter quadrant colours
cat_colors <- c(both_high = "#2166AC", study_only = "#C03B3B",
                ctrl_only = "#F4A582", both_low = "#999999")

# ---- Load PICNIC ------------------------------------------------------------
message("Reading PICNIC...")
picnic_raw <- read_csv(PICNIC_CSV, show_col_types = FALSE)
picnic_raw <- picnic_raw[picnic_raw[["Organism ID"]] == 9606, ]    # human only
message(sprintf("  PICNIC: %d human entries", nrow(picnic_raw)))

first_gene <- function(s) {
  if (is.na(s)) return(NA_character_)
  parts <- strsplit(trimws(as.character(s)), "\\s+")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) NA_character_ else parts[1]
}

picnic <- picnic_raw %>%
  transmute(gene = vapply(Genes, first_gene, character(1)),
            picnic_score = `PICNIC score`) %>%
  filter(!is.na(gene), !is.na(picnic_score)) %>%
  arrange(desc(picnic_score)) %>%             # per gene: keep highest score
  distinct(gene, .keep_all = TRUE) %>%
  mutate(high = picnic_score >= THRESH)

message(sprintf("  Per-gene PICNIC scores: %d unique HGNC genes", nrow(picnic)))
overall_high_pct <- 100 * mean(picnic$high)
message(sprintf("  Genes with score >= %.1f: %d (%.1f%%)",
                THRESH, sum(picnic$high), overall_high_pct))

score_lkp <- setNames(picnic$picnic_score, picnic$gene)
high_lkp  <- setNames(picnic$high,          picnic$gene)

# ---- Load pairs -------------------------------------------------------------
# Prefer the new matched list: gene_flags_per_gene.csv (design/pair_id/role/
# stratum/hgnc_symbol) -> matched_pairs_full.csv -> legacy matching_log_v2.csv.
load_pairs <- function() {
  if (file.exists("gene_flags_per_gene.csv")) {
    g <- read_csv("gene_flags_per_gene.csv", show_col_types = FALSE) %>%
      filter(design == "matched")
    st <- g %>% filter(role == "case")    %>% transmute(pair_id, stratum, study_gene = hgnc_symbol)
    ct <- g %>% filter(role == "control") %>% transmute(pair_id, ctrl_gene = hgnc_symbol)
    message("Pairs <- gene_flags_per_gene.csv (matched)")
    inner_join(st, ct, by = "pair_id") %>%
      transmute(study_gene, ctrl_gene,
                group = dplyr::recode(stratum, FS = "FS", SNV = "Stopgain"))
  } else if (file.exists("matched_pairs_full.csv")) {
    message("Pairs <- matched_pairs_full.csv")
    read_csv("matched_pairs_full.csv", show_col_types = FALSE) %>%
      transmute(study_gene = case_hgnc, ctrl_gene = control_hgnc,
                group = dplyr::recode(as.character(case_group), fs = "FS", snv = "Stopgain"))
  } else {
    message("Pairs <- ", PAIRS_CSV, " (legacy)")
    read_csv(PAIRS_CSV, show_col_types = FALSE)
  }
}
pairs <- load_pairs() %>%
  mutate(s_score = unname(score_lkp[study_gene]),
         c_score = unname(score_lkp[ctrl_gene]),
         s_high  = unname(high_lkp[study_gene]),
         c_high  = unname(high_lkp[ctrl_gene]))

ok <- function(g) pairs %>%
  filter(group == g, !is.na(s_score), !is.na(c_score))
fs_ok  <- ok("FS")
snv_ok <- ok("Stopgain")
message(sprintf("\nPairs with PICNIC on both sides:\n  Frameshift: %d / %d\n  Stopgain  : %d / %d",
                nrow(fs_ok), sum(pairs$group == "FS"),
                nrow(snv_ok), sum(pairs$group == "Stopgain")))

# ---- Stats ------------------------------------------------------------------
mcnemar_exact <- function(s_high, c_high) {
  s <- as.integer(s_high); c <- as.integer(c_high)
  a  <- sum(s == 1 & c == 1)   # both high
  b  <- sum(s == 1 & c == 0)   # study only
  cc <- sum(s == 0 & c == 1)   # ctrl only
  d  <- sum(s == 0 & c == 0)   # both low
  p  <- if ((b + cc) == 0) NA_real_
        else binom.test(min(b, cc), b + cc, p = 0.5,
                        alternative = "two.sided")$p.value
  list(a = a, b = b, c = cc, d = d, n = length(s), p = p,
       s_high_n = sum(s), c_high_n = sum(c))
}

paired_w <- function(s, c) {
  s <- as.numeric(s); c <- as.numeric(c)
  if (length(s) < 3 || all(s == c)) return(NA_real_)
  # match scipy.wilcoxon defaults for n>25: normal approx, no continuity correction,
  # zeros dropped (Wilcoxon method)
  suppressWarnings(
    wilcox.test(s, c, paired = TRUE, exact = FALSE, correct = FALSE)$p.value
  )
}

mc_fs  <- mcnemar_exact(fs_ok$s_high,  fs_ok$c_high)
mc_snv <- mcnemar_exact(snv_ok$s_high, snv_ok$c_high)
w_fs   <- paired_w(fs_ok$s_score,  fs_ok$c_score)
w_snv  <- paired_w(snv_ok$s_score, snv_ok$c_score)

# BH within test (2 comparisons each)
mc_padj <- p.adjust(c(mc_fs$p, mc_snv$p), method = "BH")
w_padj  <- p.adjust(c(w_fs,    w_snv),    method = "BH")
mc_fs_padj <- mc_padj[1]; mc_snv_padj <- mc_padj[2]
w_fs_padj  <- w_padj[1];  w_snv_padj  <- w_padj[2]

message("\n=== Paired McNemar (PICNIC >= 0.5) ===")
message(sprintf("  Frameshift vs FS Ctrl   n=%3d  case_high=%3d/%-3d ctrl_high=%3d/%-3d raw p=%.3g BH=%.3g",
                mc_fs$n, mc_fs$s_high_n, mc_fs$n, mc_fs$c_high_n, mc_fs$n, mc_fs$p, mc_fs_padj))
message(sprintf("  Stopgain   vs SV Ctrl   n=%3d  case_high=%3d/%-3d ctrl_high=%3d/%-3d raw p=%.3g BH=%.3g",
                mc_snv$n, mc_snv$s_high_n, mc_snv$n, mc_snv$c_high_n, mc_snv$n, mc_snv$p, mc_snv_padj))
message("=== Paired Wilcoxon (continuous PICNIC) ===")
message(sprintf("  Frameshift vs FS Ctrl   raw p=%.3g BH=%.3g", w_fs,  w_fs_padj))
message(sprintf("  Stopgain   vs SV Ctrl   raw p=%.3g BH=%.3g", w_snv, w_snv_padj))

# ---- Save stats -------------------------------------------------------------
sig <- function(p) ifelse(is.na(p), "NA",
                   ifelse(p < 0.001, "***",
                   ifelse(p < 0.01,  "**",
                   ifelse(p < 0.05,  "*", "ns"))))
stats_df <- bind_rows(
  tibble(test = "McNemar (PICNIC>=0.5)", comparison = "Frameshift vs FS Control",
         n_pairs = mc_fs$n, case_high = mc_fs$s_high_n, ctrl_high = mc_fs$c_high_n,
         concord_both_high = mc_fs$a, concord_both_low = mc_fs$d,
         discord_case_only = mc_fs$b, discord_ctrl_only = mc_fs$c,
         p_raw = mc_fs$p, p_BH = mc_fs_padj),
  tibble(test = "McNemar (PICNIC>=0.5)", comparison = "Stopgain vs Stopgain Control",
         n_pairs = mc_snv$n, case_high = mc_snv$s_high_n, ctrl_high = mc_snv$c_high_n,
         concord_both_high = mc_snv$a, concord_both_low = mc_snv$d,
         discord_case_only = mc_snv$b, discord_ctrl_only = mc_snv$c,
         p_raw = mc_snv$p, p_BH = mc_snv_padj),
  tibble(test = "Paired Wilcoxon", comparison = "Frameshift vs FS Control",
         n_pairs = mc_fs$n, p_raw = w_fs, p_BH = w_fs_padj),
  tibble(test = "Paired Wilcoxon", comparison = "Stopgain vs Stopgain Control",
         n_pairs = mc_snv$n, p_raw = w_snv, p_BH = w_snv_padj)
) %>% mutate(signif = sig(p_BH))
write_csv(stats_df, OUT_STATS)

# ---- Long frame for plotting ------------------------------------------------
long_block <- function(df, study_lbl, ctrl_lbl) {
  tibble(Group = c(rep(study_lbl, nrow(df)), rep(ctrl_lbl, nrow(df))),
         score = c(df$s_score, df$c_score),
         high  = c(df$s_high,  df$c_high))
}
long_df <- bind_rows(long_block(fs_ok,  "Frameshift", "FS Control"),
                     long_block(snv_ok, "Stopgain",   "Stopgain Control")) %>%
  mutate(Group = factor(Group, levels = group_order))

fmt_p <- function(p) {
  if (is.na(p)) return("p = n/a")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

# =============================================================================
# Shared theme
# =============================================================================
INK <- "#1a1a1a"   # near-black used for all text (no grey labels)
base_theme <- theme_classic(base_size = 17) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 21, colour = INK,
                                     margin = margin(b = 4)),
        plot.subtitle = element_text(hjust = 0.5, size = 14, colour = INK,
                                     margin = margin(b = 8)),
        axis.title    = element_text(face = "bold", size = 17, colour = INK),
        axis.text     = element_text(colour = INK, size = 14),
        axis.line     = element_line(colour = INK, linewidth = 0.5),
        axis.ticks    = element_line(colour = INK, linewidth = 0.5),
        plot.tag      = element_text(face = "bold", size = 26, colour = INK),
        plot.margin   = margin(10, 14, 10, 10),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3))

# shaded backgrounds: FS pair region (groups 1-2) blue, Stopgain region (3-4) red
shade_bg <- list(
  annotate("rect", xmin = 0.4, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = "#2166AC", alpha = 0.04),
  annotate("rect", xmin = 2.5, xmax = 4.6, ymin = -Inf, ymax = Inf,
           fill = "#B2182B", alpha = 0.04)
)

# bracket helper: returns annotate layers for a significance bracket
bracket <- function(x1, x2, y, label, tick, sig_bold) {
  fw <- if (isTRUE(sig_bold)) "bold" else "plain"
  list(
    annotate("segment", x = x1, xend = x1, y = y - tick, yend = y, linewidth = 0.6, colour = INK),
    annotate("segment", x = x2, xend = x2, y = y - tick, yend = y, linewidth = 0.6, colour = INK),
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.6, colour = INK),
    annotate("text", x = (x1 + x2) / 2, y = y, label = label, colour = INK,
             vjust = -0.3, size = 5, fontface = fw, lineheight = 0.9)
  )
}

# =============================================================================
# (A) Bar chart
# =============================================================================
bar_df <- tibble(
  Group = factor(group_order, levels = group_order),
  x     = 1:4,
  pct   = c(100 * mc_fs$s_high_n  / mc_fs$n,  100 * mc_fs$c_high_n  / mc_fs$n,
            100 * mc_snv$s_high_n / mc_snv$n, 100 * mc_snv$c_high_n / mc_snv$n),
  count = c(sprintf("n=%d/%d", mc_fs$s_high_n,  mc_fs$n),
            sprintf("n=%d/%d", mc_fs$c_high_n,  mc_fs$n),
            sprintf("n=%d/%d", mc_snv$s_high_n, mc_snv$n),
            sprintf("n=%d/%d", mc_snv$c_high_n, mc_snv$n)))
ymax_bar <- max(bar_df$pct)
yb_bar   <- min(95, ymax_bar + 12)

pA <- ggplot(bar_df, aes(x, pct, fill = Group)) +
  shade_bg +
  geom_col(width = 0.7, colour = "#333333", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.6,
            fontface = "bold", size = 5.4, colour = INK) +
  geom_text(aes(y = pmax(pct - 5, 4), label = count),
            colour = "white", fontface = "bold", size = 4.6) +
  bracket(1, 2, yb_bar, sprintf("%s\n(n=%d pairs)", fmt_p(mc_fs_padj),  mc_fs$n),
          1.5, !is.na(mc_fs_padj)  && mc_fs_padj  < 0.05) +
  bracket(3, 4, yb_bar, sprintf("%s\n(n=%d pairs)", fmt_p(mc_snv_padj), mc_snv$n),
          1.5, !is.na(mc_snv_padj) && mc_snv_padj < 0.05) +
  scale_fill_manual(values = palette, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = group_order, expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, max(95, ymax_bar + 22)),
                     breaks = seq(0, 90, 20), labels = paste0(seq(0, 90, 20), "%"),
                     expand = c(0, 0)) +
  labs(title = sprintf("Proportion of Genes with PICNIC \u2265 %.1f", THRESH),
       subtitle = sprintf("Paired McNemar's exact (BH) \u00b7 %.1f%% of all PICNIC proteins score \u2265 %.1f",
                          overall_high_pct, THRESH),
       x = NULL, y = sprintf("%% of pairs with PICNIC \u2265 %.1f", THRESH), tag = "A") +
  base_theme

# =============================================================================
# (B) Violin + box + jitter
# =============================================================================
set.seed(7)
ymax_v <- max(long_df$score, na.rm = TRUE)
yb_v   <- min(1.10, ymax_v + 0.10)
long_df$xnum <- as.integer(long_df$Group)

pB <- ggplot(long_df, aes(xnum, score, fill = Group, colour = Group, group = xnum)) +
  shade_bg +
  geom_hline(yintercept = THRESH, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_violin(trim = FALSE, alpha = 0.55, linewidth = 0.5, scale = "width", width = 0.85) +
  geom_jitter(width = 0.13, height = 0, size = 1.1, alpha = 0.55) +
  geom_boxplot(width = 0.18, fill = "white", colour = "#222222",
               outlier.shape = NA, linewidth = 0.5, alpha = 0.95,
               fatten = 2.2) +
  annotate("text", x = 3.6, y = THRESH + 0.01, label = sprintf("PICNIC = %.1f", THRESH),
           hjust = 1, vjust = 0, size = 4.2, colour = INK, fontface = "bold") +
  bracket(1, 2, yb_v, sprintf("%s\n(n=%d pairs)", fmt_p(w_fs_padj),  mc_fs$n),
          0.012, !is.na(w_fs_padj)  && w_fs_padj  < 0.05) +
  bracket(3, 4, yb_v, sprintf("%s\n(n=%d pairs)", fmt_p(w_snv_padj), mc_snv$n),
          0.012, !is.na(w_snv_padj) && w_snv_padj < 0.05) +
  scale_fill_manual(values = palette, guide = "none") +
  scale_colour_manual(values = palette, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = group_order, expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(-0.05, 1.30), breaks = seq(0, 1.25, 0.25)) +
  labs(title = "PICNIC Score Distribution",
       subtitle = sprintf("Dashed line = threshold (%.1f) \u00b7 Paired Wilcoxon signed-rank (BH)", THRESH),
       x = NULL, y = "PICNIC score", tag = "B") +
  base_theme

# =============================================================================
# (C) ECDF
# =============================================================================
medians <- long_df %>% group_by(Group) %>%
  summarise(m = median(score, na.rm = TRUE), .groups = "drop")
lty_map <- c("Frameshift" = "solid", "FS Control" = "dashed",
             "Stopgain" = "solid", "Stopgain Control" = "dashed")

pC <- ggplot(long_df, aes(score, colour = Group, linetype = Group)) +
  geom_vline(xintercept = THRESH, linetype = "dashed", colour = "grey30", linewidth = 0.45) +
  geom_vline(data = medians, aes(xintercept = m, colour = Group),
             linetype = "dotted", linewidth = 0.5, alpha = 0.7, show.legend = FALSE) +
  stat_ecdf(geom = "step", linewidth = 1.0, pad = FALSE) +
  annotate("text", x = THRESH + 0.01, y = 0.02, label = sprintf("%.1f", THRESH),
           hjust = 0, size = 4.2, colour = INK, fontface = "bold") +
  scale_colour_manual(values = palette, name = NULL) +
  scale_linetype_manual(values = lty_map, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%")) +
  labs(title = "ECDF \u2014 PICNIC Score",
       subtitle = "Dashed = threshold (0.5) \u00b7 Dotted = group medians",
       x = "PICNIC score", y = "Cumulative proportion", tag = "C") +
  base_theme +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.width = unit(1.8, "lines"),
        legend.text = element_text(size = 13, colour = INK),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))

# =============================================================================
# (D) Paired scatter with quadrant counts (FS pairs | Stopgain pairs)
# =============================================================================
classify <- function(df, facet_lbl) {
  if (nrow(df) == 0) return(tibble())
  df %>% transmute(
    facet = facet_lbl, c_score, s_score,
    cat = case_when(
      s_score >= THRESH & c_score >= THRESH ~ "both_high",
      s_score >= THRESH & c_score <  THRESH ~ "study_only",
      s_score <  THRESH & c_score >= THRESH ~ "ctrl_only",
      TRUE                                  ~ "both_low"))
}
scatter_df <- bind_rows(classify(fs_ok,  "FS pairs"),
                        classify(snv_ok, "Stopgain pairs")) %>%
  mutate(facet = factor(facet, levels = c("FS pairs", "Stopgain pairs")),
         cat   = factor(cat, levels = names(cat_colors)))

# per-facet quadrant counts placed at fixed corners
quad_pos <- tibble(
  cat   = c("both_high", "study_only", "ctrl_only", "both_low"),
  x     = c(0.78, 0.20, 0.78, 0.20),
  y     = c(0.78, 0.78, 0.20, 0.20),
  lbl0  = c("Both high", "Study high only", "Control high only", "Both low"))
counts_df <- scatter_df %>% count(facet, cat, name = "n") %>%
  tidyr::complete(facet, cat = factor(names(cat_colors), levels = names(cat_colors)),
                  fill = list(n = 0)) %>%
  left_join(quad_pos, by = "cat") %>%
  mutate(label = sprintf("%s\nn = %d", lbl0, n),
         txtcol = ifelse(cat == "ctrl_only", "#C2733E",
                  ifelse(cat == "both_low", INK, cat_colors[as.character(cat)])))

pD <- ggplot(scatter_df, aes(c_score, s_score, colour = cat)) +
  geom_hline(yintercept = THRESH, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = THRESH, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  geom_point(size = 2.1, alpha = 0.8) +
  geom_label(data = counts_df, aes(x, y, label = label),
             colour = counts_df$txtcol, inherit.aes = FALSE,
             fontface = "bold", size = 3.9, label.size = 0.3,
             label.padding = unit(0.22, "lines"), fill = "white") +
  facet_wrap(~ facet, nrow = 1) +
  scale_colour_manual(values = cat_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Matched Pair PICNIC Scores",
       subtitle = "Each point = one study-control pair \u00b7 Dashed lines = threshold (0.5)",
       x = "Control gene PICNIC score", y = "Study gene PICNIC score", tag = "D") +
  base_theme +
  theme(strip.background = element_rect(fill = "#ededed", colour = "#888888"),
        strip.text = element_text(face = "bold", size = 15, colour = INK),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        panel.spacing = unit(1.1, "lines"))

# =============================================================================
# Compose 2x2 and save
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = sprintf("PICNIC Score \u2265 %.1f \u2014 Paired Analysis Across 4 Gene Groups", THRESH),
    subtitle = sprintf("Frameshift (n=%d pairs) \u00b7 Stopgain (n=%d pairs) \u00b7 McNemar & paired Wilcoxon, BH adjusted",
                       mc_fs$n, mc_snv$n),
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 25, colour = INK),
                  plot.subtitle = element_text(hjust = 0.5, size = 15, colour = INK)))

ggsave(OUT_PNG, fig, width = 20, height = 15, dpi = 170, bg = "white")
message(sprintf("\nFigure saved -> %s\nStats saved  -> %s", OUT_PNG, OUT_STATS))
