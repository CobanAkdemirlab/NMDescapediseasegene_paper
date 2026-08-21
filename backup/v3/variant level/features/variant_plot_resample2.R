# =============================================================================
# Variant-level enrichment figure after 100 transcript-level resamplings
# Keep the actual median flag value (do NOT binarize 0.5)
#
# Workflow:
#   1. Read variant-level data
#   2. Repeat 100 times:
#        - randomly select 1 variant per transcript within each source
#   3. For each source × transcript × feature:
#        - compute median flag value across 100 resamples
#   4. For each group:
#        - compute mean median flag value × 100 as percentage
#   5. Compare groups with Wilcoxon rank-sum test
#   6. Draw barplots with:
#        - percentage labels
#        - n labels
#        - p-value brackets
# =============================================================================

rm(list = ls())
gc()

# ── 1. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("tidyverse", "patchwork", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

set.seed(123)

# ── 2. Read data ──────────────────────────────────────────────────────────────
df <- read_csv("Downloads/variants_all_with_motif_LCS_flags (1).csv", show_col_types = FALSE)

# ── 3. Basic cleaning ─────────────────────────────────────────────────────────
df <- df %>%
  filter(source %in% c("fs", "fs_control", "snv", "snv_control"))

if (!"transcript" %in% colnames(df)) {
  stop("The column 'transcript' is required but was not found.")
}

# ── 4. Define binary features to analyze ──────────────────────────────────────
feature_map <- tibble::tribble(
  ~feature,                  ~title,                                      ~ylab,                         ~panel,
  "variant_ppi_overlap",     "PPI Residue Enrichment",                    "Mean proportion (%)",         "A",
  "ptc_before_max_pfam_end", "Pfam Domain Enrichment",                    "Mean proportion (%)",         "B",
  "variant_nls_flag",        "Nuclear Localization Signal Enrichment",    "Mean proportion (%)",         "C",
  "variant_ptm_flag",        "Post-Translational Modification Enrichment","Mean proportion (%)",         "D",
  "variant_slim_flag",       "Short Linear Motif Enrichment",             "Mean proportion (%)",         "E",
  "variant_LCS_flag",        "Low Complexity Sequence Enrichment",        "Mean proportion (%)",         "F"
)

feature_map <- feature_map %>%
  filter(feature %in% colnames(df))

if (nrow(feature_map) == 0) {
  stop("None of the requested feature columns were found in the CSV.")
}

flag_cols <- feature_map$feature

# ── 5. Convert feature columns to numeric 0/1 if needed ──────────────────────
for (cc in flag_cols) {
  if (is.logical(df[[cc]])) {
    df[[cc]] <- as.numeric(df[[cc]])
  }
  
  if (is.character(df[[cc]])) {
    df[[cc]] <- case_when(
      df[[cc]] %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      df[[cc]] %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(df[[cc]]))
    )
  }
}

# ── 6. Nice group labels ──────────────────────────────────────────────────────
group_labels <- c(
  "fs" = "Frameshift",
  "fs_control" = "FS Control",
  "snv" = "Nonsense",
  "snv_control" = "Nonsense Control"
)

group_levels <- c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")

# ── 7. Resampling function ────────────────────────────────────────────────────
sample_one_variant_per_transcript <- function(data, iter_id) {
  data %>%
    group_by(source, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    mutate(iter = iter_id)
}

# ── 8. Run 100 resamplings ────────────────────────────────────────────────────
n_resamples <- 100

resampled_df <- map_dfr(seq_len(n_resamples), function(i) {
  sample_one_variant_per_transcript(df, i)
})

# ── 9. Long format ────────────────────────────────────────────────────────────
resampled_long <- resampled_df %>%
  select(iter, source, transcript, all_of(flag_cols)) %>%
  pivot_longer(
    cols = all_of(flag_cols),
    names_to = "feature",
    values_to = "value"
  )

# ── 10. Median flag across 100 resamplings (NO binarization) ─────────────────
median_df <- resampled_long %>%
  group_by(source, transcript, feature) %>%
  summarise(
    median_flag = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = recode(as.character(source), !!!group_labels),
    group = factor(group, levels = group_levels)
  )

# Save transcript-level median values
write_csv(
  median_df,
  "median_flag_by_transcript_after_100_resamples_continuous.csv"
)

# ── 11. Plot colors ───────────────────────────────────────────────────────────
group_colors <- c(
  "Frameshift"       = "#D07A3A",
  "FS Control"       = "#E2C7B5",
  "Nonsense"         = "#8ABA67",
  "Nonsense Control" = "#C9D6BE"
)

# ── 12. Helper functions ──────────────────────────────────────────────────────
format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

summarise_feature <- function(data, feat) {
  data %>%
    filter(feature == feat) %>%
    filter(!is.na(median_flag)) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_val = mean(median_flag, na.rm = TRUE),
      median_val = median(median_flag, na.rm = TRUE),
      percent = mean_val * 100,
      .groups = "drop"
    ) %>%
    right_join(
      tibble(group = factor(group_levels, levels = group_levels)),
      by = "group"
    ) %>%
    mutate(
      n = replace_na(n, 0),
      mean_val = replace_na(mean_val, 0),
      median_val = replace_na(median_val, 0),
      percent = replace_na(percent, 0)
    )
}

get_p <- function(data, feat, g1, g2) {
  dsub <- data %>%
    filter(feature == feat, group %in% c(g1, g2)) %>%
    filter(!is.na(median_flag))
  
  x <- dsub$median_flag[dsub$group == g1]
  y <- dsub$median_flag[dsub$group == g2]
  
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  
  tryCatch(
    wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

make_bracket_df <- function(summary_df, p1, p2) {
  ymax <- max(summary_df$percent, na.rm = TRUE)
  offset <- max(6, ymax * 0.12)
  
  tibble(
    xstart = c(1, 3),
    xend   = c(2, 4),
    xmid   = c(1.5, 3.5),
    y      = c(ymax + offset, ymax + offset),
    label  = c(format_p(p1), format_p(p2))
  )
}

# ── 13. Plot function ─────────────────────────────────────────────────────────
plot_bar <- function(data, feat, title, ylab, panel_letter) {
  df_sum <- summarise_feature(data, feat)
  
  p_fs  <- get_p(data, feat, "Frameshift", "FS Control")
  p_snv <- get_p(data, feat, "Nonsense", "Nonsense Control")
  
  bracket_df <- make_bracket_df(df_sum, p_fs, p_snv)
  
  ymax <- max(df_sum$percent, na.rm = TRUE)
  ylim_top <- max(100, ymax + max(15, ymax * 0.30))
  
  ggplot(df_sum, aes(x = group, y = percent, fill = group)) +
    geom_col(width = 0.62, color = "grey40", alpha = 0.95) +
    
    geom_text(
      aes(label = sprintf("%.1f%%", percent)),
      vjust = -0.45,
      size = 3.8,
      fontface = "bold",
      color = "grey15"
    ) +
    
    geom_label(
      aes(
        y = pmax(percent * 0.5, 5),
        label = paste0("n=", n)
      ),
      size = 3.0,
      label.size = 0,
      fill = "white",
      alpha = 0.78
    ) +
    
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey20"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xstart, y = y - 1.5, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey20"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xend, xend = xend, y = y - 1.5, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey20"
    ) +
    geom_text(
      data = bracket_df,
      aes(x = xmid, y = y + 2.2, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "grey15"
    ) +
    
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(
      limits = c(0, ylim_top),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      x = "Variant Category",
      y = ylab
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14.5),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 12, 12, 12)
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = panel_letter,
      hjust = -0.8, vjust = 1.4,
      size = 6,
      fontface = "bold"
    )
}

# ── 14. Generate all panels ───────────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[ feature_map$panel[i] ]] <- plot_bar(
    data = median_df,
    feat = feature_map$feature[i],
    title = feature_map$title[i],
    ylab = feature_map$ylab[i],
    panel_letter = feature_map$panel[i]
  )
}

# ── 15. Arrange combined figure ───────────────────────────────────────────────
ordered_panels <- c("A", "B", "C", "D", "E", "F")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Variant-level Enrichment of Protein Features After 100 Transcript-Level Resamplings",
    subtitle = paste(
      "100 transcript-level resamplings;",
      "1 random variant selected per transcript per resample;",
      "actual median flag values used;",
      "Wilcoxon rank-sum test for group comparisons"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0),
      plot.subtitle = element_text(size = 10.5, hjust = 0)
    )
  )

# ── 16. Summary table ─────────────────────────────────────────────────────────
summary_table <- map_dfr(feature_map$feature, function(feat_now) {
  title_now <- feature_map$title[feature_map$feature == feat_now][1]
  
  summarise_feature(median_df, feat_now) %>%
    mutate(
      feature = feat_now,
      title = title_now,
      .before = 1
    )
})

write_csv(
  summary_table,
  "Supplemental_Figure_3_variant_level_100_resamples_continuous_summary.csv"
)

# ── 17. Save outputs ──────────────────────────────────────────────────────────
ggsave(
  "Supplemental_Figure_3_variant_level_100_resamples_continuous.pdf",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

ggsave(
  "Supplemental_Figure_3_variant_level_100_resamples_continuous.png",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

# ── 18. Print summaries ───────────────────────────────────────────────────────
cat("\nCounts per source in original data:\n")
print(table(df$source))

cat("\nSummary table:\n")
print(summary_table)

cat("\nDone!\n")
cat("Output files:\n")
cat("  median_flag_by_transcript_after_100_resamples_continuous.csv\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_continuous.pdf\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_continuous.png\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_continuous_summary.csv\n")