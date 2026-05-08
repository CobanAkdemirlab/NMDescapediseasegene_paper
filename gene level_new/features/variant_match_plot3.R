# =============================================================================
# Variant-level enrichment of protein features after 100 transcript-level
# matched resamplings
#
# Y-axis:
#   mean proportion (%) = mean(transcript-level median_flag) * 100
#
# P-value:
#   paired Wilcoxon signed-rank test
#   using matched shared transcripts within each pair
#
# Matching:
#   - fs vs fs_control: shared transcripts only
#   - snv vs snv_control: shared transcripts only
#
# Resampling:
#   - 100 times
#   - 1 random variant per transcript per group per resample
# =============================================================================

rm(list = ls())
gc()

# ── 1. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("tidyverse", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

set.seed(123)

# ── 2. Read data ──────────────────────────────────────────────────────────────
variant_df <- readr::read_csv(
  "Downloads/variants_all_with_motif_LCS_flags (1).csv",
  show_col_types = FALSE
)

# ── 3. Keep groups of interest ────────────────────────────────────────────────
variant_df <- variant_df %>%
  filter(source %in% c("fs", "fs_control", "snv", "snv_control")) %>%
  filter(!is.na(transcript))

# ── 4. Features to plot ───────────────────────────────────────────────────────
feature_map <- tibble::tribble(
  ~feature,                  ~title,                                      ~panel,
  "variant_ppi_overlap",     "PPI Residue Enrichment",                    "A",
  "ptc_before_max_pfam_end", "Pfam Domain Enrichment",                    "B",
  "variant_nls_flag",        "Nuclear Localization Signal Enrichment",    "C",
  "variant_ptm_flag",        "Post-Translational Modification Enrichment","D",
  "variant_slim_flag",       "Short Linear Motif Enrichment",             "E",
  "variant_LCS_flag",        "Low Complexity Sequence Enrichment",        "F"
) %>%
  filter(feature %in% colnames(variant_df))

if (nrow(feature_map) == 0) {
  stop("None of the requested feature columns were found in the data.")
}

flag_cols <- feature_map$feature

# ── 5. Convert binary features to numeric 0/1 ────────────────────────────────
to_binary_numeric <- function(x) {
  if (is.logical(x)) return(as.numeric(x))
  if (is.character(x)) {
    return(case_when(
      x %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      x %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(x))
    ))
  }
  if (is.factor(x)) {
    x2 <- as.character(x)
    return(case_when(
      x2 %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      x2 %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(x2))
    ))
  }
  as.numeric(x)
}

for (cc in flag_cols) {
  variant_df[[cc]] <- to_binary_numeric(variant_df[[cc]])
}

# ── 6. Shared transcripts within each matched pair ───────────────────────────
get_shared_transcripts <- function(data, case_group, control_group) {
  intersect(
    unique(data$transcript[data$source == case_group]),
    unique(data$transcript[data$source == control_group])
  )
}

fs_shared  <- get_shared_transcripts(variant_df, "fs", "fs_control")
snv_shared <- get_shared_transcripts(variant_df, "snv", "snv_control")

cat("Shared transcripts:\n")
cat("  fs vs fs_control:", length(fs_shared), "\n")
cat("  snv vs snv_control:", length(snv_shared), "\n")

# ── 7. One matched resample ───────────────────────────────────────────────────
sample_once_per_group_transcript <- function(data, iter_id) {
  
  fs_part <- data %>%
    filter(source %in% c("fs", "fs_control"),
           transcript %in% fs_shared) %>%
    group_by(source, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  snv_part <- data %>%
    filter(source %in% c("snv", "snv_control"),
           transcript %in% snv_shared) %>%
    group_by(source, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  bind_rows(fs_part, snv_part) %>%
    mutate(iter = iter_id)
}

# ── 8. Repeat 100 times ───────────────────────────────────────────────────────
n_resamples <- 100

resampled_df <- purrr::map_dfr(
  seq_len(n_resamples),
  ~sample_once_per_group_transcript(variant_df, .x)
)

readr::write_csv(
  resampled_df,
  "matched_variant_level_one_variant_per_transcript_per_group_100.csv"
)

# ── 9. Long format ────────────────────────────────────────────────────────────
resampled_long <- resampled_df %>%
  select(iter, source, transcript, all_of(flag_cols)) %>%
  pivot_longer(
    cols = all_of(flag_cols),
    names_to = "feature",
    values_to = "value"
  )

# ── 10. Transcript-level median flag across 100 resamples ────────────────────
# Keep actual median value; do NOT binarize 0.5
median_flag_df <- resampled_long %>%
  group_by(source, transcript, feature) %>%
  summarise(
    median_flag = median(value, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  median_flag_df,
  "matched_variant_level_transcript_median_flags_100_resamples.csv"
)

# ── 11. Labels / colors ───────────────────────────────────────────────────────
group_labels <- c(
  "fs" = "Frameshift",
  "fs_control" = "FS Control",
  "snv" = "Nonsense",
  "snv_control" = "Nonsense Control"
)

group_levels <- c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")

group_colors <- c(
  "Frameshift"       = "#D07A3A",
  "FS Control"       = "#E2C7B5",
  "Nonsense"         = "#8ABA67",
  "Nonsense Control" = "#C9D6BE"
)

median_flag_df <- median_flag_df %>%
  mutate(
    group = recode(as.character(source), !!!group_labels),
    group = factor(group, levels = group_levels)
  )

# ── 12. Summary helper ────────────────────────────────────────────────────────
summarise_feature <- function(data, feat) {
  data %>%
    filter(feature == feat) %>%
    filter(!is.na(median_flag)) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_val = mean(median_flag, na.rm = TRUE),
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
      percent = replace_na(percent, 0)
    )
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

# ── 13. Matched/paired p-value helper ─────────────────────────────────────────
# Pair by transcript, then use paired Wilcoxon on transcript-level median_flag
get_paired_p <- function(data, feat, case_source, control_source) {
  dsub <- data %>%
    filter(feature == feat, source %in% c(case_source, control_source)) %>%
    filter(!is.na(median_flag)) %>%
    select(source, transcript, median_flag)
  
  wide <- dsub %>%
    pivot_wider(
      names_from = source,
      values_from = median_flag
    ) %>%
    drop_na()
  
  if (nrow(wide) < 2) return(NA_real_)
  
  x <- wide[[case_source]]
  y <- wide[[control_source]]
  
  tryCatch(
    wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value,
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

# ── 14. Panel plot function ───────────────────────────────────────────────────
plot_enrichment_bar <- function(data, feat, title, panel_letter) {
  df_sum <- summarise_feature(data, feat)
  
  p_fs <- get_paired_p(data, feat, "fs", "fs_control")
  p_snv <- get_paired_p(data, feat, "snv", "snv_control")
  
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
      y = "Mean proportion (%)"
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

# ── 15. Generate all panels ───────────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[ feature_map$panel[i] ]] <- plot_enrichment_bar(
    data = median_flag_df,
    feat = feature_map$feature[i],
    title = feature_map$title[i],
    panel_letter = feature_map$panel[i]
  )
}

# ── 16. Combine figure ────────────────────────────────────────────────────────
ordered_panels <- c("A", "B", "C", "D", "E", "F")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Variant-level Enrichment of Protein Features After 100 Transcript-Level Resamplings",
    subtitle = paste(
      "100 transcript-level matched resamplings;",
      "1 random variant selected per transcript per group per resample;",
      "actual median flag values used;",
      "paired Wilcoxon test for matched transcript comparisons"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0),
      plot.subtitle = element_text(size = 10.5, hjust = 0)
    )
  )

# ── 17. Summary table with paired p-values ────────────────────────────────────
summary_table <- purrr::map_dfr(feature_map$feature, function(feat_now) {
  title_now <- feature_map$title[feature_map$feature == feat_now][1]
  
  sum_df <- summarise_feature(median_flag_df, feat_now)
  
  p_fs <- get_paired_p(median_flag_df, feat_now, "fs", "fs_control")
  p_snv <- get_paired_p(median_flag_df, feat_now, "snv", "snv_control")
  
  sum_df %>%
    mutate(
      feature = feat_now,
      title = title_now,
      paired_p_fs_vs_control = p_fs,
      paired_p_snv_vs_control = p_snv,
      .before = 1
    )
})

readr::write_csv(
  summary_table,
  "matched_variant_level_enrichment_summary_100_resamples_paired_p.csv"
)

# ── 18. Save ──────────────────────────────────────────────────────────────────
ggsave(
  "matched_variant_level_enrichment_100_resamples_paired_p.pdf",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

ggsave(
  "matched_variant_level_enrichment_100_resamples_paired_p.png",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

cat("\nDone!\n")
cat("Output files:\n")
cat("  matched_variant_level_one_variant_per_transcript_per_group_100.csv\n")
cat("  matched_variant_level_transcript_median_flags_100_resamples.csv\n")
cat("  matched_variant_level_enrichment_summary_100_resamples_paired_p.csv\n")
cat("  matched_variant_level_enrichment_100_resamples_paired_p.pdf\n")
cat("  matched_variant_level_enrichment_100_resamples_paired_p.png\n")