# =============================================================================
# Plot features from variants_all_with_motif_LCS_flags CSV
# Groups: fs | fs_control | snv | snv_control
# Variant-level features
# =============================================================================

rm(list = ls())
gc()

# ── 1. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("ggplot2", "dplyr", "readr", "ggpubr", "patchwork", "scales", "tidyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(tidyr)
})

# ── 2. Read data ──────────────────────────────────────────────────────────────
df <- read_csv("Downloads/variants_all_with_motif_LCS_flags (1).csv", show_col_types = FALSE)

# ── 3. Basic cleaning ─────────────────────────────────────────────────────────
df <- df %>%
  filter(source %in% c("fs", "fs_control", "snv", "snv_control")) %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control")),
    
    # convert logical/binary columns to numeric for plotting if needed
    variant_ppi_overlap    = as.numeric(variant_ppi_overlap),
    ptc_after_max_pfam_end = as.numeric(ptc_after_max_pfam_end),
    ptc_before_max_pfam_end= as.numeric(ptc_before_max_pfam_end),
    variant_protein_flag   = as.numeric(variant_protein_flag),
    variant_domain_flag    = as.numeric(variant_domain_flag),
    variant_slim_flag      = as.numeric(variant_slim_flag),
    variant_morf_flag      = as.numeric(variant_morf_flag),
    variant_ptm_flag       = as.numeric(variant_ptm_flag),
    variant_nls_flag       = as.numeric(variant_nls_flag),
    variant_LCS_flag       = as.numeric(variant_LCS_flag)
  )

# ── 4. Colors / theme / comparisons ──────────────────────────────────────────
group_colors <- c(
  "fs" = "#2166AC",
  "fs_control" = "#92C5DE",
  "snv" = "#B2182B",
  "snv_control" = "#F4A582"
)

comparisons_list <- list(
  c("snv", "snv_control"),
  c("fs", "fs_control")
)

base_theme <- theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ── 5. Feature list ───────────────────────────────────────────────────────────
continuous_features <- tibble::tribble(
  ~feature,                                 ~ylab,                                         ~use_log10,
  "dist_to_cds_end",                        "Distance from PTC to CDS end (bp)",           TRUE,
  "cds_end",                                "CDS length / CDS end",                        TRUE,
  "cds_ptc_loc",                            "PTC location in CDS (bp)",                    TRUE,
  "variant_ppi_nearest_interface_bp",       "Nearest PPI interface position (bp)",         TRUE,
  "variant_ppi_dist_to_nearest_interface_bp","Distance to nearest PPI interface (bp)",     TRUE,
  "ptc_aa",                                 "PTC amino-acid position",                     TRUE,
  "max_pfam_end",                           "Max Pfam end",                                TRUE,
  "dist_ptc_to_max_pfam_end_aa",            "Distance from PTC to max Pfam end (aa)",      FALSE
)

binary_features <- c(
  "variant_ppi_overlap",
  "ptc_after_max_pfam_end",
  "ptc_before_max_pfam_end",
  "variant_protein_flag",
  "variant_domain_flag",
  "variant_slim_flag",
  "variant_morf_flag",
  "variant_ptm_flag",
  "variant_nls_flag",
  "variant_LCS_flag"
)

# ── 6. Helper: continuous plot ────────────────────────────────────────────────
make_continuous_plot <- function(data, feat, ylab = NULL, use_log10 = FALSE) {
  
  plot_data <- data %>%
    filter(!is.na(.data[[feat]]))
  
  p <- ggplot(plot_data, aes(x = source, y = .data[[feat]], fill = source)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(
      width = 0.14,
      outlier.size = 0.5,
      color = "black",
      fill = "white"
    ) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.format"
    ) +
    scale_fill_manual(values = group_colors) +
    labs(
      title = feat,
      x = NULL,
      y = ifelse(is.null(ylab), feat, ylab)
    ) +
    base_theme
  
  if (use_log10) {
    vals <- plot_data[[feat]]
    vals <- vals[is.finite(vals) & !is.na(vals) & vals > 0]
    
    if (length(vals) > 0) {
      p <- p + scale_y_log10()
    }
  }
  
  return(p)
}

# ── 7. Helper: binary/discrete plot ──────────────────────────────────────────
make_binary_plot <- function(data, feat) {
  
  plot_data <- data %>%
    filter(!is.na(.data[[feat]]))
  
  ggplot(plot_data, aes(x = source, y = .data[[feat]], fill = source)) +
    geom_boxplot(width = 0.45, outlier.size = 0.5, color = "black") +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.format"
    ) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(breaks = c(0, 1)) +
    labs(
      title = feat,
      x = NULL,
      y = feat
    ) +
    base_theme
}

# ── 8. Make continuous plots ─────────────────────────────────────────────────
continuous_plot_list <- list()

for (i in seq_len(nrow(continuous_features))) {
  feat <- continuous_features$feature[i]
  ylab <- continuous_features$ylab[i]
  use_log10 <- continuous_features$use_log10[i]
  
  if (feat %in% colnames(df)) {
    continuous_plot_list[[feat]] <- make_continuous_plot(
      data = df,
      feat = feat,
      ylab = ylab,
      use_log10 = use_log10
    )
  }
}

# ── 9. Make binary plots ──────────────────────────────────────────────────────
binary_plot_list <- list()

for (feat in binary_features) {
  if (feat %in% colnames(df)) {
    binary_plot_list[[feat]] <- make_binary_plot(df, feat)
  }
}

# ── 10. Save individual plots ────────────────────────────────────────────────
dir.create("variant_feature_plots", showWarnings = FALSE)

for (nm in names(continuous_plot_list)) {
  ggsave(
    filename = file.path("variant_feature_plots", paste0(nm, ".pdf")),
    plot = continuous_plot_list[[nm]],
    width = 6,
    height = 5
  )
  
  ggsave(
    filename = file.path("variant_feature_plots", paste0(nm, ".png")),
    plot = continuous_plot_list[[nm]],
    width = 6,
    height = 5,
    dpi = 300
  )
}

for (nm in names(binary_plot_list)) {
  ggsave(
    filename = file.path("variant_feature_plots", paste0(nm, ".pdf")),
    plot = binary_plot_list[[nm]],
    width = 6,
    height = 5
  )
  
  ggsave(
    filename = file.path("variant_feature_plots", paste0(nm, ".png")),
    plot = binary_plot_list[[nm]],
    width = 6,
    height = 5,
    dpi = 300
  )
}

# ── 11. Combined continuous figure ───────────────────────────────────────────
if (length(continuous_plot_list) > 0) {
  combined_continuous <- wrap_plots(continuous_plot_list, ncol = 2) +
    plot_annotation(
      title = "Variant-level continuous features across groups",
      subtitle = "Wilcoxon test: snv vs snv_control; fs vs fs_control"
    )
  
  ggsave(
    "variants_all_continuous_features.pdf",
    combined_continuous,
    width = 14,
    height = 4 * ceiling(length(continuous_plot_list) / 2)
  )
  
  ggsave(
    "variants_all_continuous_features.png",
    combined_continuous,
    width = 14,
    height = 4 * ceiling(length(continuous_plot_list) / 2),
    dpi = 300
  )
}

# ── 12. Combined binary figure ───────────────────────────────────────────────
if (length(binary_plot_list) > 0) {
  combined_binary <- wrap_plots(binary_plot_list, ncol = 2) +
    plot_annotation(
      title = "Variant-level motif / annotation flags across groups",
      subtitle = "Wilcoxon test: snv vs snv_control; fs vs fs_control"
    )
  
  ggsave(
    "variants_all_binary_features.pdf",
    combined_binary,
    width = 14,
    height = 4 * ceiling(length(binary_plot_list) / 2)
  )
  
  ggsave(
    "variants_all_binary_features.png",
    combined_binary,
    width = 14,
    height = 4 * ceiling(length(binary_plot_list) / 2),
    dpi = 300
  )
}

# ── 13. Optional: proportion barplots for binary flags ───────────────────────
binary_long <- df %>%
  select(source, all_of(binary_features)) %>%
  pivot_longer(
    cols = -source,
    names_to = "feature",
    values_to = "value"
  ) %>%
  group_by(source, feature) %>%
  summarise(
    proportion = mean(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    .groups = "drop"
  )

p_bar <- ggplot(binary_long, aes(x = source, y = proportion, fill = source)) +
  geom_col() +
  facet_wrap(~feature, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of variants with each binary feature",
    x = NULL,
    y = "Proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave("variants_all_binary_proportions.pdf", p_bar, width = 14, height = 10)
ggsave("variants_all_binary_proportions.png", p_bar, width = 14, height = 10, dpi = 300)

# ── 14. Summary stats ────────────────────────────────────────────────────────
summary_continuous <- df %>%
  group_by(source) %>%
  summarise(
    n = n(),
    dist_to_cds_end_median = median(dist_to_cds_end, na.rm = TRUE),
    cds_end_median = median(cds_end, na.rm = TRUE),
    cds_ptc_loc_median = median(cds_ptc_loc, na.rm = TRUE),
    variant_ppi_nearest_interface_bp_median =
      median(variant_ppi_nearest_interface_bp, na.rm = TRUE),
    variant_ppi_dist_to_nearest_interface_bp_median =
      median(variant_ppi_dist_to_nearest_interface_bp, na.rm = TRUE),
    ptc_aa_median = median(ptc_aa, na.rm = TRUE),
    max_pfam_end_median = median(max_pfam_end, na.rm = TRUE),
    dist_ptc_to_max_pfam_end_aa_median =
      median(dist_ptc_to_max_pfam_end_aa, na.rm = TRUE),
    .groups = "drop"
  )

summary_binary <- df %>%
  group_by(source) %>%
  summarise(
    n = n(),
    variant_ppi_overlap_prop = mean(variant_ppi_overlap, na.rm = TRUE),
    ptc_after_max_pfam_end_prop = mean(ptc_after_max_pfam_end, na.rm = TRUE),
    ptc_before_max_pfam_end_prop = mean(ptc_before_max_pfam_end, na.rm = TRUE),
    variant_protein_flag_prop = mean(variant_protein_flag, na.rm = TRUE),
    variant_domain_flag_prop = mean(variant_domain_flag, na.rm = TRUE),
    variant_slim_flag_prop = mean(variant_slim_flag, na.rm = TRUE),
    variant_morf_flag_prop = mean(variant_morf_flag, na.rm = TRUE),
    variant_ptm_flag_prop = mean(variant_ptm_flag, na.rm = TRUE),
    variant_nls_flag_prop = mean(variant_nls_flag, na.rm = TRUE),
    variant_LCS_flag_prop = mean(variant_LCS_flag, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_continuous, "variants_all_continuous_summary.csv")
write_csv(summary_binary, "variants_all_binary_summary.csv")

# ── 15. Print basic info ─────────────────────────────────────────────────────
cat("\nColumns in data:\n")
print(colnames(df))

cat("\nCounts per group:\n")
print(table(df$source))

cat("\nContinuous summary:\n")
print(summary_continuous)

cat("\nBinary summary:\n")
print(summary_binary)

cat("\nDone!\n")
cat("Output files:\n")
cat("  variant_feature_plots/                      <- individual plots\n")
cat("  variants_all_continuous_features.pdf/png\n")
cat("  variants_all_binary_features.pdf/png\n")
cat("  variants_all_binary_proportions.pdf/png\n")
cat("  variants_all_continuous_summary.csv\n")
cat("  variants_all_binary_summary.csv\n")