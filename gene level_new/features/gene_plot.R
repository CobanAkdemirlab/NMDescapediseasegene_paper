#This R script ie to plot gene level and variant level features

gene_all0407_1_ <- read_csv("Downloads/gene_all0407 (1).csv")
variants_all_with_motif_LCS_flags_1_ <- read_csv("Downloads/variants_all_with_motif_LCS_flags (1).csv")

# =============================================================================
# Plot features from gene_all0407 CSV
# Groups: fs | fs_control | snv | snv_control
# Output: individual plots + combined plots
# =============================================================================

rm(list = ls())
gc()

# ── 1. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("ggplot2", "dplyr", "readr", "ggpubr", "patchwork", "scales")
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
})

# ── 2. Read data ──────────────────────────────────────────────────────────────
df <- read_csv("Downloads/gene_all0407 (1).csv", show_col_types = FALSE)

# If needed, you can use:
# df <- read_csv("gene_all0407.csv", show_col_types = FALSE)

# ── 3. Basic cleaning ─────────────────────────────────────────────────────────
# Keep only the 4 groups of interest
df <- df %>%
  filter(group %in% c("fs", "fs_control", "snv", "snv_control")) %>%
  mutate(
    group = factor(group, levels = c("fs", "fs_control", "snv", "snv_control"))
  )

# Make sure logical columns are numeric when needed
df <- df %>%
  mutate(
    ppi_overlap = as.numeric(ppi_overlap),
    pfam_overlap_flag = as.numeric(pfam_overlap_flag),
    gene_protein_flag = as.numeric(gene_protein_flag),
    gene_domains_flag = as.numeric(gene_domains_flag),
    gene_slim_flag = as.numeric(gene_slim_flag),
    gene_morf_flag = as.numeric(gene_morf_flag),
    gene_ptm_flag = as.numeric(gene_ptm_flag),
    gene_nls_flag = as.numeric(gene_nls_flag)
  )

# ── 4. Colors and theme ───────────────────────────────────────────────────────
group_colors <- c(
  "fs" = "#2166AC",
  "fs_control" = "#92C5DE",
  "snv" = "#B2182B",
  "snv_control" = "#F4A582"
)

base_theme <- theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

comparisons_list <- list(
  c("snv", "snv_control"),
  c("fs", "fs_control")
)

# ── 5. Features to plot ───────────────────────────────────────────────────────
# Continuous features
feature_info <- tibble::tribble(
  ~feature,                        ~ylab,                              ~use_log10,
  "gc_content",                    "GC content",                       FALSE,
  "nmdesc_gc_content",             "NMDesc GC content",                FALSE,
  "repeat_fraction",               "Repeat fraction",                  FALSE,
  "nmdesc_repeat_fraction",        "NMDesc repeat fraction",           FALSE,
  "homopolymer_fraction",          "Homopolymer fraction",             FALSE,
  "nmdesc_homopolymer_fraction",   "NMDesc homopolymer fraction",      FALSE,
  "NMDesc_region_length",          "NMDesc region length (bp)",        TRUE,
  "cds_length",                    "CDS length (bp)",                  TRUE,
  "ppi_overlap",                   "PPI overlap",                      FALSE,
  "pfam_overlap_length",           "Pfam overlap length",              TRUE,
  "pfam_overlap_fraction",         "Pfam overlap fraction",            FALSE,
  "n_overlapping_pfam",            "Number of overlapping Pfam",       FALSE
)

# Optional binary gene-level flags
binary_features <- c(
  "gene_protein_flag",
  "gene_domains_flag",
  "gene_slim_flag",
  "gene_morf_flag",
  "gene_ptm_flag",
  "gene_nls_flag",
  "pfam_overlap_flag"
)

# ── 6. Helper function: violin + boxplot ─────────────────────────────────────
make_feature_plot <- function(data, feat, ylab = NULL, use_log10 = FALSE) {
  
  plot_data <- data %>%
    filter(!is.na(.data[[feat]]))
  
  p <- ggplot(plot_data, aes(x = group, y = .data[[feat]], fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA) +
    geom_boxplot(
      width = 0.14,
      outlier.size = 0.6,
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
    positive_values <- plot_data[[feat]][plot_data[[feat]] > 0]
    if (length(positive_values) > 0) {
      p <- p + scale_y_log10()
    }
  }
  
  return(p)
}

# ── 7. Helper function: boxplot for binary/discrete flags ────────────────────
make_binary_plot <- function(data, feat) {
  
  plot_data <- data %>%
    filter(!is.na(.data[[feat]]))
  
  ggplot(plot_data, aes(x = group, y = .data[[feat]], fill = group)) +
    geom_boxplot(width = 0.5, outlier.size = 0.6, color = "black") +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.format"
    ) +
    scale_fill_manual(values = group_colors) +
    labs(
      title = feat,
      x = NULL,
      y = feat
    ) +
    base_theme
}

# ── 8. Generate plots for continuous features ────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_info))) {
  feat <- feature_info$feature[i]
  ylab <- feature_info$ylab[i]
  use_log10 <- feature_info$use_log10[i]
  
  if (feat %in% colnames(df)) {
    plot_list[[feat]] <- make_feature_plot(
      data = df,
      feat = feat,
      ylab = ylab,
      use_log10 = use_log10
    )
  }
}

# ── 9. Generate plots for binary features ────────────────────────────────────
binary_plot_list <- list()

for (feat in binary_features) {
  if (feat %in% colnames(df)) {
    binary_plot_list[[feat]] <- make_binary_plot(df, feat)
  }
}

# ── 10. Save individual continuous plots ─────────────────────────────────────
dir.create("feature_plots", showWarnings = FALSE)

for (nm in names(plot_list)) {
  ggsave(
    filename = file.path("feature_plots", paste0(nm, ".pdf")),
    plot = plot_list[[nm]],
    width = 6,
    height = 5
  )
  
  ggsave(
    filename = file.path("feature_plots", paste0(nm, ".png")),
    plot = plot_list[[nm]],
    width = 6,
    height = 5,
    dpi = 300
  )
}

# ── 11. Save individual binary plots ─────────────────────────────────────────
for (nm in names(binary_plot_list)) {
  ggsave(
    filename = file.path("feature_plots", paste0(nm, ".pdf")),
    plot = binary_plot_list[[nm]],
    width = 6,
    height = 5
  )
  
  ggsave(
    filename = file.path("feature_plots", paste0(nm, ".png")),
    plot = binary_plot_list[[nm]],
    width = 6,
    height = 5,
    dpi = 300
  )
}

# ── 12. Combined plots: continuous features ──────────────────────────────────
if (length(plot_list) > 0) {
  combined_continuous <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "Feature comparison across variant groups",
      subtitle = "Wilcoxon test: snv vs snv_control; fs vs fs_control"
    )
  
  ggsave(
    "all_continuous_features.pdf",
    combined_continuous,
    width = 14,
    height = 4 * ceiling(length(plot_list) / 2)
  )
  
  ggsave(
    "all_continuous_features.png",
    combined_continuous,
    width = 14,
    height = 4 * ceiling(length(plot_list) / 2),
    dpi = 300
  )
}

# ── 13. Combined plots: binary features ──────────────────────────────────────
if (length(binary_plot_list) > 0) {
  combined_binary <- wrap_plots(binary_plot_list, ncol = 2) +
    plot_annotation(
      title = "Binary / flag features across variant groups",
      subtitle = "Wilcoxon test: snv vs snv_control; fs vs fs_control"
    )
  
  ggsave(
    "all_binary_features.pdf",
    combined_binary,
    width = 14,
    height = 4 * ceiling(length(binary_plot_list) / 2)
  )
  
  ggsave(
    "all_binary_features.png",
    combined_binary,
    width = 14,
    height = 4 * ceiling(length(binary_plot_list) / 2),
    dpi = 300
  )
}

# ── 14. Summary statistics ────────────────────────────────────────────────────
summary_stats <- df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    gc_content_median = median(gc_content, na.rm = TRUE),
    nmdesc_gc_content_median = median(nmdesc_gc_content, na.rm = TRUE),
    repeat_fraction_median = median(repeat_fraction, na.rm = TRUE),
    nmdesc_repeat_fraction_median = median(nmdesc_repeat_fraction, na.rm = TRUE),
    homopolymer_fraction_median = median(homopolymer_fraction, na.rm = TRUE),
    nmdesc_homopolymer_fraction_median = median(nmdesc_homopolymer_fraction, na.rm = TRUE),
    NMDesc_region_length_median = median(NMDesc_region_length, na.rm = TRUE),
    cds_length_median = median(cds_length, na.rm = TRUE),
    ppi_overlap_median = median(ppi_overlap, na.rm = TRUE),
    pfam_overlap_length_median = median(pfam_overlap_length, na.rm = TRUE),
    pfam_overlap_fraction_median = median(pfam_overlap_fraction, na.rm = TRUE),
    n_overlapping_pfam_median = median(n_overlapping_pfam, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats, "feature_summary_statistics.csv")

# ── 15. Print summary ─────────────────────────────────────────────────────────
print(colnames(df))
print(table(df$group))
print(summary_stats)

cat("\nDone!\n")
cat("Output files:\n")
cat("  feature_plots/                <- individual plots\n")
cat("  all_continuous_features.pdf/png\n")
cat("  all_binary_features.pdf/png\n")
cat("  feature_summary_statistics.csv\n")