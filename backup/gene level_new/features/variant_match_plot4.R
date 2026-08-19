# =============================================================================
# Matched resampling analysis for continuous variant-level features
# Variables:
#   1) dist_to_cds_end
#   2) variant_ppi_dist_to_nearest_interface_bp
#   3) dist_ptc_to_max_pfam_end_aa
#
# Matched unit: transcript
# Resampling: 100 times
# Per iteration:
#   - for each shared transcript, randomly select 1 variant in case
#   - for each shared transcript, randomly select 1 variant in control
#   - compute matched difference
# =============================================================================

rm(list = ls())
gc()

# ── 1. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "purrr", "ggplot2", "readr", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

set.seed(123)

# ── 2. Read data ──────────────────────────────────────────────────────────────
variant_df <- read_csv("variants_all_with_motif_LCS_flags (1).csv", show_col_types = FALSE)

# ── 3. Basic cleaning ─────────────────────────────────────────────────────────
variant_df <- variant_df %>%
  filter(source %in% c("fs", "fs_control", "snv", "snv_control")) %>%
  filter(!is.na(transcript)) %>%
  mutate(
    source = factor(source, levels = c("fs_control", "fs", "snv_control", "snv"))
  )

cat("Rows after filtering:", nrow(variant_df), "\n")
cat("Unique transcripts:", dplyr::n_distinct(variant_df$transcript), "\n")
print(table(variant_df$source))

# ── 4. Variables to analyze ───────────────────────────────────────────────────
vars_to_run <- c(
  "dist_to_cds_end",
  "variant_ppi_dist_to_nearest_interface_bp",
  "dist_ptc_to_max_pfam_end_aa"
)

vars_to_run <- vars_to_run[vars_to_run %in% colnames(variant_df)]

if (length(vars_to_run) == 0) {
  stop("None of the requested variables were found in the dataset.")
}

pretty_names <- c(
  "dist_to_cds_end" = "Distance to CDS end (bp)",
  "variant_ppi_dist_to_nearest_interface_bp" = "Distance to nearest PPI interface (bp)",
  "dist_ptc_to_max_pfam_end_aa" = "Distance from PTC to max Pfam end (aa)"
)

# ── 5. Shared transcripts ─────────────────────────────────────────────────────
get_shared_transcripts <- function(data, case_group, control_group, value_col) {
  dsub <- data %>%
    filter(!is.na(.data[[value_col]]))
  
  intersect(
    unique(dsub$transcript[dsub$source == case_group]),
    unique(dsub$transcript[dsub$source == control_group])
  )
}

# ── 6. One matched resample ───────────────────────────────────────────────────
sample_matched_once <- function(data, case_group, control_group, transcript_set, value_col, iter_id) {
  case_df <- data %>%
    filter(source == case_group, transcript %in% transcript_set, !is.na(.data[[value_col]])) %>%
    group_by(transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(transcript, case_value = all_of(value_col))
  
  control_df <- data %>%
    filter(source == control_group, transcript %in% transcript_set, !is.na(.data[[value_col]])) %>%
    group_by(transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(transcript, control_value = all_of(value_col))
  
  case_df %>%
    inner_join(control_df, by = "transcript") %>%
    mutate(
      iter = iter_id,
      diff = case_value - control_value
    )
}

# ── 7. Empirical p value ──────────────────────────────────────────────────────
compute_empirical_p <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  
  p_left  <- mean(x <= 0)
  p_right <- mean(x >= 0)
  
  min(2 * min(p_left, p_right), 1)
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

# ── 8. Main function for one variable ─────────────────────────────────────────
run_matched_resampling <- function(data, value_col, n_iter = 100) {
  dsub <- data %>%
    filter(!is.na(.data[[value_col]]))
  
  # positive-only for log-scale panel
  dsub_pos <- dsub %>%
    filter(.data[[value_col]] > 0)
  
  fs_tr <- get_shared_transcripts(dsub_pos, "fs", "fs_control", value_col)
  snv_tr <- get_shared_transcripts(dsub_pos, "snv", "snv_control", value_col)
  
  cat("\nVariable:", value_col, "\n")
  cat("Shared transcripts: fs vs fs_control =", length(fs_tr), "\n")
  cat("Shared transcripts: snv vs snv_control =", length(snv_tr), "\n")
  
  fs_resampled <- map_dfr(
    seq_len(n_iter),
    ~sample_matched_once(dsub_pos, "fs", "fs_control", fs_tr, value_col, .x)
  )
  
  snv_resampled <- map_dfr(
    seq_len(n_iter),
    ~sample_matched_once(dsub_pos, "snv", "snv_control", snv_tr, value_col, .x)
  )
  
  fs_summary <- fs_resampled %>%
    group_by(iter) %>%
    summarise(
      median_case = median(case_value, na.rm = TRUE),
      median_control = median(control_value, na.rm = TRUE),
      median_diff = median(diff, na.rm = TRUE),
      mean_diff = mean(diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(comparison = "Frameshift - FS Control")
  
  snv_summary <- snv_resampled %>%
    group_by(iter) %>%
    summarise(
      median_case = median(case_value, na.rm = TRUE),
      median_control = median(control_value, na.rm = TRUE),
      median_diff = median(diff, na.rm = TRUE),
      mean_diff = mean(diff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(comparison = "Nonsense - Nonsense Control")
  
  fs_p <- compute_empirical_p(fs_summary$median_diff)
  snv_p <- compute_empirical_p(snv_summary$median_diff)
  
  fs_wilcox <- wilcox.test(
    fs_summary$median_case,
    fs_summary$median_control,
    paired = TRUE,
    exact = FALSE
  )
  
  snv_wilcox <- wilcox.test(
    snv_summary$median_case,
    snv_summary$median_control,
    paired = TRUE,
    exact = FALSE
  )
  
  fs_p_label <- format_p(fs_p)
  snv_p_label <- format_p(snv_p)
  fs_w_label <- format_p(fs_wilcox$p.value)
  snv_w_label <- format_p(snv_wilcox$p.value)
  
  # Plot 1: difference distribution
  plot_diff <- bind_rows(fs_summary, snv_summary) %>%
    mutate(
      comparison = factor(
        comparison,
        levels = c("Frameshift - FS Control", "Nonsense - Nonsense Control")
      )
    )
  
  y1 <- max(plot_diff$median_diff, na.rm = TRUE)
  
  p1 <- ggplot(plot_diff, aes(x = comparison, y = median_diff, fill = comparison)) +
    geom_violin(trim = TRUE, alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.18, fill = "white", color = "black", outlier.size = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
    scale_fill_manual(values = c(
      "Frameshift - FS Control" = "#2166AC",
      "Nonsense - Nonsense Control" = "#B2182B"
    )) +
    labs(
      title = paste0(pretty_names[[value_col]], ": matched difference"),
      subtitle = paste0("Empirical p-values: FS ", fs_p_label, " | SNV ", snv_p_label),
      x = NULL,
      y = "Median difference (case - control)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = 1, y = y1 * 0.95, label = fs_p_label, size = 4) +
    annotate("text", x = 2, y = y1 * 0.95, label = snv_p_label, size = 4)
  
  # Plot 2: case/control medians
  plot_groups <- bind_rows(
    fs_summary %>%
      select(iter, case = median_case, control = median_control) %>%
      pivot_longer(cols = c(case, control), names_to = "group_type", values_to = "median_value") %>%
      mutate(group = ifelse(group_type == "case", "Frameshift", "FS Control")),
    snv_summary %>%
      select(iter, case = median_case, control = median_control) %>%
      pivot_longer(cols = c(case, control), names_to = "group_type", values_to = "median_value") %>%
      mutate(group = ifelse(group_type == "case", "Nonsense", "Nonsense Control"))
  ) %>%
    mutate(
      group = factor(group, levels = c("FS Control", "Frameshift", "Nonsense Control", "Nonsense"))
    )
  
  ymax <- max(plot_groups$median_value, na.rm = TRUE)
  
  p2 <- ggplot(plot_groups, aes(x = group, y = median_value, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.16, fill = "white", color = "black", outlier.size = 0.4) +
    scale_fill_manual(values = c(
      "FS Control" = "#92C5DE",
      "Frameshift" = "#2166AC",
      "Nonsense Control" = "#F4A582",
      "Nonsense" = "#B2182B"
    )) +
    scale_y_log10() +
    labs(
      title = paste0(pretty_names[[value_col]], ": matched resampling"),
      subtitle = paste0("Paired Wilcoxon: FS ", fs_w_label, " | SNV ", snv_w_label),
      x = NULL,
      y = pretty_names[[value_col]]
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    annotate("segment", x = 1, xend = 2, y = ymax * 1.20, yend = ymax * 1.20, linewidth = 0.6) +
    annotate("text", x = 1.5, y = ymax * 1.33, label = fs_w_label, size = 4) +
    annotate("segment", x = 3, xend = 4, y = ymax * 1.20, yend = ymax * 1.20, linewidth = 0.6) +
    annotate("text", x = 3.5, y = ymax * 1.33, label = snv_w_label, size = 4)
  
  summary_table <- bind_rows(
    fs_summary %>%
      summarise(
        comparison = "Frameshift - FS Control",
        variable = value_col,
        n_iter = n(),
        n_shared_transcripts = length(fs_tr),
        median_of_median_diff = median(median_diff, na.rm = TRUE),
        mean_of_median_diff = mean(median_diff, na.rm = TRUE),
        empirical_p = fs_p,
        paired_wilcox_p = fs_wilcox$p.value
      ),
    snv_summary %>%
      summarise(
        comparison = "Nonsense - Nonsense Control",
        variable = value_col,
        n_iter = n(),
        n_shared_transcripts = length(snv_tr),
        median_of_median_diff = median(median_diff, na.rm = TRUE),
        mean_of_median_diff = mean(median_diff, na.rm = TRUE),
        empirical_p = snv_p,
        paired_wilcox_p = snv_wilcox$p.value
      )
  )
  
  list(
    fs_resampled = fs_resampled,
    snv_resampled = snv_resampled,
    fs_summary = fs_summary,
    snv_summary = snv_summary,
    p1 = p1,
    p2 = p2,
    combined_plot = p1 / p2 + plot_annotation(
      title = paste0(pretty_names[[value_col]], " matched analysis (100 resamplings)")
    ),
    summary_table = summary_table
  )
}

# ── 9. Run for all requested variables ────────────────────────────────────────
results_list <- map(vars_to_run, ~run_matched_resampling(variant_df, .x, n_iter = 100))
names(results_list) <- vars_to_run

# ── 10. Save outputs ──────────────────────────────────────────────────────────
all_summary <- bind_rows(map(results_list, "summary_table"))
write_csv(all_summary, "matched_continuous_features_resampling_summary_100.csv")

for (v in names(results_list)) {
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", v)
  
  write_csv(results_list[[v]]$fs_resampled, paste0("matched_fs_resampled_", safe_name, "_100.csv"))
  write_csv(results_list[[v]]$snv_resampled, paste0("matched_snv_resampled_", safe_name, "_100.csv"))
  write_csv(results_list[[v]]$fs_summary, paste0("matched_fs_summary_", safe_name, "_100.csv"))
  write_csv(results_list[[v]]$snv_summary, paste0("matched_snv_summary_", safe_name, "_100.csv"))
  
  ggsave(
    paste0("matched_", safe_name, "_resampling_plot_100.pdf"),
    results_list[[v]]$combined_plot,
    width = 11,
    height = 10
  )
  
  ggsave(
    paste0("matched_", safe_name, "_resampling_plot_100.png"),
    results_list[[v]]$combined_plot,
    width = 11,
    height = 10,
    dpi = 300
  )
}

# ── 11. Optional combined main figure ─────────────────────────────────────────
main_fig <- wrap_plots(
  lapply(results_list, function(x) x$p2),
  ncol = 1
) + plot_annotation(
  title = "Matched resampling analysis for continuous variant-level features"
)

ggsave("matched_continuous_features_main_figure_100.pdf", main_fig, width = 10, height = 15)
ggsave("matched_continuous_features_main_figure_100.png", main_fig, width = 10, height = 15, dpi = 300)

cat("\nDone!\n")
cat("Summary file:\n")
cat("  matched_continuous_features_resampling_summary_100.csv\n")
cat("\nPer-variable files were saved for:\n")
print(names(results_list))