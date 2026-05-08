# =============================================================================
# Variant-level enrichment figure after 100 transcript-level resamplings
# For each resample:
#   - randomly select 1 variant per transcript within each source
# Across 100 resamples:
#   - compute median flag value for each transcript and feature
#   - convert median flag to binary
# Then:
#   - draw enrichment barplots similar to Supplemental Figure 3
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
  filter(source %in% c("fs", "fs_control", "snv", "snv_control")) %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control"))
  )

# ── 4. Define binary flag columns ─────────────────────────────────────────────
flag_cols <- c(
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

flag_cols <- flag_cols[flag_cols %in% colnames(df)]

if (length(flag_cols) == 0) {
  stop("No binary flag columns found in the input file.")
}

# convert logical/character flags to numeric 0/1
for (cc in flag_cols) {
  if (is.logical(df[[cc]])) df[[cc]] <- as.numeric(df[[cc]])
  if (is.character(df[[cc]])) {
    df[[cc]] <- case_when(
      df[[cc]] %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      df[[cc]] %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(df[[cc]]))
    )
  }
}

# ── 5. Make sure transcript column exists ─────────────────────────────────────
if (!"transcript" %in% colnames(df)) {
  stop("The column 'transcript' is required but was not found.")
}

# ── 6. Group labels ───────────────────────────────────────────────────────────
group_labels <- c(
  "fs" = "Frameshift",
  "fs_control" = "FS Control",
  "snv" = "Nonsense",
  "snv_control" = "Nonsense Control"
)

# ── 7. Resampling function ────────────────────────────────────────────────────
# For each source x transcript, randomly pick one variant
sample_one_variant_per_transcript <- function(data, iter_id) {
  data %>%
    group_by(source, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    mutate(resample_id = iter_id)
}

# ── 8. Run 100 resamplings ────────────────────────────────────────────────────
n_resamples <- 100

resampled_list <- lapply(seq_len(n_resamples), function(i) {
  sample_one_variant_per_transcript(df, i)
})

resampled_df <- bind_rows(resampled_list)

# ── 9. Compute median flag across resamplings for each transcript ────────────
# This gives one median flag per source/transcript/feature
resampled_long <- resampled_df %>%
  select(resample_id, source, transcript, all_of(flag_cols)) %>%
  pivot_longer(
    cols = all_of(flag_cols),
    names_to = "feature",
    values_to = "flag_value"
  )

median_flag_df <- resampled_long %>%
  group_by(source, transcript, feature) %>%
  summarise(
    median_flag = median(flag_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # if median is 0.5, count as positive
    median_flag_binary = ifelse(is.na(median_flag), NA_real_,
                                ifelse(median_flag >= 0.5, 1, 0)),
    group_label = recode(as.character(source), !!!group_labels),
    group_label = factor(
      group_label,
      levels = c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")
    )
  )

# optional save
write_csv(median_flag_df, "median_flag_by_transcript_after_100_resamples.csv")

# ── 10. Colors ────────────────────────────────────────────────────────────────
group_colors <- c(
  "Frameshift"       = "#D07A3A",
  "FS Control"       = "#E2C7B5",
  "Nonsense"         = "#8ABA67",
  "Nonsense Control" = "#C9D6BE"
)

# ── 11. Features to draw ──────────────────────────────────────────────────────
feature_map <- tribble(
  ~feature,                  ~title,                                      ~ylab,                              ~panel,
  "variant_ppi_overlap",     "PPI Residue Enrichment",                    "Percentage with PPI (%)",          "A",
  "ptc_before_max_pfam_end", "Pfam Domain Enrichment",                    "Percentage with Domains (%)",      "B",
  "variant_nls_flag",        "Nuclear Localization Signal Enrichment",    "Percentage with NLS (%)",          "C",
  "variant_ptm_flag",        "Post-Translational Modification Enrichment","Percentage with PTM (%)",          "D",
  "variant_slim_flag",       "Short Linear Motif Enrichment",             "Percentage with SLiM (%)",         "E",
  "variant_LCS_flag",        "Low Complexity Sequence Enrichment",        "Percentage with LCS (%)",          "F"
) %>%
  filter(feature %in% unique(median_flag_df$feature))

if (nrow(feature_map) == 0) {
  stop("None of the requested plotting features were found.")
}

# ── 12. Helper functions ──────────────────────────────────────────────────────
format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.2f", p)
}

summarise_feature <- function(data, feature_name) {
  data %>%
    filter(feature == feature_name) %>%
    filter(!is.na(median_flag_binary)) %>%
    group_by(group_label) %>%
    summarise(
      total_n = n(),
      positive_n = sum(median_flag_binary == 1, na.rm = TRUE),
      percent = 100 * positive_n / total_n,
      .groups = "drop"
    )
}

get_pairwise_p <- function(data, feature_name, g1, g2) {
  dsub <- data %>%
    filter(feature == feature_name) %>%
    filter(group_label %in% c(g1, g2)) %>%
    filter(!is.na(median_flag_binary))
  
  tab <- dsub %>%
    group_by(group_label) %>%
    summarise(
      pos = sum(median_flag_binary == 1, na.rm = TRUE),
      neg = sum(median_flag_binary == 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(tab) < 2) return(NA_real_)
  
  mat <- matrix(c(tab$pos, tab$neg), nrow = 2, byrow = FALSE)
  
  tryCatch(
    fisher.test(mat)$p.value,
    error = function(e) NA_real_
  )
}

make_bracket_df <- function(summary_df, p1, p2) {
  ymax <- max(summary_df$percent, na.rm = TRUE)
  yrange <- max(5, ymax * 0.12)
  
  tibble(
    xstart = c(1, 3),
    xend   = c(2, 4),
    xmid   = c(1.5, 3.5),
    y      = c(ymax + yrange * 0.55, ymax + yrange * 0.55),
    label  = c(format_p(p1), format_p(p2))
  )
}

plot_enrichment_bar <- function(data, feature_name, title, ylab, panel_letter) {
  
  summary_df <- summarise_feature(data, feature_name)
  
  summary_df <- tibble(
    group_label = factor(
      c("Frameshift", "FS Control", "Nonsense", "Nonsense Control"),
      levels = c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")
    )
  ) %>%
    left_join(summary_df, by = "group_label") %>%
    mutate(
      total_n = replace_na(total_n, 0),
      positive_n = replace_na(positive_n, 0),
      percent = replace_na(percent, 0)
    )
  
  p_fs  <- get_pairwise_p(data, feature_name, "Frameshift", "FS Control")
  p_snv <- get_pairwise_p(data, feature_name, "Nonsense", "Nonsense Control")
  
  bracket_df <- make_bracket_df(summary_df, p_fs, p_snv)
  
  ymax <- max(summary_df$percent, na.rm = TRUE)
  ylim_top <- max(100, ymax + max(10, ymax * 0.25))
  
  ggplot(summary_df, aes(x = group_label, y = percent, fill = group_label)) +
    geom_col(width = 0.62, color = "grey45", alpha = 0.95) +
    geom_text(
      aes(label = sprintf("%.1f%%", percent)),
      vjust = -0.55,
      size = 3.5,
      fontface = "bold",
      color = "grey20"
    ) +
    geom_label(
      aes(
        y = pmax(percent * 0.48, 6),
        label = paste0("n=", positive_n)
      ),
      size = 3.0,
      label.size = 0,
      fill = "white",
      alpha = 0.75
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xstart, y = y - 1.2, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xend, xend = xend, y = y - 1.2, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    geom_text(
      data = bracket_df,
      aes(x = xmid, y = y + 2, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "grey20"
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
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 12, 12, 12)
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = panel_letter,
      hjust = -0.8, vjust = 1.5,
      size = 20 / .pt,
      fontface = "bold"
    )
}

# ── 13. Generate plots ────────────────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[ feature_map$panel[i] ]] <- plot_enrichment_bar(
    data = median_flag_df,
    feature_name = feature_map$feature[i],
    title = feature_map$title[i],
    ylab = feature_map$ylab[i],
    panel_letter = feature_map$panel[i]
  )
}

# ── 14. Arrange combined figure ───────────────────────────────────────────────
ordered_panels <- c("A", "B", "C", "D", "E", "F")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Variant-Level Enrichment of Protein Features After 100 Transcript-Level Resamplings",
    subtitle = "100 transcript-level resamplings; 1 random variant selected per transcript per resample; median flag used for plotting",
    theme = theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0)
    )
  )

# ── 15. Save outputs ──────────────────────────────────────────────────────────
ggsave(
  "Supplemental_Figure_3_variant_level_100_resamples_median_flag.pdf",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

ggsave(
  "Supplemental_Figure_3_variant_level_100_resamples_median_flag.png",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

# ── 16. Save summary table ────────────────────────────────────────────────────
summary_table <- map_dfr(feature_map$feature, function(feat) {
  title_now <- feature_map$title[feature_map$feature == feat][1]
  
  summarise_feature(median_flag_df, feat) %>%
    mutate(feature = feat, title = title_now, .before = 1)
})

write_csv(
  summary_table,
  "Supplemental_Figure_3_variant_level_100_resamples_summary.csv"
)

cat("\nDone!\n")
cat("Output files:\n")
cat("  median_flag_by_transcript_after_100_resamples.csv\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_median_flag.pdf\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_median_flag.png\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_summary.csv\n")