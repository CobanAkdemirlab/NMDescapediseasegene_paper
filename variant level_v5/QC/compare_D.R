library(tidyverse)
library(patchwork)
library(readr)
library(ggplot2)

# =========================================================
# Compare control sets separately:
#   1) fs_control   vs fs_control2
#   2) snv_control  vs snv_control2
#
# Assumes these objects already exist in memory:
#   - variant_plp_control
#   - control_variants_all2
# =========================================================

# ── 1. Define binary features and labels ──────────────────────────────────────
feature_map <- tibble::tribble(
  ~feature,                  ~title,                                   ~panel,
  "variant_ppi_overlap",     "PPI Residue Overlap",                    "A",
  "ptc_before_max_pfam_end", "Pfam Domain Overlap",                    "B",
  "variant_nls_flag",        "Nuclear Localization Signal",            "C",
  "variant_ptm_flag",        "Post-Translational Modification",        "D",
  "variant_slim_flag",       "Short Linear Motif",                     "E",
  "variant_LCS_flag",        "Low Complexity Sequence",                "F"
)

feature_map <- feature_map %>%
  filter(
    feature %in% colnames(variant_plp_control),
    feature %in% colnames(control_variants_all2)
  )

features <- feature_map$feature

if (length(features) == 0) {
  stop("No shared binary feature columns found in the two datasets.")
}

# continuous feature
dist_col_set2 <- "dist_to_cds_end"
dist_col_set3 <- if ("dist_to_cds_end" %in% colnames(control_variants_all2)) {
  "dist_to_cds_end"
} else if ("ptc_distance_to_cds_end" %in% colnames(control_variants_all2)) {
  "ptc_distance_to_cds_end"
} else {
  NA_character_
}

if (is.na(dist_col_set3)) {
  stop("Neither 'dist_to_cds_end' nor 'ptc_distance_to_cds_end' was found in control_variants_all2.")
}

# ── 2. Helper functions ───────────────────────────────────────────────────────
to_binary <- function(x) {
  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_real_, as.numeric(x)))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    return(ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0)))
  }
  
  if (is.character(x)) {
    x2 <- dplyr::case_when(
      x %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      x %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(x))
    )
    return(ifelse(is.na(x2), NA_real_, ifelse(x2 > 0, 1, 0)))
  }
  
  return(rep(NA_real_, length(x)))
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ sprintf("p = %.3f", p)
  )
}

run_test <- function(data, feat, g1, g2) {
  dsub <- data %>%
    filter(feature == feat, group %in% c(g1, g2), !is.na(value)) %>%
    mutate(value = as.numeric(value))
  
  x1 <- dsub %>% filter(group == g1) %>% pull(value)
  x2 <- dsub %>% filter(group == g2) %>% pull(value)
  
  if (length(x1) == 0 || length(x2) == 0) {
    return(tibble(
      p = NA_real_,
      method = NA_character_,
      g1_0 = NA_real_,
      g1_1 = NA_real_,
      g2_0 = NA_real_,
      g2_1 = NA_real_
    ))
  }
  
  tab <- matrix(
    c(sum(x1 == 0), sum(x1 == 1),
      sum(x2 == 0), sum(x2 == 1)),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(tab) <- c(g1, g2)
  colnames(tab) <- c("0", "1")
  
  chisq_obj <- suppressWarnings(chisq.test(tab, correct = FALSE))
  expected_min <- min(chisq_obj$expected)
  
  if (expected_min < 5) {
    p <- fisher.test(tab)$p.value
    method <- "Fisher"
  } else {
    p <- suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
    method <- "Chi-square"
  }
  
  tibble(
    p = p,
    method = method,
    g1_0 = tab[1, "0"],
    g1_1 = tab[1, "1"],
    g2_0 = tab[2, "0"],
    g2_1 = tab[2, "1"]
  )
}

run_continuous_test <- function(data, g1, g2) {
  x1 <- data %>% filter(group == g1) %>% pull(dist_to_cds_end)
  x2 <- data %>% filter(group == g2) %>% pull(dist_to_cds_end)
  
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]
  
  if (length(x1) < 2 || length(x2) < 2) {
    return(tibble(
      p = NA_real_,
      method = NA_character_,
      n_g1 = length(x1),
      n_g2 = length(x2)
    ))
  }
  
  p <- tryCatch(
    wilcox.test(x1, x2, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    p = p,
    method = "Wilcoxon rank-sum",
    n_g1 = length(x1),
    n_g2 = length(x2)
  )
}

plot_feature <- function(feat, title, panel_letter, summary_df, p_label_df, group_levels, group_colors) {
  df_sum <- summary_df %>%
    filter(feature == feat) %>%
    mutate(group = factor(group, levels = group_levels))
  
  bracket_df <- p_label_df %>%
    filter(feature == feat)
  
  ymax <- max(df_sum$percent, na.rm = TRUE)
  offset <- max(6, ymax * 0.12)
  
  bracket_df <- bracket_df %>%
    mutate(y = ymax + offset)
  
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
      y = "Proportion (%)"
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

plot_distance <- function(dist_df, dist_p_df, group_levels, group_colors) {
  bracket_df <- dist_p_df %>%
    mutate(
      xstart = c(1, 3),
      xend   = c(2, 4),
      xmid   = c(1.5, 3.5)
    )
  
  ymax <- max(dist_df$dist_to_cds_end, na.rm = TRUE)
  ymin_pos <- min(dist_df$dist_to_cds_end[dist_df$dist_to_cds_end > 0], na.rm = TRUE)
  
  bracket_df <- bracket_df %>%
    mutate(
      y = c(ymax * 1.25, ymax * 1.25),
      label = format_p(p)
    )
  
  ggplot(
    dist_df %>% mutate(group = factor(group, levels = group_levels)),
    aes(x = group, y = dist_to_cds_end, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.14, outlier.size = 0.5, color = "black", alpha = 0.9) +
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey20"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xstart, y = y / 1.05, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey20"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xend, xend = xend, y = y / 1.05, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey20"
    ) +
    geom_text(
      data = bracket_df,
      aes(x = xmid, y = y * 1.05, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "grey15"
    ) +
    scale_fill_manual(values = group_colors) +
    scale_y_log10() +
    coord_cartesian(ylim = c(ymin_pos, ymax * 1.6)) +
    labs(
      title = "PTC Distance to CDS End",
      x = "Variant Category",
      y = "PTC distance to CDS end (bp)"
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
      label = "G",
      hjust = -0.8, vjust = 1.4,
      size = 6,
      fontface = "bold"
    )
}

# ── 3. Prepare set2 data ──────────────────────────────────────────────────────
df_set2 <- variant_plp_control %>%
  mutate(
    group = case_when(
      source == "fs_control"  ~ "FS Control",
      source == "snv_control" ~ "SNV Control",
      TRUE ~ NA_character_
    )
  ) %>%
  select(group, all_of(features)) %>%
  mutate(across(all_of(features), to_binary)) %>%
  filter(!is.na(group))

# ── 4. Prepare set3 data ──────────────────────────────────────────────────────
df_set3 <- control_variants_all2 %>%
  mutate(
    group = case_when(
      group == "fs_control2"  ~ "FS Control2",
      group == "snv_control2" ~ "SNV Control2",
      TRUE ~ NA_character_
    )
  ) %>%
  select(group, all_of(features)) %>%
  mutate(across(all_of(features), to_binary)) %>%
  filter(!is.na(group))

# ── 5. Combine and reshape binary features ────────────────────────────────────
df_all <- bind_rows(df_set2, df_set3)

group_levels <- c("FS Control", "FS Control2", "SNV Control", "SNV Control2")

df_long <- df_all %>%
  pivot_longer(
    cols = all_of(features),
    names_to = "feature",
    values_to = "value"
  ) %>%
  mutate(group = factor(group, levels = group_levels))

# ── 6. Summary table for binary plotting ──────────────────────────────────────
summary_df <- df_long %>%
  filter(!is.na(value)) %>%
  group_by(feature, group) %>%
  summarise(
    n = n(),
    n_positive = sum(value == 1, na.rm = TRUE),
    proportion = mean(value == 1, na.rm = TRUE),
    percent = proportion * 100,
    .groups = "drop"
  ) %>%
  right_join(
    expand_grid(
      feature = features,
      group = factor(group_levels, levels = group_levels)
    ),
    by = c("feature", "group")
  ) %>%
  mutate(
    n = replace_na(n, 0),
    n_positive = replace_na(n_positive, 0),
    proportion = replace_na(proportion, 0),
    percent = replace_na(percent, 0)
  ) %>%
  left_join(feature_map, by = "feature")

# ── 7. Run binary tests ───────────────────────────────────────────────────────
test_results <- purrr::map_dfr(features, function(feat) {
  bind_rows(
    run_test(df_long, feat, "FS Control", "FS Control2") %>%
      mutate(feature = feat, comparison = "FS"),
    run_test(df_long, feat, "SNV Control", "SNV Control2") %>%
      mutate(feature = feat, comparison = "SNV")
  )
}) %>%
  left_join(feature_map, by = "feature")

# ── 8. Binary p-label table ───────────────────────────────────────────────────
p_label_df <- test_results %>%
  mutate(
    xstart = ifelse(comparison == "FS", 1, 3),
    xend   = ifelse(comparison == "FS", 2, 4),
    xmid   = ifelse(comparison == "FS", 1.5, 3.5),
    label  = format_p(p)
  )

# ── 9. Prepare continuous distance data ───────────────────────────────────────
dist_set2 <- variant_plp_control %>%
  mutate(
    group = case_when(
      source == "fs_control"  ~ "FS Control",
      source == "snv_control" ~ "SNV Control",
      TRUE ~ NA_character_
    ),
    dist_to_cds_end = suppressWarnings(as.numeric(.data[[dist_col_set2]]))
  ) %>%
  select(group, dist_to_cds_end) %>%
  filter(!is.na(group), !is.na(dist_to_cds_end), dist_to_cds_end > 0)

dist_set3 <- control_variants_all2 %>%
  mutate(
    group = case_when(
      group == "fs_control2"  ~ "FS Control2",
      group == "snv_control2" ~ "SNV Control2",
      TRUE ~ NA_character_
    ),
    dist_to_cds_end = suppressWarnings(as.numeric(.data[[dist_col_set3]]))
  ) %>%
  select(group, dist_to_cds_end) %>%
  filter(!is.na(group), !is.na(dist_to_cds_end), dist_to_cds_end > 0)

dist_df <- bind_rows(dist_set2, dist_set3) %>%
  mutate(group = factor(group, levels = group_levels))

# ── 10. Run continuous tests ──────────────────────────────────────────────────
dist_p_df <- bind_rows(
  run_continuous_test(dist_df, "FS Control", "FS Control2") %>%
    mutate(comparison = "FS"),
  run_continuous_test(dist_df, "SNV Control", "SNV Control2") %>%
    mutate(comparison = "SNV")
)

dist_summary <- dist_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_dist = mean(dist_to_cds_end, na.rm = TRUE),
    median_dist = median(dist_to_cds_end, na.rm = TRUE),
    .groups = "drop"
  )

# ── 11. Plot colors ───────────────────────────────────────────────────────────
group_colors <- c(
  "FS Control"   = "#E2C7B5",
  "FS Control2"  = "#D07A3A",
  "SNV Control"  = "#C9D6BE",
  "SNV Control2" = "#8ABA67"
)

# ── 12. Generate all binary panels ────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[ feature_map$panel[i] ]] <- plot_feature(
    feat = feature_map$feature[i],
    title = feature_map$title[i],
    panel_letter = feature_map$panel[i],
    summary_df = summary_df,
    p_label_df = p_label_df,
    group_levels = group_levels,
    group_colors = group_colors
  )
}

# add continuous distance panel
plot_list[["G"]] <- plot_distance(
  dist_df = dist_df,
  dist_p_df = dist_p_df,
  group_levels = group_levels,
  group_colors = group_colors
)

# ── 13. Combine figure ────────────────────────────────────────────────────────
ordered_panels <- c("A", "B", "C", "D", "E", "F", "G")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Comparison of Control Variant Sets Across Protein Features",
    subtitle = paste(
      "FS comparison: fs_control vs fs_control2;",
      "SNV comparison: snv_control vs snv_control2;",
      "Binary features: Fisher's exact test or chi-squared test;",
      "PTC distance to CDS end: Wilcoxon rank-sum test"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0),
      plot.subtitle = element_text(size = 10.5, hjust = 0)
    )
  )

# ── 14. Output tables ─────────────────────────────────────────────────────────
results_table <- summary_df %>%
  select(feature, title, panel, group, n, n_positive, proportion, percent)

pvalue_table <- test_results %>%
  mutate(p_label = format_p(p)) %>%
  select(feature, title, panel, comparison, method, p, p_label, g1_0, g1_1, g2_0, g2_1)

dist_pvalue_table <- dist_p_df %>%
  mutate(
    feature = "dist_to_cds_end",
    title = "PTC Distance to CDS End",
    panel = "G",
    p_label = format_p(p)
  ) %>%
  select(feature, title, panel, comparison, method, p, p_label, n_g1, n_g2)

# ── 15. Print results ─────────────────────────────────────────────────────────
print(results_table)
print(pvalue_table)
print(dist_summary)
print(dist_pvalue_table)
print(combined_fig)

# ── 16. Optional save ─────────────────────────────────────────────────────────
# write_csv(results_table, "control_set_comparison_summary.csv")
# write_csv(pvalue_table, "control_set_comparison_pvalues.csv")
# write_csv(dist_summary, "control_set_distance_summary.csv")
# write_csv(dist_pvalue_table, "control_set_distance_pvalues.csv")
# ggsave("control_set_comparison_plot_with_distance.pdf", combined_fig, width = 17, height = 13, dpi = 300)
# ggsave("control_set_comparison_plot_with_distance.png", combined_fig, width = 17, height = 13, dpi = 300)