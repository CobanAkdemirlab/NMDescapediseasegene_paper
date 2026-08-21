library(tidyverse)
library(patchwork)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)

# ── 1. Read data ──────────────────────────────────────────────────────────────
df <- read_csv("Downloads/variants_all_with_motif_LCS_flags (1).csv", show_col_types = FALSE)

# ── 2. Basic cleaning ─────────────────────────────────────────────────────────
df <- df %>%
  filter(source %in% c("fs", "fs_control", "snv", "snv_control"))

if (!"transcript" %in% colnames(df)) {
  stop("The column 'transcript' is required but was not found.")
}

# ── 3. Define binary features to analyze ──────────────────────────────────────
feature_map <- tibble::tribble(
  ~feature,                  ~title,                                      ~ylab,                 ~panel,
  "variant_ppi_overlap",     "PPI Residue Enrichment",                    "Mean proportion (%)", "A",
  "ptc_before_max_pfam_end", "Pfam Domain Enrichment",                    "Mean proportion (%)", "B",
  "variant_nls_flag",        "Nuclear Localization Signal Enrichment",    "Mean proportion (%)", "C",
  "variant_ptm_flag",        "Post-Translational Modification Enrichment","Mean proportion (%)", "D",
  "variant_slim_flag",       "Short Linear Motif Enrichment",             "Mean proportion (%)", "E",
  "variant_LCS_flag",        "Low Complexity Sequence Enrichment",        "Mean proportion (%)", "F"
)

feature_map <- feature_map %>%
  filter(feature %in% colnames(df))

if (nrow(feature_map) == 0) {
  stop("None of the requested feature columns were found in the CSV.")
}

flag_cols <- feature_map$feature

# ── 4. Convert feature columns to numeric 0/1 if needed ──────────────────────
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

# force binary 0/1
for (cc in flag_cols) {
  df[[cc]] <- ifelse(is.na(df[[cc]]), NA, ifelse(df[[cc]] > 0, 1, 0))
}

# ── 5. Nice group labels ──────────────────────────────────────────────────────
group_labels <- c(
  "fs" = "Frameshift",
  "fs_control" = "FS Control",
  "snv" = "Nonsense",
  "snv_control" = "Nonsense Control"
)

group_levels <- c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")

# ── 6. Resampling function: one variant per source x transcript ───────────────
sample_one_variant_per_transcript <- function(data, iter_id) {
  data %>%
    group_by(source, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    mutate(iter = iter_id)
}

# ── 7. Run resampling ─────────────────────────────────────────────────────────
n_resamples <- 100

resampled_df <- map_dfr(seq_len(n_resamples), function(i) {
  sample_one_variant_per_transcript(df, i)
})

# ── 8. Long format ────────────────────────────────────────────────────────────
resampled_long <- resampled_df %>%
  select(iter, source, transcript, all_of(flag_cols)) %>%
  pivot_longer(
    cols = all_of(flag_cols),
    names_to = "feature",
    values_to = "value"
  ) %>%
  mutate(
    group = recode(as.character(source), !!!group_labels),
    group = factor(group, levels = group_levels)
  )

# ── 9. Helper: format p-value ────────────────────────────────────────────────
format_p <- function(p) {
  case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ sprintf("p = %.3f", p)
  )
}
# ── 10. McNemar test helper ──────────────────────────────────────────────────
# discordant pairs:
# b = study=1, control=0
# c = study=0, control=1
# if b+c <= 25 -> exact McNemar via binom.test
# else -> mcnemar.test with continuity correction
run_mcnemar_from_vectors <- function(study, control, exact_cutoff = 25) {
  keep <- !is.na(study) & !is.na(control)
  study <- study[keep]
  control <- control[keep]
  
  if (length(study) == 0) {
    return(list(
      p = NA_real_,
      test_type = NA_character_,
      total_pairs = 0,
      non_missing_pairs = 0,
      both_0 = 0,
      control_only = 0,
      study_only = 0,
      both_1 = 0,
      discordant_n = 0
    ))
  }
  
  study <- ifelse(study > 0, 1, 0)
  control <- ifelse(control > 0, 1, 0)
  
  both_0 <- sum(study == 0 & control == 0)
  control_only <- sum(study == 0 & control == 1) # c
  study_only   <- sum(study == 1 & control == 0) # b
  both_1 <- sum(study == 1 & control == 1)
  
  discordant_n <- study_only + control_only
  
  if (discordant_n == 0) {
    return(list(
      p = NA_real_,
      test_type = "No discordant pairs",
      total_pairs = length(study),
      non_missing_pairs = length(study),
      both_0 = both_0,
      control_only = control_only,
      study_only = study_only,
      both_1 = both_1,
      discordant_n = discordant_n
    ))
  }
  
  tab <- matrix(
    c(both_0, control_only,
      study_only, both_1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      study = c("0", "1"),
      control = c("0", "1")
    )
  )
  
  p <- tryCatch(
    if (discordant_n <= exact_cutoff) {
      binom.test(
        x = study_only,
        n = discordant_n,
        p = 0.5,
        alternative = "two.sided"
      )$p.value
    } else {
      mcnemar.test(tab, correct = TRUE)$p.value
    },
    error = function(e) NA_real_
  )
  
  test_type <- ifelse(
    discordant_n <= exact_cutoff,
    "Exact McNemar",
    "McNemar chi-squared"
  )
  
  list(
    p = p,
    test_type = test_type,
    total_pairs = length(study),
    non_missing_pairs = length(study),
    both_0 = both_0,
    control_only = control_only,
    study_only = study_only,
    both_1 = both_1,
    discordant_n = discordant_n
  )
}

# ── 11. Matched test for one feature in one iteration ────────────────────────
get_matched_test_one_iter <- function(data, feat, study_group, control_group, exact_cutoff = 25) {
  dsub <- data %>%
    filter(feature == feat, group %in% c(study_group, control_group)) %>%
    select(iter, transcript, group, value) %>%
    pivot_wider(names_from = group, values_from = value)
  
  if (!(study_group %in% colnames(dsub)) || !(control_group %in% colnames(dsub))) {
    return(tibble(
      p = NA_real_,
      test_type = NA_character_,
      total_pairs = 0,
      non_missing_pairs = 0,
      both_0 = 0,
      control_only = 0,
      study_only = 0,
      both_1 = 0,
      discordant_n = 0
    ))
  }
  
  out <- run_mcnemar_from_vectors(
    study = dsub[[study_group]],
    control = dsub[[control_group]],
    exact_cutoff = exact_cutoff
  )
  
  tibble(
    p = out$p,
    test_type = out$test_type,
    total_pairs = out$total_pairs,
    non_missing_pairs = out$non_missing_pairs,
    both_0 = out$both_0,
    control_only = out$control_only,
    study_only = out$study_only,
    both_1 = out$both_1,
    discordant_n = out$discordant_n
  )
}

# ── 12. Run matched tests across all resamples ───────────────────────────────
matched_results <- map_dfr(seq_len(n_resamples), function(i) {
  data_i <- resampled_long %>% filter(iter == i)
  
  map_dfr(feature_map$feature, function(feat_now) {
    bind_rows(
      get_matched_test_one_iter(data_i, feat_now, "Frameshift", "FS Control") %>%
        mutate(iter = i, feature = feat_now, match_type = "Frameshift"),
      get_matched_test_one_iter(data_i, feat_now, "Nonsense", "Nonsense Control") %>%
        mutate(iter = i, feature = feat_now, match_type = "Nonsense")
    )
  })
})

# ── 13. Summarize p-values across 100 resamples ──────────────────────────────
matched_p_summary <- matched_results %>%
  group_by(feature, match_type) %>%
  summarise(
    median_p = median(p, na.rm = TRUE),
    mean_p = mean(p, na.rm = TRUE),
    median_discordant_n = median(discordant_n, na.rm = TRUE),
    mean_discordant_n = mean(discordant_n, na.rm = TRUE),
    n_exact = sum(test_type == "Exact McNemar", na.rm = TRUE),
    n_chisq = sum(test_type == "McNemar chi-squared", na.rm = TRUE),
    n_no_discordant = sum(test_type == "No discordant pairs", na.rm = TRUE),
    .groups = "drop"
  )

# ── 14. Summarize paired counts across 100 resamples ─────────────────────────
matched_n_table <- matched_results %>%
  group_by(feature, match_type) %>%
  summarise(
    total_pairs = round(mean(total_pairs, na.rm = TRUE)),
    non_missing_pairs = round(mean(non_missing_pairs, na.rm = TRUE)),
    both_0 = round(mean(both_0, na.rm = TRUE)),
    control_only = round(mean(control_only, na.rm = TRUE)),
    study_only = round(mean(study_only, na.rm = TRUE)),
    both_1 = round(mean(both_1, na.rm = TRUE)),
    discordant_n = round(mean(discordant_n, na.rm = TRUE)),
    .groups = "drop"
  )

# ── 15. Mean proportion plot summary (y-axis still mean proportion) ──────────
summary_df <- resampled_long %>%
  group_by(feature, group, iter) %>%
  summarise(
    proportion = mean(value, na.rm = TRUE) * 100,
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(feature, group) %>%
  summarise(
    percent = mean(proportion, na.rm = TRUE),
    n = round(mean(n, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  right_join(
    expand_grid(
      feature = feature_map$feature,
      group = factor(group_levels, levels = group_levels)
    ),
    by = c("feature", "group")
  ) %>%
  mutate(
    percent = replace_na(percent, 0),
    n = replace_na(n, 0)
  )

# ── 16. Add p-value labels from matched tests ────────────────────────────────
p_label_df <- matched_p_summary %>%
  mutate(
    group1 = ifelse(match_type == "Frameshift", "Frameshift", "Nonsense"),
    group2 = ifelse(match_type == "Frameshift", "FS Control", "Nonsense Control"),
    xstart = ifelse(match_type == "Frameshift", 1, 3),
    xend   = ifelse(match_type == "Frameshift", 2, 4),
    xmid   = ifelse(match_type == "Frameshift", 1.5, 3.5),
    label  = format_p(median_p)
  )

# ── 17. Plot colors ──────────────────────────────────────────────────────────
group_colors <- c(
  "Frameshift"       = "#D07A3A",
  "FS Control"       = "#E2C7B5",
  "Nonsense"         = "#8ABA67",
  "Nonsense Control" = "#C9D6BE"
)

# ── 18. Plot function ────────────────────────────────────────────────────────
plot_bar <- function(feat, title, ylab, panel_letter) {
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

# ── 19. Generate all panels ──────────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[ feature_map$panel[i] ]] <- plot_bar(
    feat = feature_map$feature[i],
    title = feature_map$title[i],
    ylab = feature_map$ylab[i],
    panel_letter = feature_map$panel[i]
  )
}

# ── 20. Arrange combined figure ──────────────────────────────────────────────
ordered_panels <- c("A", "B", "C", "D", "E", "F")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Variant-level Enrichment of Protein Features After 100 Transcript-Level Resamplings",
    subtitle = paste(
      "100 transcript-level resamplings;",
      "1 random variant selected per transcript per group per resample;",
      "y-axis shows mean proportion (%);",
      "matched comparisons used McNemar's test;",
      "exact McNemar test when discordant pairs ≤ 25;",
      "McNemar chi-squared test when discordant pairs > 25"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0),
      plot.subtitle = element_text(size = 10.5, hjust = 0)
    )
  )

# ── 21. Output summary table for plotting ────────────────────────────────────
summary_table <- summary_df %>%
  left_join(
    feature_map %>% select(feature, title, panel),
    by = "feature"
  ) %>%
  select(feature, title, panel, group, n, percent)

# ── 22. Output p-value summary table ─────────────────────────────────────────
p_table <- matched_p_summary %>%
  left_join(
    feature_map %>% select(feature, title, panel),
    by = "feature"
  ) %>%
  select(feature, title, panel, match_type, everything())

# ── 23. Save outputs ─────────────────────────────────────────────────────────
write_csv(
  summary_table,
  "Supplemental_Figure_3_variant_level_100_resamples_matched_summary.csv"
)

write_csv(
  matched_results,
  "Supplemental_Figure_3_variant_level_100_resamples_mcnemar_per_resample.csv"
)

write_csv(
  matched_n_table,
  "Supplemental_Figure_3_variant_level_100_resamples_mcnemar_counts.csv"
)

write_csv(
  p_table,
  "Supplemental_Figure_3_variant_level_100_resamples_mcnemar_pvalues.csv"
)

ggsave(
  "Supplemental_Figure_3_variant_level_100_resamples_matched_mcnemar.pdf",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

ggsave(
  "Supplemental_Figure_3_variant_level_100_resamples_matched_mcnemar.png",
  combined_fig,
  width = 17,
  height = 10,
  dpi = 300
)

# ── 24. Print summaries ──────────────────────────────────────────────────────
cat("\nCounts per source in original data:\n")
print(table(df$source))

cat("\nMatched McNemar count table:\n")
print(matched_n_table)

cat("\nMatched McNemar p-value summary:\n")
print(p_table)

cat("\nDone!\n")
cat("Output files:\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_matched_summary.csv\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_mcnemar_per_resample.csv\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_mcnemar_counts.csv\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_mcnemar_pvalues.csv\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_matched_mcnemar.pdf\n")
cat("  Supplemental_Figure_3_variant_level_100_resamples_matched_mcnemar.png\n")