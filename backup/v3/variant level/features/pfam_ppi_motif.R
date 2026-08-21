library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(purrr)
library(biomaRt)

set.seed(123)

############################################################
## 0. input
############################################################
# required:
# variants_all2
# human_1_
# ensembl

# make sure source order is fixed
variants_all2$source <- factor(
  variants_all2$source,
  levels = c("snv", "snv_control", "fs", "fs_control")
)

group_colors2 <- c(
  "snv" = "#2ca02c",
  "snv_control" = "#98df8a",
  "fs" = "#1f77b4",
  "fs_control" = "#aec7e8"
)

############################################################
## 1. helper functions
############################################################

convert_to_c <- function(x) {
  if (is.na(x) || x == "") return(numeric(0))
  
  x <- gsub("\\[|\\]", "", x)
  x <- gsub(" ", "", x)
  
  if (x == "") return(numeric(0))
  
  vals <- unlist(strsplit(x, ","))
  vals <- suppressWarnings(as.numeric(vals))
  vals <- vals[!is.na(vals)]
  vals
}

# one variant per transcript per source
sample_one_per_transcript_source <- function(df) {
  df %>%
    filter(!is.na(transcript), transcript != "") %>%
    group_by(transcript, source) %>%
    slice_sample(n = 1) %>%
    ungroup()
}

# empirical two-sided p value from resampling distribution
empirical_p_two_sided <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  2 * min(mean(x >= 0), mean(x <= 0))
}

############################################################
## 2. annotate PPI metrics
############################################################

annotate_ppi_metrics <- function(df, human_1_) {
  
  df <- df %>%
    mutate(
      matched_uniprot = 0,
      nearest_interface_bp = NA_real_,
      dist_to_nearest_interface_bp = NA_real_
    )
  
  for (i in seq_len(nrow(df))) {
    
    uid <- df$uniprotswissprot[i]
    cds_ptc_loc <- df$cds_ptc_loc[i]
    
    if (is.na(uid) || uid == "") next
    if (is.na(cds_ptc_loc)) next
    
    re_1 <- unlist(lapply(
      human_1_$interface_residues1[human_1_$uniprot1 == uid],
      convert_to_c
    )) * 3
    
    re_2 <- unlist(lapply(
      human_1_$interface_residues2[human_1_$uniprot2 == uid],
      convert_to_c
    )) * 3
    
    all_interface_bp <- c(re_1, re_2)
    all_interface_bp <- all_interface_bp[!is.na(all_interface_bp)]
    
    if (length(all_interface_bp) == 0) next
    
    downstream_interface_bp <- all_interface_bp[all_interface_bp >= cds_ptc_loc]
    
    if (length(downstream_interface_bp) > 0) {
      df$matched_uniprot[i] <- 1
      df$nearest_interface_bp[i] <- min(downstream_interface_bp, na.rm = TRUE)
      df$dist_to_nearest_interface_bp[i] <- min(downstream_interface_bp, na.rm = TRUE) - cds_ptc_loc
    }
  }
  
  df
}

############################################################
## 3. annotate PFAM metrics
############################################################

annotate_pfam_metrics <- function(df, ensembl) {
  
  pfam_domains2 <- getBM(
    attributes = c("ensembl_transcript_id", "pfam_end"),
    filters    = "ensembl_transcript_id",
    values     = unique(df$transcript),
    mart       = ensembl
  )
  
  max_pfam_end_df <- pfam_domains2 %>%
    filter(!is.na(pfam_end)) %>%
    group_by(ensembl_transcript_id) %>%
    summarise(max_pfam_end = max(pfam_end), .groups = "drop")
  
  df %>%
    mutate(
      ptc_aa = ceiling(cds_ptc_loc / 3)
    ) %>%
    left_join(max_pfam_end_df, by = c("transcript" = "ensembl_transcript_id")) %>%
    mutate(
      dist_ptc_to_max_pfam_end_aa = ptc_aa - max_pfam_end,
      ptc_after_max_pfam_end  = ifelse(!is.na(dist_ptc_to_max_pfam_end_aa) &
                                         dist_ptc_to_max_pfam_end_aa > 0, 1, 0),
      ptc_before_max_pfam_end = ifelse(!is.na(dist_ptc_to_max_pfam_end_aa) &
                                         dist_ptc_to_max_pfam_end_aa < 0, 1, 0)
    )
}

############################################################
## 4. summarize one iteration
############################################################

summarise_one_resample <- function(df) {
  
  ppi_prop <- df %>%
    group_by(source) %>%
    summarise(
      n = n(),
      n_ppi_match = sum(matched_uniprot == 1, na.rm = TRUE),
      ppi_prop_match = mean(matched_uniprot == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  ppi_dist <- df %>%
    filter(!is.na(dist_to_nearest_interface_bp)) %>%
    group_by(source) %>%
    summarise(
      ppi_dist_median = median(dist_to_nearest_interface_bp, na.rm = TRUE),
      ppi_dist_mean   = mean(dist_to_nearest_interface_bp, na.rm = TRUE),
      .groups = "drop"
    )
  
  pfam_prop <- df %>%
    group_by(source) %>%
    summarise(
      n_pfam_before = sum(ptc_before_max_pfam_end == 1, na.rm = TRUE),
      pfam_prop_before = mean(ptc_before_max_pfam_end == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  pfam_dist <- df %>%
    filter(!is.na(dist_ptc_to_max_pfam_end_aa)) %>%
    group_by(source) %>%
    summarise(
      pfam_dist_median = median(dist_ptc_to_max_pfam_end_aa, na.rm = TRUE),
      pfam_dist_abs_median = median(abs(dist_ptc_to_max_pfam_end_aa), na.rm = TRUE),
      .groups = "drop"
    )
  
  ppi_prop %>%
    left_join(ppi_dist, by = "source") %>%
    left_join(pfam_prop, by = "source") %>%
    left_join(pfam_dist, by = "source")
}

############################################################
## 5. run 1000 resamples
############################################################

n_resample <- 1000
resample_results <- vector("list", n_resample)

for (b in seq_len(n_resample)) {
  
  df_b <- sample_one_per_transcript_source(variants_all2)
  df_b <- annotate_ppi_metrics(df_b, human_1_)
  df_b <- annotate_pfam_metrics(df_b, ensembl)
  
  resample_results[[b]] <- summarise_one_resample(df_b) %>%
    mutate(resample_id = b)
}

resample_summary_all <- bind_rows(resample_results)

############################################################
## 6. summarize across 1000 resamples
############################################################

summary_by_group <- resample_summary_all %>%
  group_by(source) %>%
  summarise(
    ppi_prop_match_median = median(ppi_prop_match, na.rm = TRUE),
    ppi_prop_match_lwr    = quantile(ppi_prop_match, 0.025, na.rm = TRUE),
    ppi_prop_match_upr    = quantile(ppi_prop_match, 0.975, na.rm = TRUE),
    
    ppi_dist_median_median = median(ppi_dist_median, na.rm = TRUE),
    ppi_dist_median_lwr    = quantile(ppi_dist_median, 0.025, na.rm = TRUE),
    ppi_dist_median_upr    = quantile(ppi_dist_median, 0.975, na.rm = TRUE),
    
    pfam_prop_before_median = median(pfam_prop_before, na.rm = TRUE),
    pfam_prop_before_lwr    = quantile(pfam_prop_before, 0.025, na.rm = TRUE),
    pfam_prop_before_upr    = quantile(pfam_prop_before, 0.975, na.rm = TRUE),
    
    pfam_dist_median_median = median(pfam_dist_median, na.rm = TRUE),
    pfam_dist_median_lwr    = quantile(pfam_dist_median, 0.025, na.rm = TRUE),
    pfam_dist_median_upr    = quantile(pfam_dist_median, 0.975, na.rm = TRUE),
    
    pfam_dist_abs_median_median = median(pfam_dist_abs_median, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_group)

############################################################
## 7. pairwise resampling-based tests
############################################################

get_pairwise_resample_diff <- function(df, value_col, g1, g2) {
  wide_df <- df %>%
    select(resample_id, source, !!sym(value_col)) %>%
    pivot_wider(names_from = source, values_from = !!sym(value_col))
  
  wide_df[[g1]] - wide_df[[g2]]
}

pairwise_tests <- bind_rows(
  tibble(
    metric = "ppi_prop_match",
    comparison = "snv - snv_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "ppi_prop_match", "snv", "snv_control"))
  ),
  tibble(
    metric = "ppi_prop_match",
    comparison = "fs - fs_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "ppi_prop_match", "fs", "fs_control"))
  ),
  tibble(
    metric = "ppi_dist_median",
    comparison = "snv - snv_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "ppi_dist_median", "snv", "snv_control"))
  ),
  tibble(
    metric = "ppi_dist_median",
    comparison = "fs - fs_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "ppi_dist_median", "fs", "fs_control"))
  ),
  tibble(
    metric = "pfam_prop_before",
    comparison = "snv - snv_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "pfam_prop_before", "snv", "snv_control"))
  ),
  tibble(
    metric = "pfam_prop_before",
    comparison = "fs - fs_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "pfam_prop_before", "fs", "fs_control"))
  ),
  tibble(
    metric = "pfam_dist_median",
    comparison = "snv - snv_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "pfam_dist_median", "snv", "snv_control"))
  ),
  tibble(
    metric = "pfam_dist_median",
    comparison = "fs - fs_control",
    diff = list(get_pairwise_resample_diff(resample_summary_all, "pfam_dist_median", "fs", "fs_control"))
  )
) %>%
  rowwise() %>%
  mutate(
    median_diff = median(diff, na.rm = TRUE),
    lwr_95 = quantile(diff, 0.025, na.rm = TRUE),
    upr_95 = quantile(diff, 0.975, na.rm = TRUE),
    emp_p = empirical_p_two_sided(unlist(diff))
  ) %>%
  ungroup() %>%
  select(metric, comparison, median_diff, lwr_95, upr_95, emp_p)

print(pairwise_tests)

############################################################
## 8. plots of 1000-resample distributions
############################################################

plot_metric_distribution <- function(df, metric_col, title_text, ylab_text) {
  ggplot(df, aes(x = source, y = .data[[metric_col]], fill = source)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.4, color = "black") +
    scale_fill_manual(values = group_colors2, guide = "none") +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = NULL,
      y = ylab_text
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

p_ppi_prop <- plot_metric_distribution(
  resample_summary_all,
  "ppi_prop_match",
  "1000 resamples: proportion with downstream PPI interface residue",
  "Proportion"
)

p_ppi_dist <- plot_metric_distribution(
  resample_summary_all,
  "ppi_dist_median",
  "1000 resamples: median distance from PTC to downstream PPI interface",
  "Median distance (bp)"
)

p_pfam_prop <- plot_metric_distribution(
  resample_summary_all,
  "pfam_prop_before",
  "1000 resamples: proportion with PTC before largest PFAM end",
  "Proportion"
)

p_pfam_dist <- plot_metric_distribution(
  resample_summary_all,
  "pfam_dist_median",
  "1000 resamples: median distance from PTC to largest PFAM end",
  "Median distance (aa)"
)

print(p_ppi_prop)
print(p_ppi_dist)
print(p_pfam_prop)
print(p_pfam_dist)

############################################################
## 9. optional: one final sampled dataset for original-style plots
############################################################

df_one <- sample_one_per_transcript_source(variants_all2)
df_one <- annotate_ppi_metrics(df_one, human_1_)
df_one <- annotate_pfam_metrics(df_one, ensembl)

# PPI proportion barplot
ppi_flag_summary <- df_one %>%
  group_by(source) %>%
  summarise(
    n = n(),
    n_match = sum(matched_uniprot == 1, na.rm = TRUE),
    prop_match = n_match / n,
    .groups = "drop"
  )

tab_snv <- matrix(c(
  ppi_flag_summary$n_match[ppi_flag_summary$source == "snv"],
  ppi_flag_summary$n[ppi_flag_summary$source == "snv"] -
    ppi_flag_summary$n_match[ppi_flag_summary$source == "snv"],
  ppi_flag_summary$n_match[ppi_flag_summary$source == "snv_control"],
  ppi_flag_summary$n[ppi_flag_summary$source == "snv_control"] -
    ppi_flag_summary$n_match[ppi_flag_summary$source == "snv_control"]
), nrow = 2, byrow = TRUE)

tab_fs <- matrix(c(
  ppi_flag_summary$n_match[ppi_flag_summary$source == "fs"],
  ppi_flag_summary$n[ppi_flag_summary$source == "fs"] -
    ppi_flag_summary$n_match[ppi_flag_summary$source == "fs"],
  ppi_flag_summary$n_match[ppi_flag_summary$source == "fs_control"],
  ppi_flag_summary$n[ppi_flag_summary$source == "fs_control"] -
    ppi_flag_summary$n_match[ppi_flag_summary$source == "fs_control"]
), nrow = 2, byrow = TRUE)

p_snv <- fisher.test(tab_snv)$p.value
p_fs  <- fisher.test(tab_fs)$p.value

ymax <- max(ppi_flag_summary$prop_match, na.rm = TRUE)

pval_df <- data.frame(
  group1 = c("snv", "fs"),
  group2 = c("snv_control", "fs_control"),
  y.position = c(ymax * 1.08, ymax * 1.20),
  label = c(
    paste0("p = ", signif(p_snv, 3)),
    paste0("p = ", signif(p_fs, 3))
  )
)

p1 <- ggplot(ppi_flag_summary, aes(x = source, y = prop_match, fill = source)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Proportion of variants with downstream PPI interface residue",
    x = NULL,
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(aes(label = paste0(n_match, "/", n)), vjust = -0.3, size = 4) +
  stat_pvalue_manual(
    pval_df,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  ) +
  expand_limits(y = ymax * 1.30)

print(p1)

# PPI distance plot
p2 <- ggplot(
  df_one %>% filter(!is.na(dist_to_nearest_interface_bp)),
  aes(x = source, y = dist_to_nearest_interface_bp, fill = source)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to nearest downstream PPI interface residue",
    x = NULL,
    y = "Distance (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(
    comparisons = list(c("snv", "snv_control"), c("fs", "fs_control")),
    method = "wilcox.test",
    label = "p.format"
  )

print(p2)

# PFAM distance plot
p3 <- ggplot(
  df_one %>% filter(!is.na(dist_ptc_to_max_pfam_end_aa)),
  aes(x = source, y = dist_ptc_to_max_pfam_end_aa, fill = source)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to the largest PFAM end",
    x = NULL,
    y = "PTC aa position - largest PFAM end (aa)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(
    comparisons = list(c("snv", "snv_control"), c("fs", "fs_control")),
    method = "wilcox.test",
    label = "p.format"
  )

print(p3)

############################################################
## 10. save
############################################################

write.csv(resample_summary_all, "resample_1000_transcript_source_summary.csv", row.names = FALSE)
write.csv(summary_by_group, "resample_1000_transcript_source_group_summary.csv", row.names = FALSE)
write.csv(pairwise_tests, "resample_1000_transcript_source_pairwise_tests.csv", row.names = FALSE)
write.csv(df_one, "one_resampled_variants_with_ppi_pfam.csv", row.names = FALSE)