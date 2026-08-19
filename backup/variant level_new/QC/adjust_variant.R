############################################################
## Clean code for:
## 1) Full/unmatched 4-group violin+box plots
##    - CDS length
##    - NMDesc region length
##    - Distance from PTC to CDS end
## 2) Matched 4-group violin+box plots
##    - same three outcomes
## 3) Unmixed and mixed models
##    - full FS, matched FS, full SNV, matched SNV
## 4) Repeated sampling (1000x) on matched variants
##    - one variant per transcript per source
##    - estimate and fold-change histograms
############################################################

############################
## 0. Packages
############################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)

############################
## 1. Prepare data
############################
variants_all2 <- variants_all %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control")),
    NMDesc_region_length = NMD_region_end - NMD_region_start + 1
  ) %>%
   mutate(
    log_dist_to_cds_end = log10(dist_to_cds_end + 1),
    log_cds_end = log10(cds_end),
    log_NMDesc_region_length = log10(NMDesc_region_length)
  ) 

############################
## 2. Shared transcripts for matched analysis
############################
fs_tr <- variants_all2 %>%
  filter(source == "fs") %>%
  pull(transcript) %>%
  unique()

fs_ctrl_tr <- variants_all2 %>%
  filter(source == "fs_control") %>%
  pull(transcript) %>%
  unique()

snv_tr <- variants_all2 %>%
  filter(source == "snv") %>%
  pull(transcript) %>%
  unique()

snv_ctrl_tr <- variants_all2 %>%
  filter(source == "snv_control") %>%
  pull(transcript) %>%
  unique()

shared_fs <- intersect(fs_tr, fs_ctrl_tr)
shared_snv <- intersect(snv_tr, snv_ctrl_tr)

############################
## 3. Datasets
############################
## full/unmatched, pairwise for models
dat_full_fs <- variants_all2 %>%
  filter(source %in% c("fs", "fs_control")) %>%
  droplevels()

dat_full_snv <- variants_all2 %>%
  filter(source %in% c("snv", "snv_control")) %>%
  droplevels()

## matched, pairwise for models
dat_matched_fs <- variants_all2 %>%
  filter(
    transcript %in% shared_fs,
    source %in% c("fs", "fs_control")
  ) %>%
  droplevels()

dat_matched_snv <- variants_all2 %>%
  filter(
    transcript %in% shared_snv,
    source %in% c("snv", "snv_control")
  ) %>%
  droplevels()

## matched, all four groups for plots
variants_matched4 <- bind_rows(
  variants_all2 %>%
    filter(
      transcript %in% shared_fs,
      source %in% c("fs", "fs_control")
    ),
  variants_all2 %>%
    filter(
      transcript %in% shared_snv,
      source %in% c("snv", "snv_control")
    )
) %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control"))
  )


############################
## 5. Plot settings
############################
group_colors <- c(
  "fs" = "#4C78A8",
  "fs_control" = "#9FB7D3",
  "snv" = "#59A14F",
  "snv_control" = "#8CD17D"
)

pair_comparisons <- list(
  c("fs", "fs_control"),
  c("snv", "snv_control")
)

plot_feature_4groups <- function(dat, yvar, ylab, title_text) {
  ggplot(dat, aes(x = source, y = .data[[yvar]], fill = source)) +
    geom_violin(trim = TRUE, alpha = 0.85, color = "black") +
    geom_boxplot(width = 0.15, outlier.size = 0.8, color = "black") +
    scale_y_log10() +
    scale_fill_manual(values = group_colors, guide = "none") +
    labs(
      title = title_text,
      x = NULL,
      y = ylab
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title.y = element_text(face = "bold")
    ) +
    stat_compare_means(
      comparisons = pair_comparisons,
      method = "wilcox.test",
      label = "p.format",
      step.increase = 0.12
    )
}

############################
## 6. Full/unmatched plots
############################
p_cds_full <- plot_feature_4groups(
  dat = variants_all2,
  yvar = "cds_end",
  ylab = "CDS length (bp, log10 scale)",
  title_text = "CDS length across variant groups"
)

p_nmdesc_full <- plot_feature_4groups(
  dat = variants_all2,
  yvar = "NMDesc_region_length",
  ylab = "NMDesc region length (bp, log10 scale)",
  title_text = "NMDesc region length across variant groups"
)

p_dist_full <- plot_feature_4groups(
  dat = variants_all2,
  yvar = "dist_to_cds_end",
  ylab = "Distance to CDS end (bp, log10 scale)",
  title_text = "Distance from PTC to CDS end across variant groups"
)

print(p_cds_full)
print(p_nmdesc_full)
print(p_dist_full)

############################
## 7. Matched plots (still four groups)
############################
p_cds_matched <- plot_feature_4groups(
  dat = variants_matched4,
  yvar = "cds_end",
  ylab = "CDS length (bp, log10 scale)",
  title_text = "Matched analysis: CDS length"
)

p_nmdesc_matched <- plot_feature_4groups(
  dat = variants_matched4,
  yvar = "NMDesc_region_length",
  ylab = "NMDesc region length (bp, log10 scale)",
  title_text = "Matched analysis: NMDesc region length"
)

p_dist_matched <- plot_feature_4groups(
  dat = variants_matched4,
  yvar = "dist_to_cds_end",
  ylab = "Distance to CDS end (bp, log10 scale)",
  title_text = "Matched analysis: distance from PTC to CDS end"
)

print(p_cds_matched)
print(p_nmdesc_matched)
print(p_dist_matched)

############################
## 8. Unmixed and mixed models
## outcome = log_dist_to_cds_end
## unmixed: log_dist_to_cds_end ~ source
## mixed:   log_dist_to_cds_end ~ source + (1 | transcript)
############################
run_models <- function(dat, dataset_name) {
  
  cat("\n========================================\n")
  cat("Dataset:", dataset_name, "\n")
  cat("========================================\n")
  
  cat("\nSource counts:\n")
  print(table(dat$source))
  
  cat("\nNumber of transcripts:\n")
  print(length(unique(dat$transcript)))
  
  ## unmixed
  fit_unmixed <- lm(
    log_dist_to_cds_end ~ source,
    data = dat
  )
  
  cat("\n----------------------------------------\n")
  cat("Unmixed model\n")
  cat("Formula: log_dist_to_cds_end ~ source\n")
  cat("----------------------------------------\n")
  print(summary(fit_unmixed))
  
  unmixed_coef <- as.data.frame(summary(fit_unmixed)$coefficients)
  cat("\nUnmixed coefficients:\n")
  print(unmixed_coef)
  
  ## mixed
  fit_mixed <- lmer(
    log_dist_to_cds_end ~ source + (1 | transcript),
    data = dat,
    REML = TRUE
  )
  
  cat("\n----------------------------------------\n")
  cat("Mixed model\n")
  cat("Formula: log_dist_to_cds_end ~ source + (1 | transcript)\n")
  cat("----------------------------------------\n")
  print(summary(fit_mixed))
  
  mixed_coef <- as.data.frame(summary(fit_mixed)$coefficients)
  cat("\nMixed fixed effects:\n")
  print(mixed_coef)
  
  cat("\nMixed random effects:\n")
  print(VarCorr(fit_mixed))
  
  ## fold change on original scale
  cat("\nFold change on original scale (10^estimate):\n")
  if ("sourcefs_control" %in% rownames(unmixed_coef)) {
    cat("Unmixed fs_control / fs =", 10^(unmixed_coef["sourcefs_control", "Estimate"]), "\n")
  }
  if ("sourcesnv_control" %in% rownames(unmixed_coef)) {
    cat("Unmixed snv_control / snv =", 10^(unmixed_coef["sourcesnv_control", "Estimate"]), "\n")
  }
  if ("sourcefs_control" %in% rownames(mixed_coef)) {
    cat("Mixed fs_control / fs =", 10^(mixed_coef["sourcefs_control", "Estimate"]), "\n")
  }
  if ("sourcesnv_control" %in% rownames(mixed_coef)) {
    cat("Mixed snv_control / snv =", 10^(mixed_coef["sourcesnv_control", "Estimate"]), "\n")
  }
  
  return(list(
    data = dat,
    unmixed = fit_unmixed,
    mixed = fit_mixed,
    unmixed_coef = unmixed_coef,
    mixed_coef = mixed_coef
  ))
}

############################
## 9. Run all four model sets
############################
res_full_fs <- run_models(dat_full_fs, "Full_FS")
res_matched_fs <- run_models(dat_matched_fs, "Matched_FS")
res_full_snv <- run_models(dat_full_snv, "Full_SNV")
res_matched_snv <- run_models(dat_matched_snv, "Matched_SNV")

############################
## 10. Summary table for key source effects
############################
extract_key_result <- function(res_obj, dataset_name, model_type, coef_name) {
  
  coef_tab <- if (model_type == "unmixed") res_obj$unmixed_coef else res_obj$mixed_coef
  
  if (!(coef_name %in% rownames(coef_tab))) return(NULL)
  
  data.frame(
    dataset = dataset_name,
    model = model_type,
    term = coef_name,
    estimate = coef_tab[coef_name, "Estimate"],
    std_error = coef_tab[coef_name, "Std. Error"],
    stat = if ("t value" %in% colnames(coef_tab)) coef_tab[coef_name, "t value"] else coef_tab[coef_name, "z value"],
    fold_change = 10^(coef_tab[coef_name, "Estimate"]),
    row.names = NULL
  )
}

results_table <- bind_rows(
  extract_key_result(res_full_fs, "Full_FS", "unmixed", "sourcefs_control"),
  extract_key_result(res_full_fs, "Full_FS", "mixed",   "sourcefs_control"),
  extract_key_result(res_matched_fs, "Matched_FS", "unmixed", "sourcefs_control"),
  extract_key_result(res_matched_fs, "Matched_FS", "mixed",   "sourcefs_control"),
  extract_key_result(res_full_snv, "Full_SNV", "unmixed", "sourcesnv_control"),
  extract_key_result(res_full_snv, "Full_SNV", "mixed",   "sourcesnv_control"),
  extract_key_result(res_matched_snv, "Matched_SNV", "unmixed", "sourcesnv_control"),
  extract_key_result(res_matched_snv, "Matched_SNV", "mixed",   "sourcesnv_control")
)

cat("\n========================================\n")
cat("Summary table of key source effects\n")
cat("========================================\n")
print(results_table)

write.csv(results_table, "model_comparison_unmixed_vs_mixed.csv", row.names = FALSE)

############################
## 11. Repeated sampling (1000x) on matched variants
## one variant per transcript per source
## model: log_dist_to_cds_end ~ source
############################
run_one_sample_unmixed <- function(dat, case_group, control_group, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dat_sub <- dat %>%
    filter(source %in% c(case_group, control_group)) %>%
    droplevels()
  
  sampled_dat <- dat_sub %>%
    group_by(source, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    droplevels()
  
  fit <- lm(log_dist_to_cds_end ~ source, data = sampled_dat)
  coef_tab <- summary(fit)$coefficients
  
  term_name <- paste0("source", control_group)
  
  data.frame(
    n_variants = nrow(sampled_dat),
    n_transcripts = length(unique(sampled_dat$transcript)),
    estimate = coef_tab[term_name, "Estimate"],
    std_error = coef_tab[term_name, "Std. Error"],
    t_value = coef_tab[term_name, "t value"],
    p_value = coef_tab[term_name, "Pr(>|t|)"],
    fold_change = 10^(coef_tab[term_name, "Estimate"]),
    stringsAsFactors = FALSE
  )
}

n_iter <- 1000

## FS repeated sampling
set.seed(1)
res_fs_repeat <- bind_rows(
  lapply(seq_len(n_iter), function(i) {
    run_one_sample_unmixed(
      dat = dat_matched_fs,
      case_group = "fs",
      control_group = "fs_control",
      seed = i
    ) %>%
      mutate(iter = i, comparison = "fs_control_vs_fs")
  })
)

## SNV repeated sampling
set.seed(1)
res_snv_repeat <- bind_rows(
  lapply(seq_len(n_iter), function(i) {
    run_one_sample_unmixed(
      dat = dat_matched_snv,
      case_group = "snv",
      control_group = "snv_control",
      seed = i
    ) %>%
      mutate(iter = i, comparison = "snv_control_vs_snv")
  })
)

############################
## 12. Summaries of repeated sampling
############################
cat("\n====================================\n")
cat("FS repeated sampling results\n")
cat("====================================\n")
print(summary(res_fs_repeat$estimate))
print(summary(res_fs_repeat$fold_change))
cat("\nProportion p < 0.05:", mean(res_fs_repeat$p_value < 0.05), "\n")
cat("Proportion estimate < 0:", mean(res_fs_repeat$estimate < 0), "\n")
cat("Median fold change:", median(res_fs_repeat$fold_change), "\n\n")

cat("====================================\n")
cat("SNV repeated sampling results\n")
cat("====================================\n")
print(summary(res_snv_repeat$estimate))
print(summary(res_snv_repeat$fold_change))
cat("\nProportion p < 0.05:", mean(res_snv_repeat$p_value < 0.05), "\n")
cat("Proportion estimate < 0:", mean(res_snv_repeat$estimate < 0), "\n")
cat("Median fold change:", median(res_snv_repeat$fold_change), "\n\n")

fs_summary_table <- data.frame(
  comparison = "fs_control_vs_fs",
  n_iter = n_iter,
  n_transcripts = unique(res_fs_repeat$n_transcripts),
  mean_estimate = mean(res_fs_repeat$estimate),
  median_estimate = median(res_fs_repeat$estimate),
  mean_fold_change = mean(res_fs_repeat$fold_change),
  median_fold_change = median(res_fs_repeat$fold_change),
  prop_p_less_0_05 = mean(res_fs_repeat$p_value < 0.05),
  prop_estimate_less_0 = mean(res_fs_repeat$estimate < 0)
)

snv_summary_table <- data.frame(
  comparison = "snv_control_vs_snv",
  n_iter = n_iter,
  n_transcripts = unique(res_snv_repeat$n_transcripts),
  mean_estimate = mean(res_snv_repeat$estimate),
  median_estimate = median(res_snv_repeat$estimate),
  mean_fold_change = mean(res_snv_repeat$fold_change),
  median_fold_change = median(res_snv_repeat$fold_change),
  prop_p_less_0_05 = mean(res_snv_repeat$p_value < 0.05),
  prop_estimate_less_0 = mean(res_snv_repeat$estimate < 0)
)

cat("====================================\n")
cat("FS summary table\n")
cat("====================================\n")
print(fs_summary_table)

cat("\n====================================\n")
cat("SNV summary table\n")
cat("====================================\n")
print(snv_summary_table)

write.csv(res_fs_repeat, "repeated_sampling_1000_matched_fs.csv", row.names = FALSE)
write.csv(res_snv_repeat, "repeated_sampling_1000_matched_snv.csv", row.names = FALSE)
write.csv(bind_rows(fs_summary_table, snv_summary_table),
          "repeated_sampling_1000_matched_summary.csv",
          row.names = FALSE)

############################
## 13. Repeated sampling plots: estimate distributions
############################
p_est_fs <- ggplot(res_fs_repeat, aes(x = estimate)) +
  geom_histogram(bins = 40, color = "black", fill = "#9FB7D3") +
  geom_vline(xintercept = 0, linetype = 2, color = "red", linewidth = 1) +
  theme_classic(base_size = 18) +
  labs(
    title = "Repeated sampling (1000x): FS estimate distribution",
    x = "Estimate for fs_control vs fs",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

p_est_snv <- ggplot(res_snv_repeat, aes(x = estimate)) +
  geom_histogram(bins = 30, color = "black", fill = "#8CD17D") +
  geom_vline(xintercept = 0, linetype = 2, color = "red", linewidth = 1) +
  theme_classic(base_size = 18) +
  labs(
    title = "Repeated sampling (1000x): SNV estimate distribution",
    x = "Estimate for snv_control vs snv",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p_est_fs)
print(p_est_snv)

############################
## 14. Repeated sampling plots: fold-change distributions
############################
p_fc_fs <- ggplot(res_fs_repeat, aes(x = fold_change)) +
  geom_histogram(bins = 40, color = "black", fill = "#9FB7D3") +
  geom_vline(xintercept = 1, linetype = 2, color = "red", linewidth = 1) +
  theme_classic(base_size = 18) +
  labs(
    title = "Frameshift variants: control vs case distance to CDS end (ratio)",
    subtitle = paste0(
      "Median = ", round(median(res_fs_repeat$fold_change), 2)
    ),
    x = "Distance ratio (fs_control / fs)",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

p_fc_snv <- ggplot(res_snv_repeat, aes(x = fold_change)) +
  geom_histogram(bins = 30, color = "black", fill = "#8CD17D") +
  geom_vline(xintercept = 1, linetype = 2, color = "red", linewidth = 1) +
  theme_classic(base_size = 18) +
  labs(
    title = "Nonsense variants: control vs case distance to CDS end (ratio)",
    subtitle = paste0(
      "Median = ", round(median(res_snv_repeat$fold_change), 2)
    ),
    x = "Distance ratio (snv_control / snv)",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

print(p_fc_snv)
print(p_fc_fs)

############################
## 15. Save plots
############################
ggsave("cds_length_full_4groups.png", p_cds_full, width = 9, height = 6, dpi = 300)
ggsave("nmdesc_region_length_full_4groups.png", p_nmdesc_full, width = 9, height = 6, dpi = 300)
ggsave("dist_to_cds_end_full_4groups.png", p_dist_full, width = 9, height = 6, dpi = 300)

ggsave("cds_length_matched_4groups.png", p_cds_matched, width = 9, height = 6, dpi = 300)
ggsave("nmdesc_region_length_matched_4groups.png", p_nmdesc_matched, width = 9, height = 6, dpi = 300)
ggsave("dist_to_cds_end_matched_4groups.png", p_dist_matched, width = 9, height = 6, dpi = 300)

ggsave("repeat1000_fs_estimate_hist.png", p_est_fs, width = 8, height = 5, dpi = 300)
ggsave("repeat1000_snv_estimate_hist.png", p_est_snv, width = 8, height = 5, dpi = 300)
ggsave("repeat1000_fs_foldchange_hist.png", p_fc_fs, width = 8, height = 5, dpi = 300)
ggsave("repeat1000_snv_foldchange_hist.png", p_fc_snv, width = 8, height = 5, dpi = 300)