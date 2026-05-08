############################
## Repeated sampling 1000 times on matched variants
## One variant per transcript per source
## Since in your dataset one gene corresponds to one transcript,
## we use transcript as the unit.
############################

############################
## 0. Packages
############################
library(dplyr)
library(ggplot2)

############################
## 1. Prepare data
############################
variants_all2 <- variants_all %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control")),
    NMDesc_region_length = NMD_region_end - NMD_region_start + 1
  ) %>%
  filter(
    !is.na(source),
    !is.na(dist_to_cds_end),
    dist_to_cds_end >= 0,
    !is.na(cds_end),
    cds_end > 0,
    !is.na(cds_ptc_loc),
    cds_ptc_loc > 0,
    !is.na(NMD_region_start),
    !is.na(NMD_region_end),
    !is.na(NMDesc_region_length),
    NMDesc_region_length > 0
  ) %>%
  mutate(
    log_dist_to_cds_end = log10(dist_to_cds_end + 1),
    log_cds_end = log10(cds_end),
    log_NMDesc_region_length = log10(NMDesc_region_length)
  ) %>%
  filter(
    is.finite(log_dist_to_cds_end),
    is.finite(log_cds_end),
    is.finite(log_NMDesc_region_length)
  )

############################
## 2. Define matched transcript sets
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
## 3. Create matched datasets
############################
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

############################
## 4. Quick check
############################
cat("Matched FS variants:", nrow(dat_matched_fs), "\n")
cat("Matched FS transcripts:", length(unique(dat_matched_fs$transcript)), "\n\n")

cat("Matched SNV variants:", nrow(dat_matched_snv), "\n")
cat("Matched SNV transcripts:", length(unique(dat_matched_snv$transcript)), "\n\n")

############################
## 5. Function: one random sample per transcript, then run unmixed model
## Model: log_dist_to_cds_end ~ source
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
  
  out <- data.frame(
    n_variants = nrow(sampled_dat),
    n_transcripts = length(unique(sampled_dat$transcript)),
    estimate = coef_tab[term_name, "Estimate"],
    std_error = coef_tab[term_name, "Std. Error"],
    t_value = coef_tab[term_name, "t value"],
    p_value = coef_tab[term_name, "Pr(>|t|)"],
    fold_change = 10^(coef_tab[term_name, "Estimate"]),
    stringsAsFactors = FALSE
  )
  
  return(out)
}

############################
## 6. Repeat 1000 times
############################
n_iter <- 1000

## FS
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

## SNV
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
## 7. Summarize results
############################
cat("====================================\n")
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

############################
## 8. Print compact summary tables
############################
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

############################
## 9. Combine results for export
############################
res_repeat_all <- bind_rows(res_fs_repeat, res_snv_repeat)
summary_repeat_all <- bind_rows(fs_summary_table, snv_summary_table)

write.csv(res_fs_repeat, "repeated_sampling_1000_matched_fs.csv", row.names = FALSE)
write.csv(res_snv_repeat, "repeated_sampling_1000_matched_snv.csv", row.names = FALSE)
write.csv(res_repeat_all, "repeated_sampling_1000_matched_all.csv", row.names = FALSE)
write.csv(summary_repeat_all, "repeated_sampling_1000_matched_summary.csv", row.names = FALSE)

############################
## 10. Plot distributions of estimates
############################
p_est_fs <- ggplot(res_fs_repeat, aes(x = estimate)) +
  geom_histogram(bins = 40, color = "black", fill = "#9FB7D3") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  theme_classic(base_size = 16) +
  labs(
    title = "Repeated sampling (1000x): FS estimate distribution",
    x = "Estimate for fs_control vs fs",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_est_snv <- ggplot(res_snv_repeat, aes(x = estimate)) +
  geom_histogram(bins = 40, color = "black", fill = "#8CD17D") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  theme_classic(base_size = 16) +
  labs(
    title = "Repeated sampling (1000x): SNV estimate distribution",
    x = "Estimate for snv_control vs snv",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_est_fs)
print(p_est_snv)

############################
## 11. Plot distributions of fold change
############################
p_fc_fs <- ggplot(res_fs_repeat, aes(x = fold_change)) +
  geom_histogram(bins = 40, color = "black", fill = "#9FB7D3") +
  geom_vline(xintercept = 1, linetype = 2, color = "red") +
  theme_classic(base_size = 16) +
  labs(
    title = "Repeated sampling (1000x): FS fold-change distribution",
    x = "Fold change (fs_control / fs)",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_fc_snv <- ggplot(res_snv_repeat, aes(x = fold_change)) +
  geom_histogram(bins = 40, color = "black", fill = "#8CD17D") +
  geom_vline(xintercept = 1, linetype = 2, color = "red") +
  theme_classic(base_size = 16) +
  labs(
    title = "Repeated sampling (1000x): SNV fold-change distribution",
    x = "Fold change (snv_control / snv)",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_fc_fs)
print(p_fc_snv)

############################
## 12. Optional: save plots
############################
ggsave("repeat1000_fs_estimate_hist.png", p_est_fs, width = 7, height = 5, dpi = 300)
ggsave("repeat1000_snv_estimate_hist.png", p_est_snv, width = 7, height = 5, dpi = 300)
ggsave("repeat1000_fs_foldchange_hist.png", p_fc_fs, width = 7, height = 5, dpi = 300)
ggsave("repeat1000_snv_foldchange_hist.png", p_fc_snv, width = 7, height = 5, dpi = 300)