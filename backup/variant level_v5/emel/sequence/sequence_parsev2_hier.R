library(tidyverse)
library(brms)
library(bayesplot)
library(posterior)

# ═════════════════════════════════════════════════════════════════════
# 1. Read data
# ═════════════════════════════════════════════════════════════════════

nmd_length_NMDPos_PARSE_v2 <- read_csv(
  "~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/parse_v2/nmd_length_NMDPos_PARSE_v2.csv",
  show_col_types = FALSE
)

# ═════════════════════════════════════════════════════════════════════
# 2. Tidy data
# ═════════════════════════════════════════════════════════════════════

parse_dat <- nmd_length_NMDPos_PARSE_v2 %>%
  rename(
    var_classifier_distance_P = `Σ classifier distance P...4`,
    var_longest_PS_IDR        = `length of longest PS IDR...5`,
    wt_classifier_distance_P  = `Σ classifier distance P...9`,
    wt_longest_PS_IDR         = `length of longest PS IDR...10`
  ) %>%
  mutate(
    source_folder = case_when(
      str_detect(id, "^fs_disease")  ~ "fs_disease",
      str_detect(id, "^fs_control")  ~ "fs_control",
      str_detect(id, "^snv_disease") ~ "snv_disease",
      str_detect(id, "^snv_control") ~ "snv_control",
      TRUE ~ NA_character_
    ),
    source_folder = factor(
      source_folder,
      levels = c(
        "fs_control",
        "fs_disease",
        "snv_control",
        "snv_disease"
      )
    ),
    uniprot = factor(uniprot)
  ) %>%
  filter(!is.na(source_folder), !is.na(uniprot))

# ═════════════════════════════════════════════════════════════════════
# 3. Compute VAR - WT
# ═════════════════════════════════════════════════════════════════════

parse_dat <- parse_dat %>%
  mutate(
    d_NMD_length = var_NMD_length_old - wt_NMD_length_old,
    d_classifier_distance_P = var_classifier_distance_P - wt_classifier_distance_P,
    d_longest_PS_IDR = var_longest_PS_IDR - wt_longest_PS_IDR
  )

d_vars <- c(
  "d_NMD_length",
  "d_classifier_distance_P",
  "d_longest_PS_IDR"
)

# ═════════════════════════════════════════════════════════════════════
# 4. Check counts per group
# ═════════════════════════════════════════════════════════════════════

print(table(parse_dat$source_folder))

parse_dat %>%
  count(source_folder, uniprot) %>%
  count(source_folder, name = "n_proteins") %>%
  print()

# ═════════════════════════════════════════════════════════════════════
# 5. prior scale function
# ═════════════════════════════════════════════════════════════════════

get_prior_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  ifelse(is.na(s) | s == 0, 1, round(s * 2, 2))
}

# ═════════════════════════════════════════════════════════════════════
# 6. Batch fit Bayesian hierarchical Gaussian models
# ═════════════════════════════════════════════════════════════════════

model_list <- list()

for (var in d_vars) {
  
  cat("\n══════════════════════════════════════\n")
  cat("Fitting variable:", var, "\n")
  cat("══════════════════════════════════════\n")
  
  model_data <- parse_dat %>%
    select(uniprot, source_folder, all_of(var)) %>%
    filter(!is.na(.data[[var]])) %>%
    rename(y = all_of(var))
  
  prior_scale <- get_prior_scale(model_data$y)
  cat("prior_scale =", prior_scale, "\n")
  
  prior_use <- c(
    set_prior("normal(0, 10)", class = "Intercept"),
    set_prior(paste0("normal(0, ", prior_scale, ")"), class = "b"),
    set_prior("exponential(1)", class = "sd"),
    set_prior("exponential(1)", class = "sigma")
  )
  
  model <- brm(
    formula = y ~ source_folder + (source_folder | uniprot),
    data    = model_data,
    family  = gaussian(),
    prior   = prior_use,
    chains  = 4,
    iter    = 4000,
    warmup  = 1000,
    cores   = 4,
    seed    = 42,
    control = list(adapt_delta = 0.95),
    silent  = 2
  )
  
  model_list[[var]] <- model
  
  cat("Done:", var, "\n")
}

# ═════════════════════════════════════════════════════════════════════
# 7. Extract two main comparisons
#    fs_disease vs fs_control
#    snv_disease vs snv_control
# ═════════════════════════════════════════════════════════════════════

results_key <- map_dfr(d_vars, function(var) {
  
  model <- model_list[[var]]
  
  fs_hyp <- hypothesis(
    model,
    "source_folderfs_disease = 0"
  )$hypothesis
  
  snv_hyp <- hypothesis(
    model,
    "source_foldersnv_disease - source_foldersnv_control = 0"
  )$hypothesis
  
  bind_rows(
    data.frame(
      variable   = var,
      comparison = "fs_disease vs fs_control",
      Estimate   = fs_hyp$Estimate,
      Est.Error  = fs_hyp$Est.Error,
      Q2.5       = fs_hyp$CI.Lower,
      Q97.5      = fs_hyp$CI.Upper
    ),
    data.frame(
      variable   = var,
      comparison = "snv_disease vs snv_control",
      Estimate   = snv_hyp$Estimate,
      Est.Error  = snv_hyp$Est.Error,
      Q2.5       = snv_hyp$CI.Lower,
      Q97.5      = snv_hyp$CI.Upper
    )
  )
}) %>%
  mutate(
    significant = ifelse(
      Q2.5 > 0 | Q97.5 < 0,
      "Significant",
      "Not significant"
    )
  )

print(results_key)

# ═════════════════════════════════════════════════════════════════════
# 8. Forest plot
# ═════════════════════════════════════════════════════════════════════

ggplot(results_key, aes(x = Estimate, y = variable, color = significant)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ comparison, scales = "free_x") +
  scale_color_manual(
    values = c(
      "Significant" = "coral",
      "Not significant" = "steelblue"
    )
  ) +
  labs(
    title = "Hierarchical model for PARSE NMD-region features",
    subtitle = "Point = posterior mean; line = 95% credible interval",
    x = "Effect size: disease - control",
    y = "Feature",
    color = ""
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ═════════════════════════════════════════════════════════════════════
# 9. Heatmap
# ═════════════════════════════════════════════════════════════════════

results_key %>%
  mutate(
    z_score = Estimate / Est.Error,
    label = paste0(
      round(Estimate, 2),
      ifelse(significant == "Significant", "*", "")
    )
  ) %>%
  ggplot(aes(x = comparison, y = variable, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "coral",
    midpoint = 0,
    name = "Estimate / SE"
  ) +
  labs(
    title = "Heatmap of PARSE NMD-region feature differences",
    subtitle = "* = 95% credible interval excludes 0",
    x = "",
    y = "Feature"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ═════════════════════════════════════════════════════════════════════
# 10. Posterior predictive check
# ═════════════════════════════════════════════════════════════════════

for (var in d_vars) {
  print(
    pp_check(model_list[[var]], ndraws = 100) +
      labs(title = paste("Posterior predictive check:", var))
  )
}

# ═════════════════════════════════════════════════════════════════════
# 11. Rhat convergence check
# ═════════════════════════════════════════════════════════════════════

rhat_summary <- map_dfr(d_vars, function(var) {
  
  rhat_vals <- rhat(model_list[[var]])
  
  data.frame(
    variable  = var,
    max_rhat  = max(rhat_vals, na.rm = TRUE),
    converged = ifelse(
      max(rhat_vals, na.rm = TRUE) < 1.01,
      "converged",
      "not converged"
    )
  )
})

print(rhat_summary)

# ═════════════════════════════════════════════════════════════════════
# 12. Save results
# ═════════════════════════════════════════════════════════════════════

write_csv(results_key, "parse_NMDPos_hierarchical_results.csv")
write_csv(rhat_summary, "parse_NMDPos_hierarchical_convergence.csv")
saveRDS(model_list, "parse_NMDPos_hierarchical_model_list.rds")

cat("\nAll done!\n")
cat("Output files:\n")
cat("parse_NMDPos_hierarchical_results.csv\n")
cat("parse_NMDPos_hierarchical_convergence.csv\n")
cat("parse_NMDPos_hierarchical_model_list.rds\n")