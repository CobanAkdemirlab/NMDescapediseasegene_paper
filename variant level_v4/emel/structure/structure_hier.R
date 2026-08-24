library(tidyverse)
library(brms)
library(bayesplot)
library(posterior)

# ═════════════════════════════════════════════════════════════════════
# 1. Select secondary structure variables to analyze
# ═════════════════════════════════════════════════════════════════════

sec_vars = c(
  "VAR_sec_struc_NMDPos_helix_percent",
  "VAR_sec_struc_NMDPos_beta_percent",
  "VAR_sec_struc_NMDPos_coil_percent",
  "WT_sec_struc_NMDPos_helix_percent",
  "WT_sec_struc_NMDPos_beta_percent",
  "WT_sec_struc_NMDPos_coil_percent",
  "d_sec_struc_NMDPos_helix_percent",
  "d_sec_struc_NMDPos_beta_percent",
  "d_sec_struc_NMDPos_coil_percent"
)

# ═════════════════════════════════════════════════════════════════════
# 2. Set source_folder reference level
# ═════════════════════════════════════════════════════════════════════

VAR_WT_structural_results = VAR_WT_structural_results %>%
  mutate(
    source_folder = factor(
      source_folder,
      levels = c(
        "fs_control",
        "fs_disease",
        "snv_control",
        "snv_disease"
      )
    ),
    uniprot_id = factor(uniprot_id)
  )

# ═════════════════════════════════════════════════════════════════════
# 3. Beta regression for VAR/WT percent
#    Gaussian for d_ difference
# ═════════════════════════════════════════════════════════════════════

to_beta_01 = function(x) {
  # Convert 0-100 percent to 0-1
  if (max(x, na.rm = TRUE) > 1) {
    x = x / 100
  }
  
  # Beta family cannot have 0 or 1, so compress slightly
  n = sum(!is.na(x))
  x_beta = (x * (n - 1) + 0.5) / n
  
  return(x_beta)
}

get_prior_scale = function(x) {
  s = sd(x, na.rm = TRUE)
  ifelse(is.na(s) | s == 0, 1, round(s * 2, 2))
}

# ═════════════════════════════════════════════════════════════════════
# 4. Batch fit hierarchical models
# ═════════════════════════════════════════════════════════════════════

model_list = list()

for (var in sec_vars) {
  
  cat("\n══════════════════════════════════════\n")
  cat("Fitting variable:", var, "\n")
  cat("══════════════════════════════════════\n")
  
  model_data = VAR_WT_structural_results %>%
    select(uniprot_id, source_folder, all_of(var)) %>%
    filter(!is.na(.data[[var]]))
  
  # d_ variables can be negative, so use Gaussian
  if (str_detect(var, "^d_")) {
    
    model_data = model_data %>%
      rename(y = all_of(var))
    
    prior_scale = get_prior_scale(model_data$y)
    
    model = brm(
      formula = y ~ source_folder + (1 | uniprot_id),
      data    = model_data,
      family  = gaussian(),
      prior   = c(
        prior(normal(0, 10), class = "Intercept"),  # mean d_var
        prior(normal(0, prior_scale), class = "b"), # group level effect
        prior(exponential(1), class = "sd"), # gene level variance
        prior(exponential(1), class = "sigma")# variant level variance
      ),
      chains  = 4,
      iter    = 4000,
      warmup  = 1000,
      cores   = 4,
      seed    = 42,
      control = list(adapt_delta = 0.95),
      silent  = 2
    )
    
  } else {
    
    # Use Beta for VAR/WT percent
    model_data = model_data %>%
      mutate(y = to_beta_01(.data[[var]]))
    
    model = brm(
      formula = y ~ source_folder + (1 | uniprot_id),
      data    = model_data,
      family  = Beta(link = "logit"),
      prior   = c(
        prior(normal(0, 2), class = "Intercept"),
        prior(normal(0, 1), class = "b"),
        prior(exponential(1), class = "sd"),
        prior(exponential(1), class = "phi")
      ),
      chains  = 4,
      iter    = 4000,
      warmup  = 1000,
      cores   = 4,
      seed    = 42,
      control = list(adapt_delta = 0.95),
      silent  = 2
    )
  }
  
  model_list[[var]] = model
  
  cat("Done:", var, "\n")
}

# ═════════════════════════════════════════════════════════════════════
# 5. Extract two main comparisons
#    fs_disease vs fs_control
#    snv_disease vs snv_control
# ═════════════════════════════════════════════════════════════════════

results_key = map_dfr(sec_vars, function(var) {
  
  model = model_list[[var]]
  
  # fs_disease vs fs_control
  fs_hyp = hypothesis(
    model,
    "source_folderfs_disease = 0"
  )$hypothesis
  
  # snv_disease vs snv_control
  snv_hyp = hypothesis(
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
    significant = ifelse(Q2.5 > 0 | Q97.5 < 0, "Significant", "Not significant")
  )

print(results_key)

# ═════════════════════════════════════════════════════════════════════
# 6. Forest plot
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
    title = "Hierarchical model: secondary structure differences",
    subtitle = "Point = posterior mean; line = 95% credible interval",
    x = "Effect size",
    y = "Variable",
    color = ""
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ═════════════════════════════════════════════════════════════════════
# 7. Rhat convergence check
# ═════════════════════════════════════════════════════════════════════

rhat_summary = map_dfr(sec_vars, function(var) {
  
  rhat_vals = rhat(model_list[[var]])
  
  data.frame(
    variable  = var,
    max_rhat  = max(rhat_vals, na.rm = TRUE),
    converged = ifelse(max(rhat_vals, na.rm = TRUE) < 1.01, "converged", "not converged")
  )
})

print(rhat_summary)

# ═════════════════════════════════════════════════════════════════════
# 8. Save results
# ═════════════════════════════════════════════════════════════════════

write_csv(results_key, "secondary_structure_hierarchical_results.csv")
write_csv(rhat_summary, "secondary_structure_hierarchical_convergence.csv")
saveRDS(model_list, "secondary_structure_hierarchical_model_list.rds")