library(tidyverse)
library(brms)
library(bayesplot)

# ══════════════════════════════════════════════════════════════════════════════
# 1. View overall mean distribution
# ══════════════════════════════════════════════════════════════════════════════
mean_dis = NMD_region_NMDPos_peptide_props %>%
  group_by(uniprot_id, source_folder) %>%
  summarise(across(all_of(d_vars), mean, na.rm = TRUE), .groups = "drop")

for (var in d_vars) {
  fs_vals = mean_dis %>%
    filter(source_folder == "fs_control") %>%
    pull(!!sym(var))
  
  hist(fs_vals, breaks = 30,
       main = paste("Distribution of mean", var, "for fs_control"),
       xlab = paste("mean", var))
  
  cat("\n── Shapiro test for", var, "──\n")
  print(shapiro.test(fs_vals))
  cat("── Summary ──\n")
  print(summary(fs_vals))
  cat("SD:", sd(fs_vals, na.rm = TRUE), "\n\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. View observation counts and distribution per uniprot_id
# ══════════════════════════════════════════════════════════════════════════════
within_patient_summary = NMD_region_NMDPos_peptide_props %>%
  group_by(uniprot_id, source_folder) %>%
  summarise(
    n_obs = n(),
    .groups = "drop"
  )

# Observation count distribution
obs_count_dist = within_patient_summary %>%
  count(n_obs) %>%
  arrange(n_obs)
print(obs_count_dist)

ggplot(within_patient_summary, aes(x = n_obs)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  facet_wrap(~ source_folder) +
  labs(title = "Distribution of observation counts per uniprot_id",
       x = "Observation count", y = "Frequency") +
  theme_bw()

# ══════════════════════════════════════════════════════════════════════════════
# 3. Automatically compute prior scale
# ══════════════════════════════════════════════════════════════════════════════
get_prior_scale = function(var, data) {
  s = sd(data[[var]], na.rm = TRUE)
  round(s * 2, 2)
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Batch fit Bayesian hierarchical models
# ══════════════════════════════════════════════════════════════════════════════
model_list = list()

for (var in d_vars) {
  cat("\n══════════════════════════════════════\n")
  cat("Fitting variable:", var, "\n")
  cat("══════════════════════════════════════\n")
  
  # Prepare data
  model_data = NMD_region_NMDPos_peptide_props %>%
    select(uniprot_id, source_folder, all_of(var)) %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(
      source_folder = factor(source_folder),
      uniprot_id    = factor(uniprot_id)
    )
  
  # Automatic prior
  prior_scale = get_prior_scale(var, model_data)
  cat("Automatic prior scale =", prior_scale, "\n")
  
  # Fit model
  model = brm(
    formula = as.formula(paste(var, "~ source_folder + (1 | uniprot_id)")),
    data    = model_data,
    family  = gaussian(),
    prior   = c(
      prior_string(paste0("normal(0, ", prior_scale * 2, ")"), class = "Intercept"),
      prior_string(paste0("normal(0, ", prior_scale, ")"),     class = "b"),
      prior_string(paste0("normal(0, ", prior_scale, ")"),     class = "sd"),
      prior_string(paste0("normal(0, ", prior_scale, ")"),     class = "sigma")
    ),
    chains  = 4,
    iter    = 4000,
    warmup  = 1000,
    cores   = 4,
    seed    = 42,
    control = list(adapt_delta = 0.95),
    silent  = 2
  )
  
  model_list[[var]] = model
  cat("✅ Done:", var, "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Extract results for two key comparisons
# ══════════════════════════════════════════════════════════════════════════════
results_key = map_dfr(d_vars, function(var) {
  
  # Comparison 1: fs_disease vs fs_control (read directly from fixef)
  fe = fixef(model_list[[var]])
  fs_row = data.frame(
    variable   = var,
    comparison = "fs_disease vs fs_control",
    Estimate   = fe["source_folderfs_disease", "Estimate"],
    Est.Error  = fe["source_folderfs_disease", "Est.Error"],
    Q2.5       = fe["source_folderfs_disease", "Q2.5"],
    Q97.5      = fe["source_folderfs_disease", "Q97.5"]
  )
  
  # Comparison 2: snv_disease vs snv_control (difference via hypothesis)
  snv_hyp = hypothesis(
    model_list[[var]],
    "source_foldersnv_disease - source_foldersnv_control = 0"
  )$hypothesis
  
  snv_row = data.frame(
    variable   = var,
    comparison = "snv_disease vs snv_control",
    Estimate   = snv_hyp$Estimate,
    Est.Error  = snv_hyp$Est.Error,
    Q2.5       = snv_hyp$CI.Lower,
    Q97.5      = snv_hyp$CI.Upper
  )
  
  bind_rows(fs_row, snv_row)
}) %>%
  mutate(
    significant = ifelse(Q2.5 > 0 | Q97.5 < 0, "✅ Significant", "❌ Not significant")
  )

print(results_key)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Visualization
# ══════════════════════════════════════════════════════════════════════════════

# ── 6a. Forest plot ───────────────────────────────────────────────────────────
ggplot(results_key, aes(x = Estimate, y = variable, color = significant)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ comparison, scales = "free_x") +
  scale_color_manual(values = c("Significant" = "coral", "Not significant" = "steelblue")) +
  labs(
    title    = "Effect size of disease vs control per variable",
    subtitle = "Point = posterior mean, line = 95% credible interval",
    x        = "Effect size (disease - control)",
    y        = "Variable",
    color    = ""
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 6b. Heatmap (z-score normalized for cross-variable comparison) ─────────────
results_key %>%
  mutate(
    z_score = Estimate / Est.Error,
    label   = paste0(round(Estimate, 2), ifelse(significant == "Significant", "*", ""))
  ) %>%
  ggplot(aes(x = comparison, y = variable, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "coral",
    midpoint = 0,
    name     = "Effect / SE"
  ) +
  labs(
    title    = "Heatmap of group differences across variables",
    subtitle = "* = 95% credible interval excludes 0 (significant)\nReference: respective control group",
    x        = "",
    y        = "Variable"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ── 6c. Posterior predictive check (per model) ──────────────────────────────────────────────
for (var in d_vars) {
  print(
    pp_check(model_list[[var]], ndraws = 100) +
      labs(title = paste("Posterior predictive check:", var))
  )
}

# ── 6d. Convergence diagnostics (Rhat overview) ────────────────────────────────────────────────
rhat_summary = map_dfr(d_vars, function(var) {
  rhat_vals = rhat(model_list[[var]])
  data.frame(
    variable = var,
    max_rhat = max(rhat_vals, na.rm = TRUE),
    converged = ifelse(max(rhat_vals, na.rm = TRUE) < 1.01, "✅ Converged", "⚠️ Not converged")
  )
})
print(rhat_summary)

# ══════════════════════════════════════════════════════════════════════════════
# 7. Save results
# ══════════════════════════════════════════════════════════════════════════════
write_csv(results_key,    "hierarchical_model_key_comparisons.csv")
write_csv(rhat_summary,   "hierarchical_model_convergence.csv")

# Save model object for later reload without rerunning
saveRDS(model_list, "hierarchical_model_list.rds")

cat("\n✅ All done!\n")
cat("── Output files ──\n")
cat("  hierarchical_model_key_comparisons.csv  ← main results\n")
cat("  hierarchical_model_convergence.csv      ← convergence diagnostics\n")
cat("  hierarchical_model_list.rds             ← model object\n")
cat("── Access individual model ──\n")
cat("  model_list[['d_boman']]                 ← example\n")