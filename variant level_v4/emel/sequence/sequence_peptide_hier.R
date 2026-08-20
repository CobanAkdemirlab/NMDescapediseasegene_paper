library(tidyverse)
library(brms)
library(bayesplot)

# ══════════════════════════════════════════════════════════════════════════════
# 1. 查看总体均值分布
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
# 2. 查看每个 uniprot_id 内部的观测数与分布
# ══════════════════════════════════════════════════════════════════════════════
within_patient_summary = NMD_region_NMDPos_peptide_props %>%
  group_by(uniprot_id, source_folder) %>%
  summarise(
    n_obs = n(),
    .groups = "drop"
  )

# 观测数分布
obs_count_dist = within_patient_summary %>%
  count(n_obs) %>%
  arrange(n_obs)
print(obs_count_dist)

ggplot(within_patient_summary, aes(x = n_obs)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  facet_wrap(~ source_folder) +
  labs(title = "每个 uniprot_id 的观测数分布",
       x = "观测数", y = "频数") +
  theme_bw()

# ══════════════════════════════════════════════════════════════════════════════
# 3. 自动计算先验 scale
# ══════════════════════════════════════════════════════════════════════════════
get_prior_scale = function(var, data) {
  s = sd(data[[var]], na.rm = TRUE)
  round(s * 2, 2)
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. 批量拟合贝叶斯层次模型
# ══════════════════════════════════════════════════════════════════════════════
model_list = list()

for (var in d_vars) {
  cat("\n══════════════════════════════════════\n")
  cat("正在拟合变量：", var, "\n")
  cat("══════════════════════════════════════\n")
  
  # 准备数据
  model_data = NMD_region_NMDPos_peptide_props %>%
    select(uniprot_id, source_folder, all_of(var)) %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(
      source_folder = factor(source_folder),
      uniprot_id    = factor(uniprot_id)
    )
  
  # 自动先验
  prior_scale = get_prior_scale(var, model_data)
  cat("自动先验 scale =", prior_scale, "\n")
  
  # 拟合模型
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
  cat("✅ 完成：", var, "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. 提取两个关键对比的结果
# ══════════════════════════════════════════════════════════════════════════════
results_key = map_dfr(d_vars, function(var) {
  
  # 对比1：fs_disease vs fs_control（直接从 fixef 读取）
  fe = fixef(model_list[[var]])
  fs_row = data.frame(
    variable   = var,
    comparison = "fs_disease vs fs_control",
    Estimate   = fe["source_folderfs_disease", "Estimate"],
    Est.Error  = fe["source_folderfs_disease", "Est.Error"],
    Q2.5       = fe["source_folderfs_disease", "Q2.5"],
    Q97.5      = fe["source_folderfs_disease", "Q97.5"]
  )
  
  # 对比2：snv_disease vs snv_control（hypothesis 做差）
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
    significant = ifelse(Q2.5 > 0 | Q97.5 < 0, "✅ 显著", "❌ 不显著")
  )

print(results_key)

# ══════════════════════════════════════════════════════════════════════════════
# 6. 可视化
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

# ── 6c. 后验预测检验（每个模型）──────────────────────────────────────────────
for (var in d_vars) {
  print(
    pp_check(model_list[[var]], ndraws = 100) +
      labs(title = paste("后验预测检验：", var))
  )
}

# ── 6d. 收敛诊断（Rhat 总览）────────────────────────────────────────────────
rhat_summary = map_dfr(d_vars, function(var) {
  rhat_vals = rhat(model_list[[var]])
  data.frame(
    variable = var,
    max_rhat = max(rhat_vals, na.rm = TRUE),
    converged = ifelse(max(rhat_vals, na.rm = TRUE) < 1.01, "✅ 收敛", "⚠️ 未收敛")
  )
})
print(rhat_summary)

# ══════════════════════════════════════════════════════════════════════════════
# 7. 保存结果
# ══════════════════════════════════════════════════════════════════════════════
write_csv(results_key,    "hierarchical_model_key_comparisons.csv")
write_csv(rhat_summary,   "hierarchical_model_convergence.csv")

# 保存模型对象（方便以后重新加载，不用重跑）
saveRDS(model_list, "hierarchical_model_list.rds")
# 重新加载：model_list = readRDS("hierarchical_model_list.rds")

cat("\n✅ 全部完成！\n")
cat("── 输出文件 ──\n")
cat("  hierarchical_model_key_comparisons.csv  ← 主要结果\n")
cat("  hierarchical_model_convergence.csv      ← 收敛诊断\n")
cat("  hierarchical_model_list.rds             ← 模型对象\n")
cat("── 访问单个模型 ──\n")
cat("  model_list[['d_boman']]                 ← 示例\n")