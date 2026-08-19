library(tidyverse)
library(lme4)
library(lmerTest)
library(gtsummary)
library(gt)
library(ggpubr)
library(patchwork)

# ══════════════════════════════════════════════════════════════
# 0. 参数
# ══════════════════════════════════════════════════════════════

MIN_PEPTIDE_LEN <- 10   # 最小肽段长度（aa）；短肽的理化指标不可靠

# ══════════════════════════════════════════════════════════════
# 1. 读入
# ══════════════════════════════════════════════════════════════

base_dir <- "~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/peptides"

div_raw <- read_csv(file.path(base_dir, "NMD_region_DivergentPos_peptide_props.csv"),
                    show_col_types = FALSE)
nmd_raw <- read_csv(file.path(base_dir, "NMD_region_NMDPos_peptide_props.csv"),
                    show_col_types = FALSE)


# ══════════════════════════════════════════════════════════════
# 2. file_name -> key，并检查与 all_variants 的匹配率
# ══════════════════════════════════════════════════════════════
# file_name: "1_196745951_C_T"  ->  key: "1:196745951|C|T"
# 注意染色体位置可能有前导零（如 17_067975880_GT_G），as.numeric 会去掉。
# 若匹配率很低，说明转换规则需要调整 —— 先看下面的诊断输出。

make_key <- function(fn) {
  parts <- str_split_fixed(fn, "_", 4)
  paste0(parts[, 1], ":", as.numeric(parts[, 2]), "|", parts[, 3], "|", parts[, 4])
}

check_key_match <- function(d, nm) {
  d2 <- d %>% mutate(key = make_key(file_name))
  dis <- d2 %>% filter(str_detect(source_folder, "disease"))
  tibble(
    file            = nm,
    n_disease_rows  = nrow(dis),
    n_matched       = sum(dis$key %in% all_variants$key),
    pct_matched     = round(100 * mean(dis$key %in% all_variants$key), 1)
  )
}

cat("\n=== key 匹配诊断（疾病组行）===\n")
print(as.data.frame(bind_rows(
  check_key_match(div_raw, "DivergentPos"),
  check_key_match(nmd_raw, "NMDPos")
)))
cat("pct_matched 应与结构分析中的疾病组保留比例接近；过低说明 key 转换有误。\n\n")


# ══════════════════════════════════════════════════════════════
# 3. 清洗：分组、队列筛选、去除报错行、长度阈值
# ══════════════════════════════════════════════════════════════

prep_region <- function(d, region_name, len_var, len_wt) {
  
  d <- d %>%
    mutate(
      key   = make_key(file_name),
      group = case_when(
        source_folder == "snv_disease" ~ "SNV Disease",
        source_folder == "snv_control" ~ "SNV Control",
        source_folder == "fs_disease"  ~ "FS Disease",
        source_folder == "fs_control"  ~ "FS Control",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group))
  
  # 队列筛选（与结构分析一致）
  d <- bind_rows(
    d %>% filter(key %in% all_variants$key),
    d %>% filter(source_folder %in% c("snv_control", "fs_control"))
  )
  
  n0 <- nrow(d)
  
  # 去除计算失败的行
  d <- d %>% filter(is.na(var_error))
  n1 <- nrow(d)
  
  # 肽段长度：Divergent 文件有长度列，NMD 文件需从序列计算
  d <- d %>%
    mutate(
      var_len = if (!is.null(len_var) && len_var %in% names(.)) .data[[len_var]]
      else nchar(.data[["var_sequence_after_NMD_old"]]),
      wt_len  = if (!is.null(len_wt) && len_wt %in% names(.)) .data[[len_wt]]
      else nchar(.data[["wt_sequence_after_NMD_old"]])
    )
  
  d <- d %>% filter(!is.na(var_len), var_len >= MIN_PEPTIDE_LEN)
  n2 <- nrow(d)
  
  cat(sprintf("[%s] 筛选后 %d 行 -> 去除 var_error %d 行 -> 长度 >= %d aa 后 %d 行\n",
              region_name, n0, n0 - n1, MIN_PEPTIDE_LEN, n2))
  
  d %>%
    mutate(
      region = region_name,
      group  = factor(group, levels = c("SNV Disease", "SNV Control",
                                        "FS Disease", "FS Control"))
    )
}

cat("=== 清洗记录 ===\n")
df_div <- prep_region(div_raw, "Divergent", "var_NMD_length", "wt_NMD_length")
df_nmd <- prep_region(nmd_raw, "NMD", NULL, NULL)
cat("\n")


# ══════════════════════════════════════════════════════════════
# 4. 变量与标签
# ══════════════════════════════════════════════════════════════
# length_proxy = TRUE 的变量与肽段长度近似成正比，
# 组间差异很可能只是截断程度的重述，解读时需谨慎。

var_info <- tribble(
  ~var,                 ~label,                ~length_proxy,
  "aliphatic_index",    "Aliphatic Index",      FALSE,
  "boman",              "Boman Index",          FALSE,
  "charge",             "Charge",               TRUE,
  "entropy",            "Entropy",              TRUE,
  "instability_index",  "Instability Index",    FALSE,
  "isoelectric_point",  "Isoelectric Point",    FALSE,
  "molecular_weight",   "Molecular Weight",     TRUE,
  "mz",                 "m/z",                  TRUE,
  "longest_run",        "Longest Run",          TRUE,
  "max_frequency",      "Max Frequency",        TRUE,
  "energy_cost",        "Energy Cost",          TRUE,
  "nutrient_cost",      "Nutrient Cost",        TRUE,
  "mass_shift",         "Mass Shift",           TRUE,
  "hydrophobicity",     "Hydrophobicity",       FALSE,
  "hydrophobic_moment", "Hydrophobic Moment",   FALSE
)

numeric_vars <- var_info$var


# ══════════════════════════════════════════════════════════════
# 5. 原始数据：描述性表格（中位数 IQR）
# ══════════════════════════════════════════════════════════════
# 展示 VAR、WT、Δ 三套值。WT 为变异体水平展示（描述队列现状），
# 推断检验见第 6 节。

make_desc_table <- function(d, region_name) {
  
  d_tbl <- d %>%
    mutate(across(all_of(paste0("var_", numeric_vars)), identity))
  
  for (v in numeric_vars) {
    d_tbl[[paste0("delta_", v)]] <- d_tbl[[paste0("var_", v)]] - d_tbl[[paste0("wt_", v)]]
  }
  
  lbls <- c(
    setNames(paste0("VAR ", var_info$label), paste0("var_", var_info$var)),
    setNames(paste0("WT ",  var_info$label), paste0("wt_",  var_info$var)),
    setNames(paste0("Δ ",   var_info$label), paste0("delta_", var_info$var)),
    var_len = "VAR Peptide Length (aa)",
    wt_len  = "WT Peptide Length (aa)"
  )
  
  grp <- c(
    setNames(rep("Peptide Length", 2), c("var_len", "wt_len")),
    setNames(rep("VAR Properties", length(numeric_vars)), paste0("var_", numeric_vars)),
    setNames(rep("WT Properties",  length(numeric_vars)), paste0("wt_",  numeric_vars)),
    setNames(rep("Δ (VAR − WT)",   length(numeric_vars)), paste0("delta_", numeric_vars))
  )
  
  sel <- c("group", "var_len", "wt_len",
           paste0("var_", numeric_vars),
           paste0("wt_",  numeric_vars),
           paste0("delta_", numeric_vars))
  
  d_tbl %>%
    select(all_of(sel)) %>%
    tbl_summary(
      by        = group,
      statistic = all_continuous() ~ "{median} ({p25}, {p75})",
      digits    = all_continuous() ~ 2,
      missing   = "no",
      label     = as.list(setNames(
        lapply(names(lbls), function(x) lbls[[x]]), names(lbls)
      ))
    ) %>%
    add_overall() %>%
    modify_table_body(~ .x %>% mutate(groupname_col = grp[variable])) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_caption(paste0("**Peptide properties — ", region_name,
                          " region (descriptive, median [IQR])**")) %>%
    bold_labels()
}

desc_div <- make_desc_table(df_div, "Divergent")
desc_nmd <- make_desc_table(df_nmd, "NMD")

desc_div
desc_nmd

desc_div %>% as_gt() %>% gtsave("Peptide_descriptive_Divergent.html")
desc_nmd %>% as_gt() %>% gtsave("Peptide_descriptive_NMD.html")
desc_div %>% as_flex_table() %>% flextable::save_as_docx(
  path = "Peptide_descriptive_Divergent.docx")
desc_nmd %>% as_flex_table() %>% flextable::save_as_docx(
  path = "Peptide_descriptive_NMD.docx")


# ══════════════════════════════════════════════════════════════
# 6. 混合模型：value ~ group + (1 | uniprot_id)
# ══════════════════════════════════════════════════════════════
# 只对 VAR 建模。原因：WT 在基因内恒定，被随机截距完全吸收，
# 因此 Δ 模型与 VAR 模型的组效应估计和 p 值在数学上完全相同，
# 同时报告两套等于重复计数并加倍 FDR 惩罚。

fit_lmm <- function(data, outcome, g_dis, g_con, family_label) {
  
  d <- data %>%
    select(uniprot_id, group, value = all_of(outcome)) %>%
    filter(group %in% c(g_dis, g_con), !is.na(value)) %>%
    mutate(group = factor(group, levels = c(g_con, g_dis)))
  
  m <- tryCatch(
    lmerTest::lmer(value ~ group + (1 | uniprot_id), data = d),
    error = function(e) NULL
  )
  
  if (is.null(m)) {
    return(tibble(family = family_label, outcome = outcome,
                  estimate = NA_real_, se = NA_real_, ci_lo = NA_real_,
                  ci_hi = NA_real_, p_raw = NA_real_, sd_value = NA_real_,
                  n_obs = nrow(d), n_genes = n_distinct(d$uniprot_id)))
  }
  
  cf  <- coef(summary(m))
  row <- grep("^group", rownames(cf))[1]
  est <- cf[row, "Estimate"]; se <- cf[row, "Std. Error"]
  
  tibble(
    family = family_label, outcome = outcome,
    estimate = est, se = se,
    ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
    p_raw = cf[row, "Pr(>|t|)"],
    sd_value = sd(d$value),
    n_obs = nrow(d), n_genes = n_distinct(d$uniprot_id)
  )
}

run_region_lmm <- function(d, region_name) {
  
  outcomes <- c("var_len", paste0("var_", numeric_vars))
  
  bind_rows(
    map_dfr(outcomes, ~ fit_lmm(d, .x, "SNV Disease", "SNV Control", "SNV")),
    map_dfr(outcomes, ~ fit_lmm(d, .x, "FS Disease",  "FS Control",  "FS"))
  ) %>%
    mutate(region = region_name)
}

lmm_all <- bind_rows(
  run_region_lmm(df_div, "Divergent"),
  run_region_lmm(df_nmd, "NMD")
) %>%
  group_by(family) %>%                       # 每个比较内校正（跨两个 region）
  mutate(q_fdr = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    var    = str_remove(outcome, "^var_"),
    label  = if_else(var == "len", "Peptide Length (aa)",
                     var_info$label[match(var, var_info$var)]),
    length_proxy = if_else(var == "len", TRUE,
                           var_info$length_proxy[match(var, var_info$var)]),
    label_disp = if_else(length_proxy, paste0(label, ""), label),
    family = factor(family, levels = c("SNV", "FS")),
    region = factor(region, levels = c("Divergent", "NMD")),
    std_est = estimate / sd_value,
    std_lo  = ci_lo    / sd_value,
    std_hi  = ci_hi    / sd_value,
    sig = case_when(
      is.na(q_fdr)  ~ "",
      q_fdr < 0.001 ~ "***",
      q_fdr < 0.01  ~ "**",
      q_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

write_csv(lmm_all, "Peptide_LMM_results.csv")

cat("=== 混合模型结果 ===\n")
lmm_all %>%
  select(region, family, label, estimate, p_raw, q_fdr, sig, n_obs, n_genes) %>%
  mutate(across(c(estimate, p_raw, q_fdr), ~ signif(.x, 3))) %>%
  arrange(region, family, q_fdr) %>%
  print(n = Inf)


# ══════════════════════════════════════════════════════════════
# 7. 混合模型表格：SNV 与 FS 并排
# ══════════════════════════════════════════════════════════════

fmt_p <- function(x) {
  ifelse(is.na(x), "—",
         ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                formatC(x, format = "f", digits = 3)))
}

make_lmm_table <- function(res, region_name) {
  
  n_info <- res %>%
    filter(region == region_name) %>%
    group_by(family) %>%
    summarise(n_obs = max(n_obs), n_genes = max(n_genes), .groups = "drop")
  
  res %>%
    filter(region == region_name) %>%
    mutate(ci = sprintf("(%s, %s)", signif(ci_lo, 3), signif(ci_hi, 3))) %>%
    select(label_disp, family, estimate, ci, p_raw, q_fdr, sig) %>%
    pivot_wider(names_from = family,
                values_from = c(estimate, ci, p_raw, q_fdr, sig),
                names_sep = "_") %>%
    select(label_disp,
           estimate_SNV, ci_SNV, p_raw_SNV, q_fdr_SNV, sig_SNV,
           estimate_FS,  ci_FS,  p_raw_FS,  q_fdr_FS,  sig_FS) %>%
    gt() %>%
    fmt_number(columns = c(estimate_SNV, estimate_FS), n_sigfig = 3) %>%
    fmt(columns = c(p_raw_SNV, q_fdr_SNV, p_raw_FS, q_fdr_FS), fns = fmt_p) %>%
    cols_label(
      label_disp   = "Variable",
      estimate_SNV = "Estimate", ci_SNV = "95% CI",
      p_raw_SNV = "p", q_fdr_SNV = "q (FDR)", sig_SNV = "",
      estimate_FS  = "Estimate", ci_FS  = "95% CI",
      p_raw_FS  = "p", q_fdr_FS  = "q (FDR)", sig_FS  = ""
    ) %>%
    tab_spanner(
      label = md(sprintf("**SNV** (n = %s, %s genes)",
                         n_info$n_obs[n_info$family == "SNV"],
                         n_info$n_genes[n_info$family == "SNV"])),
      columns = c(estimate_SNV, ci_SNV, p_raw_SNV, q_fdr_SNV, sig_SNV)
    ) %>%
    tab_spanner(
      label = md(sprintf("**FS** (n = %s, %s genes)",
                         n_info$n_obs[n_info$family == "FS"],
                         n_info$n_genes[n_info$family == "FS"])),
      columns = c(estimate_FS, ci_FS, p_raw_FS, q_fdr_FS, sig_FS)
    ) %>%
    tab_header(
      title = paste0("Peptide properties — ", region_name, " region"),
      subtitle = "Mixed-effects models; Estimate = Disease − Control, original units"
    ) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_body(
                columns = c(estimate_SNV, ci_SNV, p_raw_SNV, q_fdr_SNV, sig_SNV),
                rows = sig_SNV != "")) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_body(
                columns = c(estimate_FS, ci_FS, p_raw_FS, q_fdr_FS, sig_FS),
                rows = sig_FS != "")) %>%
    tab_style(style = cell_borders(sides = "left", color = "gray70",
                                   weight = px(1.5)),
              locations = cells_body(columns = estimate_FS)) %>%
    tab_source_note(md(paste0(
      "Linear mixed-effects model: value ~ group + (1 | gene). Random intercept ",
      "for gene (UniProt ID) accounts for multiple variants per gene. ",
      "P-values: Satterthwaite approximation (lmerTest); q: Benjamini-Hochberg ",
      "FDR applied within each comparison across both regions. ",
      "\\* q<0.05, \\*\\* q<0.01, \\*\\*\\* q<0.001. ",
      "Peptides shorter than ", MIN_PEPTIDE_LEN, " aa were excluded. ",
      "† scales with peptide length — differences may reflect truncation extent ",
      "rather than composition. WT properties are constant within gene and were ",
      "not modeled; Δ (VAR−WT) models are mathematically identical to VAR models ",
      "under this random-intercept specification."
    )))
}

tbl_lmm_div <- make_lmm_table(lmm_all, "Divergent")
tbl_lmm_nmd <- make_lmm_table(lmm_all, "NMD")

tbl_lmm_div
tbl_lmm_nmd

gtsave(tbl_lmm_div, "Peptide_LMM_table_Divergent.html")
gtsave(tbl_lmm_nmd, "Peptide_LMM_table_NMD.html")
try(gtsave(tbl_lmm_div, "Peptide_LMM_table_Divergent.docx"), silent = TRUE)
try(gtsave(tbl_lmm_nmd, "Peptide_LMM_table_NMD.docx"),       silent = TRUE)


# ══════════════════════════════════════════════════════════════
# 8. 森林图：SNV 与 FS 同图，两个 region 分面
# ══════════════════════════════════════════════════════════════

plot_dat <- lmm_all %>%
  filter(
    !is.na(estimate),
  #  region != "NMD"
  ) %>%
  mutate(
    label_disp = fct_rev(
      factor(label_disp, levels = unique(label_disp))
    )
  )
p_forest <- ggplot(plot_dat, aes(x = std_est, y = label_disp, color = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = std_lo, xmax = std_hi), height = 0,
                 linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = q_fdr < 0.05), size = 2.2,
             position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig), size = 3, hjust = -0.5, vjust = 0.8,
            position = position_dodge(width = 0.6), show.legend = FALSE) +
  facet_wrap(~ region, nrow = 1) +
  scale_color_manual(values = c(SNV = "#2C7FB8", FS = "#C0392B"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  labs(
    title    = "Peptide properties: Disease vs Control",
    subtitle = "Mixed-effects estimates standardized by outcome SD; bars = 95% CI",
    x = "Standardized estimate (Disease − Control)", y = NULL,
    caption = paste0(
      "Linear mixed-effects model: value ~ group + (1 | gene). Filled points BH-adjusted p < 0.05; ",
      "* p<0.05  ** p<0.01  *** p<0.001."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "gray95", color = NA),
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "top",
    plot.title         = element_text(face = "bold"),
    plot.caption       = element_text(hjust = 0, size = 8, color = "gray40")
  )

print(p_forest)

ggsave("Peptide_LMM_forest.pdf", p_forest, width = 8, height = 7,
       dpi = 300, device = cairo_pdf)
ggsave("Peptide_LMM_forest.png", p_forest, width = 8, height = 7, dpi = 300)

