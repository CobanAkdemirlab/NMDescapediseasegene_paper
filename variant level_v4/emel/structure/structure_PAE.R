library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ══════════════════════════════════════════════════════════════
# 1. 队列筛选（与 Table 1 / Table 2 / pLDDT 图一致）
# ══════════════════════════════════════════════════════════════

df_all <- VAR_WT_structural_results %>%
  mutate(
    group = case_when(
      source_folder == "snv_disease" ~ "SNV\nDisease",
      source_folder == "snv_control" ~ "SNV\nControl",
      source_folder == "fs_disease"  ~ "FS\nDisease",
      source_folder == "fs_control"  ~ "FS\nControl",
      TRUE ~ NA_character_
    ),
    group = factor(
      group,
      levels = c("SNV\nDisease", "SNV\nControl", "FS\nDisease", "FS\nControl")
    )
  ) %>%
  filter(!is.na(group))

df_disease <- df_all %>% filter(key %in% all_variants$key)
df_control <- df_all %>% filter(source_folder %in% c("snv_control", "fs_control"))
df_all     <- bind_rows(df_disease, df_control)


# ══════════════════════════════════════════════════════════════
# 2. 数据准备（WT 做蛋白水平去重）
# ══════════════════════════════════════════════════════════════

# WT 的 PAE 对同一蛋白的所有变异体完全相同，不去重会造成伪重复
prepare_pae_data <- function(data, pae_col, panel_type) {
  
  d <- data %>%
    select(group, uniprot_id, PAE = all_of(pae_col)) %>%
    filter(!is.na(PAE))
  
  if (panel_type == "WT") {
    d <- d %>% distinct(group, uniprot_id, .keep_all = TRUE)
  }
  
  d
}


# ══════════════════════════════════════════════════════════════
# 3. p 值与 BH-FDR 校正
# ══════════════════════════════════════════════════════════════

panel_spec <- tibble::tribble(
  ~panel_id,  ~panel_type, ~col,                         ~title,
  "var_full", "VAR",       "VAR_pae_full_mean",          "VAR PAE — Full",
  "var_nmd",  "VAR",       "VAR_pae_NMDPos_mean",        "VAR PAE — NMD",
  "var_div",  "VAR",       "VAR_pae_DivergentPos_mean",  "VAR PAE — Divergent",
  "wt_full",  "WT",        "WT_pae_full_mean",           "WT PAE — Full",
  "wt_nmd",   "WT",        "WT_pae_NMDPos_mean",         "WT PAE — NMD",
  "wt_div",   "WT",        "WT_pae_DivergentPos_mean",   "WT PAE — Divergent"
)

p_raw <- pmap_dfr(
  list(panel_spec$panel_id, panel_spec$col, panel_spec$panel_type),
  function(pid, col, ptype) {
    
    d <- prepare_pae_data(df_all, col, ptype)
    
    test_pair <- function(g1, g2) {
      dd <- d %>%
        filter(group %in% c(g1, g2)) %>%
        mutate(group = droplevels(group))
      tryCatch(
        wilcox.test(PAE ~ group, data = dd)$p.value,
        error = function(e) NA_real_
      )
    }
    
    tibble(
      panel_id = pid,
      family   = c("SNV", "FS"),
      group1   = c("SNV\nDisease", "FS\nDisease"),
      group2   = c("SNV\nControl", "FS\nControl"),
      p_raw    = c(test_pair("SNV\nDisease", "SNV\nControl"),
                   test_pair("FS\nDisease",  "FS\nControl"))
    )
  }
)

# 默认在本图的 12 个检验内校正，SNV 与 FS 家族分开。
# 若要与 Table 1 完全一致，见文件末尾的说明。
p_annot <- p_raw %>%
  group_by(family) %>%
  mutate(q_val = p.adjust(p_raw, method = "BH")) %>%
  ungroup()

q_stars <- function(q) {
  case_when(
    is.na(q)   ~ "NA",
    q < 0.0001 ~ "****",
    q < 0.001  ~ "***",
    q < 0.01   ~ "**",
    q < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

p_annot <- p_annot %>% mutate(label = q_stars(q_val))


# ══════════════════════════════════════════════════════════════
# 4. 单图绘制函数
# ══════════════════════════════════════════════════════════════

fill_colors <- c("SNV\nDisease" = "#E07B7B",
                 "SNV\nControl" = "#7BA7D4",
                 "FS\nDisease"  = "#E8B49A",
                 "FS\nControl"  = "#A8C8E8")

plot_pae <- function(df, title, annot, y_label = "PAE (Å, mean)") {
  
  y_min <- min(df$PAE, na.rm = TRUE)
  y_max <- max(df$PAE, na.rm = TRUE)
  y_rng <- y_max - y_min
  
  # 括号位置 + 上方留白
  annot <- annot %>% mutate(y.position = y_max + 0.08 * y_rng)
  
  # PAE 解读参考线：仅在落入本图数值范围内时才画
  # (<5 Å 相对位置可信；>10 Å 结构域间位置基本不可确定)
  ref_lines <- tibble(yint = c(5, 10),
                      col  = c("#4dac26", "#f1a340")) %>%
    filter(yint >= y_min, yint <= y_max)
  
  ggplot(df, aes(x = group, y = PAE, fill = group)) +
    
    geom_violin(alpha = 0.6, trim = FALSE, color = NA) +
    
    geom_boxplot(width = 0.15,
                 outlier.size  = 0.5,
                 outlier.alpha = 0.3,
                 color = "gray30",
                 fill  = "white",
                 alpha = 0.8) +
    
    # 显著性括号：FDR 校正后的 q 值
    stat_pvalue_manual(
      annot,
      label        = "label",
      tip.length   = 0.01,
      bracket.size = 0.4,
      size         = 3.5
    ) +
    
    geom_vline(xintercept = 2.5,
               linetype  = "dashed",
               color     = "gray50",
               linewidth = 0.5) +
    
    # 原脚本这里画的是「四组混合的整体中位数」，位置取决于各组样本量，
    # 容易被误读为基准线，已改为固定的 PAE 解读阈值。
    geom_hline(data = ref_lines,
               aes(yintercept = yint),
               color     = ref_lines$col,
               linetype  = "dotted",
               linewidth = 0.4,
               inherit.aes = FALSE) +
    
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    
    labs(title = title, x = NULL, y = y_label) +
    
    theme_classic(base_size = 11) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x        = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray92")
    )
}


# ══════════════════════════════════════════════════════════════
# 5. 生成六张图
# ══════════════════════════════════════════════════════════════

plots <- pmap(
  list(panel_spec$panel_id, panel_spec$col,
       panel_spec$panel_type, panel_spec$title),
  function(pid, col, ptype, ttl) {
    plot_pae(
      df    = prepare_pae_data(df_all, col, ptype),
      title = ttl,
      annot = p_annot %>%
        filter(panel_id == pid) %>%
        select(group1, group2, label)
    )
  }
) %>% set_names(panel_spec$panel_id)


# ══════════════════════════════════════════════════════════════
# 6. 拼图与导出
# ══════════════════════════════════════════════════════════════

combined <-
  (plots$var_full | plots$var_nmd | plots$var_div) /
  (plots$wt_full  | plots$wt_nmd  | plots$wt_div) +
  plot_annotation(
    title    = "PAE Comparison: SNV Disease vs Control  |  FS Disease vs Control",
    subtitle = "Top row: VAR protein (variant-level)  |  Bottom row: WT protein (protein-level, deduplicated)",
    caption  = paste0(
      "PAE: lower = more confident inter-residue positioning. ",
      "Dotted reference lines at 5 Å (green) and 10 Å (orange), shown only where in range.\n",
      "Significance: Wilcoxon rank-sum test, Benjamini-Hochberg FDR-adjusted across the ",
      "12 tests in this figure (SNV and FS families corrected separately).\n",
      "* q<0.05  ** q<0.01  *** q<0.001  **** q<0.0001  ns = not significant\n",
      "WT panels are deduplicated to one value per protein per group; VAR panels are variant-level. ",
      "Y-axes are panel-specific — do not compare heights across panels.\n",
      "Note: mean PAE over a whole protein scales with protein length and domain count; ",
      "Full-region differences may partly reflect length rather than structural disruption."
    ),
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      plot.caption  = element_text(hjust = 0, size = 8, color = "gray50")
    )
  )

print(combined)

ggsave("PAE_comparison_fdr.pdf", plot = combined,
       width = 14, height = 10, dpi = 300)
ggsave("PAE_comparison_fdr.png", plot = combined,
       width = 14, height = 10, dpi = 300)


# ══════════════════════════════════════════════════════════════
# 7. 核对输出
# ══════════════════════════════════════════════════════════════

p_annot %>%
  select(panel_id, family, p_raw, q_val, label) %>%
  mutate(across(c(p_raw, q_val), ~ signif(.x, 3))) %>%
  print(n = Inf)

pmap_dfr(
  list(panel_spec$panel_id, panel_spec$col, panel_spec$panel_type),
  function(pid, col, ptype) {
    prepare_pae_data(df_all, col, ptype) %>%
      count(group) %>%
      mutate(panel_id = pid, .before = 1)
  }
) %>%
  pivot_wider(names_from = group, values_from = n) %>%
  print()


# ── 与 Table 1 统一校正（可选）────────────────────────────────
# 本图默认在自己的 12 个检验内校正，因此 q 值与 Table 1（54 个检验）
# 不同。若要完全对齐，先运行 Table 1 脚本保留其 p_all，然后：
#
# t1_lookup <- tibble::tribble(
#   ~panel_id,  ~t1_variable,
#   "var_full", "VAR_pae_full_mean",
#   "var_nmd",  "VAR_pae_NMDPos_mean",
#   "var_div",  "VAR_pae_DivergentPos_mean",
#   "wt_full",  "WT_pae_full_mean",
#   "wt_nmd",   "WT_pae_NMDPos_mean",
#   "wt_div",   "WT_pae_DivergentPos_mean"
# )
# p_annot <- p_raw %>%
#   left_join(t1_lookup, by = "panel_id") %>%
#   left_join(p_all, by = c("t1_variable" = "variable")) %>%
#   mutate(q_val = if_else(family == "SNV", snv_q_raw, fs_q_raw),
#          label = q_stars(q_val))
# 然后重新运行第 5、6 节。