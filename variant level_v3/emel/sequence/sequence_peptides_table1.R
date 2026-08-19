
library(tidyverse)
library(gtsummary)
library(gt)
library(flextable)
library(officer)

# ══════════════════════════════════════════════════════════
# 通用函数
# ══════════════════════════════════════════════════════════

make_group <- function(source_folder) {
  factor(
    case_when(
      source_folder == "snv_disease" ~ "SNV Disease",
      source_folder == "snv_control" ~ "SNV Control",
      source_folder == "fs_disease"  ~ "FS Disease",
      source_folder == "fs_control"  ~ "FS Control"
    ),
    levels = c("SNV Disease", "SNV Control", "FS Disease", "FS Control")
  )
}

# ══════════════════════════════════════════════════════════
# 计算 VAR - WT 差值
# ══════════════════════════════════════════════════════════

compute_diff <- function(data, has_length = TRUE) {
  
  df <- data %>%
    mutate(
      # 连续变量差值
      diff_aliphatic_index    = var_aliphatic_index    - wt_aliphatic_index,
      diff_boman              = var_boman               - wt_boman,
      diff_charge             = var_charge              - wt_charge,
      diff_entropy            = var_entropy             - wt_entropy,
      diff_instability_index  = var_instability_index  - wt_instability_index,
      diff_isoelectric_point  = var_isoelectric_point  - wt_isoelectric_point,
      diff_molecular_weight   = var_molecular_weight   - wt_molecular_weight,
      diff_mz                 = var_mz                 - wt_mz,
      diff_max_frequency      = var_max_frequency      - wt_max_frequency,
      diff_energy_cost        = var_energy_cost        - wt_energy_cost,
      diff_nutrient_cost      = var_nutrient_cost      - wt_nutrient_cost,
      diff_mass_shift         = var_mass_shift         - wt_mass_shift,
      diff_hydrophobicity     = var_hydrophobicity     - wt_hydrophobicity,
      diff_hydrophobic_moment = var_hydrophobic_moment - wt_hydrophobic_moment,
      
      # 结构类别是否改变
      structural_class_changed = factor(
        ifelse(var_structural_class == wt_structural_class,
               "Unchanged", "Changed")
      )
    )
  
  if (has_length) {
    df <- df %>%
      mutate(diff_length = var_length - wt_length)
  }
  
  return(df)
}

# ══════════════════════════════════════════════════════════
# 建表函数
# ══════════════════════════════════════════════════════════

make_diff_merged_table <- function(data,
                                   has_length = TRUE,
                                   caption) {
  
  # 选列
  base_vars <- c(
    "group",
    "diff_aliphatic_index",
    "diff_boman",
    "diff_charge",
    "diff_entropy",
    "diff_instability_index",
    "diff_isoelectric_point",
    "diff_molecular_weight",
    "diff_mz",
    "diff_max_frequency",
    "diff_energy_cost",
    "diff_nutrient_cost",
    "diff_mass_shift",
    "diff_hydrophobicity",
    "diff_hydrophobic_moment",
    "structural_class_changed"
  )
  
  if (has_length) {
    all_vars <- c("group", "diff_length", base_vars[-1])
  } else {
    all_vars <- base_vars
  }
  
  df_sel <- data %>% select(all_of(all_vars))
  
  # labels
  base_labels <- list(
    diff_aliphatic_index    ~ "Δ Aliphatic Index (VAR−WT)",
    diff_boman              ~ "Δ Boman Index (VAR−WT)",
    diff_charge             ~ "Δ Net Charge (VAR−WT)",
    diff_entropy            ~ "Δ Sequence Entropy (VAR−WT)",
    diff_instability_index  ~ "Δ Instability Index (VAR−WT)",
    diff_isoelectric_point  ~ "Δ Isoelectric Point (VAR−WT)",
    diff_molecular_weight   ~ "Δ Molecular Weight (VAR−WT, Da)",
    diff_mz                 ~ "Δ m/z (VAR−WT)",
    diff_max_frequency      ~ "Δ Max Residue Frequency (VAR−WT)",
    diff_energy_cost        ~ "Δ Energy Cost (VAR−WT)",
    diff_nutrient_cost      ~ "Δ Nutrient Cost (VAR−WT)",
    diff_mass_shift         ~ "Δ Mass Shift (VAR−WT)",
    diff_hydrophobicity     ~ "Δ Hydrophobicity (VAR−WT)",
    diff_hydrophobic_moment ~ "Δ Hydrophobic Moment (VAR−WT)",
    structural_class_changed ~ "Structural Class Changed"
  )
  
  if (has_length) {
    all_labels <- c(
      list(diff_length ~ "Δ Sequence Length (VAR−WT, aa)"),
      base_labels
    )
  } else {
    all_labels <- base_labels
  }
  
  # group_map
  base_map <- c(
    diff_aliphatic_index     = "VAR − WT Peptide Property Changes",
    diff_boman               = "VAR − WT Peptide Property Changes",
    diff_charge              = "VAR − WT Peptide Property Changes",
    diff_entropy             = "VAR − WT Peptide Property Changes",
    diff_instability_index   = "VAR − WT Peptide Property Changes",
    diff_isoelectric_point   = "VAR − WT Peptide Property Changes",
    diff_molecular_weight    = "VAR − WT Peptide Property Changes",
    diff_mz                  = "VAR − WT Peptide Property Changes",
    diff_max_frequency       = "VAR − WT Peptide Property Changes",
    diff_energy_cost         = "VAR − WT Peptide Property Changes",
    diff_nutrient_cost       = "VAR − WT Peptide Property Changes",
    diff_mass_shift          = "VAR − WT Peptide Property Changes",
    diff_hydrophobicity      = "VAR − WT Peptide Property Changes",
    diff_hydrophobic_moment  = "VAR − WT Peptide Property Changes",
    structural_class_changed = "Structural Class"
  )
  
  if (has_length) {
    group_map <- c(
      diff_length = "Sequence Length Change",
      base_map
    )
  } else {
    group_map <- base_map
  }
  
  # ── SNV 子表 ──────────────────────────────────────────
  tbl_snv <- df_sel %>%
    filter(group %in% c("SNV Disease", "SNV Control")) %>%
    mutate(group = droplevels(group)) %>%
    tbl_summary(
      by        = group,
      statistic = list(
        all_continuous()  ~ "{median} ({p25}, {p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits  = list(
        all_continuous()  ~ 2,
        all_categorical() ~ c(0, 1)
      ),
      missing = "no",
      label   = all_labels
    ) %>%
    add_p(
      test = list(
        all_continuous()  ~ "wilcox.test",
        all_categorical() ~ "chisq.test"
      ),
      pvalue_fun = ~style_pvalue(.x, digits = 3)
    ) %>%
    modify_header(
      label           ~ "**Variable**",
      all_stat_cols() ~ "**{level}**\nN = {n}"
    ) %>%
    bold_labels()
  
  # ── FS 子表 ───────────────────────────────────────────
  tbl_fs <- df_sel %>%
    filter(group %in% c("FS Disease", "FS Control")) %>%
    mutate(group = droplevels(group)) %>%
    tbl_summary(
      by        = group,
      statistic = list(
        all_continuous()  ~ "{median} ({p25}, {p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits  = list(
        all_continuous()  ~ 2,
        all_categorical() ~ c(0, 1)
      ),
      missing = "no",
      label   = all_labels
    ) %>%
    add_p(
      test = list(
        all_continuous()  ~ "wilcox.test",
        all_categorical() ~ "chisq.test"
      ),
      pvalue_fun = ~style_pvalue(.x, digits = 3)
    ) %>%
    modify_header(
      all_stat_cols() ~ "**{level}**\nN = {n}"
    ) %>%
    bold_labels()
  
  # ── 合并 ──────────────────────────────────────────────
  tbl_merged <- tbl_merge(
    tbls        = list(tbl_snv, tbl_fs),
    tab_spanner = c("**SNV Disease vs Control**",
                    "**FS Disease vs Control**")
  ) %>%
    modify_table_body(
      ~ .x %>% mutate(groupname_col = group_map[variable])
    ) %>%
    modify_caption(paste0("**", caption, "**"))
  
  return(tbl_merged)
}

# ══════════════════════════════════════════════════════════
# 准备数据
# ══════════════════════════════════════════════════════════

df_div_raw <- NMD_region_DivergentPos_peptide_props %>%
  mutate(
    group      = make_group(source_folder),
    wt_length  = wt_NMD_length,
    var_length = var_NMD_length,
    wt_structural_class  = factor(wt_structural_class),
    var_structural_class = factor(var_structural_class)
  )

df_nmd_raw <- NMD_region_NMDPos_peptide_props %>%
  mutate(
    group = make_group(source_folder),
    wt_structural_class  = factor(wt_structural_class),
    var_structural_class = factor(var_structural_class)
  )

# 计算差值
df_div_diff <- compute_diff(df_div_raw, has_length = TRUE)
df_nmd_diff <- compute_diff(df_nmd_raw, has_length = FALSE)

# ══════════════════════════════════════════════════════════
# 生成表格
# ══════════════════════════════════════════════════════════

tbl_div_merged <- make_diff_merged_table(
  data       = df_div_diff,
  has_length = TRUE,
  caption    = "Table 3a: VAR vs WT Peptide Property Changes — Divergent Region"
)

tbl_nmd_merged <- make_diff_merged_table(
  data       = df_nmd_diff,
  has_length = FALSE,
  caption    = "Table 3b: VAR vs WT Peptide Property Changes — NMD Region"
)

# ── 打印 ──────────────────────────────────────────────────

tbl_div_merged
tbl_nmd_merged

# ── 导出 HTML ─────────────────────────────────────────────

tbl_div_merged %>% as_gt() %>% gt::gtsave("Table3a_Divergent_diff.html")
tbl_nmd_merged %>% as_gt() %>% gt::gtsave("Table3b_NMD_diff.html")

# ── 导出 RTF ─────────────────────────────────────────────

tbl_div_merged %>% as_gt() %>% gt::gtsave("Table3a_Divergent_diff.rtf")
tbl_nmd_merged %>% as_gt() %>% gt::gtsave("Table3b_NMD_diff.rtf")

# ── 导出 Word ─────────────────────────────────────────────

doc <- read_docx() %>%
  body_add_par(
    "Table 3a: VAR vs WT Changes — Divergent Region",
    style = "heading 1"
  ) %>%
  body_add_par("SNV Disease vs SNV Control", style = "heading 2") %>%
  body_add_flextable(
    tbl_div_merged$tbls[[1]] %>% as_flex_table()
  ) %>%
  body_add_break() %>%
  body_add_par("FS Disease vs FS Control", style = "heading 2") %>%
  body_add_flextable(
    tbl_div_merged$tbls[[2]] %>% as_flex_table()
  ) %>%
  body_add_break() %>%
  body_add_par(
    "Table 3b: VAR vs WT Changes — NMD Region",
    style = "heading 1"
  ) %>%
  body_add_par("SNV Disease vs SNV Control", style = "heading 2") %>%
  body_add_flextable(
    tbl_nmd_merged$tbls[[1]] %>% as_flex_table()
  ) %>%
  body_add_break() %>%
  body_add_par("FS Disease vs FS Control", style = "heading 2") %>%
  body_add_flextable(
    tbl_nmd_merged$tbls[[2]] %>% as_flex_table()
  )

print(doc, target = "Table3_VAR_WT_diff.docx")