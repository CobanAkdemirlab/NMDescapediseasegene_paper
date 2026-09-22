library(tidyverse)
library(gtsummary)
library(gt)
library(flextable)

# ── Grouping function ───────────────────────────────────────

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

# ── P-value calculation function ────────────────────────────

calc_pairwise_p <- function(data, variables) {
  
  map_dfr(variables, function(v) {
    
    snv_data <- data %>%
      filter(group %in% c("SNV Disease", "SNV Control")) %>%
      select(group, value = all_of(v)) %>%
      drop_na()
    
    fs_data <- data %>%
      filter(group %in% c("FS Disease", "FS Control")) %>%
      select(group, value = all_of(v)) %>%
      drop_na()
    
    snv_p <- tryCatch(
      wilcox.test(value ~ group, data = snv_data)$p.value,
      error = function(e) NA_real_
    )
    
    fs_p <- tryCatch(
      wilcox.test(value ~ group, data = fs_data)$p.value,
      error = function(e) NA_real_
    )
    
    tibble(
      variable = v,
      snv_p = style_pvalue(snv_p, digits = 3),
      fs_p  = style_pvalue(fs_p, digits = 3)
    )
  })
}

# ── Table building function ─────────────────────────────────

make_cider_table <- function(data, title) {
  
  cider_vars <- c(
    "wt_length", "var_length",
    
    "wt_isoelectric_point", "wt_molecular_weight", "wt_fcr_7",
    "wt_count_neg", "wt_count_pos",
    "wt_fraction_negative", "wt_fraction_positive",
    "wt_fraction_expanding_7", "wt_fraction_disorder_promoting",
    "wt_mean_hydropathy", "wt_uversky_hydropathy",
    
    "wt_kappa", "wt_Omega", "wt_delta", "wt_deltaMax",
    
    "var_isoelectric_point", "var_molecular_weight", "var_fcr_7",
    "var_count_neg", "var_count_pos",
    "var_fraction_negative", "var_fraction_positive",
    "var_fraction_expanding_7", "var_fraction_disorder_promoting",
    "var_mean_hydropathy", "var_uversky_hydropathy",
    
    "var_kappa", "var_Omega", "var_delta", "var_deltaMax"
  )
  
  group_map <- c(
    wt_length = "Sequence Length",
    var_length = "Sequence Length",
    
    wt_isoelectric_point = "WT Physicochemical Properties",
    wt_molecular_weight = "WT Physicochemical Properties",
    wt_fcr_7 = "WT Physicochemical Properties",
    wt_count_neg = "WT Physicochemical Properties",
    wt_count_pos = "WT Physicochemical Properties",
    wt_fraction_negative = "WT Physicochemical Properties",
    wt_fraction_positive = "WT Physicochemical Properties",
    wt_fraction_expanding_7 = "WT Physicochemical Properties",
    wt_fraction_disorder_promoting = "WT Physicochemical Properties",
    wt_mean_hydropathy = "WT Physicochemical Properties",
    wt_uversky_hydropathy = "WT Physicochemical Properties",
    
    wt_kappa = "WT Charge Patterning",
    wt_Omega = "WT Charge Patterning",
    wt_delta = "WT Charge Patterning",
    wt_deltaMax = "WT Charge Patterning",
    
    var_isoelectric_point = "VAR Physicochemical Properties",
    var_molecular_weight = "VAR Physicochemical Properties",
    var_fcr_7 = "VAR Physicochemical Properties",
    var_count_neg = "VAR Physicochemical Properties",
    var_count_pos = "VAR Physicochemical Properties",
    var_fraction_negative = "VAR Physicochemical Properties",
    var_fraction_positive = "VAR Physicochemical Properties",
    var_fraction_expanding_7 = "VAR Physicochemical Properties",
    var_fraction_disorder_promoting = "VAR Physicochemical Properties",
    var_mean_hydropathy = "VAR Physicochemical Properties",
    var_uversky_hydropathy = "VAR Physicochemical Properties",
    
    var_kappa = "VAR Charge Patterning",
    var_Omega = "VAR Charge Patterning",
    var_delta = "VAR Charge Patterning",
    var_deltaMax = "VAR Charge Patterning"
  )
  
  p_tbl <- calc_pairwise_p(data, cider_vars)
  
  tbl <- data %>%
    select(group, all_of(cider_vars)) %>%
    tbl_summary(
      by = group,
      statistic = all_continuous() ~ "{median} ({p25}, {p75})",
      digits = all_continuous() ~ 3,
      missing = "no",
      label = list(
        wt_length ~ "WT Sequence Length (aa)",
        var_length ~ "VAR Sequence Length (aa)",
        
        wt_isoelectric_point ~ "WT Isoelectric Point",
        wt_molecular_weight ~ "WT Molecular Weight (Da)",
        wt_fcr_7 ~ "WT Fraction Charged Residues (pH 7)",
        wt_count_neg ~ "WT Count Negative Residues",
        wt_count_pos ~ "WT Count Positive Residues",
        wt_fraction_negative ~ "WT Fraction Negative",
        wt_fraction_positive ~ "WT Fraction Positive",
        wt_fraction_expanding_7 ~ "WT Fraction Expanding (pH 7)",
        wt_fraction_disorder_promoting ~ "WT Fraction Disorder Promoting",
        wt_mean_hydropathy ~ "WT Mean Hydropathy",
        wt_uversky_hydropathy ~ "WT Uversky Hydropathy",
        
        wt_kappa ~ "WT κ (Charge Segregation)",
        wt_Omega ~ "WT Ω (Charge-Hydropathy Mixing)",
        wt_delta ~ "WT δ (Charge Asymmetry)",
        wt_deltaMax ~ "WT δ Max",
        
        var_isoelectric_point ~ "VAR Isoelectric Point",
        var_molecular_weight ~ "VAR Molecular Weight (Da)",
        var_fcr_7 ~ "VAR Fraction Charged Residues (pH 7)",
        var_count_neg ~ "VAR Count Negative Residues",
        var_count_pos ~ "VAR Count Positive Residues",
        var_fraction_negative ~ "VAR Fraction Negative",
        var_fraction_positive ~ "VAR Fraction Positive",
        var_fraction_expanding_7 ~ "VAR Fraction Expanding (pH 7)",
        var_fraction_disorder_promoting ~ "VAR Fraction Disorder Promoting",
        var_mean_hydropathy ~ "VAR Mean Hydropathy",
        var_uversky_hydropathy ~ "VAR Uversky Hydropathy",
        
        var_kappa ~ "VAR κ (Charge Segregation)",
        var_Omega ~ "VAR Ω (Charge-Hydropathy Mixing)",
        var_delta ~ "VAR δ (Charge Asymmetry)",
        var_deltaMax ~ "VAR δ Max"
      )
    ) %>%
    add_overall() %>%
    modify_table_body(
      ~ .x %>%
        left_join(p_tbl, by = "variable") %>%
        mutate(groupname_col = group_map[variable])
    ) %>%
    modify_header(
      label ~ "**Variable**",
      snv_p ~ "**SNV Disease vs Control**",
      fs_p ~ "**FS Disease vs Control**"
    ) %>%
    modify_caption(paste0("**", title, "**")) %>%
    bold_labels()
  
  return(tbl)
}

# ══════════════════════════════════════════════════════════
# Table 1: Divergent Region
# ══════════════════════════════════════════════════════════

df_div <- NMD_region_DivergentPos_CIDER %>%
  mutate(group = make_group(source_folder))

tbl_div <- make_cider_table(
  df_div,
  "Table 2a: CIDER Features — Divergent Region"
)

# ══════════════════════════════════════════════════════════
# Table 2: NMD Region
# ══════════════════════════════════════════════════════════

df_nmd <- NMD_region_NMDPos_CIDER %>%
  mutate(group = make_group(source_folder))

tbl_nmd <- make_cider_table(
  df_nmd,
  "Table 2b: CIDER Features — NMD Region"
)

# ── Print ────────────────────────────────────────────────

tbl_div
tbl_nmd

# ── Export Word ──────────────────────────────────────────

tbl_div %>%
  as_flex_table() %>%
  save_as_docx(path = "Table2a_CIDER_Divergent.docx")

tbl_nmd %>%
  as_flex_table() %>%
  save_as_docx(path = "Table2b_CIDER_NMD.docx")

# ── Export HTML ──────────────────────────────────────────

tbl_div %>%
  as_gt() %>%
  gt::gtsave("Table2a_CIDER_Divergent.html")

tbl_nmd %>%
  as_gt() %>%
  gt::gtsave("Table2b_CIDER_NMD.html")