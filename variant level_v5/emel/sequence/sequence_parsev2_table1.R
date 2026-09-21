library(tidyverse)
library(gtsummary)
library(gt)
library(flextable)
library(glue)

extract_group <- function(id_col) {
  factor(
    case_when(
      str_starts(id_col, "snv_disease") ~ "SNV Disease",
      str_starts(id_col, "snv_control") ~ "SNV Control",
      str_starts(id_col, "fs_disease")  ~ "FS Disease",
      str_starts(id_col, "fs_control")  ~ "FS Control",
      TRUE ~ NA_character_
    ),
    levels = c("SNV Disease", "SNV Control", "FS Disease", "FS Control")
  )
}

calc_pairwise_p <- function(data, variables) {
  variables <- variables[variables %in% colnames(data)]
  
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

make_parse_table <- function(data, vars, labels, group_labels, title) {
  vars <- vars[vars %in% colnames(data)]
  p_tbl <- calc_pairwise_p(data, vars)
  
  data %>%
    select(group, all_of(vars)) %>%
    tbl_summary(
      by = group,
      statistic = all_continuous() ~ "{median} ({p25}, {p75})",
      digits = all_continuous() ~ 2,
      missing = "no",
      label = labels
    ) %>%
    add_overall() %>%
    bold_labels() %>%
    modify_table_body(
      ~ .x %>%
        left_join(p_tbl, by = "variable") %>%
        mutate(groupname_col = group_labels[variable])
    ) %>%
    modify_caption(paste0("**", title, "**")) %>%
    modify_header(
      label ~ "**Variable**",
      snv_p ~ "**SNV Disease vs Control**",
      fs_p  ~ "**FS Disease vs Control**"
    )
}

# ══════════════════════════════════════════════════════════
# Table 1a: Full length
# ══════════════════════════════════════════════════════════

df_full <- full_length_PARSE_v2 %>%
  rename(
    var_classifier_dist = `Σ classifier distance P...4`,
    var_longest_PS_IDR  = `length of longest PS IDR...5`,
    wt_classifier_dist  = `Σ classifier distance P...9`,
    wt_longest_PS_IDR   = `length of longest PS IDR...10`
  ) %>%
  mutate(
    group = extract_group(id),
    truncation_length = wt_length - var_length,
    truncation_pct = (wt_length - var_length) / wt_length * 100,
    diff_classifier_dist = var_classifier_dist - wt_classifier_dist,
    diff_longest_PS_IDR = var_longest_PS_IDR - wt_longest_PS_IDR
  )

full_vars <- c(
  "wt_length", "var_length",
  "truncation_length", "truncation_pct",
  "var_classifier_dist", "wt_classifier_dist", "diff_classifier_dist",
  "var_longest_PS_IDR", "wt_longest_PS_IDR", "diff_longest_PS_IDR"
)

full_labels <- list(
  wt_length ~ "WT Protein Length (aa)",
  var_length ~ "VAR Protein Length (aa)",
  truncation_length ~ "Truncation Length (aa)",
  truncation_pct ~ "Truncation (%)",
  var_classifier_dist ~ "VAR Σ Classifier Distance P",
  wt_classifier_dist ~ "WT Σ Classifier Distance P",
  diff_classifier_dist ~ "VAR - WT Σ Classifier Distance P",
  var_longest_PS_IDR ~ "VAR Length of Longest PS IDR",
  wt_longest_PS_IDR ~ "WT Length of Longest PS IDR",
  diff_longest_PS_IDR ~ "VAR - WT Length of Longest PS IDR"
)

full_group_labels <- c(
  wt_length = "Protein Length",
  var_length = "Protein Length",
  truncation_length = "Protein Length",
  truncation_pct = "Protein Length",
  var_classifier_dist = "PARSE Features — Full",
  wt_classifier_dist = "PARSE Features — Full",
  diff_classifier_dist = "PARSE Features — Full",
  var_longest_PS_IDR = "PARSE Features — Full",
  wt_longest_PS_IDR = "PARSE Features — Full",
  diff_longest_PS_IDR = "PARSE Features — Full"
)

tbl_full <- make_parse_table(
  df_full,
  full_vars,
  full_labels,
  full_group_labels,
  "Table 1a: PARSE Features — Full Length"
)

# ══════════════════════════════════════════════════════════
# Table 1b: Divergent region
# ══════════════════════════════════════════════════════════

df_div <- nmd_length_DivergentPos_PARSE_v2 %>%
  rename(
    var_classifier_dist = `Σ classifier distance P...4`,
    var_longest_PS_IDR  = `length of longest PS IDR...5`,
    wt_classifier_dist  = `Σ classifier distance P...9`,
    wt_longest_PS_IDR   = `length of longest PS IDR...10`
  ) %>%
  mutate(
    group = extract_group(id),
    truncation_length = wt_NMD_length - var_NMD_length,
    truncation_pct = (wt_NMD_length - var_NMD_length) / wt_NMD_length * 100,
    diff_classifier_dist = var_classifier_dist - wt_classifier_dist,
    diff_longest_PS_IDR = var_longest_PS_IDR - wt_longest_PS_IDR
  )

div_vars <- c(
  "wt_NMD_length", "var_NMD_length",
  "truncation_length", "truncation_pct",
  "var_classifier_dist", "wt_classifier_dist", "diff_classifier_dist",
  "var_longest_PS_IDR", "wt_longest_PS_IDR", "diff_longest_PS_IDR"
)

div_labels <- list(
  wt_NMD_length ~ "WT NMD Region Length (aa)",
  var_NMD_length ~ "VAR NMD Region Length (aa)",
  truncation_length ~ "Truncation Length (aa)",
  truncation_pct ~ "Truncation (%)",
  var_classifier_dist ~ "VAR Σ Classifier Distance P",
  wt_classifier_dist ~ "WT Σ Classifier Distance P",
  diff_classifier_dist ~ "VAR - WT Σ Classifier Distance P",
  var_longest_PS_IDR ~ "VAR Length of Longest PS IDR",
  wt_longest_PS_IDR ~ "WT Length of Longest PS IDR",
  diff_longest_PS_IDR ~ "VAR - WT Length of Longest PS IDR"
)

div_group_labels <- c(
  wt_NMD_length = "Sequence Length",
  var_NMD_length = "Sequence Length",
  truncation_length = "Sequence Length",
  truncation_pct = "Sequence Length",
  var_classifier_dist = "PARSE Features — Divergent Region",
  wt_classifier_dist = "PARSE Features — Divergent Region",
  diff_classifier_dist = "PARSE Features — Divergent Region",
  var_longest_PS_IDR = "PARSE Features — Divergent Region",
  wt_longest_PS_IDR = "PARSE Features — Divergent Region",
  diff_longest_PS_IDR = "PARSE Features — Divergent Region"
)

tbl_div <- make_parse_table(
  df_div,
  div_vars,
  div_labels,
  div_group_labels,
  "Table 1b: PARSE Features — Divergent Region"
)

# ══════════════════════════════════════════════════════════
# Table 1c: NMD region
# ══════════════════════════════════════════════════════════

df_nmd <- nmd_length_NMDPos_PARSE_v2 %>%
  rename(
    var_classifier_dist = `Σ classifier distance P...4`,
    var_longest_PS_IDR  = `length of longest PS IDR...5`,
    wt_classifier_dist  = `Σ classifier distance P...9`,
    wt_longest_PS_IDR   = `length of longest PS IDR...10`
  ) %>%
  rename(
    wt_NMD_length = any_of("wt_NMD_length_old"),
    var_NMD_length = any_of("var_NMD_length_old")
  ) %>%
  mutate(
    group = extract_group(id),
    truncation_length = wt_NMD_length - var_NMD_length,
    truncation_pct = (wt_NMD_length - var_NMD_length) / wt_NMD_length * 100,
    diff_classifier_dist = var_classifier_dist - wt_classifier_dist,
    diff_longest_PS_IDR = var_longest_PS_IDR - wt_longest_PS_IDR
  )

nmd_vars <- c(
  "wt_NMD_length", "var_NMD_length",
  "truncation_length", "truncation_pct",
  "var_classifier_dist", "wt_classifier_dist", "diff_classifier_dist",
  "var_longest_PS_IDR", "wt_longest_PS_IDR", "diff_longest_PS_IDR"
)

nmd_labels <- list(
  wt_NMD_length ~ "WT NMD Region Length (aa)",
  var_NMD_length ~ "VAR NMD Region Length (aa)",
  truncation_length ~ "Truncation Length (aa)",
  truncation_pct ~ "Truncation (%)",
  var_classifier_dist ~ "VAR Σ Classifier Distance P",
  wt_classifier_dist ~ "WT Σ Classifier Distance P",
  diff_classifier_dist ~ "VAR - WT Σ Classifier Distance P",
  var_longest_PS_IDR ~ "VAR Length of Longest PS IDR",
  wt_longest_PS_IDR ~ "WT Length of Longest PS IDR",
  diff_longest_PS_IDR ~ "VAR - WT Length of Longest PS IDR"
)

nmd_group_labels <- c(
  wt_NMD_length = "Sequence Length",
  var_NMD_length = "Sequence Length",
  truncation_length = "Sequence Length",
  truncation_pct = "Sequence Length",
  var_classifier_dist = "PARSE Features — NMD Region",
  wt_classifier_dist = "PARSE Features — NMD Region",
  diff_classifier_dist = "PARSE Features — NMD Region",
  var_longest_PS_IDR = "PARSE Features — NMD Region",
  wt_longest_PS_IDR = "PARSE Features — NMD Region",
  diff_longest_PS_IDR = "PARSE Features — NMD Region"
)

tbl_nmd <- make_parse_table(
  df_nmd,
  nmd_vars,
  nmd_labels,
  nmd_group_labels,
  "Table 1c: PARSE Features — NMD Region"
)

# ── print ────────────────────────────────────────────────

tbl_full
tbl_div
tbl_nmd

# ── export Word ──────────────────────────────────────────

tbl_full %>%
  as_flex_table() %>%
  save_as_docx(path = "Table1a_PARSE_Full.docx")

tbl_div %>%
  as_flex_table() %>%
  save_as_docx(path = "Table1b_PARSE_Divergent.docx")

tbl_nmd %>%
  as_flex_table() %>%
  save_as_docx(path = "Table1c_PARSE_NMD.docx")

# ── export HTML ──────────────────────────────────────────

tbl_full %>%
  as_gt() %>%
  gt::gtsave("Table1a_PARSE_Full.html")

tbl_div %>%
  as_gt() %>%
  gt::gtsave("Table1b_PARSE_Divergent.html")

tbl_nmd %>%
  as_gt() %>%
  gt::gtsave("Table1c_PARSE_NMD.html")