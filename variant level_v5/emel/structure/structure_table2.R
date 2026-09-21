library(tidyverse)
library(gtsummary)
library(gt)
library(flextable)
library(glue)

# ── Prepare data ──────────────────────────────────────────────

df <- VAR_WT_structural_results %>%
  mutate(
    group = case_when(
      source_folder == "snv_disease" ~ "SNV Disease",
      source_folder == "snv_control" ~ "SNV Control",
      source_folder == "fs_disease"  ~ "FS Disease",
      source_folder == "fs_control"  ~ "FS Control",
      TRUE ~ NA_character_
    ),
    group = factor(
      group,
      levels = c("SNV Disease", "SNV Control", "FS Disease", "FS Control")
    )
  ) %>%
  filter(!is.na(group))

# Same cohort filter as Table 1: disease variants restricted to the
# curated set in all_variants, controls taken as-is.
df_disease <- df %>% filter(key %in% all_variants$key)
df_control <- df %>% filter(source_folder %in% c("snv_control", "fs_control"))
df <- bind_rows(df_disease, df_control)


# ── Delta variable construction ───────────────────────────────

delta_defs <- function(region) {
  tibble::tribble(
    ~delta,              ~var_col,                                          ~wt_col,
    "delta_plddt",       glue("VAR_plddt_{region}_mean"),                   glue("WT_plddt_{region}_mean"),
    "delta_pae",         glue("VAR_pae_{region}_mean"),                     glue("WT_pae_{region}_mean"),
    "delta_rel_asa",     glue("VAR_rel_asa_{region}_mean"),                 glue("WT_rel_asa_{region}_mean"),
    "delta_abs_sasa",    glue("VAR_abs_sasa_dssp_{region}_mean"),           glue("WT_abs_sasa_dssp_{region}_mean"),
    "delta_total_sasa",  glue("VAR_SASA_total_{region}"),                   glue("WT_SASA_total_{region}"),
    "delta_helix",       glue("VAR_sec_struc_{region}_helix_percent"),      glue("WT_sec_struc_{region}_helix_percent"),
    "delta_beta",        glue("VAR_sec_struc_{region}_beta_percent"),       glue("WT_sec_struc_{region}_beta_percent"),
    "delta_coil",        glue("VAR_sec_struc_{region}_coil_percent"),       glue("WT_sec_struc_{region}_coil_percent")
  )
}

add_deltas <- function(data, region) {
  d <- delta_defs(region)
  out <- data
  for (i in seq_len(nrow(d))) {
    out[[d$delta[i]]] <- data[[d$var_col[i]]] - data[[d$wt_col[i]]]
  }
  out
}

regions <- c(full = "Full", NMDPos = "NMD", DivergentPos = "Divergent")
delta_vars <- delta_defs("full")$delta   # names are region-independent


# ── Raw p-values across ALL regions, then BH-FDR ──────────────
# Delta variables are variant-level (one independent value per variant),
# no protein-level deduplication needed here
# features in Table 1.

get_raw_p <- function(data, vars) {
  map_dfr(vars, function(v) {
    d <- data %>%
      select(group, value = all_of(v)) %>%
      filter(!is.na(value))
    
    test_pair <- function(g1, g2) {
      dd <- d %>%
        filter(group %in% c(g1, g2)) %>%
        mutate(group = droplevels(group))
      tryCatch(
        wilcox.test(value ~ group, data = dd)$p.value,
        error = function(e) NA_real_
      )
    }
    
    tibble(
      variable  = v,
      snv_p_raw = test_pair("SNV Disease", "SNV Control"),
      fs_p_raw  = test_pair("FS Disease",  "FS Control")
    )
  })
}

# 8 deltas x 3 regions = 24 unique tests per comparison family.
p_all <- map_dfr(names(regions), function(rg) {
  add_deltas(df, rg) %>%
    get_raw_p(delta_vars) %>%
    mutate(region = rg)
}) %>%
  mutate(
    snv_q_raw = p.adjust(snv_p_raw, method = "BH"),
    fs_q_raw  = p.adjust(fs_p_raw,  method = "BH")
  )

n_tests <- nrow(p_all)


# ── Main table function ───────────────────────────────────────

make_delta_table <- function(data, region, title, p_all) {
  
  suffix <- regions[[region]]
  df_tbl <- add_deltas(data, region) %>%
    select(group, all_of(delta_vars))
  
  p_tbl <- p_all %>%
    filter(region == !!region) %>%
    transmute(
      variable,
      snv_p = style_pvalue(snv_p_raw, digits = 3),
      snv_q = style_pvalue(snv_q_raw, digits = 3),
      fs_p  = style_pvalue(fs_p_raw,  digits = 3),
      fs_q  = style_pvalue(fs_q_raw,  digits = 3)
    )
  
  group_map <- setNames(
    rep(paste0("Δ Structural Features — ", suffix), length(delta_vars)),
    delta_vars
  )
  
  df_tbl %>%
    tbl_summary(
      by = group,
      statistic = all_continuous() ~ "{median} ({p25}, {p75})",
      digits = all_continuous() ~ 2,
      missing = "no",
      label = list(
        delta_plddt      ~ glue("Δ pLDDT — {suffix}"),
        delta_pae        ~ glue("Δ PAE — {suffix}"),
        delta_rel_asa    ~ glue("Δ Relative ASA — {suffix}"),
        delta_abs_sasa   ~ glue("Δ Absolute SASA — {suffix}"),
        delta_total_sasa ~ glue("Δ Total SASA — {suffix}"),
        delta_helix      ~ glue("Δ Helix % — {suffix}"),
        delta_beta       ~ glue("Δ Beta % — {suffix}"),
        delta_coil       ~ glue("Δ Coil % — {suffix}")
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
      snv_p ~ "**SNV p**",
      snv_q ~ "**SNV q (FDR)**",
      fs_p  ~ "**FS p**",
      fs_q  ~ "**FS q (FDR)**"
    ) %>%
    modify_footnote(
      c(snv_q, fs_q) ~ glue(
        "Benjamini-Hochberg FDR-adjusted across all {n_tests} Wilcoxon ",
        "tests (8 delta variables x 3 regions), separately within the ",
        "SNV and FS comparison families."
      )
    ) %>%
    modify_caption(paste0("**", title, "**")) %>%
    bold_labels()
}


# ── Generate tables ───────────────────────────────────────────

tbl2_full <- make_delta_table(df, "full",
                              "Table 2a: VAR - WT Structural Differences — Full Region", p_all)
tbl2_nmd  <- make_delta_table(df, "NMDPos",
                              "Table 2b: VAR - WT Structural Differences — NMD Region", p_all)
tbl2_div  <- make_delta_table(df, "DivergentPos",
                              "Table 2c: VAR - WT Structural Differences — Divergent Region", p_all)

tbl2_full
tbl2_nmd
tbl2_div


# ── Export Word ───────────────────────────────────────────────

tbl2_full %>% as_flex_table() %>% save_as_docx(path = "Table2a_Delta_Full_fdr.docx")
tbl2_nmd  %>% as_flex_table() %>% save_as_docx(path = "Table2b_Delta_NMD_fdr.docx")
tbl2_div  %>% as_flex_table() %>% save_as_docx(path = "Table2c_Delta_Divergent_fdr.docx")


# ── Export HTML ───────────────────────────────────────────────

tbl2_full %>% as_gt() %>% gt::gtsave("Table2a_Delta_Full_fdr.html")
tbl2_nmd  %>% as_gt() %>% gt::gtsave("Table2b_Delta_NMD_fdr.html")
tbl2_div  %>% as_gt() %>% gt::gtsave("Table2c_Delta_Divergent_fdr.html")


# ── Optional: joint correction across Table 1 + Table 2 ───────
# If both tables appear in the same paper, a reviewer may expect one
# FDR family. Run the Table 1 script first (it creates its own p_all),
# save it as p_all_t1, then:
#
# p_joint <- bind_rows(
#   p_all_t1 %>% mutate(table = "T1"),
#   p_all    %>% mutate(table = "T2")
# ) %>%
#   mutate(
#     snv_q_raw = p.adjust(snv_p_raw, method = "BH"),
#     fs_q_raw  = p.adjust(fs_p_raw,  method = "BH")
#   )
#
# then rebuild both sets of tables from the p_joint subsets.