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
    ),
    truncation_length = wt_length - var_length,
    truncation_pct = (wt_length - var_length) / wt_length * 100
  ) %>%
  filter(!is.na(group))

df_disease <- df %>% filter(key %in% all_variants$key)
df_control <- df %>% filter(source_folder %in% c("snv_control", "fs_control"))
df <- bind_rows(df_disease, df_control)


# ── Variable maps ─────────────────────────────────────────────

# Length/position variables shared by all three tables (tested ONCE)
shared_vars <- c(
  "wt_length", "var_length", "NMD_region_start_prot",
  "Divergence_start", "truncation_length", "truncation_pct"
)

# Region-specific structural columns: display name -> source column
region_cols <- function(region) {
  c(
    VAR_plddt      = glue("VAR_plddt_{region}_mean"),
    VAR_pae        = glue("VAR_pae_{region}_mean"),
    VAR_rel_asa    = glue("VAR_rel_asa_{region}_mean"),
    VAR_abs_sasa   = glue("VAR_abs_sasa_dssp_{region}_mean"),
    VAR_SASA_total = glue("VAR_SASA_total_{region}"),
    VAR_helix      = glue("VAR_sec_struc_{region}_helix_percent"),
    VAR_beta       = glue("VAR_sec_struc_{region}_beta_percent"),
    VAR_coil       = glue("VAR_sec_struc_{region}_coil_percent"),
    WT_plddt       = glue("WT_plddt_{region}_mean"),
    WT_pae         = glue("WT_pae_{region}_mean"),
    WT_rel_asa     = glue("WT_rel_asa_{region}_mean"),
    WT_abs_sasa    = glue("WT_abs_sasa_dssp_{region}_mean"),
    WT_SASA_total  = glue("WT_SASA_total_{region}"),
    WT_helix       = glue("WT_sec_struc_{region}_helix_percent"),
    WT_beta        = glue("WT_sec_struc_{region}_beta_percent"),
    WT_coil        = glue("WT_sec_struc_{region}_coil_percent")
  )
}

all_test_cols <- unique(c(
  shared_vars,
  unname(region_cols("full")),
  unname(region_cols("NMDPos")),
  unname(region_cols("DivergentPos"))
))

# Protein-level variables: identical for every variant of the same protein.
# Testing them on the variant-level rows pseudo-replicates (one gene with
# many variants contributes many copies of the same value), so these are
# deduplicated to one value per protein per group before testing.
protein_level_vars <- c(
  "wt_length",
  grep("^WT_", all_test_cols, value = TRUE)
  # If these are constant within a protein in your data (they appear to be),
  # uncomment to treat them as protein-level too:
   , "NMD_region_start_prot", "Divergence_start"
)


# ── Raw p-values (numeric) ────────────────────────────────────

get_raw_p <- function(data, vars, protein_level_vars = character()) {
  map_dfr(vars, function(v) {
    d <- data %>%
      select(group, uniprot_id, value = all_of(v)) %>%
      filter(!is.na(value))
    
    if (v %in% protein_level_vars) {
      d <- d %>% distinct(group, uniprot_id, .keep_all = TRUE)
    }
    
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

# One test per unique variable, then BH-FDR across the full family
# (54 tests per comparison: 6 shared + 16 structural x 3 regions).
p_all <- get_raw_p(df, all_test_cols, protein_level_vars) %>%
  mutate(
    snv_q_raw = p.adjust(snv_p_raw, method = "BH"),
    fs_q_raw  = p.adjust(fs_p_raw,  method = "BH")
  )

# Build the per-table p/q lookup keyed on the display variable names
make_p_tbl <- function(p_all, region) {
  cols <- region_cols(region)
  tibble(
    variable = c(shared_vars, names(cols)),
    src      = c(shared_vars, unname(cols))
  ) %>%
    left_join(p_all, by = c("src" = "variable")) %>%
    transmute(
      variable,
      snv_p = style_pvalue(snv_p_raw, digits = 3),
      snv_q = style_pvalue(snv_q_raw, digits = 3),
      fs_p  = style_pvalue(fs_p_raw,  digits = 3),
      fs_q  = style_pvalue(fs_q_raw,  digits = 3)
    )
}


# ── Main table function ───────────────────────────────────────

make_table1 <- function(data, region, title, p_all) {
  
  suffix <- c(full = "Full", NMDPos = "NMD", DivergentPos = "Divergent")[[region]]
  cols   <- region_cols(region)
  
  df_tbl <- data %>%
    mutate(!!!syms(setNames(unname(cols), names(cols)))) %>%
    select(group, all_of(shared_vars), all_of(names(cols)))
  
  p_tbl <- make_p_tbl(p_all, region)
  
  group_map <- c(
    setNames(rep("Protein Length & Position", length(shared_vars)), shared_vars),
    setNames(rep(paste0("VAR Structural Features — ", suffix), 8),
             names(cols)[1:8]),
    setNames(rep(paste0("WT Structural Features — ", suffix), 8),
             names(cols)[9:16])
  )
  
  df_tbl %>%
    tbl_summary(
      by = group,
      statistic = all_continuous() ~ "{median} ({p25}, {p75})",
      digits = all_continuous() ~ 2,
      missing = "no",
      label = list(
        wt_length             ~ "WT Protein Length (aa)",
        var_length            ~ "VAR Protein Length (aa)",
        NMD_region_start_prot ~ "NMD Region Start (aa)",
        Divergence_start      ~ "Divergence Start (aa)",
        truncation_length     ~ "Truncation Length (aa)",
        truncation_pct        ~ "Truncation (%)",
        VAR_plddt             ~ glue("VAR pLDDT — {suffix} (mean)"),
        VAR_pae               ~ glue("VAR PAE — {suffix} (mean)"),
        VAR_rel_asa           ~ glue("VAR Relative ASA — {suffix} (mean)"),
        VAR_abs_sasa          ~ glue("VAR Absolute SASA — {suffix} (mean)"),
        VAR_SASA_total        ~ glue("VAR Total SASA — {suffix}"),
        VAR_helix             ~ glue("VAR Helix % — {suffix}"),
        VAR_beta              ~ glue("VAR Beta % — {suffix}"),
        VAR_coil              ~ glue("VAR Coil % — {suffix}"),
        WT_plddt              ~ glue("WT pLDDT — {suffix} (mean)"),
        WT_pae                ~ glue("WT PAE — {suffix} (mean)"),
        WT_rel_asa            ~ glue("WT Relative ASA — {suffix} (mean)"),
        WT_abs_sasa           ~ glue("WT Absolute SASA — {suffix} (mean)"),
        WT_SASA_total         ~ glue("WT Total SASA — {suffix}"),
        WT_helix              ~ glue("WT Helix % — {suffix}"),
        WT_beta               ~ glue("WT Beta % — {suffix}"),
        WT_coil               ~ glue("WT Coil % — {suffix}")
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
      c(snv_q, fs_q) ~ paste0(
        "Benjamini–Hochberg FDR-adjusted across all ",
        nrow(p_all),
        " unique Wilcoxon tests per comparison (shared length/position ",
        "variables counted once across the three tables). ",
        "WT features and WT length tested at the protein level ",
        "(one value per protein per group)."
      )
    ) %>%
    modify_caption(paste0("**", title, "**")) %>%
    bold_labels()
}


# ── Generate three tables ─────────────────────────────────────

tbl_full <- make_table1(df, "full",         "Table 1a: Structural Features — Full Region",      p_all)
tbl_nmd  <- make_table1(df, "NMDPos",       "Table 1b: Structural Features — NMD Region",       p_all)
tbl_div  <- make_table1(df, "DivergentPos", "Table 1c: Structural Features — Divergent Region", p_all)

tbl_full
tbl_nmd
tbl_div


# ── Export Word ───────────────────────────────────────────────

tbl_full %>% as_flex_table() %>% save_as_docx(path = "Table1a_Full_p_fdr.docx")
tbl_nmd  %>% as_flex_table() %>% save_as_docx(path = "Table1b_NMD_p_fdr.docx")
tbl_div  %>% as_flex_table() %>% save_as_docx(path = "Table1c_Divergent_p_fdr.docx")


# ── Export HTML ───────────────────────────────────────────────

tbl_full %>% as_gt() %>% gt::gtsave("Table1a_Full_p_fdr.html")
tbl_nmd  %>% as_gt() %>% gt::gtsave("Table1b_NMD_p_fdr.html")
tbl_div  %>% as_gt() %>% gt::gtsave("Table1c_Divergent_p_fdr.html")