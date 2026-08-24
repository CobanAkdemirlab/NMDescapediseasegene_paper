library(tidyverse)
library(lme4)
library(lmerTest)
library(gt)

# ══════════════════════════════════════════════════════════════
# 1. Data preparation (consistent with Table 1 / Table 2)
# ══════════════════════════════════════════════════════════════

df_all <- VAR_WT_structural_results %>%
  mutate(
    group = case_when(
      source_folder == "snv_disease" ~ "SNV Disease",
      source_folder == "snv_control" ~ "SNV Control",
      source_folder == "fs_disease"  ~ "FS Disease",
      source_folder == "fs_control"  ~ "FS Control",
      TRUE ~ NA_character_
    ),
    truncation_length = wt_length - var_length,
    truncation_pct    = (wt_length - var_length) / wt_length * 100
  ) %>%
  filter(!is.na(group))

df_disease <- df_all %>% filter(key %in% all_variants$key)
df_control <- df_all %>% filter(source_folder %in% c("snv_control", "fs_control"))
df_all     <- bind_rows(df_disease, df_control)

for (rg in c("full", "NMDPos", "DivergentPos")) {
  df_all[[paste0("delta_plddt_",      rg)]] <- df_all[[paste0("VAR_plddt_", rg, "_mean")]]              - df_all[[paste0("WT_plddt_", rg, "_mean")]]
  df_all[[paste0("delta_pae_",        rg)]] <- df_all[[paste0("VAR_pae_", rg, "_mean")]]                - df_all[[paste0("WT_pae_", rg, "_mean")]]
  df_all[[paste0("delta_rel_asa_",    rg)]] <- df_all[[paste0("VAR_rel_asa_", rg, "_mean")]]            - df_all[[paste0("WT_rel_asa_", rg, "_mean")]]
  df_all[[paste0("delta_abs_sasa_",   rg)]] <- df_all[[paste0("VAR_abs_sasa_dssp_", rg, "_mean")]]      - df_all[[paste0("WT_abs_sasa_dssp_", rg, "_mean")]]
  df_all[[paste0("delta_total_sasa_", rg)]] <- df_all[[paste0("VAR_SASA_total_", rg)]]                  - df_all[[paste0("WT_SASA_total_", rg)]]
  df_all[[paste0("delta_helix_",      rg)]] <- df_all[[paste0("VAR_sec_struc_", rg, "_helix_percent")]] - df_all[[paste0("WT_sec_struc_", rg, "_helix_percent")]]
  df_all[[paste0("delta_beta_",       rg)]] <- df_all[[paste0("VAR_sec_struc_", rg, "_beta_percent")]]  - df_all[[paste0("WT_sec_struc_", rg, "_beta_percent")]]
  df_all[[paste0("delta_coil_",       rg)]] <- df_all[[paste0("VAR_sec_struc_", rg, "_coil_percent")]]  - df_all[[paste0("WT_sec_struc_", rg, "_coil_percent")]]
}


# ══════════════════════════════════════════════════════════════
# 2. Outcome variables + readable labels
# ══════════════════════════════════════════════════════════════
# Excludes WT_* / wt_length (constant within gene, no group effect identifiable)

region_label <- c(full = "Full", NMDPos = "NMD", DivergentPos = "Divergent")

feature_label <- c(
  plddt      = "pLDDT",
  pae        = "PAE",
  rel_asa    = "Relative ASA",
  abs_sasa   = "Absolute SASA",
  total_sasa = "Total SASA",
  helix      = "Helix %",
  beta       = "Beta %",
  coil       = "Coil %"
)

var_map <- bind_rows(
  tibble(
    outcome  = c("var_length", "truncation_length", "truncation_pct"),
    label    = c("VAR Length (aa)", "Truncation Length (aa)", "Truncation (%)"),
    category = "Length / Truncation",
    ord      = 1:3
  ),
  expand_grid(rg = names(region_label), ft = names(feature_label)) %>%
    mutate(
      outcome = case_when(
        ft == "total_sasa" ~ paste0("VAR_SASA_total_", rg),
        ft == "abs_sasa"   ~ paste0("VAR_abs_sasa_dssp_", rg, "_mean"),
        ft %in% c("helix", "beta", "coil") ~
          paste0("VAR_sec_struc_", rg, "_", ft, "_percent"),
        TRUE ~ paste0("VAR_", ft, "_", rg, "_mean")
      ),
      label    = feature_label[ft],
      category = paste0("VAR — ", region_label[rg]),
      ord      = match(ft, names(feature_label))
    ) %>%
    select(outcome, label, category, ord),
  expand_grid(rg = names(region_label), ft = names(feature_label)) %>%
    mutate(
      outcome  = paste0("delta_", ft, "_", rg),
      label    = paste0("Δ ", feature_label[ft]),
      category = paste0("Δ (VAR−WT) — ", region_label[rg]),
      ord      = match(ft, names(feature_label))
    ) %>%
    select(outcome, label, category, ord)
)

cat_levels <- c(
  "Length / Truncation",
  paste0("VAR — ",        c("Full", "NMD", "Divergent")),
  paste0("Δ (VAR−WT) — ", c("Full", "NMD", "Divergent"))
)


# ══════════════════════════════════════════════════════════════
# 3. Mixed model: value ~ group + (1 | uniprot_id)
# ══════════════════════════════════════════════════════════════

fit_lmm <- function(data, outcome, g_dis, g_con, family_label) {
  
  d <- data %>%
    filter(group %in% c(g_dis, g_con)) %>%
    select(uniprot_id, group, value = all_of(outcome)) %>%
    filter(!is.na(value)) %>%
    mutate(group = factor(group, levels = c(g_con, g_dis)))
  
  m <- tryCatch(
    lmerTest::lmer(value ~ group + (1 | uniprot_id), data = d),
    error = function(e) NULL
  )
  
  if (is.null(m)) {
    return(tibble(family = family_label, outcome = outcome,
                  estimate = NA_real_, se = NA_real_,
                  ci_lo = NA_real_, ci_hi = NA_real_, p_raw = NA_real_,
                  sd_value = NA_real_, n_obs = nrow(d),
                  n_proteins = n_distinct(d$uniprot_id)))
  }
  
  cf  <- coef(summary(m))
  row <- grep("^group", rownames(cf))[1]
  est <- cf[row, "Estimate"]
  se  <- cf[row, "Std. Error"]
  
  tibble(
    family     = family_label,
    outcome    = outcome,
    estimate   = est,
    se         = se,
    ci_lo      = est - 1.96 * se,
    ci_hi      = est + 1.96 * se,
    p_raw      = cf[row, "Pr(>|t|)"],
    sd_value   = sd(d$value),
    n_obs      = nrow(d),
    n_proteins = n_distinct(d$uniprot_id)
  )
}

all_outcomes <- var_map$outcome

lmm_results <- bind_rows(
  map_dfr(all_outcomes, ~ fit_lmm(df_all, .x, "SNV Disease", "SNV Control", "SNV")),
  map_dfr(all_outcomes, ~ fit_lmm(df_all, .x, "FS Disease",  "FS Control",  "FS"))
) %>%
  group_by(family) %>%
  mutate(q_fdr = p.adjust(p_raw, method = "BH")) %>%   # correct independently within each comparison
  ungroup() %>%
  left_join(var_map, by = "outcome") %>%
  mutate(
    category = factor(category, levels = cat_levels),
    family   = factor(family, levels = c("SNV", "FS")),
    std_est  = estimate / sd_value,
    std_lo   = ci_lo    / sd_value,
    std_hi   = ci_hi    / sd_value,
    sig      = case_when(
      is.na(q_fdr)  ~ "",
      q_fdr < 0.001 ~ "***",
      q_fdr < 0.01  ~ "**",
      q_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

write_csv(lmm_results, "LMM_results_all.csv")


# ══════════════════════════════════════════════════════════════
# 4. Combined table: SNV and FS side by side
# ══════════════════════════════════════════════════════════════

tbl_wide <- lmm_results %>%
  mutate(ci = sprintf("(%s, %s)", signif(ci_lo, 3), signif(ci_hi, 3))) %>%
  select(category, ord, label, family, estimate, ci, p_raw, q_fdr, sig) %>%
  pivot_wider(
    names_from  = family,
    values_from = c(estimate, ci, p_raw, q_fdr, sig),
    names_sep   = "_"
  ) %>%
  arrange(category, ord) %>%
  select(
    category, label,
    estimate_SNV, ci_SNV, p_raw_SNV, q_fdr_SNV, sig_SNV,
    estimate_FS,  ci_FS,  p_raw_FS,  q_fdr_FS,  sig_FS
  )

# Sample sizes go into footnote
n_info <- lmm_results %>%
  group_by(family) %>%
  summarise(n_obs = max(n_obs), n_prot = max(n_proteins), .groups = "drop")

fmt_p <- function(x) {
  ifelse(is.na(x), "—",
         ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                formatC(x, format = "f", digits = 3)))
}

gt_combined <- tbl_wide %>%
  gt(groupname_col = "category") %>%
  fmt_number(columns = c(estimate_SNV, estimate_FS), n_sigfig = 3) %>%
  fmt(columns = c(p_raw_SNV, q_fdr_SNV, p_raw_FS, q_fdr_FS), fns = fmt_p) %>%
  cols_label(
    label        = "Variable",
    estimate_SNV = "Estimate", ci_SNV = "95% CI",
    p_raw_SNV    = "p",        q_fdr_SNV = "q (FDR)", sig_SNV = "",
    estimate_FS  = "Estimate", ci_FS  = "95% CI",
    p_raw_FS     = "p",        q_fdr_FS  = "q (FDR)", sig_FS  = ""
  ) %>%
  tab_spanner(
    label   = md(sprintf("**SNV** (n = %s variants, %s genes)",
                         n_info$n_obs[n_info$family == "SNV"],
                         n_info$n_prot[n_info$family == "SNV"])),
    columns = c(estimate_SNV, ci_SNV, p_raw_SNV, q_fdr_SNV, sig_SNV)
  ) %>%
  tab_spanner(
    label   = md(sprintf("**FS** (n = %s variants, %s genes)",
                         n_info$n_obs[n_info$family == "FS"],
                         n_info$n_prot[n_info$family == "FS"])),
    columns = c(estimate_FS, ci_FS, p_raw_FS, q_fdr_FS, sig_FS)
  ) %>%
  tab_header(
    title    = "Mixed-effects models: Disease vs Control",
    subtitle = "Estimate = Disease − Control, original units"
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = c(estimate_SNV, ci_SNV, p_raw_SNV,
                                       q_fdr_SNV, sig_SNV),
                           rows = sig_SNV != "")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = c(estimate_FS, ci_FS, p_raw_FS,
                                       q_fdr_FS, sig_FS),
                           rows = sig_FS != "")
  ) %>%
  tab_style(
    style     = cell_borders(sides = "left", color = "gray70", weight = px(1.5)),
    locations = cells_body(columns = estimate_FS)
  ) %>%
  tab_source_note(
    "Linear mixed-effects model: value ~ group + (1 | gene). Random intercept
     for gene (UniProt ID) accounts for multiple variants per gene. P-values:
     Satterthwaite approximation (lmerTest); q: Benjamini-Hochberg FDR applied
     separately within the SNV and FS comparisons. * q<0.05, ** q<0.01,
     *** q<0.001. WT features are constant within gene and were not modeled."
  )

gt_combined
gtsave(gt_combined, "LMM_table_combined.html")
try(gtsave(gt_combined, "LMM_table_combined.docx"), silent = TRUE)

# ══════════════════════════════════════════════════════════════
# 5. Forest plot: Delta (VAR-WT) only, Full and Divergent only
# ══════════════════════════════════════════════════════════════

plot_dat <- lmm_results %>%
  filter(
    !is.na(estimate),
    str_starts(outcome, "delta_"),                      # VAR-WT only
    str_detect(outcome, "_(full|DivergentPos)$")        # Full / Divergent only
  ) %>%
  mutate(category = droplevels(category)) %>%           # drop empty facets
  arrange(category, desc(ord)) %>%
  mutate(y_key = factor(paste(category, label, sep = "___"),
                        levels = unique(paste(category, label, sep = "___"))))

p_combined <- ggplot(plot_dat,
                     aes(x = std_est, y = y_key, color = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = std_lo, xmax = std_hi),
                 height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = q_fdr < 0.05), size = 2.2,
             position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig), size = 3,
            position = position_dodge(width = 0.6),
            hjust = -0.5, vjust = 0.8, show.legend = FALSE) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_color_manual(values = c(SNV = "#2C7FB8", FS = "#C0392B"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  labs(
    title    = "Disease vs Control: variant - wild type structural changes",
    subtitle = "Estimate = Disease − Control, standardized by outcome SD; bars = 95% CI",
    x        = "Standardized estimate (Disease − Control)",
    y        = NULL,
    caption  = paste0(
      "Linear mixed-effects model: value ~ group + (1 | gene). ",
      "Filled points BH-adjusted p < 0.05; * p<0.05  ** p<0.01  *** p<0.001 "
      
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, face = "bold", size = 9),
    strip.background   = element_rect(fill = "gray95", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "top",
    legend.box         = "horizontal",
    plot.title         = element_text(face = "bold"),
    plot.caption       = element_text(hjust = 0, size = 8, color = "gray40")
  )

print(p_combined)

ggsave("LMM_forest_deltaVARWT_full_divergent.pdf", p_combined,
       width = 9, height = 5.5, dpi = 300)
ggsave("LMM_forest_deltaVARWT_full_divergent.png", p_combined,
       width = 9, height = 5.5, dpi = 300)