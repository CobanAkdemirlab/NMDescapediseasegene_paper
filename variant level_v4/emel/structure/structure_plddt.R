library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ══════════════════════════════════════════════════════════════
# 1. Cohort filtering (matches Table 1 / Table 2)
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
# 2. Plot data prep (WT deduplicated at protein level)
# ══════════════════════════════════════════════════════════════

# panel_type = "VAR" -> variant level, one row per variant
# panel_type = "WT"  -> protein level, one row per protein per group
#   WT features are identical across all variants of a protein; skipping dedup causes pseudoreplication:
#     a gene with 200 variants would be counted 200 times in the distribution

prepare_plddt_data <- function(data, plddt_col, panel_type) {
  
  d <- data %>%
    select(group, uniprot_id, pLDDT = all_of(plddt_col)) %>%
    filter(!is.na(pLDDT))
  
  if (panel_type == "WT") {
    d <- d %>% distinct(group, uniprot_id, .keep_all = TRUE)
  }
  
  d
}


# ══════════════════════════════════════════════════════════════
# 3. p-values and BH-FDR correction
# ══════════════════════════════════════════════════════════════

panel_spec <- tibble::tribble(
  ~panel_id,   ~panel_type, ~region,          ~col,                           ~title,
  "var_full",  "VAR",       "Full",           "VAR_plddt_full_mean",          "VAR pLDDT — Full",
  "var_nmd",   "VAR",       "NMD",            "VAR_plddt_NMDPos_mean",        "VAR pLDDT — NMD",
  "var_div",   "VAR",       "Divergent",      "VAR_plddt_DivergentPos_mean",  "VAR pLDDT — Divergent",
  "wt_full",   "WT",        "Full",           "WT_plddt_full_mean",           "WT pLDDT — Full",
  "wt_nmd",    "WT",        "NMD",            "WT_plddt_NMDPos_mean",         "WT pLDDT — NMD",
  "wt_div",    "WT",        "Divergent",      "WT_plddt_DivergentPos_mean",   "WT pLDDT — Divergent"
)

# Two Wilcoxon tests per panel (SNV pair, FS pair)
p_raw <- pmap_dfr(
  list(panel_spec$panel_id, panel_spec$col, panel_spec$panel_type),
  function(pid, col, ptype) {
    
    d <- prepare_plddt_data(df_all, col, ptype)
    
    test_pair <- function(g1, g2) {
      dd <- d %>%
        filter(group %in% c(g1, g2)) %>%
        mutate(group = droplevels(group))
      tryCatch(
        wilcox.test(pLDDT ~ group, data = dd)$p.value,
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

# ── FDR correction ──────────────────────────────────────────────────
# Default: correct within this figure's 12 tests (SNV and FS families separate)
# Keeps figure and table significance consistent (avoids stars vs non-significant q mismatch)

p_annot <- p_raw %>%
  group_by(family) %>%
  mutate(q_val = p.adjust(p_raw, method = "BH")) %>%
  ungroup()

# --- Reuse Table 1's q-values (optional) ---------------------------------
# t1_lookup <- tibble::tribble(
#   ~panel_id,  ~t1_variable,
#   "var_full", "VAR_plddt_full_mean",
#   "var_nmd",  "VAR_plddt_NMDPos_mean",
#   "var_div",  "VAR_plddt_DivergentPos_mean",
#   "wt_full",  "WT_plddt_full_mean",
#   "wt_nmd",   "WT_plddt_NMDPos_mean",
#   "wt_div",   "WT_plddt_DivergentPos_mean"
# )
# p_annot <- p_raw %>%
#   left_join(t1_lookup, by = "panel_id") %>%
#   left_join(p_all, by = c("t1_variable" = "variable")) %>%
#   mutate(q_val = if_else(family == "SNV", snv_q_raw, fs_q_raw))
# ----------------------------------------------------------------

# Convert q-value to significance stars
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
# 4. Single-panel plotting function
# ══════════════════════════════════════════════════════════════

fill_colors <- c("SNV\nDisease" = "#E07B7B",
                 "SNV\nControl" = "#7BA7D4",
                 "FS\nDisease"  = "#E8B49A",
                 "FS\nControl"  = "#A8C8E8")

plot_plddt <- function(df, title, annot, y_label = "pLDDT score") {
  
  # Bracket height: above data, both brackets equal height
  y_max <- max(df$pLDDT, na.rm = TRUE)
  annot <- annot %>% mutate(y.position = y_max + 6)
  
  ggplot(df, aes(x = group, y = pLDDT, fill = group)) +
    
    geom_violin(alpha = 0.6, trim = FALSE, color = NA) +
    
    geom_boxplot(width = 0.15,
                 outlier.size  = 0.5,
                 outlier.alpha = 0.3,
                 color = "gray30",
                 fill  = "white",
                 alpha = 0.8) +
    
    # Significance brackets: using FDR-corrected q-values
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
    
    geom_hline(yintercept = 90, linetype = "dotted",
               color = "#2166ac", linewidth = 0.4) +
    geom_hline(yintercept = 70, linetype = "dotted",
               color = "#4dac26", linewidth = 0.4) +
    geom_hline(yintercept = 50, linetype = "dotted",
               color = "#f1a340", linewidth = 0.4) +
    
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(limits = c(0, 108),
                       breaks = c(0, 25, 50, 70, 90, 100)) +
    
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
# 5. Generate six plots
# ══════════════════════════════════════════════════════════════

plots <- pmap(
  list(panel_spec$panel_id, panel_spec$col,
       panel_spec$panel_type, panel_spec$title),
  function(pid, col, ptype, ttl) {
    plot_plddt(
      df    = prepare_plddt_data(df_all, col, ptype),
      title = ttl,
      annot = p_annot %>%
        filter(panel_id == pid) %>%
        select(group1, group2, label)
    )
  }
) %>% set_names(panel_spec$panel_id)


# ══════════════════════════════════════════════════════════════
# 6. Combine plots and export
# ══════════════════════════════════════════════════════════════

n_wt_proteins <- df_all %>%
  distinct(group, uniprot_id) %>%
  nrow()

combined <-
  (plots$var_full | plots$var_nmd | plots$var_div) /
  (plots$wt_full  | plots$wt_nmd  | plots$wt_div) +
  plot_annotation(
    title    = "pLDDT Comparison: SNV Disease vs Control  |  FS Disease vs Control",
    subtitle = "Top row: VAR protein (variant-level)  |  Bottom row: WT protein (protein-level, deduplicated)",
    caption  = paste0(
      "pLDDT confidence thresholds: >90 (blue dotted) very high  |  ",
      "70-90 (green dotted) confident  |  50-70 (orange dotted) low\n",
      "Significance: Wilcoxon rank-sum test, Benjamini-Hochberg FDR-adjusted ",
      "across the 12 tests in this figure (SNV and FS families corrected separately).\n",
      "* q<0.05  ** q<0.01  *** q<0.001  **** q<0.0001  ns = not significant\n",
      "WT panels are deduplicated to one value per protein per group; ",
      "VAR panels are variant-level. Dashed line separates SNV and FS groups."
    ),
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      plot.caption  = element_text(hjust = 0, size = 8, color = "gray50")
    )
  )

print(combined)

ggsave("pLDDT_comparison_fdr.pdf", plot = combined,
       width = 14, height = 9.5, dpi = 300)
ggsave("pLDDT_comparison_fdr.png", plot = combined,
       width = 14, height = 9.5, dpi = 300)


# ── Check: print all p / q values for comparison with table ──────────────────

p_annot %>%
  select(panel_id, family, p_raw, q_val, label) %>%
  mutate(across(c(p_raw, q_val), ~ signif(.x, 3))) %>%
  print(n = Inf)

# ── Check: actual sample size per panel ──────────────────────────────

pmap_dfr(
  list(panel_spec$panel_id, panel_spec$col, panel_spec$panel_type),
  function(pid, col, ptype) {
    prepare_plddt_data(df_all, col, ptype) %>%
      count(group) %>%
      mutate(panel_id = pid, .before = 1)
  }
) %>%
  pivot_wider(names_from = group, values_from = n) %>%
  print()