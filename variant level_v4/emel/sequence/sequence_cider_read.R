library(tidyverse)
library(ggpubr)
library(readr)

NMD_region_DivergentPos_CIDER <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/CIDER/NMD_region_DivergentPos_CIDER.csv")
str(NMD_region_DivergentPos_CIDER)
NMD_region_NMDPos_CIDER <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/CIDER/NMD_region_NMDPos_CIDER.csv")
str(NMD_region_NMDPos_CIDER)



# ── Colors & comparisons ──────────────────────────────────────────────────────
comparisons <- list(
  c("snv_control", "snv_disease"),
  c("fs_control",  "fs_disease")
)

group_colors <- c(
  snv_control = "#7BAFD4",
  snv_disease = "#E07B7B",
  fs_control  = "#85C985",
  fs_disease  = "#E0A857"
)

# ── Add group column ──────────────────────────────────────────────────────────
NMD_region_DivergentPos_CIDER <- NMD_region_DivergentPos_CIDER %>%
  mutate(group = source_folder) %>%
  filter(group %in% c("snv_control", "snv_disease", "fs_control", "fs_disease"))

NMD_region_NMDPos_CIDER <- NMD_region_NMDPos_CIDER %>%
  mutate(group = source_folder) %>%
  filter(group %in% c("snv_control", "snv_disease", "fs_control", "fs_disease"))

# ── Variables & labels ────────────────────────────────────────────────────────
cider_vars <- c(
  "fcr_7",
  "isoelectric_point",
  "molecular_weight",
  "count_neg",
  "count_pos",
  "count_neut",
  "fraction_negative",
  "fraction_positive",
  "fraction_expanding_7",
  "fraction_disorder_promoting",
  "mean_hydropathy",
  "uversky_hydropathy",
  "kappa",
  "Omega",
  "delta",
  "deltaMax"
)

labels <- c(
  fcr_7                       = "FCR (pH 7)",
  isoelectric_point           = "Isoelectric Point",
  molecular_weight            = "Molecular Weight",
  count_neg                   = "Count Negative",
  count_pos                   = "Count Positive",
  count_neut                  = "Count Neutral",
  fraction_negative           = "Fraction Negative",
  fraction_positive           = "Fraction Positive",
  fraction_expanding_7        = "Fraction Expanding (pH 7)",
  fraction_disorder_promoting = "Fraction Disorder Promoting",
  mean_hydropathy             = "Mean Hydropathy",
  uversky_hydropathy          = "Uversky Hydropathy",
  kappa                       = "Kappa (κ)",
  Omega                       = "Omega (Ω)",
  delta                       = "Delta (δ)",
  deltaMax                    = "Delta Max (δmax)"
)

# ── Panel functions ───────────────────────────────────────────────────────────
make_panel <- function(df, y_var, y_label) {
  df %>%
    filter(!is.na(.data[[y_var]])) %>%
    ggplot(aes(x = group, y = .data[[y_var]], fill = group)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", tip.length = 0.01, size = 3) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = y_label, title = y_label) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 35, hjust = 1),
          plot.title      = element_text(size = 9, face = "bold"))
}

make_delta_panel <- function(df, v, y_label) {
  df %>%
    mutate(delta_val = .data[[paste0("var_", v)]] - .data[[paste0("wt_", v)]]) %>%
    filter(!is.na(delta_val)) %>%
    ggplot(aes(x = group, y = delta_val, fill = group)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", tip.length = 0.01, size = 3) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = paste0("Δ ", y_label), title = paste0("Δ ", y_label)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 35, hjust = 1),
          plot.title      = element_text(size = 9, face = "bold"))
}

# ══════════════════════════════════════════════════════════════════════════════
# DIVERGENT REGION — VAR properties
# ══════════════════════════════════════════════════════════════════════════════
panels_div_var <- lapply(cider_vars, function(v) {
  make_panel(NMD_region_DivergentPos_CIDER, paste0("var_", v), labels[v])
})

ggarrange(plotlist = panels_div_var, ncol = 4, nrow = 4) %>%
  annotate_figure(top = text_grob(
    "Divergent Region — Variant CIDER Properties",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# DIVERGENT REGION — DELTA (Var − WT)
# ══════════════════════════════════════════════════════════════════════════════
panels_div_delta <- lapply(cider_vars, function(v) {
  make_delta_panel(NMD_region_DivergentPos_CIDER, v, labels[v])
})

ggarrange(plotlist = panels_div_delta, ncol = 4, nrow = 4) %>%
  annotate_figure(top = text_grob(
    "Divergent Region — Δ CIDER Properties (Var − WT)",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# NMD REGION — VAR properties
# ══════════════════════════════════════════════════════════════════════════════
panels_nmd_var <- lapply(cider_vars, function(v) {
  make_panel(NMD_region_NMDPos_CIDER, paste0("var_", v), labels[v])
})

ggarrange(plotlist = panels_nmd_var, ncol = 4, nrow = 4) %>%
  annotate_figure(top = text_grob(
    "NMD Region — Variant CIDER Properties",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# NMD REGION — DELTA (Var − WT)
# ══════════════════════════════════════════════════════════════════════════════
panels_nmd_delta <- lapply(cider_vars, function(v) {
  make_delta_panel(NMD_region_NMDPos_CIDER, v, labels[v])
})

ggarrange(plotlist = panels_nmd_delta, ncol = 4, nrow = 4) %>%
  annotate_figure(top = text_grob(
    "NMD Region — Δ CIDER Properties (Var − WT)",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# OPTIONAL: Save all 4 plots to PDF
# ══════════════════════════════════════════════════════════════════════════════
# pdf("CIDER_combined.pdf", width = 20, height = 18)
#
# ggarrange(plotlist = panels_div_var, ncol = 4, nrow = 4) %>%
#   annotate_figure(top = text_grob("Divergent Region — Variant CIDER Properties",
#                                    face = "bold", size = 14)) %>% print()
#
# ggarrange(plotlist = panels_div_delta, ncol = 4, nrow = 4) %>%
#   annotate_figure(top = text_grob("Divergent Region — Δ CIDER Properties (Var − WT)",
#                                    face = "bold", size = 14)) %>% print()
#
# ggarrange(plotlist = panels_nmd_var, ncol = 4, nrow = 4) %>%
#   annotate_figure(top = text_grob("NMD Region — Variant CIDER Properties",
#                                    face = "bold", size = 14)) %>% print()
#
# ggarrange(plotlist = panels_nmd_delta, ncol = 4, nrow = 4) %>%
#   annotate_figure(top = text_grob("NMD Region — Δ CIDER Properties (Var − WT)",
#                                    face = "bold", size = 14)) %>% print()
#
# dev.off()