NMD_region_DivergentPos_peptide_props <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/peptides/NMD_region_DivergentPos_peptide_props.csv")
str(NMD_region_DivergentPos_peptide_props)
NMD_region_NMDPos_peptide_props <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/peptides/NMD_region_NMDPos_peptide_props.csv")
str(NMD_region_NMDPos_peptide_props)

library(tidyverse)
library(ggpubr)

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
NMD_region_DivergentPos_peptide_props <- NMD_region_DivergentPos_peptide_props %>%
  mutate(group = source_folder) %>%
  filter(group %in% c("snv_control", "snv_disease", "fs_control", "fs_disease"))

NMD_region_NMDPos_peptide_props <- NMD_region_NMDPos_peptide_props %>%
  mutate(group = source_folder) %>%
  filter(group %in% c("snv_control", "snv_disease", "fs_control", "fs_disease"))

# ── Variables & labels ────────────────────────────────────────────────────────
numeric_vars <- c(
  "aliphatic_index", "boman", "charge", "entropy",
  "instability_index", "isoelectric_point", "molecular_weight",
  "mz", "longest_run", "max_frequency", "energy_cost",
  "nutrient_cost", "mass_shift", "hydrophobicity", "hydrophobic_moment"
)

labels <- c(
  aliphatic_index    = "Aliphatic Index",
  boman              = "Boman Index",
  charge             = "Charge",
  entropy            = "Entropy",
  instability_index  = "Instability Index",
  isoelectric_point  = "Isoelectric Point",
  molecular_weight   = "Molecular Weight",
  mz                 = "m/z",
  longest_run        = "Longest Run",
  max_frequency      = "Max Frequency",
  energy_cost        = "Energy Cost",
  nutrient_cost      = "Nutrient Cost",
  mass_shift         = "Mass Shift",
  hydrophobicity     = "Hydrophobicity",
  hydrophobic_moment = "Hydrophobic Moment"
)

# ── Single panel function ─────────────────────────────────────────────────────
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

# ── Delta panel function ──────────────────────────────────────────────────────
make_delta_panel <- function(df, v, y_label) {
  df %>%
    mutate(delta = .data[[paste0("var_", v)]] - .data[[paste0("wt_", v)]]) %>%
    filter(!is.na(delta)) %>%
    ggplot(aes(x = group, y = delta, fill = group)) +
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
# DIVERGENT REGION — VAR properties (one combined plot)
# ══════════════════════════════════════════════════════════════════════════════
panels <- lapply(numeric_vars, function(v) {
  make_panel(NMD_region_DivergentPos_peptide_props, paste0("var_", v), labels[v])
})

ggarrange(plotlist = panels, ncol = 5, nrow = 3) %>%
  annotate_figure(top = text_grob(
    "Divergent Region — Variant Peptide Properties",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# DIVERGENT REGION — DELTA (Var − WT) (one combined plot)
# ══════════════════════════════════════════════════════════════════════════════
panels_delta <- lapply(numeric_vars, function(v) {
  make_delta_panel(NMD_region_DivergentPos_peptide_props, v, labels[v])
})

ggarrange(plotlist = panels_delta, ncol = 5, nrow = 3) %>%
  annotate_figure(top = text_grob(
    "Divergent Region — Δ Peptide Properties (Var − WT)",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# NMD REGION — VAR properties (one combined plot)
# ══════════════════════════════════════════════════════════════════════════════
panels_nmd <- lapply(numeric_vars, function(v) {
  make_panel(NMD_region_NMDPos_peptide_props, paste0("var_", v), labels[v])
})

ggarrange(plotlist = panels_nmd, ncol = 5, nrow = 3) %>%
  annotate_figure(top = text_grob(
    "NMD Region — Variant Peptide Properties",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# NMD REGION — DELTA (Var − WT) (one combined plot)
# ══════════════════════════════════════════════════════════════════════════════
panels_nmd_delta <- lapply(numeric_vars, function(v) {
  make_delta_panel(NMD_region_NMDPos_peptide_props, v, labels[v])
})

ggarrange(plotlist = panels_nmd_delta, ncol = 5, nrow = 3) %>%
  annotate_figure(top = text_grob(
    "NMD Region — Δ Peptide Properties (Var − WT)",
    face = "bold", size = 14))

# ══════════════════════════════════════════════════════════════════════════════
# OPTIONAL: Save all 4 combined plots to PDF
# ══════════════════════════════════════════════════════════════════════════════
# pdf("peptide_props_all.pdf", width = 20, height = 14)
#
# ggarrange(plotlist = panels, ncol = 5, nrow = 3) %>%
#   annotate_figure(top = text_grob("Divergent Region — Variant Peptide Properties",
#                                    face = "bold", size = 14)) %>% print()
#
# ggarrange(plotlist = panels_delta, ncol = 5, nrow = 3) %>%
#   annotate_figure(top = text_grob("Divergent Region — Δ Peptide Properties (Var − WT)",
#                                    face = "bold", size = 14)) %>% print()
#
# ggarrange(plotlist = panels_nmd, ncol = 5, nrow = 3) %>%
#   annotate_figure(top = text_grob("NMD Region — Variant Peptide Properties",
#                                    face = "bold", size = 14)) %>% print()
#
# ggarrange(plotlist = panels_nmd_delta, ncol = 5, nrow = 3) %>%
#   annotate_figure(top = text_grob("NMD Region — Δ Peptide Properties (Var − WT)",
#                                    face = "bold", size = 14)) %>% print()
#
# dev.off()