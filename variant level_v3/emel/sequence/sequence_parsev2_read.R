#parse v2
full_length_PARSE_v2 <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/parse_v2/full_length_PARSE_v2.csv")
str(full_length_PARSE_v2)

nmd_length_DivergentPos_PARSE_v2 <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/parse_v2/nmd_length_DivergentPos_PARSE_v2.csv")
str(nmd_length_DivergentPos_PARSE_v2)

nmd_length_NMDPos_PARSE_v2 <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/parse_v2/nmd_length_NMDPos_PARSE_v2.csv")
str(nmd_length_NMDPos_PARSE_v2)

library(tidyverse)
library(ggpubr)

# ── 1. Clean column names & add group ────────────────────────────────────────

full_length_PARSE_v2 <- full_length_PARSE_v2 %>%
  rename(
    sigma_P_var = `Σ classifier distance P...4`,
    longest_PS_IDR_var = `length of longest PS IDR...5`,
    sigma_P_wt = `Σ classifier distance P...9`,
    longest_PS_IDR_wt = `length of longest PS IDR...10`
  ) %>%
  mutate(group = str_extract(id, "^(snv_control|snv_disease|fs_control|fs_disease|fs|snv)")) %>%
  filter(!is.na(group))

nmd_length_DivergentPos_PARSE_v2 <- nmd_length_DivergentPos_PARSE_v2 %>%
  rename(
    sigma_P_var = `Σ classifier distance P...4`,
    longest_PS_IDR_var = `length of longest PS IDR...5`,
    sigma_P_wt = `Σ classifier distance P...9`,
    longest_PS_IDR_wt = `length of longest PS IDR...10`
  ) %>%
  mutate(group = str_extract(id, "^(snv_control|snv_disease|fs_control|fs_disease|fs|snv)")) %>%
  filter(!is.na(group))

nmd_length_NMDPos_PARSE_v2 <- nmd_length_NMDPos_PARSE_v2 %>%
  rename(
    sigma_P_var = `Σ classifier distance P...4`,
    longest_PS_IDR_var = `length of longest PS IDR...5`,
    sigma_P_wt = `Σ classifier distance P...9`,
    longest_PS_IDR_wt = `length of longest PS IDR...10`
  ) %>%
  mutate(group = str_extract(id, "^(snv_control|snv_disease|fs_control|fs_disease|fs|snv)")) %>%
  filter(!is.na(group))

# Check groups
table(full_length_PARSE_v2$group)
table(nmd_length_DivergentPos_PARSE_v2$group)
table(nmd_length_NMDPos_PARSE_v2$group)

# ── Define comparisons ───────────────────────────────────────────────────────
comparisons <- list(
  c("snv_control", "snv_disease"),
  c("fs_control", "fs_disease")
)

# ── Helper plot function ─────────────────────────────────────────────────────
plot_parse <- function(df, y_var, y_label, title) {
  df %>%
    filter(group %in% c("snv_control", "snv_disease", "fs_control", "fs_disease")) %>%
    ggplot(aes(x = group, y = .data[[y_var]], fill = group)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", tip.length = 0.01) +
    stat_compare_means(method = "kruskal.test", label.y.npc = 0.05,
                       label = "p.format") +
    scale_fill_manual(values = c(
      snv_control  = "#7BAFD4",
      snv_disease  = "#E07B7B",
      fs_control   = "#85C985",
      fs_disease   = "#E0A857"
    )) +
    labs(title = title, x = "Group", y = y_label) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 1: full_length_PARSE_v2
# ══════════════════════════════════════════════════════════════════════════════

# 1a. Σ classifier distance P — variant
plot_parse(full_length_PARSE_v2,
           y_var   = "sigma_P_var",
           y_label = "Σ Classifier Distance P (variant)",
           title   = "Full Length — Σ Classifier Distance P (Variant)")

# 1b. Σ classifier distance P — wildtype
plot_parse(full_length_PARSE_v2,
           y_var   = "sigma_P_wt",
           y_label = "Σ Classifier Distance P (wildtype)",
           title   = "Full Length — Σ Classifier Distance P (Wildtype)")

# 1c. Length of longest PS IDR — variant
plot_parse(full_length_PARSE_v2,
           y_var   = "longest_PS_IDR_var",
           y_label = "Length of Longest PS IDR (variant)",
           title   = "Full Length — Longest PS IDR (Variant)")

# 1d. Length of longest PS IDR — wildtype
plot_parse(full_length_PARSE_v2,
           y_var   = "longest_PS_IDR_wt",
           y_label = "Length of Longest PS IDR (wildtype)",
           title   = "Full Length — Longest PS IDR (Wildtype)")

# 1e. Variant length
plot_parse(full_length_PARSE_v2,
           y_var   = "var_length",
           y_label = "Variant Sequence Length",
           title   = "Full Length — Variant Sequence Length")

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 2: nmd_length_DivergentPos_PARSE_v2
# ══════════════════════════════════════════════════════════════════════════════

# 2a. Σ classifier distance P — variant
plot_parse(nmd_length_DivergentPos_PARSE_v2,
           y_var   = "sigma_P_var",
           y_label = "Σ Classifier Distance P (variant)",
           title   = "Divergent Region — Σ Classifier Distance P (Variant)")

# 2b. Σ classifier distance P — wildtype
plot_parse(nmd_length_DivergentPos_PARSE_v2,
           y_var   = "sigma_P_wt",
           y_label = "Σ Classifier Distance P (wildtype)",
           title   = "Divergent Region — Σ Classifier Distance P (Wildtype)")

# 2c. Longest PS IDR — variant
plot_parse(nmd_length_DivergentPos_PARSE_v2,
           y_var   = "longest_PS_IDR_var",
           y_label = "Length of Longest PS IDR (variant)",
           title   = "Divergent Region — Longest PS IDR (Variant)")

# 2d. Longest PS IDR — wildtype
plot_parse(nmd_length_DivergentPos_PARSE_v2,
           y_var   = "longest_PS_IDR_wt",
           y_label = "Length of Longest PS IDR (wildtype)",
           title   = "Divergent Region — Longest PS IDR (Wildtype)")

# 2e. Variant NMD length
plot_parse(nmd_length_DivergentPos_PARSE_v2,
           y_var   = "var_NMD_length",
           y_label = "Variant NMD Length",
           title   = "Divergent Region — Variant NMD Length")

# ══════════════════════════════════════════════════════════════════════════════
# DATASET 3: nmd_length_NMDPos_PARSE_v2
# ══════════════════════════════════════════════════════════════════════════════

# 3a. Σ classifier distance P — variant
plot_parse(nmd_length_NMDPos_PARSE_v2,
           y_var   = "sigma_P_var",
           y_label = "Σ Classifier Distance P (variant)",
           title   = "NMD Region — Σ Classifier Distance P (Variant)")

# 3b. Σ classifier distance P — wildtype
plot_parse(nmd_length_NMDPos_PARSE_v2,
           y_var   = "sigma_P_wt",
           y_label = "Σ Classifier Distance P (wildtype)",
           title   = "NMD Region — Σ Classifier Distance P (Wildtype)")

# 3c. Longest PS IDR — variant
plot_parse(nmd_length_NMDPos_PARSE_v2,
           y_var   = "longest_PS_IDR_var",
           y_label = "Length of Longest PS IDR (variant)",
           title   = "NMD Region — Longest PS IDR (Variant)")

# 3d. Longest PS IDR — wildtype
plot_parse(nmd_length_NMDPos_PARSE_v2,
           y_var   = "longest_PS_IDR_wt",
           y_label = "Length of Longest PS IDR (wildtype)",
           title   = "NMD Region — Longest PS IDR (Wildtype)")

# 3e. Variant NMD length
plot_parse(nmd_length_NMDPos_PARSE_v2,
           y_var   = "var_NMD_length_old",
           y_label = "Variant NMD Length",
           title   = "NMD Region — Variant NMD Length")