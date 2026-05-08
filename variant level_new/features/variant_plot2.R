# =============================================================================
# Variant-level enrichment barplot figure
# Similar style to the example figure:
#   - percentage bars
#   - n labels inside bars
#   - p-value brackets above pairs
#   - 6-panel supplemental-style figure
# =============================================================================

rm(list = ls())
gc()

# ── 1. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("tidyverse", "patchwork", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# ── 2. Read data ──────────────────────────────────────────────────────────────
df <- read_csv("Downloads/variants_all_with_motif_LCS_flags (1).csv", show_col_types = FALSE)

# ── 3. Basic cleaning ─────────────────────────────────────────────────────────
df <- df %>%
  filter(source %in% c("fs", "fs_control", "snv", "snv_control")) %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control"))
  )

# Convert logical/character binary columns to numeric 0/1 if needed
flag_cols <- c(
  "variant_ppi_overlap",
  "ptc_after_max_pfam_end",
  "ptc_before_max_pfam_end",
  "variant_protein_flag",
  "variant_domain_flag",
  "variant_slim_flag",
  "variant_morf_flag",
  "variant_ptm_flag",
  "variant_nls_flag",
  "variant_LCS_flag"
)

for (cc in flag_cols) {
  if (cc %in% colnames(df)) {
    if (is.logical(df[[cc]])) df[[cc]] <- as.numeric(df[[cc]])
    if (is.character(df[[cc]])) {
      df[[cc]] <- case_when(
        df[[cc]] %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
        df[[cc]] %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
        TRUE ~ suppressWarnings(as.numeric(df[[cc]]))
      )
    }
  }
}

# ── 4. Rename groups for nicer labels ─────────────────────────────────────────
group_labels <- c(
  "fs" = "Frameshift",
  "fs_control" = "FS Control",
  "snv" = "Nonsense",
  "snv_control" = "Nonsense Control"
)

df <- df %>%
  mutate(
    group_label = recode(as.character(source), !!!group_labels),
    group_label = factor(
      group_label,
      levels = c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")
    )
  )

# ── 5. Colors ─────────────────────────────────────────────────────────────────
group_colors <- c(
  "Frameshift"       = "#D07A3A",
  "FS Control"       = "#E2C7B5",
  "Nonsense"         = "#8ABA67",
  "Nonsense Control" = "#C9D6BE"
)

# Optional alternative palette:
# group_colors <- c(
#   "Frameshift"       = "#2166AC",
#   "FS Control"       = "#92C5DE",
#   "Nonsense"         = "#B2182B",
#   "Nonsense Control" = "#F4A582"
# )

# ── 6. Features to draw ───────────────────────────────────────────────────────
# Panel B: choose ONE Pfam-related feature.
# Here I use "ptc_before_max_pfam_end" as a domain-related enrichment measure.
# You can switch to "ptc_after_max_pfam_end" if that fits your intended figure better.

feature_map <- tribble(
  ~feature,                  ~title,                                      ~ylab,                              ~panel,
  "variant_ppi_overlap",     "PPI Residue Enrichment",                    "Percentage with PPI (%)",          "A",
  "ptc_before_max_pfam_end", "Pfam Domain Enrichment",                    "Percentage with Domains (%)",      "B",
  "variant_nls_flag",        "Nuclear Localization Signal Enrichment",    "Percentage with NLS (%)",          "C",
  "variant_ptm_flag",        "Post-Translational Modification Enrichment","Percentage with PTM (%)",          "D",
  "variant_slim_flag",       "Short Linear Motif Enrichment",             "Percentage with SLiM (%)",         "E",
  "variant_LCS_flag",        "Low Complexity Sequence Enrichment",        "Percentage with LCS (%)",          "F"
)

feature_map <- feature_map %>%
  filter(feature %in% colnames(df))

if (nrow(feature_map) == 0) {
  stop("None of the requested feature columns were found in the CSV.")
}

# ── 7. Helper functions ───────────────────────────────────────────────────────

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.2f", p)
}

summarise_feature <- function(data, feature) {
  data %>%
    filter(!is.na(.data[[feature]])) %>%
    group_by(group_label) %>%
    summarise(
      total_n = n(),
      positive_n = sum(.data[[feature]] == 1, na.rm = TRUE),
      percent = 100 * positive_n / total_n,
      .groups = "drop"
    )
}

get_pairwise_p <- function(data, feature, g1, g2) {
  dsub <- data %>%
    filter(group_label %in% c(g1, g2)) %>%
    filter(!is.na(.data[[feature]]))
  
  tab <- dsub %>%
    group_by(group_label) %>%
    summarise(
      pos = sum(.data[[feature]] == 1, na.rm = TRUE),
      neg = sum(.data[[feature]] == 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(tab) < 2) return(NA_real_)
  
  mat <- matrix(c(tab$pos, tab$neg), nrow = 2, byrow = FALSE)
  
  p <- tryCatch(
    fisher.test(mat)$p.value,
    error = function(e) NA_real_
  )
  
  p
}

make_bracket_df <- function(summary_df, p1, p2) {
  ymax <- max(summary_df$percent, na.rm = TRUE)
  yrange <- max(5, ymax * 0.12)
  
  tibble(
    xstart = c(1, 3),
    xend   = c(2, 4),
    xmid   = c(1.5, 3.5),
    y      = c(ymax + yrange * 0.55, ymax + yrange * 0.55),
    label  = c(format_p(p1), format_p(p2))
  )
}

plot_enrichment_bar <- function(data, feature, title, ylab, panel_letter) {
  
  summary_df <- summarise_feature(data, feature)
  
  summary_df <- tibble(
    group_label = factor(
      c("Frameshift", "FS Control", "Nonsense", "Nonsense Control"),
      levels = c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")
    )
  ) %>%
    left_join(summary_df, by = "group_label") %>%
    mutate(
      total_n = replace_na(total_n, 0),
      positive_n = replace_na(positive_n, 0),
      percent = replace_na(percent, 0)
    )
  
  p_fs  <- get_pairwise_p(data, feature, "Frameshift", "FS Control")
  p_snv <- get_pairwise_p(data, feature, "Nonsense", "Nonsense Control")
  
  bracket_df <- make_bracket_df(summary_df, p_fs, p_snv)
  
  ymax <- max(summary_df$percent, na.rm = TRUE)
  ylim_top <- max(100, ymax + max(10, ymax * 0.25))
  
  ggplot(summary_df, aes(x = group_label, y = percent, fill = group_label)) +
    geom_col(width = 0.62, color = "grey45", alpha = 0.95) +
    
    geom_text(
      aes(label = sprintf("%.1f%%", percent)),
      vjust = -0.55,
      size = 3.5,
      fontface = "bold",
      color = "grey20"
    ) +
    
    geom_label(
      aes(
        y = pmax(percent * 0.48, 6),
        label = paste0("n=", positive_n)
      ),
      size = 3.0,
      label.size = 0,
      fill = "white",
      alpha = 0.75
    ) +
    
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xstart, xend = xstart, y = y - 1.2, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = xend, xend = xend, y = y - 1.2, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "grey25"
    ) +
    geom_text(
      data = bracket_df,
      aes(x = xmid, y = y + 2, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "grey20"
    ) +
    
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(
      limits = c(0, ylim_top),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      x = "Variant Category",
      y = ylab
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 12, 12, 12)
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = panel_letter,
      hjust = -0.8, vjust = 1.5,
      size = 20 / .pt,
      fontface = "bold"
    )
}

# ── 8. Generate all panels ────────────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[ feature_map$panel[i] ]] <- plot_enrichment_bar(
    data = df,
    feature = feature_map$feature[i],
    title = feature_map$title[i],
    ylab = feature_map$ylab[i],
    panel_letter = feature_map$panel[i]
  )
}

# ── 9. Arrange panels ─────────────────────────────────────────────────────────
ordered_panels <- c("A", "B", "C", "D", "E", "F")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Variant-Level Feature Enrichment",
    theme = theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0)
    )
  )

# ── 10. Save outputs ──────────────────────────────────────────────────────────
ggsave("Supplemental_Figure_3_variant_level.pdf", combined_fig,
       width = 17, height = 10, dpi = 300)

ggsave("Supplemental_Figure_3_variant_level.png", combined_fig,
       width = 17, height = 10, dpi = 300)

# ── 11. Optional summary table ────────────────────────────────────────────────
summary_table <- map_dfr(feature_map$feature, function(feat) {
  title_now <- feature_map$title[feature_map$feature == feat][1]
  
  summarise_feature(df, feat) %>%
    mutate(feature = feat, title = title_now, .before = 1)
})

write_csv(summary_table, "Supplemental_Figure_3_variant_level_summary.csv")

cat("\nDone!\n")
cat("Output files:\n")
cat("  Supplemental_Figure_3_variant_level.pdf\n")
cat("  Supplemental_Figure_3_variant_level.png\n")
cat("  Supplemental_Figure_3_variant_level_summary.csv\n")