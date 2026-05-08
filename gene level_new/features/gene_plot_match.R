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

# ── 2. Read input files ───────────────────────────────────────────────────────
gene_all <- read_csv("gene_all0407 (1).csv", show_col_types = FALSE)
match_log <- read_csv("Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/gene/matching_log_v2.csv")


# ── 3. Check required columns ─────────────────────────────────────────────────
req_match_cols <- c("study_gene", "ctrl_gene", "group")
miss_match_cols <- setdiff(req_match_cols, colnames(match_log))
if (length(miss_match_cols) > 0) {
  stop("matching_log.csv is missing columns: ",
       paste(miss_match_cols, collapse = ", "))
}

if (!"hgnc_symbol" %in% colnames(gene_all)) {
  stop("gene_all0407 (1).csv must contain column: hgnc_symbol")
}

# ── 4. Standardize match type ─────────────────────────────────────────────────
# FS = frameshift
# STOPGAIN = nonsense / snv
match_log <- match_log %>%
  mutate(
    group = toupper(group),
    match_type = case_when(
      group == "FS" ~ "Frameshift",
      group == "STOPGAIN" ~ "Nonsense",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_type))

# ── 5. Feature columns ────────────────────────────────────────────────────────
flag_cols <- c(
  "ppi_overlap",
  "pfam_overlap_flag",
  "gene_nls_flag",
  "gene_ptm_flag",
  "gene_slim_flag",
  "gene_morf_flag",
  "gene_protein_flag",
  "gene_domains_flag",
  "gene_LCS_flag"
)

flag_cols <- flag_cols[flag_cols %in% colnames(gene_all)]

if (length(flag_cols) == 0) {
  stop("None of the requested flag columns were found in gene_all0407 (1).csv")
}

# ── 6. Convert feature columns to numeric 0/1 ────────────────────────────────
for (cc in flag_cols) {
  if (is.logical(gene_all[[cc]])) {
    gene_all[[cc]] <- as.numeric(gene_all[[cc]])
  } else if (is.character(gene_all[[cc]])) {
    gene_all[[cc]] <- case_when(
      gene_all[[cc]] %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      gene_all[[cc]] %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(gene_all[[cc]]))
    )
  } else {
    gene_all[[cc]] <- suppressWarnings(as.numeric(gene_all[[cc]]))
  }
}

# ── 7. Collapse to one row per gene ───────────────────────────────────────────
# If a gene appears multiple times, use max() across rows for each binary flag
gene_flags <- gene_all %>%
  select(hgnc_symbol, all_of(flag_cols)) %>%
  group_by(hgnc_symbol) %>%
  summarise(
    across(
      all_of(flag_cols),
      ~ {
        x <- suppressWarnings(as.numeric(.x))
        if (all(is.na(x))) {
          NA_real_
        } else {
          max(x, na.rm = TRUE)
        }
      }
    ),
    .groups = "drop"
  )

# ── 8. Join study/control genes to feature table ─────────────────────────────
study_flags <- gene_flags %>%
  rename(study_gene = hgnc_symbol) %>%
  rename_with(~ paste0(.x, "_study"), -study_gene)

ctrl_flags <- gene_flags %>%
  rename(ctrl_gene = hgnc_symbol) %>%
  rename_with(~ paste0(.x, "_ctrl"), -ctrl_gene)

matched_df <- match_log %>%
  left_join(study_flags, by = "study_gene") %>%
  left_join(ctrl_flags, by = "ctrl_gene")

sapply(matched_df,table,useNA = 'ifany')
# ── 9. Feature map for plotting ───────────────────────────────────────────────
feature_map <- tribble(
  ~feature,              ~title,                                      ~panel,
  "ppi_overlap",         "PPI Residue Enrichment",                    "A",
  "pfam_overlap_flag",   "Pfam Domain Enrichment",                    "B",
  "gene_nls_flag",       "Nuclear Localization Signal Enrichment",    "C",
  "gene_ptm_flag",       "Post-Translational Modification Enrichment","D",
  "gene_slim_flag",      "Short Linear Motif Enrichment",             "E",
  "gene_LCS_flag",       "Low Complexity Sequence Enrichment",        "F"
) %>%
  filter(feature %in% flag_cols)

if (nrow(feature_map) == 0) {
  stop("No features in feature_map were found in gene_all0407 (1).csv")
}

# ── 10. Helper functions ──────────────────────────────────────────────────────
format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

summarise_matched_feature <- function(data, feature, match_type_now) {
  study_col <- paste0(feature, "_study")
  ctrl_col  <- paste0(feature, "_ctrl")
  
  d <- data %>%
    filter(match_type == match_type_now) %>%
    filter(!is.na(.data[[study_col]]), !is.na(.data[[ctrl_col]])) %>%
    mutate(
      study = ifelse(.data[[study_col]] >= 1, 1, 0),
      ctrl  = ifelse(.data[[ctrl_col]]  >= 1, 1, 0)
    )
  
  ctrl_label <- ifelse(match_type_now == "Frameshift", "FS Control", "Nonsense Control")
  
  tibble(
    group_label = c(match_type_now, ctrl_label),
    positive_n = c(sum(d$study == 1), sum(d$ctrl == 1)),
    total_n = c(nrow(d), nrow(d)),
    percent = c(mean(d$study == 1) * 100, mean(d$ctrl == 1) * 100)
  )
}

get_mcnemar_p <- function(data, feature, match_type_now) {
  study_col <- paste0(feature, "_study")
  ctrl_col  <- paste0(feature, "_ctrl")
  
  d <- data %>%
    filter(match_type == match_type_now) %>%
    filter(!is.na(.data[[study_col]]), !is.na(.data[[ctrl_col]])) %>%
    mutate(
      study = ifelse(.data[[study_col]] >= 1, 1, 0),
      ctrl  = ifelse(.data[[ctrl_col]]  >= 1, 1, 0)
    )
  
  if (nrow(d) == 0) return(NA_real_)
  
  tab <- table(
    factor(d$study, levels = c(0, 1)),
    factor(d$ctrl,  levels = c(0, 1))
  )
  
  tryCatch(
    mcnemar.test(tab, correct = TRUE)$p.value,
    error = function(e) NA_real_
  )
}

# ── 11. Plot function ─────────────────────────────────────────────────────────
make_plot <- function(data, feature, title, panel_letter) {
  fs_sum  <- summarise_matched_feature(data, feature, "Frameshift")
  snv_sum <- summarise_matched_feature(data, feature, "Nonsense")
  
  summary_df <- bind_rows(fs_sum, snv_sum) %>%
    mutate(
      group_label = factor(
        group_label,
        levels = c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")
      )
    )
  
  p_fs  <- get_mcnemar_p(data, feature, "Frameshift")
  p_snv <- get_mcnemar_p(data, feature, "Nonsense")
  
  ymax <- max(summary_df$percent, na.rm = TRUE)
  ylim_top <- max(100, ymax + 18)
  
  group_colors <- c(
    "Frameshift"       = "#D07A3A",
    "FS Control"       = "#E2C7B5",
    "Nonsense"         = "#8ABA67",
    "Nonsense Control" = "#C9D6BE"
  )
  
  ggplot(summary_df, aes(x = group_label, y = percent, fill = group_label)) +
    geom_col(width = 0.62, color = "grey40", alpha = 0.95) +
    geom_text(
      aes(label = sprintf("%.1f%%", percent)),
      vjust = -0.55,
      size = 3.6,
      fontface = "bold"
    ) +
    geom_label(
      aes(
        y = pmax(percent * 0.5, 6),
        label = paste0("n=", positive_n, "/", total_n)
      ),
      size = 3.0,
      label.size = 0,
      fill = "white",
      alpha = 0.8
    ) +
    annotate("segment", x = 1, xend = 2, y = ymax + 8, yend = ymax + 8, linewidth = 0.6) +
    annotate("segment", x = 1, xend = 1, y = ymax + 6.5, yend = ymax + 8, linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2, y = ymax + 6.5, yend = ymax + 8, linewidth = 0.6) +
    annotate("text", x = 1.5, y = ymax + 10.5, label = format_p(p_fs), fontface = "bold", size = 3.4) +
    annotate("segment", x = 3, xend = 4, y = ymax + 8, yend = ymax + 8, linewidth = 0.6) +
    annotate("segment", x = 3, xend = 3, y = ymax + 6.5, yend = ymax + 8, linewidth = 0.6) +
    annotate("segment", x = 4, xend = 4, y = ymax + 6.5, yend = ymax + 8, linewidth = 0.6) +
    annotate("text", x = 3.5, y = ymax + 10.5, label = format_p(p_snv), fontface = "bold", size = 3.4) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(
      limits = c(0, ylim_top),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      x = "Matched Gene Category",
      y = "Percentage (%)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(size = 11),
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

# ── 12. Generate all panels ───────────────────────────────────────────────────
plot_list <- list()

for (i in seq_len(nrow(feature_map))) {
  plot_list[[feature_map$panel[i]]] <- make_plot(
    data = matched_df,
    feature = feature_map$feature[i],
    title = feature_map$title[i],
    panel_letter = feature_map$panel[i]
  )
}

ordered_panels <- c("A", "B", "C", "D", "E", "F")
ordered_panels <- ordered_panels[ordered_panels %in% names(plot_list)]

combined_fig <- wrap_plots(plot_list[ordered_panels], ncol = 3) +
  plot_annotation(
    title = "Matched Gene-Level Feature Enrichment",
    theme = theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0)
    )
  )

# ── 13. Summary table ─────────────────────────────────────────────────────────
summary_table <- map_dfr(feature_map$feature, function(feat) {
  this_title <- feature_map$title[match(feat, feature_map$feature)]
  
  bind_rows(
    summarise_matched_feature(matched_df, feat, "Frameshift") %>%
      mutate(comparison = "FS vs matched control"),
    summarise_matched_feature(matched_df, feat, "Nonsense") %>%
      mutate(comparison = "STOPGAIN vs matched control")
  ) %>%
    mutate(
      feature = feat,
      title = this_title,
      .before = 1
    )
})

p_table <- tibble(
  feature = feature_map$feature,
  title = feature_map$title,
  p_fs = map_dbl(feature_map$feature, ~ get_mcnemar_p(matched_df, .x, "Frameshift")),
  p_stopgain = map_dbl(feature_map$feature, ~ get_mcnemar_p(matched_df, .x, "Nonsense"))
)

summary_table <- summary_table %>%
  left_join(p_table, by = c("feature", "title"))

# ── 14. Save outputs ──────────────────────────────────────────────────────────
ggsave("Matched_Gene_Level_Enrichment.pdf", combined_fig, width = 17, height = 10, dpi = 300)
ggsave("Matched_Gene_Level_Enrichment.png", combined_fig, width = 17, height = 10, dpi = 300)
write_csv(summary_table, "Matched_Gene_Level_Enrichment_summary.csv")

cat("\nDone!\n")
cat("Output files:\n")
cat("  Matched_Gene_Level_Enrichment.pdf\n")
cat("  Matched_Gene_Level_Enrichment.png\n")
cat("  Matched_Gene_Level_Enrichment_summary.csv\n")