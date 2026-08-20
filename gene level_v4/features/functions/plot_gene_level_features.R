# 需要的包（source() 本文件即检查；缺包立即报错，不会跑到一半才失败）
for (.pkg in c(
  "biomaRt",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "readr"
)) if (!requireNamespace(.pkg, quietly = TRUE))
  stop(sprintf("缺少包 %s —— install.packages(\"%s\") 或 BiocManager::install(\"%s\")", .pkg, .pkg, .pkg))
for (.pkg in c(
  "biomaRt",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "readr"
)) suppressPackageStartupMessages(library(.pkg, character.only = TRUE))
rm(.pkg)

plot_gene_level_features <- function(
    gene_all,
    lof_metrics_path,
    ensembl     = NULL,
    out_dir     = ".",
    prefix      = "gene_level",
    comparisons = CONFIG$comparisons
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  group_colors <- c("fs" = "#ff7f0e", "fs_control" = "#ffbb78",
                    "snv" = "#2ca02c", "snv_control" = "#98df8a")
  
  # 1. Load and merge LOF metrics; bin pLI / LOEUF
  Lof_metrics <- read.delim(lof_metrics_path)
  pli_all <- merge(gene_all, Lof_metrics, by.x = "hgnc_symbol", by.y = "gene") %>%
    dplyr::select(hgnc_symbol, pLI, oe_lof_upper, group) %>%
    distinct() %>%
    mutate(
      pli_cat   = cut(pLI,          breaks = c(0, 0.35, 0.66, 1), labels = c("Low", "Medium", "High"), include.lowest = TRUE),
      loeuf_cat = cut(oe_lof_upper, breaks = c(0, 0.2, 0.6, 2),   labels = c("Low", "Medium", "High"), include.lowest = TRUE)
    )
  write_csv(pli_all, file.path(out_dir, paste0(prefix, "_pli_loeuf_category.csv")))
  
  base_box <- function(data, y, title, ylab, log_y = FALSE) {
    g <- ggplot(data, aes(x = group, y = .data[[y]], fill = group)) +
      geom_boxplot(width = 0.1, color = "black") +
      theme_minimal() +
      labs(title = title, y = ylab, x = "Group") +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            plot.title   = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none") +
      scale_fill_manual(values = group_colors) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test")
    if (log_y) g <- g + scale_y_log10()
    g
  }
  
  # 2-4. pLI / LOEUF / CDS length
  pli_plot   <- base_box(pli_all,  "pLI",          "pLI distribution by gene group",   "pLI")
  loeuf_plot <- base_box(pli_all,  "oe_lof_upper", "LOEUF distribution by gene group", "LOEUF")
  cds_plot   <- base_box(gene_all, "cds_length",   "CDS length distribution by gene group", "CDS length, bp", log_y = TRUE)
  
  ggsave(file.path(out_dir, paste0(prefix, "_pli_distribution.pdf")),        plot = pli_plot,   width = 8, height = 5)
  ggsave(file.path(out_dir, paste0(prefix, "_loeuf_distribution.pdf")),      plot = loeuf_plot, width = 8, height = 5)
  ggsave(file.path(out_dir, paste0(prefix, "_cds_length_distribution.pdf")), plot = cds_plot,   width = 8, height = 5)
  
  # 5. Exon count (optional, needs ensembl)
  exon_plot <- NULL; exon_counts <- NULL
  if (!is.null(ensembl)) {
    exons <- getBM(
      attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "rank", "strand"),
      filters    = "ensembl_transcript_id",
      values     = unique(gene_all$ensembl_transcript_id),
      mart       = ensembl
    )
    exon_counts <- exons %>%
      inner_join(gene_all, by = "ensembl_transcript_id") %>%
      group_by(ensembl_transcript_id, hgnc_symbol, strand, group) %>%
      summarise(exon_num = max(rank, na.rm = TRUE), .groups = "drop") %>%
      distinct()
    write_csv(exon_counts, file.path(out_dir, paste0(prefix, "_exon_counts.csv")))
    
    exon_plot <- base_box(exon_counts, "exon_num", "Exon count by gene group", "Exon number", log_y = TRUE)
    ggsave(file.path(out_dir, paste0(prefix, "_exon_count_distribution.pdf")), plot = exon_plot, width = 8, height = 5)
  }
  
  list(pli_all = pli_all, exon_counts = exon_counts,
       pli_plot = pli_plot, loeuf_plot = loeuf_plot,
       cds_plot = cds_plot, exon_plot = exon_plot)
}
