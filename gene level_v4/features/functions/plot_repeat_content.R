# 需要的包（source() 本文件即检查；缺包立即报错，不会跑到一半才失败）
for (.pkg in c(
  "ggplot2",
  "ggpubr",
  "stringr"
)) if (!requireNamespace(.pkg, quietly = TRUE))
  stop(sprintf("缺少包 %s —— install.packages(\"%s\") 或 BiocManager::install(\"%s\")", .pkg, .pkg, .pkg))
for (.pkg in c(
  "ggplot2",
  "ggpubr",
  "stringr"
)) suppressPackageStartupMessages(library(.pkg, character.only = TRUE))
rm(.pkg)

plot_repeat_content <- function(
    gene_all,
    comparisons = CONFIG$comparisons,
    output_csv  = "repeat_content.csv",
    output_fig  = "repeat_content.png"
) {
  detect_repeats <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(0)
    matches <- str_extract_all(sequence, "([ATGC]{1,6})\\1+")[[1]]
    if (length(matches) == 0) return(0)
    sum(nchar(matches))
  }
  repeat_fraction <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    detect_repeats(sequence) / nchar(sequence)
  }
  detect_homopolymer <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(0)
    matches <- str_extract_all(sequence, "([ATGC])\\1{3,}")[[1]]
    if (length(matches) == 0) return(0)
    sum(nchar(matches))
  }
  homopolymer_fraction <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    detect_homopolymer(sequence) / nchar(sequence)
  }
  
  gene_all$repeat_fraction             <- sapply(gene_all$coding,     repeat_fraction)
  gene_all$nmdesc_repeat_fraction      <- sapply(gene_all$nmdesc_cds, repeat_fraction)
  gene_all$homopolymer_fraction        <- sapply(gene_all$coding,     homopolymer_fraction)
  gene_all$nmdesc_homopolymer_fraction <- sapply(gene_all$nmdesc_cds, homopolymer_fraction)
  
  make_plot <- function(y_var, y_label, title) {
    ggplot(gene_all, aes(x = group, y = .data[[y_var]], fill = group)) +
      geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_fill_manual(values = CONFIG$group_colors) +
      labs(x = "", y = y_label, title = title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  }
  
  p1 <- make_plot("repeat_fraction",             "Repeat fraction",      "Overall Repeat Content")
  p2 <- make_plot("nmdesc_repeat_fraction",      "Repeat fraction",      "NMD-escape Repeat Content")
  p3 <- make_plot("homopolymer_fraction",        "Homopolymer fraction", "Overall Homopolymer Content")
  p4 <- make_plot("nmdesc_homopolymer_fraction", "Homopolymer fraction", "NMD-escape Homopolymer Content")
  combined <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  out_cols <- c("hgnc_symbol", "ensembl_transcript_id", "group",
                "repeat_fraction", "nmdesc_repeat_fraction",
                "homopolymer_fraction", "nmdesc_homopolymer_fraction")
  write.csv(gene_all[, out_cols], output_csv, row.names = FALSE)
  message("CSV saved to: ", output_csv)
  ggsave(output_fig, plot = combined, width = 12, height = 10, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = combined, data = gene_all))
}


# ── Usage ─────────────────────────────────────────────────────────────────

# ===========================================================================
# 以下为原文件末尾的示例调用，已注释。
# 原来 source() 本文件会立刻执行它们并去读数据文件 —— 函数库不应有副作用。
# 需要跑的话，复制到你的脚本里，并把路径改成 data_file("文件名")。
# ===========================================================================
# result <- plot_repeat_content(
#   gene_all   = gene_all,
#   output_csv = "repeat_content.csv",
#   output_fig = "repeat_content.png"
# )
