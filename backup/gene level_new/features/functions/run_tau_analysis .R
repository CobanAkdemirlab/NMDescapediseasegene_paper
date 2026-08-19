run_tau_analysis <- function(gene_all,
                             gtex_path,
                             output_prefix = "tau") {
  
  # 1. 加载必要包
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(readr)
  
  # 2. 读取GTEx数据
  gtex <- read_tsv(gtex_path, skip = 2)
  
  # 3. 清理Ensembl ID，按基因symbol聚合
  gtex$Name <- sub("\\..*", "", gtex$Name)
  
  gene_expr <- gtex %>%
    dplyr::select(-tidyselect::any_of("Name")) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::rename(Gene = Description)
  
  # 4. 保存基因表达矩阵
  write.csv(gene_expr, paste0(output_prefix, "_gene_matrix.csv"), row.names = FALSE)
  
  # 5. 准备矩阵
  gene.df     <- as.data.frame(gene_expr)
  gene.matrix <- as.matrix(gene.df[, -1])
  rownames(gene.matrix) <- gene.df$Gene
  
  # 6. 定义tau函数
  tau <- function(x) {
    if (any(is.na(x)))  stop("NA values found.")
    if (any(x < 0))     stop("Negative values found. Maybe data is log-transformed?")
    sum(1 - x / max(x)) / (length(x) - 1)
  }
  
  # 7. 计算每个基因的tau值
  tau_scores <- apply(gene.matrix, 1, tau)
  
  # 8. 过滤并添加tau值
  tau_all <- gene_all %>%
    filter(hgnc_symbol %in% names(tau_scores)) %>%
    mutate(tau = tau_scores[hgnc_symbol])
  
  # 9. 设置分组顺序和颜色
  tau_all$group <- factor(
    tau_all$group,
    levels = c("fs", "fs_control", "snv", "snv_control")
  )
  
  group_colors <- c(
    "fs"          = "#1f77b4",
    "fs_control"  = "#aec7e8",
    "snv"         = "#2ca02c",
    "snv_control" = "#98df8a"
  )
  
  comparisons <- list(
    c("fs",  "fs_control"),
    c("snv", "snv_control")
  )
  
  # 10. 绘图
  p <- ggplot(tau_all, aes(x = group, y = tau, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons,
                       method = "wilcox.test",
                       label  = "p.signif") +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = c(
      "fs"          = "fs",
      "fs_control"  = "fs_Control",
      "snv"         = "Nonsense",
      "snv_control" = "Nonsense_Control"
    )) +
    theme_minimal() +
    labs(
      title = "Tissue specificity (Tau) across gene categories",
      y     = "Tissue Specificity (Tau)",
      x     = "Gene group"
    ) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title.x    = element_text(size = 12, face = "bold"),
      axis.title.y    = element_text(size = 12, face = "bold"),
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      panel.border    = element_rect(colour = "black", fill = NA)
    )
  
  print(p)
  
  # 11. 保存结果
  write.csv(tau_all, paste0(output_prefix, "_all.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_violin_plot.png"), p, width = 8, height = 6, dpi = 300)
  
  # 12. 返回结果
  invisible(list(
    tau_all    = tau_all,
    tau_scores = tau_scores,
    plot       = p
  ))
}