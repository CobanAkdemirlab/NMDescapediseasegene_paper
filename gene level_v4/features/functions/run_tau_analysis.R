# 需要的包（source() 本文件即检查；缺包立即报错，不会跑到一半才失败）
for (.pkg in c(
  "dplyr",
  "ggplot2",
  "ggpubr",
  "readr"
)) if (!requireNamespace(.pkg, quietly = TRUE))
  stop(sprintf("缺少包 %s —— install.packages(\"%s\") 或 BiocManager::install(\"%s\")", .pkg, .pkg, .pkg))
for (.pkg in c(
  "dplyr",
  "ggplot2",
  "ggpubr",
  "readr"
)) suppressPackageStartupMessages(library(.pkg, character.only = TRUE))
rm(.pkg)

run_tau_analysis <- function(gene_all, gtex_path, output_prefix = "tau") {

  # GTEx 是【可选】输入：缺失时把 tau 填 NA 并返回，不中断整条流程。
  # 与 gene_main_dbh.R 里 annotate_tau() 的行为一致。
  # 注意：gene_all_matched.csv / gene_all_random.csv / tau_all.csv 里
  # 【已经有算好的 tau 列】——只有要用新版 GTEx 重算时才需要这个文件。
  if (is.null(gtex_path) || length(gtex_path) != 1 || is.na(gtex_path) ||
      !file.exists(gtex_path)) {
    message("  跳过组织特异性 tau —— GTEx 文件不可用（tau 置为 NA）")
    # 返回结构必须与正常路径一致：list(tau_all, tau_scores, plot)，
    # 否则调用方取 $tau_all 会拿到 NULL。
    return(invisible(list(
      tau_all    = gene_all %>% dplyr::mutate(tau = NA_real_),
      tau_scores = NULL,
      plot       = NULL)))
  }

  # 1. GTEx median TPM
  gtex <- read_tsv(gtex_path, skip = 2)
  gtex$Name <- sub("\\..*", "", gtex$Name)
  
  gene_expr <- gtex %>%
    dplyr::select(-tidyselect::any_of("Name")) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::rename(Gene = Description)
  
  write.csv(gene_expr, paste0(output_prefix, "_gene_matrix.csv"), row.names = FALSE)
  
  # 2. Build expression matrix
  gene.df     <- as.data.frame(gene_expr)
  gene.matrix <- as.matrix(gene.df[, -1])
  rownames(gene.matrix) <- gene.df$Gene
  
  # 3. Tau (tissue specificity)
  tau <- function(x) {
    if (any(is.na(x))) stop("NA values found.")
    if (any(x < 0))    stop("Negative values found. Maybe data is log-transformed?")
    sum(1 - x / max(x)) / (length(x) - 1)
  }
  tau_scores <- apply(gene.matrix, 1, tau)
  
  # 4. Attach tau to gene_all
  tau_all <- gene_all %>%
    filter(hgnc_symbol %in% names(tau_scores)) %>%
    mutate(tau = tau_scores[hgnc_symbol])
  
  # 5. Group order / colors
  tau_all$group <- factor(tau_all$group, levels = c("fs", "fs_control", "snv", "snv_control"))
  group_colors <- c("fs" = "#1f77b4", "fs_control" = "#aec7e8",
                    "snv" = "#2ca02c", "snv_control" = "#98df8a")
  comparisons <- list(c("fs", "fs_control"), c("snv", "snv_control"))
  
  # 6. Violin plot
  p <- ggplot(tau_all, aes(x = group, y = tau, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = c("fs" = "fs", "fs_control" = "fs_Control",
                                "snv" = "Nonsense", "snv_control" = "Nonsense_Control")) +
    theme_minimal() +
    labs(title = "Tissue specificity (Tau) across gene categories",
         y = "Tissue Specificity (Tau)", x = "Gene group") +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none",
          plot.title   = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA))
  print(p)
  
  # 7. Save
  write.csv(tau_all, paste0(output_prefix, "_all.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_violin_plot.png"), p, width = 8, height = 6, dpi = 300)
  
  invisible(list(tau_all = tau_all, tau_scores = tau_scores, plot = p))
}
