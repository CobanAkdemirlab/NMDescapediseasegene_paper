run_ppi_overlap_analysis <- function(gene_all, ppi_file_path, output_prefix = "ppi_overlap") {
  
  # 1. 读取PPI网络数据
  human_ppi <- fread(ppi_file_path, data.table = FALSE)
  
  # 2. 字符串转数值向量的辅助函数
  convert_to_c <- function(x) {
    x <- gsub("\\[|\\]", "", x)
    if (x == "") return(c())
    as.numeric(trimws(unlist(strsplit(x, ","))))
  }
  
  # 3. 计算每个基因是否有PPI overlap
  gene_all$matched_uniprot <- 0
  
  for (i in seq_len(nrow(gene_all))) {
    uid         <- gene_all$uniprot[i]
    aft_NMD_ind <- gene_all$NMD_region_start[i]
    
    if (is.na(uid) || uid == "") next
    if (is.na(aft_NMD_ind)) next
    
    re_1 <- unlist(lapply(
      human_ppi$interface_residues1[human_ppi$uniprot1 == uid],
      convert_to_c
    )) * 3 - 2
    
    re_2 <- unlist(lapply(
      human_ppi$interface_residues2[human_ppi$uniprot2 == uid],
      convert_to_c
    )) * 3 - 2
    
    if (length(re_1) == 0 && length(re_2) == 0) next
    
    if (any(re_1 >= aft_NMD_ind, na.rm = TRUE) ||
        any(re_2 >= aft_NMD_ind, na.rm = TRUE)) {
      gene_all$matched_uniprot[i] <- 1
    }
  }
  
  # 4. 重命名为ppi_overlap
  gene_all <- gene_all %>%
    rename(ppi_overlap = matched_uniprot)
  
  # 5. 汇总统计
  ppi_summary <- gene_all %>%
    group_by(group) %>%
    summarise(
      total_genes   = n(),
      matched_genes = sum(ppi_overlap),
      percentage    = matched_genes / total_genes * 100,
      .groups       = "drop"
    )
  print(ppi_summary)
  
  # 6. 设置分组顺序和颜色
  gene_all$group <- factor(
    gene_all$group,
    levels = c("snv", "snv_control", "fs", "fs_control")
  )
  
  group_colors <- c(
    "snv"         = "#2ca02c",
    "snv_control" = "#98df8a",
    "fs"          = "#1f77b4",
    "fs_control"  = "#aec7e8"
  )
  
  ppi_summary$group <- factor(
    ppi_summary$group,
    levels = c("snv", "snv_control", "fs", "fs_control")
  )
  
  # 7. 绘图：各组PPI overlap百分比
  p <- ggplot(ppi_summary, aes(x = group, y = percentage, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, guide = "none") +
    labs(
      title = "Percentage of genes with PPI network overlap",
      x     = NULL,
      y     = "Percentage (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title  = element_text(hjust = 0.5, face = "bold")
    ) +
    geom_text(aes(label = paste0(matched_genes, "/", total_genes)),
              vjust = -0.3, size = 4)
  
  print(p)
  
  # 8. Fisher检验
  for (grp_pair in list(c("snv", "snv_control"), c("fs", "fs_control"))) {
    cat("\n--- Fisher test:", grp_pair[1], "vs", grp_pair[2], "---\n")
    tab <- table(
      gene_all$group[gene_all$group %in% grp_pair],
      gene_all$ppi_overlap[gene_all$group %in% grp_pair]
    )
    print(tab)
    print(fisher.test(tab))
  }
  
  # 9. 保存结果
  write.csv(gene_all,   paste0(output_prefix, "_gene_all.csv"),  row.names = FALSE)
  write.csv(ppi_summary, paste0(output_prefix, "_summary.csv"),  row.names = FALSE)
  ggsave(paste0(output_prefix, "_barplot.png"), p, width = 8, height = 6, dpi = 300)
  
  # 10. 返回结果
  invisible(list(
    gene_all    = gene_all,
    ppi_summary = ppi_summary,
    plot        = p
  ))
}