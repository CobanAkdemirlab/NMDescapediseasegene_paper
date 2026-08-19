run_pfam_overlap_analysis <- function(gene_all, ensembl, output_prefix = "pfam_overlap") {
  
  # 1. 获取PFAM信息
  pfam_domains <- getBM(
    attributes = c("ensembl_transcript_id", "pfam", "pfam_start", "pfam_end"),
    filters    = "ensembl_transcript_id",
    values     = unique(gene_all$ensembl_transcript_id),
    mart       = ensembl
  )
  
  # 2. 合并PFAM信息，转换坐标
  gene_all_pfam <- gene_all %>%
    left_join(pfam_domains, by = "ensembl_transcript_id") %>%
    mutate(
      pfam_start_bp = pfam_start * 3 - 2,
      pfam_end_bp   = pfam_end * 3
    )
  
  # 3. 计算每条记录的overlap
  gene_all_pfam <- gene_all_pfam %>%
    mutate(
      overlap_start      = pmax(NMD_region_start, pfam_start_bp),
      overlap_end        = pmin(NMD_region_end, pfam_end_bp),
      overlap_valid      = !is.na(overlap_start) & !is.na(overlap_end) & (overlap_start <= overlap_end),
      overlap_length_raw = ifelse(overlap_valid, overlap_end - overlap_start + 1, 0)
    )
  
  # 4. 对每个row_id做区间reduce，避免重复计数
  idx_list <- split(seq_len(nrow(gene_all_pfam)), gene_all_pfam$row_id)
  
  pfam_overlap_summary <- lapply(idx_list, function(idx) {
    df <- gene_all_pfam[idx, , drop = FALSE]
    nm_len      <- df$NMDesc_region_length[1]
    this_row_id <- df$row_id[1]
    df_valid    <- df[df$overlap_valid %in% TRUE, , drop = FALSE]
    
    if (nrow(df_valid) == 0 || is.na(nm_len) || nm_len <= 0) {
      return(data.frame(
        row_id               = this_row_id,
        pfam_overlap_length  = 0,
        pfam_overlap_flag    = 0,
        pfam_overlap_fraction = 0,
        n_overlapping_pfam   = 0,
        stringsAsFactors     = FALSE
      ))
    }
    
    gr <- GenomicRanges::GRanges(
      seqnames = rep("x", nrow(df_valid)),
      ranges   = IRanges::IRanges(
        start = as.numeric(df_valid$overlap_start),
        end   = as.numeric(df_valid$overlap_end)
      )
    )
    gr_red        <- GenomicRanges::reduce(gr)
    total_overlap <- sum(IRanges::width(gr_red))
    
    data.frame(
      row_id                = this_row_id,
      pfam_overlap_length   = total_overlap,
      pfam_overlap_flag     = as.integer(total_overlap > 20),
      pfam_overlap_fraction = total_overlap / nm_len,
      n_overlapping_pfam    = nrow(df_valid),
      stringsAsFactors      = FALSE
    )
  })
  
  pfam_overlap_summary <- dplyr::bind_rows(pfam_overlap_summary)
  
  # 5. 合并回原始gene_all
  gene_all <- gene_all %>%
    left_join(pfam_overlap_summary, by = "row_id")
  
  # 6. 检查结果
  print(table(gene_all$pfam_overlap_flag, useNA = "ifany"))
  print(summary(gene_all$pfam_overlap_fraction))
  print(summary(gene_all$pfam_overlap_length))
  
  # 7. 分组设置
  gene_all$group <- factor(
    gene_all$group,
    levels = c("snv", "snv_control", "fs", "fs_control")
  )
  
  group_colors2 <- c(
    "snv"         = "#2ca02c",
    "snv_control" = "#98df8a",
    "fs"          = "#1f77b4",
    "fs_control"  = "#aec7e8"
  )
  
  comparisons <- list(
    c("snv", "snv_control"),
    c("fs", "fs_control")
  )
  
  # 8. 图1：PFAM overlap fraction
  p1 <- ggplot(gene_all, aes(x = group, y = pfam_overlap_fraction, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors2, guide = "none") +
    labs(title = "PFAM overlap fraction across variant groups",
         x = NULL, y = "PFAM overlap fraction") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    stat_compare_means(comparisons = comparisons,
                       method = "wilcox.test", label = "p.format")
  
  # 9. 图2：PFAM overlap length
  p2 <- ggplot(gene_all, aes(x = group, y = pfam_overlap_length, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors2, guide = "none") +
    labs(title = "PFAM overlap length across variant groups",
         x = NULL, y = "PFAM overlap length (bp)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    stat_compare_means(comparisons = comparisons,
                       method = "wilcox.test", label = "p.format")
  
  # 10. 图3：PFAM overlap flag 比例
  pfam_flag_summary <- gene_all %>%
    group_by(group) %>%
    summarise(
      n            = n(),
      n_overlap    = sum(pfam_overlap_flag == 1, na.rm = TRUE),
      prop_overlap = n_overlap / n,
      .groups      = "drop"
    )
  
  p3 <- ggplot(pfam_flag_summary, aes(x = group, y = prop_overlap, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors2, guide = "none") +
    labs(title = "Proportion of variants with PFAM overlap > 20 bp",
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    geom_text(aes(label = paste0(n_overlap, "/", n)),
              vjust = -0.3, size = 4)
  
  print(p1); print(p2); print(p3)
  
  # 11. Fisher检验
  for (grp_pair in list(c("snv", "snv_control"), c("fs", "fs_control"))) {
    cat("\n--- Fisher test:", grp_pair[1], "vs", grp_pair[2], "---\n")
    tab <- table(
      gene_all$group[gene_all$group %in% grp_pair],
      gene_all$pfam_overlap_flag[gene_all$group %in% grp_pair]
    )
    print(tab)
    print(fisher.test(tab))
  }
  
  # 12. 保存结果和图片
  write.csv(gene_all, paste0(output_prefix, "_gene_all.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_fraction_plot.png"), p1, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_length_plot.png"),   p2, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_flag_barplot.png"),  p3, width = 8, height = 6, dpi = 300)
  
  # 返回所有结果
  invisible(list(
    gene_all             = gene_all,
    pfam_overlap_summary = pfam_overlap_summary,
    plots                = list(p1 = p1, p2 = p2, p3 = p3)
  ))
}