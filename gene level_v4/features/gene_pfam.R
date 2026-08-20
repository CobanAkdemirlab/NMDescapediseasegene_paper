library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)


#get pfam info by gene
pfam_domains <- getBM(
  attributes = c("ensembl_transcript_id", "pfam", "pfam_start", "pfam_end"),
  filters    = "ensembl_transcript_id",
  values     = unique(gene_all$ensembl_transcript_id),
  mart       = ensembl
)

#keep all pfam records
gene_all_pfam <- gene_all %>%
  left_join(pfam_domains, by = "ensembl_transcript_id") %>%
  mutate(
    pfam_start_bp = pfam_start * 3 - 2,   # aa坐标转bp坐标，起始更严格些
    pfam_end_bp   = pfam_end * 3
  )

#get overlap for each record
gene_all_pfam <- gene_all_pfam %>%
  mutate(
    overlap_start = pmax(NMD_region_start, pfam_start_bp),
    overlap_end   = pmin(NMD_region_end, pfam_end_bp),
    overlap_valid = !is.na(overlap_start) & !is.na(overlap_end) & (overlap_start <= overlap_end),
    overlap_length_raw = ifelse(overlap_valid, overlap_end - overlap_start + 1, 0)
  )

###--------------------------------------------------
### 4. 对每个 row_id，把多个 PFAM overlap 区间 reduce 后求 union 长度
### 避免多个 PFAM 彼此重叠时重复计数
###--------------------------------------------------
idx_list <- split(seq_len(nrow(gene_all_pfam)), gene_all_pfam$row_id)

pfam_overlap_summary <- lapply(idx_list, function(idx) {
  df <- gene_all_pfam[idx, , drop = FALSE]
  
  nm_len <- df$NMDesc_region_length[1]
  this_row_id <- df$row_id[1]
  
  df_valid <- df[df$overlap_valid %in% TRUE, , drop = FALSE]
  
  if (nrow(df_valid) == 0 || is.na(nm_len) || nm_len <= 0) {
    return(data.frame(
      row_id = this_row_id,
      pfam_overlap_length = 0,
      pfam_overlap_flag = 0,
      pfam_overlap_fraction = 0,
      n_overlapping_pfam = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = rep("x", nrow(df_valid)),
    ranges = IRanges::IRanges(
      start = as.numeric(df_valid$overlap_start),
      end   = as.numeric(df_valid$overlap_end)
    )
  )
  
  gr_red <- GenomicRanges::reduce(gr)
  total_overlap <- sum(IRanges::width(gr_red))
  
  data.frame(
    row_id = this_row_id,
    pfam_overlap_length = total_overlap,
    pfam_overlap_flag = as.integer(total_overlap > 20),
    pfam_overlap_fraction = total_overlap / nm_len,
    n_overlapping_pfam = nrow(df_valid),
    stringsAsFactors = FALSE
  )
})

pfam_overlap_summary <- dplyr::bind_rows(pfam_overlap_summary)

###--------------------------------------------------
### 5. 合并回原始 gene_all
###--------------------------------------------------
gene_all <- gene_all %>%
  left_join(pfam_overlap_summary, by = "row_id")

###--------------------------------------------------
### 6. 检查结果
###--------------------------------------------------
table(gene_all$pfam_overlap_flag, useNA = "ifany")
summary(gene_all$pfam_overlap_fraction)
summary(gene_all$pfam_overlap_length)

###--------------------------------------------------
### 7. 设置分组顺序和颜色
### 你可以按自己的分组名字改
###--------------------------------------------------
gene_all$group <- factor(
  gene_all$group,
  levels = c("snv", "snv_control", "fs", "fs_control")
)

group_colors2 <- c(
  "snv" = "#2ca02c",
  "snv_control" = "#98df8a",
  "fs" = "#1f77b4",
  "fs_control" = "#aec7e8"
)

comparisons <- list(
  c("snv", "snv_control"),
  c("fs", "fs_control")
)

###--------------------------------------------------
### 8. 图1：PFAM overlap fraction 分布图
###--------------------------------------------------
p1 <- ggplot(gene_all, aes(x = group, y = pfam_overlap_fraction, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "PFAM overlap fraction across variant groups",
    x = NULL,
    y = "PFAM overlap fraction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format"
  )

print(p1)

###--------------------------------------------------
### 9. 图2：PFAM overlap length 分布图
### 如果想看 log10，可以加 scale_y_log10()
###--------------------------------------------------
p2 <- ggplot(gene_all, aes(x = source, y = pfam_overlap_length, fill = source)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "PFAM overlap length across variant groups",
    x = NULL,
    y = "PFAM overlap length (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format"
  )

print(p2)

###--------------------------------------------------
### 10. 图3：PFAM overlap flag 比例图
###--------------------------------------------------
pfam_flag_summary <- gene_all %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_overlap = sum(pfam_overlap_flag == 1, na.rm = TRUE),
    prop_overlap = n_overlap / n
  )

p3 <- ggplot(pfam_flag_summary, aes(x = group, y = prop_overlap, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Proportion of variants with PFAM overlap > 20 bp",
    x = NULL,
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(aes(label = paste0(n_overlap, "/", n)),
            vjust = -0.3, size = 4)

print(p3)

###--------------------------------------------------
### 11. 如果你想做组间比例检验（Fisher）
###--------------------------------------------------
# snv vs snv_control
tab_snv <- table(
  gene_all$source[gene_all$source %in% c("snv", "snv_control")],
  gene_all$pfam_overlap_flag[gene_all$source %in% c("snv", "snv_control")]
)
print(tab_snv)
print(fisher.test(tab_snv))

# fs vs fs_control
tab_fs <- table(
  gene_all$source[gene_all$source %in% c("fs", "fs_control")],
  gene_all$pfam_overlap_flag[gene_all$source %in% c("fs", "fs_control")]
)
print(tab_fs)
print(fisher.test(tab_fs))

###--------------------------------------------------
### 12. 保存结果和图片
###--------------------------------------------------
write.csv(gene_all, "gene_all_with_pfam_overlap.csv", row.names = FALSE)

ggsave("pfam_overlap_fraction_plot.png", p1, width = 8, height = 6, dpi = 300)
ggsave("pfam_overlap_length_plot.png", p2, width = 8, height = 6, dpi = 300)
ggsave("pfam_overlap_flag_barplot.png", p3, width = 8, height = 6, dpi = 300)