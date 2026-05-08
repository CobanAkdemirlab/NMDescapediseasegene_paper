library(readr)
 minus1_control_vep_uniprot_intact <- read_delim("Downloads/vep0512 2/minus1_control_vep_uniprot_intact.txt", 
                                                       delim = "\t", escape_double = FALSE, 
                                                       trim_ws = TRUE, skip = 94)
 minus1_vep_uniprot_intact <- read_delim("Downloads/vep0512 2/minus1_vep_uniprot_intact.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE, skip = 94)
 plus1_control_vep_uniprot_intact <- read_delim("vep0512 2/plus1_control_vep_uniprot_intact.txt", 
                                                +     delim = "\t", escape_double = FALSE, 
                                                +     trim_ws = TRUE, skip = 94)
 plus1_vep_uniprot_intact <- read_delim("vep0512 2/plus1_vep_uniprot_intact.txt",
                                             delim = "\t", escape_double = FALSE, 
                                             trim_ws = TRUE, skip = 94)
 snv_control_vep_uniprot_intact <- read_delim("vep0512 2/snv_control_vep_uniprot_intact.txt", 
                                             delim = "\t", escape_double = FALSE, 
                                             trim_ws = TRUE, skip = 94)
 snv_vep_uniprot_intact <- read_delim("vep0512 2/snv_vep_uniprot_intact.txt",
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE, skip = 94)
 
 #only keep the rows where variant and gene both match minus1_variants
minus1_vep = minus1_vep_uniprot_intact[(which(minus1_vep_uniprot_intact$Feature %in% minus1_variants$transcript)),]
minus1_control_vep = minus1_control_vep_uniprot_intact[(which(minus1_control_vep_uniprot_intact$Feature %in% minus1_control_variants$transcript)),]

plus1_vep = plus1_vep_uniprot_intact[(which(plus1_vep_uniprot_intact$Feature %in% plus1_variants$transcript)),]
plus1_control_vep = plus1_control_vep_uniprot_intact[(which(plus1_control_vep_uniprot_intact$Feature %in% plus1_control_variants$transcript)),]
snv_vep = snv_vep_uniprot_intact[(which(snv_vep_uniprot_intact$Feature %in% snv_variants$transcript)),]
snv_control_vep = snv_control_vep_uniprot_intact[(which(snv_control_vep_uniprot_intact$Feature %in% snv_control_variants$transcript)),]

table(minus1_vep$IMPACT)
length(which(minus1_vep$IMPACT == "MODIFIER"))/length(minus1_vep$IMPACT)
table(minus1_control_vep$IMPACT)
length(which(minus1_control_vep$IMPACT == "MODIFIER"))/length(minus1_control_vep$IMPACT)
table(plus1_vep$IMPACT)
length(which(plus1_vep$IMPACT == "MODIFIER"))/length(plus1_vep$IMPACT)
table(plus1_control_vep$IMPACT)
length(which(plus1_control_vep$IMPACT == "MODIFIER"))/length(plus1_control_vep$IMPACT)
table(snv_vep$IMPACT)
length(which(snv_vep$IMPACT == "MODIFIER"))/length(snv_vep$IMPACT)
table(snv_control_vep$IMPACT)
length(which(snv_control_vep$IMPACT == "MODIFIER"))/length(snv_control_vep$IMPACT)

#add cds lenth to each transcript
minus1_idr_unique <- minus1_idr %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # 或按特定规则排序后取最新值

merged_minus1 <- left_join(
  minus1_vep,
  minus1_idr_unique[, c("Transcript_ID", "NMDesc.end")],
  by = c("Feature" = "Transcript_ID")
)
merged_minus1$start_cds = as.numeric(sub("-.*", "", merged_minus1$CDS_position))
merged_minus1$vdis2end = merged_minus1$NMDesc.end - merged_minus1$start_cds

#assume the no match is NA
merged_rn_minus1 = merged_minus1[which(!is.na(merged_minus1$NMDesc.end)),]

minus1_control_idr_unique <- minus1_control_idr %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # 或按特定规则排序后取最新值

merged_minus1_control <- left_join(
  minus1_control_vep,
  minus1_control_idr_unique[, c("Transcript_ID", "NMDesc.end")],
  by = c("Feature" = "Transcript_ID")
)
merged_minus1_control$start_cds = as.numeric(sub("-.*", "", merged_minus1_control$CDS_position))
merged_minus1_control$vdis2end = merged_minus1_control$NMDesc.end - merged_minus1_control$start_cds
merged_rn_minus1_control = merged_minus1_control[which(!is.na(merged_minus1_control$NMDesc.end)),]

plus1_idr_unique <- plus1_idr %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # 或按特定规则排序后取最新值

merged_plus1 <- left_join(
  plus1_vep,
  plus1_idr_unique[, c("Transcript_ID", "NMDesc.end")],
  by = c("Feature" = "Transcript_ID")
)
merged_plus1$start_cds = as.numeric(sub("-.*", "", merged_plus1$CDS_position))
merged_plus1$vdis2end = merged_plus1$NMDesc.end - merged_plus1$start_cds

#assume the no match is NA
merged_rn_plus1 = merged_plus1[which(!is.na(merged_plus1$NMDesc.end)),]

plus1_control_idr_unique <- plus1_control_idr %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # 或按特定规则排序后取最新值

merged_plus1_control <- left_join(
  plus1_control_vep,
  plus1_control_idr_unique[, c("Transcript_ID", "NMDesc.end")],
  by = c("Feature" = "Transcript_ID")
)
merged_plus1_control$start_cds = as.numeric(sub("-.*", "", merged_plus1_control$CDS_position))
merged_plus1_control$vdis2end = merged_plus1_control$NMDesc.end - merged_plus1_control$start_cds
merged_rn_plus1_control = merged_plus1_control[which(!is.na(merged_plus1_control$NMDesc.end)),]
snv_idr_unique <- snv_idr %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # 或按特定规则排序后取最新值

merged_snv <- left_join(
  snv_vep,
  snv_idr_unique[, c("Transcript_ID", "NMDesc.end")],
  by = c("Feature" = "Transcript_ID")
)
merged_snv$start_cds = as.numeric(sub("-.*", "", merged_snv$CDS_position))
merged_snv$vdis2end = merged_snv$NMDesc.end - merged_snv$start_cds

#assume the no match is NA
merged_rn_snv = merged_snv[which(!is.na(merged_snv$NMDesc.end)),]

snv_control_idr_unique <- snv_control_idr %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # 或按特定规则排序后取最新值

merged_snv_control <- left_join(
  snv_control_vep,
  snv_control_idr_unique[, c("Transcript_ID", "NMDesc.end")],
  by = c("Feature" = "Transcript_ID")
)
merged_snv_control$start_cds = as.numeric(sub("-.*", "", merged_snv_control$CDS_position))
merged_snv_control$vdis2end = merged_snv_control$NMDesc.end - merged_snv_control$start_cds
merged_rn_snv_control = merged_snv_control[which(!is.na(merged_snv_control$NMDesc.end)),]


# 为每个组别添加标识列
merged_rn_minus1$Group <- "minus1"
merged_rn_minus1_control$Group <- "minus1_control"
merged_rn_plus1$Group <- "plus1"
merged_rn_plus1_control$Group <- "plus1_control"
merged_rn_snv$Group <- "snv"
merged_rn_snv_control$Group <- "snv_control"

# 合并所有数据
combined_data <- rbind(
  merged_rn_minus1,
  merged_rn_minus1_control,
  merged_rn_plus1,
  merged_rn_plus1_control,
  merged_rn_snv,
  merged_rn_snv_control
)

# 提取变异类型和是否为对照
combined_data <- combined_data %>%
  mutate(
    Variant_Type = case_when(
      grepl("minus1", Group) ~ "minus1",
      grepl("plus1", Group) ~ "plus1",
      grepl("snv", Group) ~ "snv"
    ),
    Is_Control = grepl("control", Group)
  )


library(ggplot2)

ggplot(combined_data, aes(x = Variant_Type, y = vdis2end, fill = Is_Control)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  scale_fill_manual(
    name = "Group",
    labels = c("P/LP", "Benign"),
    values = c("#E69F00", "#56B4E9")
  ) +
  labs(
    title = "Comparison of variant distance to cds end Across Groups",
    x = "Variant Type",
    y = "vdis2end (distance from variant to CDS end)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

library(rstatix)

# 比较每个变异类型的实验组 vs 对照组
stat_results <- combined_data %>%
  group_by(Variant_Type) %>%
  wilcox_test(vdis2end ~ Is_Control) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# 输出结果
print(stat_results)

# 创建分面数据标签
label_data <- data.frame(
  Variant_Type = c("minus1", "plus1", "snv"),
  label = paste("p =", c("1.3e-11", "1.3e-09", "1.1e-31")),  # 从 stat_results 提取
  x = 1.5,  # 分面内 x 轴居中位置
  y = max(combined_data$vdis2end, na.rm = TRUE) * 1.1
)

# 分面箱线图
ggplot(combined_data, aes(x = Is_Control, y = vdis2end)) +
  geom_boxplot(aes(fill = Is_Control)) +
  scale_y_log10() + 
  geom_text(
    data = label_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE
  ) +
  facet_wrap(~Variant_Type) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9"),
    labels = c("P/LP", "Benign")
  ) +
  theme_bw()
