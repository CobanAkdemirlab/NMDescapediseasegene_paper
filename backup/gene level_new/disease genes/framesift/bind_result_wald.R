
PTC_info = read.csv('PTC_info20260201_region.csv')
fs_NMD_result = read.csv('fs_NMD_result20260201.csv')

#add number of variants for plus1 and plus2
count_fs = fs_NMD_result %>%
  # 只保留 canonical NMDesc
  filter(NMDesc_can == TRUE) %>%
  # 按 transcript + type 统计
  group_by(transcript_id, type) %>%
  summarise(count = n(), .groups = "drop") %>%
  # 转成宽格式
  pivot_wider(
    names_from = type,
    values_from = count,
    values_fill = 0
  ) %>%
  # 重命名列（如果 type 是 plus1 / plus2）
  rename(
    plus1_count = plus1,
    plus2_count = plus2
  )

transcript_object <- list()
#transcript_set3 <- unique(PTC_info$transcript)
for (transcript in transcript_set3) {
  # For each transcript, organize the data by plus
  PTC_plus1_idx <- which(PTC_info$transcript == transcript & PTC_info$type == "plus1")
  PTC_plus2_idx <- which(PTC_info$transcript == transcript & PTC_info$type == "plus2") 
  va_plus1_idx <- which(fs_NMD_result$transcript_id == transcript & fs_NMD_result$type == "plus1")
  va_plus2_idx <- which(fs_NMD_result$transcript_id == transcript & fs_NMD_result$type == "plus2")
  
  # Add data for each plus
  transcript_object[[transcript]] <- list(
    plus1 = list(
      PTC_loc = PTC_info$PTC_loc[PTC_plus1_idx],
      PTC_NMDesc = PTC_info$NMDesc[PTC_plus1_idx],
      can_region = PTC_info$can_region[PTC_plus1_idx],
      css_region = PTC_info$css_region[PTC_plus1_idx],
      long_region = PTC_info$long_region[PTC_plus1_idx],
      NMDesc_can = fs_NMD_result$NMDesc_can[va_plus1_idx],
      NMDesc_css = fs_NMD_result$NMDesc_css[va_plus1_idx],
      NMDesc_long = fs_NMD_result$NMDesc_long[va_plus1_idx]
    ),
    plus2 = list(
      PTC_loc = PTC_info$PTC_loc[PTC_plus2_idx],
      PTC_NMDesc = PTC_info$NMDesc[PTC_plus2_idx],
      can_region = PTC_info$can_region[PTC_plus2_idx],
      css_region = PTC_info$css_region[PTC_plus2_idx],
      long_region = PTC_info$long_region[PTC_plus2_idx],
      NMDesc_can = fs_NMD_result$NMDesc_can[va_plus2_idx],
      NMDesc_css = fs_NMD_result$NMDesc_css[va_plus2_idx],
      NMDesc_long = fs_NMD_result$NMDesc_long[va_plus2_idx]
    )
  )
}

plus1_NMD_result <- lapply(transcript_object, function(transcript) {
  if (!is.null(transcript$plus1)) {
    transcript$plus1
  } else {
    NULL  # 如果没有 plus1，则返回 NULL
  }
})

#6a. get plus1 enrichment p value
plus1_results_df <- rbind(plus1_results_df, data.frame(
  transcript       = transcript,
  NMDesc_can_count = NMDesc_can_count,
  can_region       = ifelse(is.null(plus1_NMD_result[[i]][["can_region"]]), NA_integer_, plus1_NMD_result[[i]][["can_region"]]),
  all_variant      = all.variant,
  stringsAsFactors = FALSE
))
cds.info = getBM(attributes = c("ensembl_transcript_id", "chromosome_name",
                                "cds_start",'cds_end'),
                 filters = "ensembl_transcript_id",
                 values = names(plus1_NMD_result),
                 mart = ensembl)
for(i in seq_along(plus1_NMD_result)){
  print(i)
  transcript       = names(plus1_NMD_result)[i]
  fs_NMD_can       = plus1_NMD_result[[i]][["NMDesc_can"]]
  NMDesc_can_count = length(which(fs_NMD_can == TRUE))
  non_NMDesc_can_count = length(which(fs_NMD_can == FALSE))
  all.variant      = NMDesc_can_count + non_NMDesc_can_count
  
  # 默认所有结果为NA
  plus1_NMD_result[[i]][["can_p"]]           <- NA
  plus1_NMD_result[[i]][["can_wald_logOR"]]  <- NA_real_
  plus1_NMD_result[[i]][["can_wald_se"]]     <- NA_real_
  plus1_NMD_result[[i]][["can_wald_p"]]      <- NA_real_
  plus1_NMD_result[[i]][["can_wald_CI_low"]] <- NA_real_
  plus1_NMD_result[[i]][["can_wald_CI_up"]]  <- NA_real_
  can_region_val <- NA_integer_
  
  # 保护1：cds.info中找不到该转录本
  cds.rows = which(cds.info$ensembl_transcript_id == transcript)
  if(length(cds.rows) == 0){
    plus1_results_df <- rbind(plus1_results_df, data.frame(
      transcript = transcript, NMDesc_can_count = NMDesc_can_count,
      can_region = can_region_val, all_variant = all.variant,
      stringsAsFactors = FALSE))
    next
  }
  
  cds.start = 1
  cds.end   = max(cds.info$cds_end[cds.rows])
  chrom     = cds.info$chromosome_name[cds.rows][1]
  
  # 保护2：PTC_info中找不到该转录本的plus1记录
  ptc.rows = which(PTC_info$transcript == transcript & PTC_info$type == "plus1")
  if(length(ptc.rows) == 0){
    plus1_results_df <- rbind(plus1_results_df, data.frame(
      transcript = transcript, NMDesc_can_count = NMDesc_can_count,
      can_region = can_region_val, all_variant = all.variant,
      stringsAsFactors = FALSE))
    next
  }
  
  NMDesc.start = PTC_info$can_region_start[ptc.rows]
  NMDesc.end   = PTC_info$can_region_end[ptc.rows]
  NMDesc.count = get_syn_count(chrom, NMDesc.start, NMDesc.end)
  n_can        = NMDesc.count
  n_all        = get_syn_count(chrom, cds.start, cds.end)
  n_rest       = n_all - n_can
  rest.PTC     = all.variant - NMDesc_can_count
  
  # 记录can_region值（用于结果df）
  can_region_val = if(is.null(plus1_NMD_result[[i]][["can_region"]])) NA_integer_ else plus1_NMD_result[[i]][["can_region"]]
  
  # 保护3：关键数值无效，跳过检验
  if(is.na(n_can) || n_can <= 0 || is.na(n_all) || n_all <= 0 || all.variant == 0 ||
     is.null(can_region_val) || is.na(can_region_val) || can_region_val <= 0 ||
     NMDesc_can_count <= 0){
    plus1_results_df <- rbind(plus1_results_df, data.frame(
      transcript = transcript, NMDesc_can_count = NMDesc_can_count,
      can_region = can_region_val, all_variant = all.variant,
      stringsAsFactors = FALSE))
    next
  }
  
  # 二项检验（单侧，富集）
  plus1_NMD_result[[i]][["can_p"]] = binom.test(
    NMDesc_can_count, n_can, all.variant / n_all, alternative = 'greater')$p.value
  
  # Wald检验（单侧，富集）
  if(!is.na(n_rest) && n_rest > 0 && rest.PTC >= 0 && rest.PTC <= n_rest){
    a <- NMDesc_can_count + 0.5
    b <- n_can - NMDesc_can_count + 0.5
    c <- rest.PTC + 0.5
    d <- n_rest - rest.PTC + 0.5
    logOR    <- log(a * d / (b * c))
    se_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
    z_stat   <- logOR / se_logOR
    plus1_NMD_result[[i]][["can_wald_logOR"]]  <- logOR
    plus1_NMD_result[[i]][["can_wald_se"]]     <- se_logOR
    plus1_NMD_result[[i]][["can_wald_p"]]      <- pnorm(z_stat, lower.tail = FALSE)
    plus1_NMD_result[[i]][["can_wald_CI_low"]] <- logOR - 1.96 * se_logOR
    plus1_NMD_result[[i]][["can_wald_CI_up"]]  <- logOR + 1.96 * se_logOR
  }
  
  plus1_results_df <- rbind(plus1_results_df, data.frame(
    transcript       = transcript,
    NMDesc_can_count = NMDesc_can_count,
    can_region       = can_region_val,
    all_variant      = all.variant,
    stringsAsFactors = FALSE
  ))
}
saveRDS(plus1_NMD_result, file = "plus1_NMD_wald_result20260201.rds")
plus1_NMD_result = readRDS("plus1_NMD_result20260201.rds")

plus2_NMD_result <- lapply(transcript_object, function(transcript) {
  if (!is.null(transcript$plus2)) transcript$plus2 else NULL
})

#6b. get plus2 enrichment p value
plus2_results_df <- data.frame(
  transcript = character(),
  NMDesc_can_count = integer(),
  can_region = integer(),
  all_variant = integer(),
  stringsAsFactors = FALSE
)
cds.info = getBM(attributes = c("ensembl_transcript_id", "chromosome_name",
                                "cds_start",'cds_end'),
                 filters = "ensembl_transcript_id",
                 values = names(plus2_NMD_result),
                 mart = ensembl)
for(i in seq_along(plus2_NMD_result)){
  print(i)
  transcript       = names(plus2_NMD_result)[i]
  fs_NMD_can       = plus2_NMD_result[[i]][["NMDesc_can"]]
  NMDesc_can_count = length(which(fs_NMD_can == TRUE))
  non_NMDesc_can_count = length(which(fs_NMD_can == FALSE))
  all.variant      = NMDesc_can_count + non_NMDesc_can_count
  
  # 默认所有结果为NA
  plus2_NMD_result[[i]][["can_p"]]           <- NA
  plus2_NMD_result[[i]][["can_wald_logOR"]]  <- NA_real_
  plus2_NMD_result[[i]][["can_wald_se"]]     <- NA_real_
  plus2_NMD_result[[i]][["can_wald_p"]]      <- NA_real_
  plus2_NMD_result[[i]][["can_wald_CI_low"]] <- NA_real_
  plus2_NMD_result[[i]][["can_wald_CI_up"]]  <- NA_real_
  can_region_val <- NA_integer_
  
  # 保护1：cds.info中找不到该转录本
  cds.rows = which(cds.info$ensembl_transcript_id == transcript)
  if(length(cds.rows) == 0){
    plus2_results_df <- rbind(plus2_results_df, data.frame(
      transcript = transcript, NMDesc_can_count = NMDesc_can_count,
      can_region = can_region_val, all_variant = all.variant,
      stringsAsFactors = FALSE))
    next
  }
  
  cds.start = 1
  cds.end   = max(cds.info$cds_end[cds.rows])
  chrom     = cds.info$chromosome_name[cds.rows][1]
  
  # 保护2：PTC_info中找不到该转录本的plus2记录
  ptc.rows = which(PTC_info$transcript == transcript & PTC_info$type == "plus2")
  if(length(ptc.rows) == 0){
    plus2_results_df <- rbind(plus2_results_df, data.frame(
      transcript = transcript, NMDesc_can_count = NMDesc_can_count,
      can_region = can_region_val, all_variant = all.variant,
      stringsAsFactors = FALSE))
    next
  }
  
  NMDesc.start = PTC_info$can_region_start[ptc.rows]
  NMDesc.end   = PTC_info$can_region_end[ptc.rows]
  NMDesc.count = get_syn_count(chrom, NMDesc.start, NMDesc.end)
  n_can        = NMDesc.count
  n_all        = get_syn_count(chrom, cds.start, cds.end)
  n_rest       = n_all - n_can
  rest.PTC     = all.variant - NMDesc_can_count
  
  # 记录can_region值（用于结果df）
  can_region_val = if(is.null(plus2_NMD_result[[i]][["can_region"]])) NA_integer_ else plus2_NMD_result[[i]][["can_region"]]
  
  # 保护3：关键数值无效，跳过检验
  if(is.na(n_can) || n_can <= 0 || is.na(n_all) || n_all <= 0 || all.variant == 0 ||
     is.null(can_region_val) || is.na(can_region_val) || can_region_val <= 0 ||
     NMDesc_can_count <= 0){
    plus2_results_df <- rbind(plus2_results_df, data.frame(
      transcript = transcript, NMDesc_can_count = NMDesc_can_count,
      can_region = can_region_val, all_variant = all.variant,
      stringsAsFactors = FALSE))
    next
  }
  
  # 二项检验（单侧，富集）
  plus2_NMD_result[[i]][["can_p"]] = binom.test(
    NMDesc_can_count, n_can, all.variant / n_all, alternative = 'greater')$p.value
  
  # Wald检验（单侧，富集）
  if(!is.na(n_rest) && n_rest > 0 && rest.PTC >= 0 && rest.PTC <= n_rest){
    a <- NMDesc_can_count + 0.5
    b <- n_can - NMDesc_can_count + 0.5
    c <- rest.PTC + 0.5
    d <- n_rest - rest.PTC + 0.5
    logOR    <- log(a * d / (b * c))
    se_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
    z_stat   <- logOR / se_logOR
    plus2_NMD_result[[i]][["can_wald_logOR"]]  <- logOR
    plus2_NMD_result[[i]][["can_wald_se"]]     <- se_logOR
    plus2_NMD_result[[i]][["can_wald_p"]]      <- pnorm(z_stat, lower.tail = FALSE)
    plus2_NMD_result[[i]][["can_wald_CI_low"]] <- logOR - 1.96 * se_logOR
    plus2_NMD_result[[i]][["can_wald_CI_up"]]  <- logOR + 1.96 * se_logOR
  }
  
  plus2_results_df <- rbind(plus2_results_df, data.frame(
    transcript       = transcript,
    NMDesc_can_count = NMDesc_can_count,
    can_region       = can_region_val,
    all_variant      = all.variant,
    stringsAsFactors = FALSE
  ))
}
saveRDS(plus2_NMD_result, file = "plus2_NMD_result20260201.rds")
plus2_NMD_result = readRDS("plus2_NMD_result20260201.rds")
#7. for each transcript, use fisher's combined probability to combine plus1 and plus2 p value, then use p_combined to get enrichment gene list for frameshift variants
#combine p value
combine_fisher = function(pvals) {
  pvals = pvals[!is.na(pvals) & is.finite(pvals)]
  if (length(pvals) == 0) return(NA_real_)
  # 防止 p=0 导致 log(0)
  pvals[pvals == 0] <- .Machine$double.xmin
  stat = -2 * sum(log(pvals))
  df = 2 * length(pvals)
  pchisq(stat, df = df, lower.tail = FALSE)
}

tx_all = unique(c(names(plus1_NMD_result), names(plus2_NMD_result)))

p1 = sapply(tx_all, function(tx) {
  if (!is.null(plus1_NMD_result[[tx]]) && !is.null(plus1_NMD_result[[tx]][["can_p"]])) {
    plus1_NMD_result[[tx]][["can_p"]]
  } else {
    NA_real_
  }
})
p2 = sapply(tx_all, function(tx) {
  if (!is.null(plus2_NMD_result[[tx]]) && !is.null(plus2_NMD_result[[tx]][["can_p"]])) {
    plus2_NMD_result[[tx]][["can_p"]]
  } else {
    NA_real_
  }
})
combined_df = data.frame(
  transcript = tx_all,
  can_p_plus1 = as.numeric(p1),
  can_p_plus2 = as.numeric(p2),
  stringsAsFactors = FALSE
)
combined_df$fisher_p <- apply(
  combined_df[, c("can_p_plus1","can_p_plus2")],
  1,
  combine_fisher
)

# ✅ 新增：BH FDR 校正
combined_df$fisher_p_adj <- p.adjust(combined_df$fisher_p, method = "BH")

cat('Fisher合并p值 FDR < 0.05 的转录本数:',
    sum(combined_df$fisher_p_adj < 0.05, na.rm = TRUE), '\n')

# ✅ 改用校正后的 p 值筛选（原来用 fisher_p <= 0.05）
can_list_fisher <- combined_df$transcript[
  !is.na(combined_df$fisher_p_adj) & combined_df$fisher_p_adj < 0.05
]

can_gene_fisher <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = can_list_fisher,
  mart = ensembl
)

writeLines(
  c("hgnc_symbol", unique(can_gene_fisher$hgnc_symbol)),
  "fs_can_gene20260201_syn_FDR0.05.txt"
)


fs_can_gene = read.table('fs_can_gene20260201_syn.txt', header = TRUE, stringsAsFactors = FALSE)
#filter for AD genes
omim <- readLines("~/Desktop/backup/autism/data/genemap2.txt")
omim_hg38_data <- read.delim("~/Desktop/backup/autism/data/genemap2.txt", skip = 3, header = TRUE, nrows=length(omim) - 80)
final_omim_hg38_data <- dplyr::select(omim_hg38_data, X..Chromosome, Genomic.Position.Start, Genomic.Position.End, MIM.Number, Approved.Gene.Symbol, Gene.Locus.And.Other.Related.Symbols, Gene.Name, Entrez.Gene.ID, Ensembl.Gene.ID, Phenotypes, Comments)
#select pathogenic/lp
final_omim_hg38_data = final_omim_hg38_data[which(nchar(final_omim_hg38_data[,'Phenotypes'])>1),]
colnames(final_omim_hg38_data) <- c("chr", "start", "end", "Gene_stable_ID", "approved_gene_symbol", "other_gene_symbol", "gene_name", "entrez_gene_id", "Ensembl_gene_ID", "Phenotypes", "Comments")
final_omim_hg38_data$Inheritance_pattern <- ""
final_omim_hg38_data$Inheritance_pattern[grep("Autosomal dominant", final_omim_hg38_data$Phenotypes)] <- "Autosomal dominant"
final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)] <- paste0(final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)], ", ", "Autosomal recessive")
final_omim_hg38_data$Inheritance_pattern[grep("Digenic recessive", final_omim_hg38_data$Phenotypes)] <- "Digenic recessive"
final_omim_hg38_data$Inheritance_pattern[grep("X-linked recessive", final_omim_hg38_data$Phenotypes)] <- "X-linked recessive"
final_omim_hg38_data$Inheritance_pattern[grep("X-linked dominant", final_omim_hg38_data$Phenotypes)] <- "X-linked dominant"
final_omim_hg38_data$Inheritance_pattern[grep("Y-linked", final_omim_hg38_data$Phenotypes)] <- "Y-linked"
final_omim_hg38_data$Inheritance_pattern[which(final_omim_hg38_data$Inheritance_pattern  == ", Autosomal recessive")] = "Autosomal recessive"
final_omim_hg38_data$Inheritance_pattern[which(final_omim_hg38_data$Inheritance_pattern  == "")] = "Other"
table(final_omim_hg38_data$Inheritance_pattern)
omim2 = final_omim_hg38_data %>% dplyr::select(Ensembl_gene_ID, Inheritance_pattern)
genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","phenotype_description"), filters = "ensembl_gene_id", values = omim2$Ensembl_gene_ID, mart = ensembl)
omim2 = merge(omim2, genes, by.x = "Ensembl_gene_ID", by.y = "ensembl_gene_id")
omim_AD_symbols <- omim2 %>%
  filter(str_detect(Inheritance_pattern, fixed("Autosomal dominant", ignore_case = TRUE))) %>%
  mutate(hgnc_symbol = str_trim(hgnc_symbol)) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol) %>%
  pull(hgnc_symbol)
write.csv(omim_AD_symbols, "omim_AD_symbols.csv", row.names = FALSE)
omim_AD_symbols = read.csv("omim_AD_symbols.csv", header = TRUE, stringsAsFactors = FALSE)
fs_can_AD_gene = fs_can_gene %>% filter(hgnc_symbol %in% omim_AD_symbols)
fs_can_FDR0.05_gene =  data.frame(hgnc_symbol = unlist(unique(can_gene_fisher$hgnc_symbol)))
fs_can_AD_FDR0.05_gene = fs_can_FDR0.05_gene %>% filter(hgnc_symbol %in% omim_AD_symbols$x)
write.csv(fs_can_FDR0.05_gene, "fs_can_FDR0.05_gene.csv", row.names = FALSE)
write.csv(fs_can_AD_FDR0.05_gene, "fs_can_AD_FDR0.05_wald_gene.csv", row.names = FALSE)
length(unique(fs_can_AD_gene$hgnc_symbol))
writeLines(
  c("hgnc_symbol", unique(fs_can_AD_gene$hgnc_symbol)),
  "fs_can_syn_AD_gene_filtered_greater_20260201.txt"
)
fs_can_AD_gene_filtered_greater = read.csv('fs_can_syn_AD_gene_filtered_greater_20260201.txt', header = TRUE, stringsAsFactors = FALSE)
#compare gene list from synonymous variants with gene list from region length
ppi_AD_genes <- read_csv("~/Desktop/frameshift/ppi_AD_genes.csv")
old_fs_AD = ppi_AD_genes %>% filter(group %in% c('Minus1','Plus1')) %>% dplyr::select(hgnc_symbol) 
length(unique(old_fs_AD$hgnc_symbol))

length(fs_gene)
length(fs_can_AD_FDR0.05_gene$hgnc_symbol)
length(intersect(fs_gene,fs_can_AD_FDR0.05_gene$hgnc_symbol))

length(intersect(old_fs_AD$hgnc_symbol, fs_can_AD_gene$hgnc_symbol))