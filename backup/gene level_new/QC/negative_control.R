#This R script is to find plp/benign variants for the contorl gene list, 
#perform same analysis as for the case gene list, and compare results

#1. input control gene list
gene_all <- read_csv("Downloads/gene_all0407 (1).csv", show_col_types = FALSE)
snv_control = gene_all %>% filter(group == "snv_control") 
fs_control = gene_all %>% filter(group == "fs_control")

#2. from clinvar, find plp ptc NMDesc variants, the number is 0 due to its defination
snv_res_can = readRDS('/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/snv_plp_ptc_nmdesc_can_filtered20260201.rds')
#get NMDesc variants 
snv_control_clinvar_variants = data.frame(transcript = snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_control$ensembl_transcript_id)],
                                          key = snv_res_can@elementMetadata@listData[["key"]][which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_control$ensembl_transcript_id)])

fs_res = read_rds('/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/fs.rds')
fs_res_can = fs_res[which(fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
fs_control_clinvar_variants = data.frame(transcript = fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_control$ensembl_transcript_id)],
                                         key = fs_res_can@elementMetadata@listData[["key"]][which(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_control$ensembl_transcript_id)])


#3. from gnomad, find (benign) ptc NMDesc variants
ptc_can_NMD_df <- read_csv("Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/gnomad/snv_fs/ptc_can_NMD_df.csv")
snv_control_gnomad = ptc_can_NMD_df[which(ptc_can_NMD_df$transcript %in% snv_control$ensembl_transcript_id & ptc_can_NMD_df$type == 'snv'),]
fs_control_gnomad = ptc_can_NMD_df[which(ptc_can_NMD_df$transcript %in% fs_control$ensembl_transcript_id & ptc_can_NMD_df$type != 'snv'),]
#remove inframe frameshift variants based on the key
# 从 id 列提取 ref 和 alt，计算长度差，排除3的倍数（inframe）
fs_control_gnomad$ref = sub(".*\\|(.+)\\|.*", "\\1",fs_control_gnomad$id)
fs_control_gnomad$alt = sub(".*\\|.*\\|(.*)", "\\1", fs_control_gnomad$id)
fs_control_gnomad$len_diff = abs(nchar(fs_control_gnomad$ref) - nchar(fs_control_gnomad$alt))
# 保留长度差不是3的倍数的（真正的frameshift）
fs_control_gnomad_filtered = fs_control_gnomad[fs_control_gnomad$len_diff %% 3 != 0, ]

##4.1 compare ptc distance to cds end
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 108)
BM.info4 = getBM(
  attributes = c("ensembl_transcript_id","rank",'cds_start','cds_end','exon_chrom_start','exon_chrom_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  snv_control$ensembl_transcript_id,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
BM.info4 = BM.info4[order(BM.info4$rank),]
get_re_loc <- function(location, transcript_match){
  df <- BM.info4[BM.info4$ensembl_transcript_id == transcript_match, ]
  
  # 通过坐标判断 strand：rank1 比 rank2 大 → 负链（对 MN1 是这样）
  is_minus <- nrow(df) >= 2 && df$exon_chrom_start[1] > df$exon_chrom_start[2]
  
  for (i in seq_len(nrow(df))) {
    start <- df$exon_chrom_start[i]
    end   <- df$exon_chrom_end[i]
    # location 是否落在 exon 区间（不管 start/end 谁大都能判断）
    if (location >= min(start, end) && location <= max(start, end)) {
      if (is_minus) {
        offset <- max(start, end) - location
      } else {
        offset <- location - min(start, end)
      }
      return(df$cds_start[i] + offset)
    }
  }
  
  return(NA)
}
snv_control_gnomad <- snv_control_gnomad %>%
  mutate(
    mutation_loc = as.numeric(sub(".*:(\\d+)\\|.*", "\\1", id))
  )
snv_control_gnomad$cds_mutation_loc = mapply(get_re_loc, snv_control_gnomad$mutation_loc, snv_control_gnomad$transcript)
snv_control_gnomad$ptc_loc = snv_control_gnomad$cds_mutation_loc + 1
  #cds length should be max(cds_end)
snv_control_gnomad <- snv_control_gnomad %>%
  left_join(
    BM.info4 %>%
      group_by(ensembl_transcript_id) %>%
      summarise(cds_length = max(cds_end, na.rm = TRUE), .groups = "drop"),
    by = c("transcript" = "ensembl_transcript_id")
  )
snv_control_gnomad$ptc_distance_to_cds_end = snv_control_gnomad$cds_length - snv_control_gnomad$ptc_loc

BM.info4 = getBM(
  attributes = c("ensembl_transcript_id","rank",'cds_start','cds_end','exon_chrom_start','exon_chrom_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  fs_control$ensembl_transcript_id,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
BM.info4 = BM.info4[order(BM.info4$rank),]
fs_control_gnomad <- fs_control_gnomad %>%
  mutate(
    mutation_loc = as.numeric(sub(".*:(\\d+)\\|.*", "\\1", id))
  )
fs_control_gnomad$cds_mutation_loc = mapply(get_re_loc, fs_control_gnomad$mutation_loc, fs_control_gnomad$transcript)
#cds length
fs_control_gnomad <- fs_control_gnomad %>%
  left_join(
    BM.info4 %>%
      group_by(ensembl_transcript_id) %>%
      summarise(cds_length = max(cds_end, na.rm = TRUE), .groups = "drop"),
    by = c("transcript" = "ensembl_transcript_id")
  )

library(dplyr)
library(stringr)
library(tidyr)
library(biomaRt)

bp <- 50

transcript_set_fs <- unique(fs_control_gnomad_filtered$transcript)

can.info <- getBM(
  attributes = c("ensembl_transcript_id", "transcript_is_canonical"),
  filters = "ensembl_transcript_id",
  values = transcript_set_fs,
  mart = ensembl
)

## 只保留 canonical transcript
transcript_set2_fs <- can.info %>%
  filter(transcript_is_canonical == 1) %>%
  pull(ensembl_transcript_id) %>%
  unique()

BM.infoo_fs <- getBM(
  attributes = c(
    "ensembl_transcript_id", "rank", "cds_start", "cds_end",
    "exon_chrom_start", "exon_chrom_end"
  ),
  filters = "ensembl_transcript_id",
  values = transcript_set2_fs,
  mart = ensembl
)

## 只保留 coding exon 数 > 1 的 transcript
exon.num_fs <- BM.infoo_fs %>%
  group_by(ensembl_transcript_id) %>%
  summarise(max_rank = sum(!is.na(cds_start)), .groups = "drop")

transcript_set3_fs <- exon.num_fs %>%
  filter(max_rank > 1) %>%
  pull(ensembl_transcript_id) %>%
  unique()

## coding sequence
cds.info_fs <- getBM(
  attributes = c("ensembl_transcript_id", "coding"),
  filters = "ensembl_transcript_id",
  values = transcript_set3_fs,
  mart = ensembl
)

## 只保留这些 transcript 的 exon 信息
BM.info2_fs <- BM.infoo_fs %>%
  filter(ensembl_transcript_id %in% transcript_set3_fs)

BM.info2_fs$exon_length <- BM.info2_fs$cds_end - BM.info2_fs$cds_start + 1

## 只保留真正 coding 的 exon
BM.info3_fs <- BM.info2_fs %>%
  filter(!is.na(exon_length), exon_length > 0)

## 注意：这里 cds_length 应该用 sum(exon_length)，不是 max(cds_end)
BM.info4_fs_len <- BM.info3_fs %>%
  group_by(ensembl_transcript_id) %>%
  summarise(cds_length = sum(exon_length), .groups = "drop")

BM.info4_fs <- BM.info3_fs %>%
  left_join(BM.info4_fs_len, by = "ensembl_transcript_id") %>%
  arrange(ensembl_transcript_id, rank)
get_PTC_plus_ind <- function(input_seqs, type = "plus1") {
  para <- ifelse(type == "plus1", 2, 0)
  a <- unlist(input_seqs)
  ind <- unlist(gregexpr("TAG|TAA|TGA", as.character(a)))
  plus.ind <- ind[which(ind %% 3 == para)]
  paste(plus.ind, collapse = ",")
}

get_NMDesc_PTC <- function(input_seqs, bp, variant_transcript, type = "plus1") {
  para <- ifelse(type == "plus1", 2, 0)
  a <- unlist(input_seqs)
  ind <- unlist(gregexpr("TAG|TAA|TGA", as.character(a)))
  plus.ind <- ind[which(ind %% 3 == para)]
  
  df <- BM.info4_fs[BM.info4_fs$ensembl_transcript_id == variant_transcript, ]
  df <- df[!is.na(df$cds_start), ]
  
  exon.length <- (df$cds_end - df$cds_start) + 1
  
  pen.length <- exon.length[length(exon.length) - 1]
  if (pen.length < bp) {
    aft.NMD.ind <- df$cds_start[length(exon.length) - 1]
  } else {
    aft.NMD.ind <- sum(exon.length[1:(length(exon.length) - 1)]) - bp
  }
  
  NMD_status_can <- plus.ind >= aft.NMD.ind
  NMD_status_css <- plus.ind <= 150
  
  exon.long <- which(exon.length > 407)
  exon.long <- exon.long[exon.long != length(exon.length)]
  exon.long <- exon.long[exon.long != 1]
  
  long.start <- df$cds_start[exon.long]
  long.end   <- df$cds_end[exon.long]
  
  if (length(exon.long) == 0) {
    long.start <- 0
    long.end <- 0
  }
  
  NMD_status_long <- rep(FALSE, length(plus.ind))
  for (k in seq_along(plus.ind)) {
    NMD_status_long[k] <- any(plus.ind[k] >= long.start & plus.ind[k] <= long.end)
  }
  
  NMD_status <- c(NMD_status_can, NMD_status_css, NMD_status_long)
  paste(NMD_status, collapse = ",")
}

get_NMDesc_PTC_region <- function(PTC_ind, variant_transcript, PTC_status) {
  PTC_status_vector <- as.logical(strsplit(PTC_status, ",")[[1]])
  PTC_ind_vector <- as.numeric(unlist(strsplit(PTC_ind, ",")))
  
  jump <- length(PTC_status_vector) / 3
  
  can_PTC_status_vector  <- PTC_status_vector[1:jump]
  css_PTC_status_vector  <- PTC_status_vector[(jump + 1):(2 * jump)]
  long_PTC_status_vector <- PTC_status_vector[(2 * jump + 1):(3 * jump)]
  
  can_change_points <- which(diff(can_PTC_status_vector) != 0)
  css_change_points <- which(diff(css_PTC_status_vector) != 0)
  long_change_points <- which(diff(long_PTC_status_vector) != 0)
  
  can_change_points_loc <- PTC_ind_vector[can_change_points]
  css_change_points_loc <- PTC_ind_vector[css_change_points]
  long_change_points_loc <- PTC_ind_vector[long_change_points]
  
  long_FT_change_points <- long_change_points_loc[long_PTC_status_vector[long_change_points] == FALSE]
  long_TF_change_points <- long_change_points_loc[long_PTC_status_vector[long_change_points] == TRUE]
  
  can_region  <- max((PTC_ind_vector[length(PTC_ind_vector)] - (can_change_points_loc + 1)), 0)
  css_region  <- max((css_change_points_loc - 1), 0)
  long_region <- max(long_TF_change_points - (long_FT_change_points + 1), 0)
  
  df <- BM.info4_fs[BM.info4_fs$ensembl_transcript_id == variant_transcript, ]
  df <- df[!is.na(df$cds_start), ]
  
  exon.length <- (df$cds_end - df$cds_start) + 1
  pen.length <- exon.length[length(exon.length) - 1]
  
  if (pen.length > 407) {
    long_region <- long_region - 150
  }
  
  if (length(long_region) == 0 || is.na(long_region)) long_region <- 0
  if (length(css_region) == 0 || is.na(css_region)) css_region <- 0
  if (length(can_region) == 0 || is.na(can_region)) can_region <- 0
  
  list(
    can_region = can_region,
    css_region = css_region,
    long_region = long_region,
    can_region_start = max(can_change_points_loc, 0) + 1,
    can_region_end = max(PTC_ind_vector[length(PTC_ind_vector)], 0)
  )
}
PTC_info_fs <- data.frame(
  PTC_loc = rep(NA, 2 * length(transcript_set3_fs)),
  type = rep(NA, 2 * length(transcript_set3_fs)),
  transcript = rep(NA, 2 * length(transcript_set3_fs)),
  PTC_NMDesc = rep(NA, 2 * length(transcript_set3_fs)),
  can_region = rep(NA, 2 * length(transcript_set3_fs)),
  css_region = rep(NA, 2 * length(transcript_set3_fs)),
  long_region = rep(NA, 2 * length(transcript_set3_fs)),
  can_region_start = rep(NA, 2 * length(transcript_set3_fs)),
  can_region_end = rep(NA, 2 * length(transcript_set3_fs)),
  stringsAsFactors = FALSE
)

counter <- 1

for (i in seq_along(transcript_set3_fs)) {
  transcript <- transcript_set3_fs[i]
  
  for (j in c("plus1", "plus2")) {
    input_seqs <- cds.info_fs[cds.info_fs$ensembl_transcript_id == transcript, "coding"]
    
    PTC_ind <- get_PTC_plus_ind(input_seqs, type = j)
    PTC_status <- get_NMDesc_PTC(input_seqs, bp, transcript, type = j)
    region_set <- get_NMDesc_PTC_region(PTC_ind, transcript, PTC_status)
    
    PTC_info_fs[counter, "PTC_loc"] <- PTC_ind
    PTC_info_fs[counter, "type"] <- j
    PTC_info_fs[counter, "transcript"] <- transcript
    PTC_info_fs[counter, "PTC_NMDesc"] <- PTC_status
    PTC_info_fs[counter, "can_region"] <- region_set$can_region
    PTC_info_fs[counter, "css_region"] <- region_set$css_region
    PTC_info_fs[counter, "long_region"] <- region_set$long_region
    PTC_info_fs[counter, "can_region_start"] <- region_set$can_region_start
    PTC_info_fs[counter, "can_region_end"] <- region_set$can_region_end
    
    counter <- counter + 1
  }
}

PTC_info_fs <- PTC_info_fs %>%
  filter(!is.na(transcript), !is.na(type))
get_re_loc_fs <- function(location, transcript_match) {
  df <- BM.info4_fs[BM.info4_fs$ensembl_transcript_id == transcript_match, ]
  
  is_minus <- nrow(df) >= 2 && df$exon_chrom_start[1] > df$exon_chrom_start[2]
  
  for (i in seq_len(nrow(df))) {
    start <- df$exon_chrom_start[i]
    end   <- df$exon_chrom_end[i]
    
    if (location >= min(start, end) && location <= max(start, end)) {
      if (is_minus) {
        offset <- max(start, end) - location
      } else {
        offset <- location - min(start, end)
      }
      return(df$cds_start[i] + offset)
    }
  }
  
  return(NA_real_)
}

get_frameshift_type_from_id <- function(id) {
  parts <- str_split(id, "\\|", simplify = TRUE)
  ref <- parts[, 2]
  alt <- parts[, 3]
  yu <- abs(nchar(ref) - nchar(alt)) %% 3
  
  ifelse(
    (nchar(ref) > nchar(alt) & yu == 1) | (nchar(ref) < nchar(alt) & yu == 2),
    "plus1",
    ifelse(
      (nchar(ref) > nchar(alt) & yu == 2) | (nchar(ref) < nchar(alt) & yu == 1),
      "plus2",
      NA_character_
    )
  )
}

get_fs_ptc_loc <- function(variant_loc_rel, variant_transcript, fs_type, PTC_info_obj) {
  ptc_string <- PTC_info_obj[
    PTC_info_obj$transcript == variant_transcript & PTC_info_obj$type == fs_type,
    "PTC_loc"
  ]
  
  if (length(ptc_string) == 0 || is.na(ptc_string[1]) || ptc_string[1] == "") {
    return(NA_real_)
  }
  
  ptc_vec <- as.numeric(unlist(strsplit(ptc_string[1], ",")))
  ptc_vec <- ptc_vec[!is.na(ptc_vec)]
  
  if (length(ptc_vec) == 0 || is.na(variant_loc_rel)) {
    return(NA_real_)
  }
  
  downstream_ptc <- ptc_vec[ptc_vec > variant_loc_rel]
  
  if (length(downstream_ptc) == 0) {
    return(NA_real_)
  }
  
  min(downstream_ptc)
}

fs_control_gnomad_filtered2 <- fs_control_gnomad_filtered %>%
  filter(transcript %in% transcript_set3_fs) %>%
  mutate(
    mutation_loc = as.numeric(sub(".*:(\\d+)\\|.*", "\\1", id)),
    fs_type = get_frameshift_type_from_id(id)
  )

fs_control_gnomad_filtered2$cds_mutation_loc <- mapply(
  get_re_loc_fs,
  fs_control_gnomad_filtered2$mutation_loc,
  fs_control_gnomad_filtered2$transcript
)

fs_control_gnomad_filtered2 <- fs_control_gnomad_filtered2 %>%
  left_join(
    BM.info4_fs_len,
    by = c("transcript" = "ensembl_transcript_id")
  )

fs_control_gnomad_filtered2$ptc_loc <- mapply(
  get_fs_ptc_loc,
  variant_loc_rel = fs_control_gnomad_filtered2$cds_mutation_loc,
  variant_transcript = fs_control_gnomad_filtered2$transcript,
  fs_type = fs_control_gnomad_filtered2$fs_type,
  MoreArgs = list(PTC_info_obj = PTC_info_fs)
)

fs_control_gnomad_filtered2$ptc_distance_to_cds_end <-
  fs_control_gnomad_filtered2$cds_length - fs_control_gnomad_filtered2$ptc_loc

library(dplyr)

set.seed(123)

n_resample <- 100

resampling_ptc_dist <- lapply(seq_len(n_resample), function(i) {
  
  sampled_df <- control_variants_all %>%
    group_by(group, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  sampled_df %>%
    group_by(group) %>%
    summarise(
      median_dist = median(ptc_distance_to_cds_end, na.rm = TRUE),
      mean_dist   = mean(ptc_distance_to_cds_end, na.rm = TRUE),
      n           = n(),
      n_na        = sum(is.na(ptc_distance_to_cds_end)),
      .groups = "drop"
    ) %>%
    mutate(resample_id = i)
})

resampling_ptc_dist_df <- bind_rows(resampling_ptc_dist)

resampling_ptc_dist_df

ptc_dist_summary <- resampling_ptc_dist_df %>%
  group_by(group) %>%
  summarise(
    mean_of_mean_dist = mean(mean_dist, na.rm = TRUE),
    sd_of_mean_dist = sd(mean_dist, na.rm = TRUE),
    median_of_median_dist = median(median_dist, na.rm = TRUE),
    mean_of_median_dist = mean(median_dist, na.rm = TRUE),
    sd_of_median_dist = sd(median_dist, na.rm = TRUE),
    mean_n = mean(n, na.rm = TRUE),
    mean_n_na = mean(n_na, na.rm = TRUE),
    .groups = "drop"
  )

ptc_dist_summary

##4.2 compare the motif/LCS features of the control gene variants with the case gene variants, and see if there is a significant difference
fs_control_gnomad_filtered3 = fs_control_gnomad_filtered2[,c(1,2,3,7,9,10,11,12)]


control_variants_all = rbind(fs_control_gnomad_filtered3, snv_control_gnomad)
control_variants_all$group = ifelse(control_variants_all$type == 'snv', 'snv_control2', 'fs_control2')
#add uniprot id by transcript
uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = control_variants_all$transcript,
  mart = ensembl
)
control_variants_all <- control_variants_all %>%
  left_join(uniprot_mapping, by = c("transcript" = "ensembl_transcript_id"))
#rename it to uniprot
control_variants_all <- control_variants_all %>%
  rename(uniprot = uniprotswissprot)
merge(motif_max,control_variants_all, by = "uniprot") -> motif_max3
LCS_max3 = merge(LCS_max, control_variants_all, by = "uniprot")

motif_max3$variant_protein_feature_flag = motif_max3$`Protein Features`*3 >= motif_max3$cds_mutation_loc
motif_max3$variant_domains_flag = motif_max3$`Domains`*3 >= motif_max3$cds_mutation_loc
motif_max3$variant_slim_flag = motif_max3$`SLiMs`*3 >= motif_max3$cds_mutation_loc #variant location on cds scale
motif_max3$variant_morf_flag = motif_max3$MORFs*3 >= motif_max3$cds_mutation_loc
motif_max3$variant_ptm_flag = motif_max3$`PTMs`*3 >= motif_max3$cds_mutation_loc
motif_max3$variant_nls_flag = motif_max3$`NLSs`*3 >= motif_max3$cds_mutation_loc
LCS_max3$variant_LCS_flag = LCS_max3$`LCSs`*3 >= LCS_max3$cds_mutation_loc

control_variants_all$variant_protein_flag = motif_max3$variant_protein_feature_flag[match(control_variants_all$uniprot, motif_max3$uniprot)]
control_variants_all$variant_domain_flag = motif_max3$variant_domains_flag[match(control_variants_all$uniprot, motif_max3$uniprot)]
control_variants_all$variant_slim_flag = motif_max3$variant_slim_flag[match(control_variants_all$uniprot, motif_max3$uniprot)]
control_variants_all$variant_morf_flag = motif_max3$variant_morf_flag[match(control_variants_all$uniprot, motif_max3$uniprot)]
control_variants_all$variant_ptm_flag = motif_max3$variant_ptm_flag[match(control_variants_all$uniprot, motif_max3$uniprot)]
control_variants_all$variant_nls_flag = motif_max3$variant_nls_flag[match(control_variants_all$uniprot, motif_max3$uniprot)]
control_variants_all$variant_LCS_flag = LCS_max3$variant_LCS_flag[match(control_variants_all$uniprot, LCS_max3$uniprot)]

write.csv(control_variants_all, "control_variants_all_with_motif_LCS_flags.csv", row.names = FALSE)

##4.3 compare PPI and pfam
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(biomaRt)

##################################################
## 0. input assumptions
##################################################
## Required objects already in environment:
## 1) control_variants_all
## 2) human_1_
## 3) ensembl
##
## Required columns in control_variants_all:
## transcript, id, group, uniprot, ptc_loc, ptc_distance_to_cds_end
##
## human_1_ should contain:
## uniprot1, uniprot2, interface_residues1, interface_residues2

##################################################
## 1. helper function: convert interface residue string to numeric vector
##################################################
convert_to_c <- function(x) {
  if (is.na(x) || x == "") return(numeric(0))
  
  x <- gsub("\\[|\\]", "", x)
  x <- gsub(" ", "", x)
  
  if (x == "") return(numeric(0))
  
  vals <- unlist(strsplit(x, ","))
  vals <- suppressWarnings(as.numeric(vals))
  vals <- vals[!is.na(vals)]
  
  return(vals)
}

##################################################
## 2. copy input data
##################################################
control_variants_all2 <- control_variants_all %>%
  mutate(
    variant_ppi_overlap = 0,
    variant_ppi_nearest_interface_bp = NA_real_,
    variant_ppi_dist_to_nearest_interface_bp = NA_real_
  )

##################################################
## 3. generate PPI data
## logic follows your script:
## compare interface residues (aa -> bp) with ptc_loc
##################################################
for (i in seq_len(nrow(control_variants_all2))) {
  
  uid <- control_variants_all2$uniprot[i]
  cds_ptc_loc <- control_variants_all2$ptc_loc[i]
  
  if (is.na(uid) || uid == "") next
  if (is.na(cds_ptc_loc)) next
  
  re_1 <- unlist(lapply(
    human_1_$interface_residues1[human_1_$uniprot1 == uid],
    convert_to_c
  )) * 3
  
  re_2 <- unlist(lapply(
    human_1_$interface_residues2[human_1_$uniprot2 == uid],
    convert_to_c
  )) * 3
  
  all_interface_bp <- c(re_1, re_2)
  all_interface_bp <- all_interface_bp[!is.na(all_interface_bp)]
  
  if (length(all_interface_bp) == 0) next
  
  downstream_interface_bp <- all_interface_bp[all_interface_bp >= cds_ptc_loc]
  
  if (length(downstream_interface_bp) > 0) {
    control_variants_all2$variant_ppi_overlap[i] <- 1
    control_variants_all2$variant_ppi_nearest_interface_bp[i] <- min(downstream_interface_bp, na.rm = TRUE)
    control_variants_all2$variant_ppi_dist_to_nearest_interface_bp[i] <-
      min(downstream_interface_bp, na.rm = TRUE) - cds_ptc_loc
  }
}

##################################################
## 4. basic check for PPI
##################################################
table(control_variants_all2$variant_ppi_overlap, useNA = "ifany")
summary(control_variants_all2$variant_ppi_dist_to_nearest_interface_bp)

##################################################
## 5. generate PFAM data
## logic follows your script:
## ptc_aa = ceiling(ptc_loc / 3)
## then compare with max pfam_end per transcript
##################################################
pfam_domains2 <- getBM(
  attributes = c("ensembl_transcript_id", "pfam_end"),
  filters = "ensembl_transcript_id",
  values = unique(control_variants_all2$transcript),
  mart = ensembl
)

max_pfam_end_df <- pfam_domains2 %>%
  filter(!is.na(pfam_end)) %>%
  group_by(ensembl_transcript_id) %>%
  summarise(max_pfam_end = max(pfam_end), .groups = "drop")

control_variants_all2 <- control_variants_all2 %>%
  mutate(
    ptc_aa = ceiling(ptc_loc / 3)
  ) %>%
  left_join(
    max_pfam_end_df,
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  mutate(
    dist_ptc_to_max_pfam_end_aa = ptc_aa - max_pfam_end,
    ptc_after_max_pfam_end = ifelse(
      !is.na(dist_ptc_to_max_pfam_end_aa) & dist_ptc_to_max_pfam_end_aa > 0, 1, 0
    ),
    ptc_before_max_pfam_end = ifelse(
      !is.na(dist_ptc_to_max_pfam_end_aa) & dist_ptc_to_max_pfam_end_aa < 0, 1, 0
    )
  )

##################################################
## 6. basic check for PFAM
##################################################
summary(control_variants_all2$dist_ptc_to_max_pfam_end_aa)
table(control_variants_all2$ptc_after_max_pfam_end, useNA = "ifany")
table(control_variants_all2$ptc_before_max_pfam_end, useNA = "ifany")

##################################################
## 7. save full data with PPI + PFAM columns
##################################################
write.csv(
  control_variants_all2,
  "control_variants_all_with_ppi_pfam_data.csv",
  row.names = FALSE
)

##################################################
## 8. group-level summary without resampling
##################################################
ppi_summary <- control_variants_all2 %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_ppi_overlap = sum(variant_ppi_overlap == 1, na.rm = TRUE),
    prop_ppi_overlap = mean(variant_ppi_overlap == 1, na.rm = TRUE),
    mean_ppi_dist = mean(variant_ppi_dist_to_nearest_interface_bp, na.rm = TRUE),
    median_ppi_dist = median(variant_ppi_dist_to_nearest_interface_bp, na.rm = TRUE),
    n_ppi_dist_na = sum(is.na(variant_ppi_dist_to_nearest_interface_bp)),
    .groups = "drop"
  )

pfam_summary <- control_variants_all2 %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_before_pfam_end = sum(ptc_before_max_pfam_end == 1, na.rm = TRUE),
    prop_before_pfam_end = mean(ptc_before_max_pfam_end == 1, na.rm = TRUE),
    n_after_pfam_end = sum(ptc_after_max_pfam_end == 1, na.rm = TRUE),
    prop_after_pfam_end = mean(ptc_after_max_pfam_end == 1, na.rm = TRUE),
    mean_pfam_dist = mean(dist_ptc_to_max_pfam_end_aa, na.rm = TRUE),
    median_pfam_dist = median(dist_ptc_to_max_pfam_end_aa, na.rm = TRUE),
    n_pfam_dist_na = sum(is.na(dist_ptc_to_max_pfam_end_aa)),
    .groups = "drop"
  )

print(ppi_summary)
print(pfam_summary)

##################################################
## 9. 100 resampling:
## choose 1 variant per transcript per group
##################################################
set.seed(123)
n_resample <- 100

resampled_list <- lapply(seq_len(n_resample), function(i) {
  control_variants_all2 %>%
    group_by(group, transcript) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    mutate(resample_id = i)
})

resampled_all <- bind_rows(resampled_list)

write.csv(
  resampled_all,
  "control_variants_all_resampled_100.csv",
  row.names = FALSE
)

##################################################
## 10. PPI resampling summary
##################################################
ppi_resampling_results <- resampled_all %>%
  group_by(resample_id, group) %>%
  summarise(
    n = n(),
    n_ppi_overlap = sum(variant_ppi_overlap == 1, na.rm = TRUE),
    prop_ppi_overlap = mean(variant_ppi_overlap == 1, na.rm = TRUE),
    mean_ppi_dist = mean(variant_ppi_dist_to_nearest_interface_bp, na.rm = TRUE),
    median_ppi_dist = median(variant_ppi_dist_to_nearest_interface_bp, na.rm = TRUE),
    n_ppi_dist_na = sum(is.na(variant_ppi_dist_to_nearest_interface_bp)),
    .groups = "drop"
  )

ppi_resampling_summary_100 <- ppi_resampling_results %>%
  group_by(group) %>%
  summarise(
    mean_prop_ppi_overlap = mean(prop_ppi_overlap, na.rm = TRUE),
    sd_prop_ppi_overlap = sd(prop_ppi_overlap, na.rm = TRUE),
    mean_mean_ppi_dist = mean(mean_ppi_dist, na.rm = TRUE),
    sd_mean_ppi_dist = sd(mean_ppi_dist, na.rm = TRUE),
    mean_median_ppi_dist = mean(median_ppi_dist, na.rm = TRUE),
    sd_median_ppi_dist = sd(median_ppi_dist, na.rm = TRUE),
    mean_n = mean(n, na.rm = TRUE),
    mean_n_ppi_dist_na = mean(n_ppi_dist_na, na.rm = TRUE),
    .groups = "drop"
  )

print(ppi_resampling_results)
print(ppi_resampling_summary_100)

write.csv(
  ppi_resampling_results,
  "control_variants_all_ppi_resampling_results.csv",
  row.names = FALSE
)

write.csv(
  ppi_resampling_summary_100,
  "control_variants_all_ppi_resampling_summary_100.csv",
  row.names = FALSE
)

##################################################
## 11. PFAM resampling summary
##################################################
pfam_resampling_results <- resampled_all %>%
  group_by(resample_id, group) %>%
  summarise(
    n = n(),
    n_before_pfam_end = sum(ptc_before_max_pfam_end == 1, na.rm = TRUE),
    prop_before_pfam_end = mean(ptc_before_max_pfam_end == 1, na.rm = TRUE),
    n_after_pfam_end = sum(ptc_after_max_pfam_end == 1, na.rm = TRUE),
    prop_after_pfam_end = mean(ptc_after_max_pfam_end == 1, na.rm = TRUE),
    mean_pfam_dist = mean(dist_ptc_to_max_pfam_end_aa, na.rm = TRUE),
    median_pfam_dist = median(dist_ptc_to_max_pfam_end_aa, na.rm = TRUE),
    n_pfam_dist_na = sum(is.na(dist_ptc_to_max_pfam_end_aa)),
    .groups = "drop"
  )

pfam_resampling_summary_100 <- pfam_resampling_results %>%
  group_by(group) %>%
  summarise(
    mean_prop_before_pfam_end = mean(prop_before_pfam_end, na.rm = TRUE),
    sd_prop_before_pfam_end = sd(prop_before_pfam_end, na.rm = TRUE),
    mean_prop_after_pfam_end = mean(prop_after_pfam_end, na.rm = TRUE),
    sd_prop_after_pfam_end = sd(prop_after_pfam_end, na.rm = TRUE),
    mean_mean_pfam_dist = mean(mean_pfam_dist, na.rm = TRUE),
    sd_mean_pfam_dist = sd(mean_pfam_dist, na.rm = TRUE),
    mean_median_pfam_dist = mean(median_pfam_dist, na.rm = TRUE),
    sd_median_pfam_dist = sd(median_pfam_dist, na.rm = TRUE),
    mean_n = mean(n, na.rm = TRUE),
    mean_n_pfam_dist_na = mean(n_pfam_dist_na, na.rm = TRUE),
    .groups = "drop"
  )

print(pfam_resampling_results)
print(pfam_resampling_summary_100)

write.csv(
  pfam_resampling_results,
  "control_variants_all_pfam_resampling_results.csv",
  row.names = FALSE
)

write.csv(
  pfam_resampling_summary_100,
  "control_variants_all_pfam_resampling_summary_100.csv",
  row.names = FALSE
)

##################################################
## 12. combine PPI + PFAM resampling summary
##################################################
combined_resampling_summary_100 <- ppi_resampling_summary_100 %>%
  full_join(pfam_resampling_summary_100, by = c("group", "mean_n"))

print(combined_resampling_summary_100)

write.csv(
  combined_resampling_summary_100,
  "control_variants_all_combined_ppi_pfam_resampling_summary_100.csv",
  row.names = FALSE
)

##################################################
## 13. example plots
##################################################

## PPI overlap proportion across 100 resamplings
p_ppi_prop <- ggplot(
  ppi_resampling_results,
  aes(x = group, y = prop_ppi_overlap, fill = group)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, color = "black") +
  theme_classic() +
  labs(
    title = "PPI overlap proportion across 100 resamplings",
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(p_ppi_prop)

## PPI distance across 100 resamplings
p_ppi_dist <- ggplot(
  ppi_resampling_results,
  aes(x = group, y = mean_ppi_dist, fill = group)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, color = "black") +
  theme_classic() +
  labs(
    title = "Mean distance from PTC to nearest downstream PPI interface",
    x = NULL,
    y = "Mean distance (bp)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(p_ppi_dist)

## PFAM before proportion across 100 resamplings
p_pfam_prop <- ggplot(
  pfam_resampling_results,
  aes(x = group, y = prop_before_pfam_end, fill = group)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, color = "black") +
  theme_classic() +
  labs(
    title = "Proportion of variants with PTC before max PFAM end",
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(p_pfam_prop)

## PFAM distance across 100 resamplings
p_pfam_dist <- ggplot(
  pfam_resampling_results,
  aes(x = group, y = mean_pfam_dist, fill = group)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, color = "black") +
  theme_classic() +
  labs(
    title = "Mean distance from PTC to max PFAM end",
    x = NULL,
    y = "Mean distance (aa)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(p_pfam_dist)

##################################################
## 14. optional quick view
##################################################
control_variants_all2 %>%
  select(
    transcript, id, group, uniprot, ptc_loc, ptc_distance_to_cds_end,
    variant_ppi_overlap,
    variant_ppi_nearest_interface_bp,
    variant_ppi_dist_to_nearest_interface_bp,
    ptc_aa,
    max_pfam_end,
    dist_ptc_to_max_pfam_end_aa,
    ptc_after_max_pfam_end,
    ptc_before_max_pfam_end
  ) %>%
  head(20)

write.csv(
  control_variants_all2,
  "control_variants_all2_with_ppi_pfam_data.csv",
  row.names = FALSE
)
