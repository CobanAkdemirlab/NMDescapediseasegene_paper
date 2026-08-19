transcript_object <- list()

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


----------
plus1_NMD_result <- lapply(transcript_object, function(transcript) {
    if (!is.null(transcript$plus1)) {
      transcript$plus1
    } else {
      NULL  # 如果没有 plus1，则返回 NULL
    }
  })

#6. get the length of NMDesc variant region, count of variants that is NMDesc, use binomial test to get the p-value
#calculate p value
results_df <- data.frame(
  transcript = character(),
  NMDesc_can_count = integer(),
  can_region = integer(),
  all_variant = integer(),
  tr_length = integer(),
  stringsAsFactors = FALSE
)
for(i in seq_along(plus1_NMD_result)){
  transcript = names(plus1_NMD_result)[i]
  tr_length = BM.info4[which(BM.info4$ensembl_transcript_id == transcript),'cds_length']
  #count of NMDesc variants
  fs_NMD_can = plus1_NMD_result[[i]][["NMDesc_can"]]
  NMDesc_can_count = length(which(fs_NMD_can==TRUE))
  #count of non_NMDesc variants
  non_NMDesc_can_count = length(which(fs_NMD_can==FALSE))
  #all variants
  all.variant = NMDesc_can_count + non_NMDesc_can_count
  #binomial test
  if(!is.null(plus1_NMD_result[[i]][["can_region"]]) && length(tr_length)>0 && length(plus1_NMD_result[[i]][["can_region"]])>0 && !is.na(plus1_NMD_result[[i]][["can_region"]]) && NMDesc_can_count >0 && plus1_NMD_result[[i]][["can_region"]]>0){
    plus1_NMD_result[[i]][["can_p"]] = binom.test(NMDesc_can_count,plus1_NMD_result[[i]][["can_region"]],all.variant/tr_length[1],alternative='less')$p.value
  }else(plus1_NMD_result[[i]][["can_p"]] = NA)
  
  results_df <- rbind(results_df, data.frame(
    transcript = transcript,
    NMDesc_can_count = NMDesc_can_count,
    can_region = plus1_NMD_result[[i]][["can_region"]],
    all_variant = all.variant,
    tr_length = tr_length[1],
    stringsAsFactors = FALSE
  ))
}


#7. use 90% as cutoff, get the clinvar frameshift gene list that is NMDesc variants enriched
can_list = c()

for (i in 1:length(plus1_NMD_result)) {
  if (!is.na(plus1_NMD_result[[i]][["can_p"]]) && plus1_NMD_result[[i]][["can_p"]] >= 0.9) {
    can_list <- c(can_list, names(plus1_NMD_result)[i])
  }
}

##transform transcript to gene
can_gene = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values = can_list,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
writeLines(c("hgnc_symbol", unique(can_gene$hgnc_symbol)), 'plus1_can_gene0217.txt')


---------
  plus2_NMD_result <- lapply(transcript_object, function(transcript) {
    if (!is.null(transcript$plus2)) {
      transcript$plus2
    } else {
      NULL  # 如果没有 plus2，则返回 NULL
    }
  })

#6. get the length of NMDesc variant region, count of variants that is NMDesc, use binomial test to get the p-value
#calculate p value
for(i in seq_along(plus2_NMD_result)){
  transcript = names(plus2_NMD_result)[i]
  tr_length = BM.info4[which(BM.info4$ensembl_transcript_id == transcript),'cds_length']
  #count of NMDesc variants
  fs_NMD_can = plus2_NMD_result[[i]][["NMDesc_can"]]
  NMDesc_can_count = length(which(fs_NMD_can==TRUE))
  #count of non_NMDesc variants
  non_NMDesc_can_count = length(which(fs_NMD_can==FALSE))
  #all variants
  all.variant = NMDesc_can_count + non_NMDesc_can_count
  #binomial test
  if(!is.null(plus2_NMD_result[[i]][["can_region"]]) && length(tr_length)>0 && length(plus2_NMD_result[[i]][["can_region"]])>0 && !is.na(plus2_NMD_result[[i]][["can_region"]]) && NMDesc_can_count >0 && plus2_NMD_result[[i]][["can_region"]]>0){
    plus2_NMD_result[[i]][["can_p"]] = binom.test(NMDesc_can_count,plus2_NMD_result[[i]][["can_region"]],all.variant/tr_length[1],alternative='less')$p.value
  }else(plus2_NMD_result[[i]][["can_p"]] = NA)
 }


#7. use 90% as cutoff, get the clinvar frameshift gene list that is NMDesc variants enriched
can_list = c()

for (i in 1:length(plus2_NMD_result)) {
  if (!is.na(plus2_NMD_result[[i]][["can_p"]]) && plus2_NMD_result[[i]][["can_p"]] >= 0.9) {
    can_list <- c(can_list, names(plus2_NMD_result)[i])
  }
}

##transform transcript to gene
can_gene = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values = can_list,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
writeLines(c("hgnc_symbol", unique(can_gene$hgnc_symbol)), 'plus2_can_gene0217.txt')