get_snv_idr_match2 = function(input_path="/Users/jxu14/Desktop/idr/output/snv_disorder_predictions.csv"){
  idr_loc = read.csv(input_path)
  idr_loc = idr_loc %>%
    filter(!is.na(Disorder_Domains) & Disorder_Domains != "[]") %>%  # Remove empty disorder domains
    mutate(C_end_IDR = sapply(Disorder_Domains, get_largest_number),
           C_start_IDR = sapply(Disorder_Domains, get_second_largest_number))
  key_to_transcript_snv = data.frame(
    Variant_Key = snv@elementMetadata@listData[["id"]],
    Transcript_ID = snv@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]],
    stringsAsFactors = FALSE
  )
  key_to_transcript_snv = data.frame(
    Variant_Key = snv_NMD_result2$id,
    Transcript_ID = snv_NMD_result2$transcript,
    stringsAsFactors = FALSE
  )
      
  # use key get transcript id
  idr_loc2 <- merge(idr_loc, key_to_transcript_snv, by = "Variant_Key", all.x = TRUE)
  #from transcript id, get cds information
  tr2cds = getBM(attributes = c('ensembl_transcript_id','cds_start','cds_end'),
                 filters = 'ensembl_transcript_id',
                 values = idr_loc2$Transcript_ID,
                 mart = ensembl)
  tr2cds <- tr2cds[!is.na(tr2cds$cds_start),]
  
  aft_NMD_results <- tr2cds %>%
    group_by(ensembl_transcript_id) %>%
    group_split() %>%
    lapply(calculate_snv_aft_NMD, 50) %>%
    bind_rows() %>%
    left_join(tr2cds, by = "ensembl_transcript_id")  %>%
    dplyr::select(ensembl_transcript_id, aft.NMD.ind, cds_length) %>%
    distinct(ensembl_transcript_id, .keep_all = TRUE) 
  
  #merge with idr_loc
  idr_NMD = merge(aft_NMD_results,idr_loc2,by.x = 'ensembl_transcript_id',by.y = 'Transcript_ID')
  #match NMD and idr
  idr_NMD$match_start = (idr_NMD$aft.NMD.ind - idr_NMD$C_start_IDR*3) >= 20*3
  idr_NMD$match_end = idr_NMD$C_end_IDR*3 > idr_NMD$aft.NMD.ind
  idr_NMD$match = idr_NMD$match_start & idr_NMD$match_end
  #return the count of matched idr and total gene count
  return(data.frame(matched_idr = sum(idr_NMD$match,na.rm=T), total_gene = nrow(idr_NMD)))
  }
