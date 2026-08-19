get_fs_idr_match = function(input_path="~/Downloads/plus1_disorder_predictions.csv",type='plus2'){
  idr_loc = read.csv(input_path)
  idr_loc = idr_loc %>%
    filter(!is.na(Disorder_Domains) & Disorder_Domains != "[]") %>%  # Remove empty disorder domains
    mutate(C_end_IDR = sapply(Disorder_Domains, get_largest_number),
           C_start_IDR = sapply(Disorder_Domains, get_second_largest_number))
  
  # use key get transcript id
  idr_loc2 <- merge(idr_loc, key_to_transcript, by = "Variant_Key", all.x = TRUE)

  # Apply function to each transcript group
  idr_NMD <-  idr_loc2 %>%
    group_by(Transcript_ID) %>%
    group_split() %>%
    lapply(function(df) calculate_fs_aft_NMD(df$Transcript_ID[1], type=type)) %>%
    bind_rows() %>%
    left_join(idr_loc2, by = "Transcript_ID") 
  
   #keep one row for each gene
  idr_NMD$match_start = pmax(idr_NMD$NMDesc.start,idr_NMD$C_start_IDR*3)
  idr_NMD$match_end = pmin(idr_NMD$NMDesc.end,idr_NMD$C_end_IDR*3)
  idr_NMD$match = (idr_NMD$match_end - idr_NMD$match_start) >= 20*3

  return(data.frame(matched_idr = sum(idr_NMD$match,na.rm=T), total_variants = nrow(idr_NMD)))
}
