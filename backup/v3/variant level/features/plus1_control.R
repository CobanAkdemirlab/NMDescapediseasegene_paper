# 1. Read + map
plus1_control_disorder_predictions <- read_csv("~/Downloads/idr/plus1_control_disorder_predictions.csv")

plus1_control_key_to_transcript <- data.frame(
  Variant_Key = plus1_gnomAD_variants$id,
  Transcript_ID = plus1_gnomAD_variants$transcript,
  stringsAsFactors = FALSE
)

plus1_control_key_to_transcript <- plus1_control_key_to_transcript %>%
  filter(Transcript_ID %in% plus1_control_tr$ensembl_transcript_id)

plus1_control_disorder_predictions_mapped <- merge(
  plus1_control_disorder_predictions,
  plus1_control_key_to_transcript,
  by = "Variant_Key",
  all.x = TRUE
)

# 2. Get exon info
plus1_control_tx_list <- unique(plus1_control_key_to_transcript$Transcript_ID)

plus1_control_exon_info <- getBM(
  attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "rank", "strand"),
  filters = "ensembl_transcript_id",
  values = plus1_control_tx_list,
  mart = ensembl
)

# 3. Get CDS sequences
plus1_control_cds_info <- getBM(
  attributes = c("ensembl_transcript_id", "coding"),
  filters = "ensembl_transcript_id",
  values = plus1_control_tx_list,
  mart = ensembl
)

# 4. Loop through transcripts to compute NMDesc region
plus1_control_nmdesc_coords <- data.frame(
  Transcript_ID = character(),
  NMDesc.start = numeric(),
  NMDesc.end = numeric(),
  stringsAsFactors = FALSE
)

for (tx in plus1_control_tx_list) {
  cds_entry <- plus1_control_cds_info %>%
    filter(ensembl_transcript_id == tx)
  
  exon_df <- plus1_control_exon_info %>%
    filter(ensembl_transcript_id == tx) %>%
    filter(!is.na(cds_start) & !is.na(cds_end)) %>%
    arrange(rank)
  
  if (nrow(cds_entry) == 0 || nrow(exon_df) < 2) {
    next
  }
  
  input_seq <- cds_entry$coding
  
  tryCatch({
    PTC_ind <- get_PTC_plus_ind(input_seq, type = "plus1")
    PTC_status <- get_NMDesc_PTC(input_seq, bp = 50, variant_transcript = tx, type = "plus1")
    region_set <- get_NMDesc_PTC_region(PTC_ind, tx, PTC_status)
    
    cds_length <- sum((exon_df$cds_end - exon_df$cds_start) + 1)
    
    nmdesc_row <- data.frame(
      Transcript_ID = tx,
      NMDesc.start = cds_length - as.numeric(region_set$can_region) + 1,
      NMDesc.end = cds_length
    )
    
    plus1_control_nmdesc_coords <- rbind(plus1_control_nmdesc_coords, nmdesc_row)
  }, error = function(e) {
    message(sprintf("Skipped %s due to error: %s", tx, e$message))
  })
}

# 5. Join into disorder predictions
plus1_control_idr <- plus1_control_disorder_predictions_mapped %>%
  left_join(plus1_control_nmdesc_coords, by = "Transcript_ID")

# 6. Process overlaps
plus1_control_idr <- plus1_control_idr %>%
  filter(!is.na(Disorder_Domains) & Disorder_Domains != "[]") %>%
  mutate(
    C_end_IDR = sapply(Disorder_Domains, get_largest_number),
    C_start_IDR = sapply(Disorder_Domains, get_second_largest_number)
  )

plus1_control_idr$match_start <- pmax(plus1_control_idr$NMDesc.start, plus1_control_idr$C_start_IDR * 3)
plus1_control_idr$match_end <- pmin(plus1_control_idr$NMDesc.end, plus1_control_idr$C_end_IDR * 3)
plus1_control_idr$match <- (plus1_control_idr$match_end - plus1_control_idr$match_start) >= 20 * 3
plus1_control_idr %>% group_by(Transcript_ID) %>%
  summarise(
    match = any(match),
    .groups = "drop"
  ) -> plus1_control_idr3

# add ptc_loc ,variant_loc, wild_type_idr
#change percentage figure, 
#add ptc loc
ptc_loc = PTC_info %>%
  filter(type == "plus1") %>%
  dplyr::select(Transcript_ID = transcript, PTC_loc)

#add variant_loc
plus1_control_idr <- plus1_control_idr %>%
  left_join(PTC_info, by = c("Transcript_ID" = "transcript")) %>%
  mutate(
    ptc_loc = as.numeric(str_extract(PTC_loc, "\\d+")),
    variant_loc = as.numeric(str_extract(Variant_Loc, "\\d+")),
    wild_type_idr = ifelse(match, "Yes", "No")
  )

# View results
table(plus1_control_idr$match, useNA = "always")
