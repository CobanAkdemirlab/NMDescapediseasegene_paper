# Assign control gene list
snv_control_gene <- snv_control_pli$gene

# Get canonical transcripts for the control gene list using biomaRt
snv_control_tr <- getBM(
  attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'transcript_is_canonical'),
  filters = 'hgnc_symbol',
  values = snv_control_gene,
  mart = ensembl
) %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::pull(ensembl_transcript_id)

# Subset snv_gnomAD_variants to those matching the canonical control transcripts
snv_control_variants <- snv_gnomAD_variants[snv_gnomAD_variants$transcript %in% snv_control_tr, ]
#add uniprot id
snv_control_tr <- getBM(
  attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot', 'transcript_is_canonical'),
  filters = 'hgnc_symbol',
  values = snv_control_gene,
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)

# Pull just the transcript IDs for filtering
snv_control_tr_ids <- snv_control_tr$ensembl_transcript_id

# Subset snv_gnomAD_variants to those matching the canonical control transcripts
snv_control_variants <- snv_gnomAD_variants %>%
  filter(transcript %in% snv_control_tr_ids)

# Merge UniProt IDs into the variant table
snv_control_variants <- snv_control_variants %>%
  left_join(
    snv_control_tr %>% select(ensembl_transcript_id, uniprotswissprot),
    by = c("transcript" = "ensembl_transcript_id")
  )

# Write the filtered variant dataframe to CSV
write.csv(snv_control_variants, 'snv_gnomAD_variants0406.csv', row.names = FALSE)
