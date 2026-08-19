snv_control_tr2gene = getBM(attributes = c('hgnc_symbol','ensembl_transcript_id','transcript_is_canonical'),
                    filters = 'ensembl_transcript_id',
                    values = snv_control_idr$Transcript_ID,
                    mart = ensembl)
snv_control_idr <- snv_control_idr %>%
  select(-gene) %>%  # Remove the old gene column
  left_join(snv_control_tr2gene, by = c("Transcript_ID" = "ensembl_transcript_id")) %>%
  rename(gene = hgnc_symbol)
#from snv_control_uniprot_ids_out(Uniprot_ID) get the gene 
snv_control_uniprot_ids_out$gene = getBM(attributes = c('hgnc_symbol','uniprotswissprot'),
                                 filters = 'uniprotswissprot',
                                 values = snv_control_uniprot_ids_out$UniProt_ID,
                                 mart = ensembl)$hgnc_symbol
# join snv_control_uniprot_ids_out to snv_control_idr by gene, the Disordered_Domain_Boundaries in snv_control_uniprot_ids_out
# will be used to be the wildtype_Disordered_Domain_Boundaries in snv_control_idr
snv_control_idr$wild_Disordered_Domain_Boundaries = snv_control_uniprot_ids_out$Disordered_Domain_Boundaries[match(snv_control_idr$gene, snv_control_uniprot_ids_out$gene)]
snv_control_idr <- snv_control_idr %>%
  mutate(
    wild_c_idr_end = sapply(wild_Disordered_Domain_Boundaries, get_largest_number),
    wild_c_idr_start = sapply(wild_Disordered_Domain_Boundaries, get_second_largest_number)
  )

snv_control_cds_lengths$Transcript_ID = snv_control_cds_lengths$ensembl_transcript_id

snv_control_idr <- snv_control_idr %>%
  left_join(snv_control_cds_lengths %>% select(Transcript_ID, cds_length),
            by = "Transcript_ID")

snv_control_idr <- snv_control_idr %>%
  left_join(snv_control_nmdesc_lengths %>% select(gene, nmdesc_length),
            by = "gene") 