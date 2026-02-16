library(readr)
plus1_control_uniprot_ids_out <- read_delim("~/Downloads/idr0512/plus1_control_uniprot_ids_out.txt", 
                                                   delim = "\t", escape_double = FALSE, 
                                                   trim_ws = TRUE)
minus1_control_uniprot_ids_out <- read_delim("~/Downloads/idr0512/minus1_control_uniprot_ids_out.txt", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)
plus1_uniprot_ids_out <- read_delim("~/Downloads/idr0512/plus1_uniprot_ids_out.txt", 
                                                   delim = "\t", escape_double = FALSE, 
                                                   trim_ws = TRUE)
minus1_uniprot_ids_out <- read_delim("~/Downloads/idr0512/minus1_uniprot_ids_out.txt",
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)
snv_control_uniprot_ids_out <- read_delim("~/Downloads/idr0512/snv_control_uniprot_ids_out.txt",
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)
snv_uniprot_ids_out <- read_delim("~/Downloads/idr0512/snv_uniprot_ids_out.txt",
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)
#merge snv_uniprot_ids_out(Uniprot_ID) with snv_idr(Transcript ID) by snv_can_uni(hgnc_symbol), the purpose is to 
#add wildtype idr_loc to snv_idr
snv_tr2gene = getBM(attributes = c('hgnc_symbol','ensembl_transcript_id','transcript_is_canonical'),
                                    filters = 'ensembl_transcript_id',
                                    values = snv_idr$Transcript_ID,
                                    mart = ensembl)
snv_idr <- snv_idr %>%
  select(-gene) %>%  # Remove the old gene column
  left_join(snv_tr2gene, by = c("Transcript_ID" = "ensembl_transcript_id")) %>%
  rename(gene = hgnc_symbol)
#from snv_uniprot_ids_out(Uniprot_ID) get the gene 
snv_uniprot_ids_out$gene = getBM(attributes = c('hgnc_symbol','uniprotswissprot'),
                                    filters = 'uniprotswissprot',
                                    values = snv_uniprot_ids_out$UniProt_ID,
                                    mart = ensembl)$hgnc_symbol
#mjoin snv_uniprot_ids_out to snv_idr by gene, the Disordered_Domain_Boundaries in snv_uniprot_ids_out
# will be used to be the wildtype_Disordered_Domain_Boundaries in snv_idr
snv_idr$wild_Disordered_Domain_Boundaries = snv_uniprot_ids_out$Disordered_Domain_Boundaries[match(snv_idr$gene, snv_uniprot_ids_out$gene)]
snv_idr <- snv_idr %>%
  mutate(
    wild_c_idr_end = sapply(wild_Disordered_Domain_Boundaries, get_largest_number),
    wild_c_idr_start = sapply(wild_Disordered_Domain_Boundaries, get_second_largest_number)
  )





