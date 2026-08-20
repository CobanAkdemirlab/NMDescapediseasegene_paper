library(biomaRt)
library(dplyr)

snv_gnomAD_variants <- fs_gnomAD_variants %>%
  filter(fs_type == 'snv')

# snv CONTROL
snv_control_gene <- snv_control_pli$gene

snv_control_tr <- getBM(
  attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot', 'transcript_is_canonical'),
  filters = 'hgnc_symbol',
  values = snv_control_gene,
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)

snv_control_tr_ids <- snv_control_tr$ensembl_transcript_id

snv_control_variants <- snv_gnomAD_variants %>%
  filter(transcript %in% snv_control_tr_ids) %>%
  left_join(
    snv_control_tr %>% select(ensembl_transcript_id, uniprotswissprot),
    by = c("transcript" = "ensembl_transcript_id")
  )

write.csv(snv_control_variants, 'snv_control_gnomAD_variants0406.csv', row.names = FALSE)

